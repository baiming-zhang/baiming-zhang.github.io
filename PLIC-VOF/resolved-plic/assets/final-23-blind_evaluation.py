"""New blinded accuracy only, generated strictly after deployment selection.

Default entry is called by autofinish_resolved.py after the deployment lock.
--check-gate-only never generates data or initializes CUDA. --wait may wait for
that lock; without --wait an absent lock returns a truthful waiting status.
"""
from pathlib import Path
import os,sys
sys.dont_write_bytecode=True
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[name]='1'
import argparse,csv,ctypes,datetime,fcntl,gc,hashlib,importlib.util,json,time,traceback
import numpy as np
R=Path(__file__).resolve().parent;W=R.parent;B=R/'blind_test'
if R.name=='blind_test':R=R.parent;W=R.parent;B=R/'blind_test'
sys.path.insert(0,str(R));sys.path.insert(0,str(W/'operator-ablation-20261122/shared'))
from gpu_lease import gpu_lease
CIRCLE_SEED=202611270398;PHYSICAL_SEED=202611270399;CIRCLE_COUNT=32768
NS=[32,64,128,256,512,1024,2048,4096,8192]
SEEDS=[4101,4102,4103];START=time.monotonic()

def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def save(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(obj,indent=2,allow_nan=False))
    for j in range(40):
        try:temp.replace(path);return
        except PermissionError:
            if j==39:raise
            time.sleep(.05)
def status(stage,**extra):
    record=dict(state=stage,pid=os.getpid(),updated_utc=now(),elapsed_seconds=time.monotonic()-START,formal_timing=False);record.update(extra);save(B/'status.json',record)
def pathof(p):
    p=Path(p);return p if p.is_absolute() else W/p
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def atomic_npz(path,**arrays):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.tmp.npz');np.savez_compressed(tmp,**arrays);tmp.replace(path)
def csvout(path,rows):
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8-sig') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
def assert_gate(identity=None):
    lock=R/'deployment_selection_locked.json'
    if not lock.exists():raise RuntimeError('Blind geometry forbidden before deployment_selection_locked.json')
    if identity and sha(lock)!=identity['deployment_lock_sha256']:raise RuntimeError('Deployment lock changed after blind registration; preserve all results, do not generate a new test.')

def locked_plan():
    assert_gate();lock=read(R/'deployment_selection_locked.json');manifest_path=pathof(lock.get('model_manifest_path',str((R/'inference/model_manifest.json').relative_to(W))));manifest=read(manifest_path);models=manifest['models']
    assert manifest.get('complete') and len(models)==24 and all(m['status']=='ready' for m in models)
    assert set(lock['model_ids'])=={m['model_id'] for m in models} and lock['all_models_retained']
    assert set(lock['model_runtimes'])==set(lock['model_ids'])
    for m in models:
        kind=lock['model_runtimes'][m['model_id']];assert kind in ['onnx_cuda','torchscript_cuda'],kind
        key='onnx' if kind=='onnx_cuda' else 'torchscript';assert m[key+'_status']=='ready'
        assert sha(pathof(m[key+'_path']))==m[key+'_sha256'],m['model_id']
    identity=dict(deployment_lock_sha256=sha(R/'deployment_selection_locked.json'),export_manifest_path=str(manifest_path.relative_to(W)),export_manifest_sha256=sha(manifest_path),circle_seed=CIRCLE_SEED,circle_count=CIRCLE_COUNT,physical_seed=PHYSICAL_SEED,resolutions=NS,
        models=[dict(model_id=m['model_id'],runtime=lock['model_runtimes'][m['model_id']],artifact_sha256=m[('onnx' if lock['model_runtimes'][m['model_id']]=='onnx_cuda' else 'torchscript')+'_sha256']) for m in models])
    p=B/'registration.json'
    if p.exists():assert read(p)['identity']==identity,'Blind registration identity changed; do not silently reuse/resample.'
    else:save(p,dict(identity=identity,registered_utc=now(),script_sha256=sha(__file__),no_prelock_generation=True,selection_source='Validation accuracy/scout latency only;all24/allseeds retained;no blind selection',truth='Complete geometry cell/strip area;raw7f32VOF(+CFL) only as operator input',speed_evidence=False))
    return lock,models,identity

def radius_bins(radius):
    return [('all',np.ones(len(radius),bool)),('4_to_8',(radius>=4)&(radius<8)),('8_to_32',(radius>=8)&(radius<32)),('32_to_128',(radius>=32)&(radius<128)),('128_to_512',(radius>=128)&(radius<=512)),('above_512',radius>512)]

def circle_data(identity):
    assert_gate(identity);dest=B/'datasets/circle.npz';record=B/'datasets/circle.json'
    if dest.exists() and record.exists():assert sha(dest)==read(record)['sha256'];return read(record)
    import torch
    G=module('_blind_circle_generator',R/'data_generation.py');G.DATA=B/'circle';G.DATA.mkdir(parents=True,exist_ok=True)
    G.status=lambda **kw:status('circle_'+str(kw.pop('state','progress')),**kw)
    save(B/'circle/generation_protocol.json',dict(source_generator_sha256=sha(R/'data_generation.py'),count=CIRCLE_COUNT,seed=CIRCLE_SEED,id_start=3000000000,radius_cells=[4,512],new_after_deployment_lock=True,source_config=G.CONFIG))
    with gpu_lease('resolved-blind-circle-generation',formal_timing=False):
        report=G.generate('test',CIRCLE_COUNT,CIRCLE_SEED,3000000000);torch.cuda.synchronize()
    source=G.DATA/'test.npz';assert sha(source)==report['sha256'];d=dict(np.load(source));keys=G.geometry_keys(d)
    assert len(np.unique(keys))==CIRCLE_COUNT
    overlaps={}
    for split in ['train','validation']:
        with np.load(R/'data'/f'{split}.npz') as z:other={k:z[k] for k in ['cx','cy','radius_cells']}
        overlaps[split]=int(np.intersect1d(keys,G.geometry_keys(other)).size);assert overlaps[split]==0,split
    n=len(d['flux']);d.update(q=np.where(d['complement'],-1.,1.)/d['radius_cells'],flux_ratio=d['swept_volume']/d['center_fraction'],N=np.full(n,-1,dtype='i8'),h=np.ones(n),family=np.zeros(n,dtype='i4'),shape=np.full(n,'circle',dtype='U10'),geometry_id=d['base_id'].copy(),query_id=np.arange(n,dtype='i8'),sample_id=np.arange(n,dtype='i8'),transform=2*d['complement'].astype('i4'))
    assert np.all(abs(d['q'])<=.25+1e-14)
    atomic_npz(dest,**d)
    rec=dict(dataset='circle',path=str(dest.relative_to(W)),sha256=sha(dest),samples=n,new_seed=CIRCLE_SEED,parent_geometry_count=CIRCLE_COUNT,split_overlap=overlaps,source_generation=report,source_raw_path=str(source.relative_to(W)),N_convention='N=-1,h=1:cell-normalized radius sweep,not a physical refinement curve',radius_filter_removed=0,created_utc=now())
    save(record,rec);return rec

def geometry_helpers():
    folder=W/'operator-ablation-20261122/01-convergence';sys.path.insert(0,str(folder))
    I=module('_blind_physical_geometry_helpers',folder/'independent.py');sys.path.insert(0,str(R))
    return I

def full_area(G,I,p,xlo,xhi,ylo,yhi,h,sx,sy):
    """Complete closed ellipse/circle; sine is an unbounded single graph."""
    import torch
    n=len(p);a=p[:,1];b=p[:,2];t=p[:,3];cx=-a*torch.cos(t)/h-sx;cy=-b*torch.sin(t)/h-sy;A=a/h;BB=b/h
    xl=torch.as_tensor(xlo,device=p.device,dtype=p.dtype).expand(n,-1);xh=torch.as_tensor(xhi,device=p.device,dtype=p.dtype).expand(n,-1)
    yl=torch.as_tensor(ylo,device=p.device,dtype=p.dtype).expand(n,-1);yh=torch.as_tensor(yhi,device=p.device,dtype=p.dtype).expand(n,-1)
    zero=torch.zeros(n,device=p.device,dtype=p.dtype);one=torch.ones_like(zero)
    ellipse=G.disk_rectangle(zero,zero,one,(xl-cx[:,None])/A[:,None],(xh-cx[:,None])/A[:,None],(yl-cy[:,None])/BB[:,None],(yh-cy[:,None])/BB[:,None])*A[:,None]*BB[:,None]
    sine=I.area(p,xl,xh,yl,yh,h,sx,sy)
    return torch.where((p[:,0]<2)[:,None],ellipse,sine).clamp_min(0).minimum((xh-xl)*(yh-yl))

def physical_data(identity,resolved_policy='common_across_N'):
    assert_gate(identity);dest=B/'datasets/physical.npz';report=B/'datasets/physical.json'
    if dest.exists() and report.exists():assert sha(dest)==read(report)['sha256'];return read(report)
    import torch
    G=module('_blind_full_disk_helpers',R/'data_generation.py');I=geometry_helpers();rng=np.random.default_rng(PHYSICAL_SEED);P=[]
    for family in range(3):
        for _ in range(16):
            if family==0:a=float(rng.uniform(.24,.8));b=a;t=float(np.pi/2+rng.uniform(-.25,.25));linear=A=w=phi=0.
            elif family==1:a=float(rng.uniform(.45,.9));b=float(rng.uniform(.15,.4));t=float(np.pi/2+rng.uniform(-.38,.38));linear=A=w=phi=0.
            else:a=b=1.;t=np.pi/2;linear=float(rng.choice([-1,1])*rng.uniform(.1,.8));A=float(rng.uniform(.0003,.0012));w=float(rng.uniform(3,6));phi=float(rng.uniform(.3,2.8))
            P.append([family,a,b,t,linear,A,w,phi])
    P=np.asarray(P);offset=np.concatenate((np.zeros((1,2)),rng.uniform(-.4,.4,(8,2))))
    specs=np.asarray([(g,o,t) for g in range(48) for o in range(9) for t in range(4)]);g,oi,tr=specs.T;n=len(g)
    # Same Courant for phase/transpose variants and at every N; fixed before area integration.
    cf=rng.random(48*9).astype('f4');mode=rng.random(len(cf));low=mode<.10;cf[low]=np.exp(rng.uniform(np.log(1e-4),np.log(.05),int(low.sum()))).astype('f4');cf[(mode>=.10)&(mode<.12)]=0;cf[(mode>=.12)&(mode<.14)]=1;cf=np.repeat(cf,4)
    old_paths=[W/'curvature-neural-response-20261126/test/physical.npz',W/'operator-ablation-20261122/01-convergence/independent_fixtures.npz'];overlap={}
    pk=np.ascontiguousarray(P).view(np.dtype((np.void,64))).reshape(-1)
    for old in old_paths:
        if old.exists():
            pp=np.load(old)['physical_parameters'];ok=np.ascontiguousarray(pp).view(np.dtype((np.void,64))).reshape(-1);overlap[str(old.relative_to(W))]=int(np.intersect1d(pk,ok).size);assert not overlap[str(old.relative_to(W))]
    yy,xx=np.meshgrid(np.arange(3,-4,-1),np.arange(-3,4),indexing='ij');xl=xx.reshape(1,-1)-.5;yl=yy.reshape(1,-1)-.5;blocks=[]
    for ni,N in enumerate(NS):
        part=B/'physical_parts'/f'N{N}.npz';meta=part.with_suffix('.json')
        if part.exists() and meta.exists():assert read(meta)['seed']==PHYSICAL_SEED and sha(part)==read(meta)['sha256'];blocks.append(dict(np.load(part)));continue
        xs=[];qs=[];fv=[];vx=[];vy=[]
        status('generating_physical',N=N,completed_resolutions=ni,total_resolutions=len(NS))
        with gpu_lease('resolved-blind-physical-'+str(N),formal_timing=False):
            for start in range(0,n,128):
                sl=slice(start,start+128);p=torch.as_tensor(P[g[sl]],device='cuda',dtype=torch.float64);h=torch.full((len(p),),1/N,device='cuda',dtype=torch.float64);sx=torch.as_tensor(offset[oi[sl],0],device='cuda');sy=torch.as_tensor(offset[oi[sl],1],device='cuda');c=torch.as_tensor(cf[sl].astype('f8'),device='cuda')
                x=full_area(G,I,p,xl,xl+1,yl,yl+1,h,sx,sy).reshape(-1,7,7)
                f=full_area(G,I,p,[[-.5]],[[.5]],[[-.5]],[[.5]],h,sx,sy)[:,0]
                right=full_area(G,I,p,(.5-c)[:,None],torch.full_like(c,.5)[:,None],torch.full_like(c,-.5)[:,None],torch.full_like(c,.5)[:,None],h,sx,sy)[:,0]
                top=full_area(G,I,p,torch.full_like(c,-.5)[:,None],torch.full_like(c,.5)[:,None],(.5-c)[:,None],torch.full_like(c,.5)[:,None],h,sx,sy)[:,0]
                q=I.curvature(p,h,sx,sy);q=torch.where(p[:,0]==0,h/p[:,1],q)
                xs.append(x.cpu().numpy());fv.append(f.cpu().numpy());vx.append(right.cpu().numpy());vy.append(top.cpu().numpy());qs.append(q.cpu().numpy())
            torch.cuda.synchronize()
        x=np.concatenate(xs);f=np.concatenate(fv);v=np.concatenate(vx);top=np.concatenate(vy);q=np.concatenate(qs);transpose=(tr&1)>0;comp=(tr&2)>0
        x[transpose]=x[transpose].transpose(0,2,1)[:,::-1,::-1];v[transpose]=top[transpose];x[comp]=1-x[comp];f[comp]=1-f[comp];v[comp]=cf[comp].astype('f8')-v[comp];q[comp]*=-1
        # Exact geometry limits only in truth, never in numerical/native predictions.
        v=np.where(cf==0,0.,np.where(cf==1,f,v));assert np.isfinite(x).all() and (f>0).all() and (f<1).all()
        block=dict(stencils=x.astype('f4'),stencils64=x,courants=cf.copy(),flux=(v/f).astype('f4'),curvature=q.astype('f4'),q=q,flux_ratio=v/f,center_fraction=f,swept_volume=v,N=np.full(n,N,dtype='i8'),h=np.full(n,1/N),family=P[g,0].astype('i4'),shape=np.array(['circle','ellipse','weak_sine'],dtype='U10')[P[g,0].astype('i4')],geometry_id=4000000000+g,query_id=ni*n+np.arange(n),sample_id=np.arange(n),physical_id=g,offset_id=oi,transform=tr,radius_cells=np.divide(1,abs(q),out=np.full(n,np.inf),where=q!=0))
        atomic_npz(part,**block);save(meta,dict(seed=PHYSICAL_SEED,N=N,samples=n,sha256=sha(part)));blocks.append(block)
    d={k:np.concatenate([z[k] for z in blocks]) for k in blocks[0]};resolved=(abs(d['q']).reshape(len(NS),n)<=.25+1e-14).all(0)
    if resolved_policy=='common_across_N':keep=np.tile(resolved,len(NS))
    elif resolved_policy=='preserve_all':keep=np.ones(len(d['q']),dtype=bool)
    else:raise ValueError(resolved_policy)
    removed=int((~keep).sum());d={k:v[keep] for k,v in d.items()}
    d.update(physical_parameters=P,offsets=offset)
    # Independent CPU full-ellipse areas / adaptive graph quadrature and nearest-distance curvature.
    from neuralvof.geometry import EllipseGeometry,geometry_rectangle_area
    V=module('_blind_cpu_nearest_oracle',W/'operator-ablation-20261122/01-convergence/verify_curvature.py');indices=[]
    for N in NS:
        for fam in range(3):
            for transform in range(4):indices.append(int(np.flatnonzero((d['N']==N)&(d['family']==fam)&(d['transform']==transform))[0]))
    max_cell=max_ratio=max_q=0.
    for checked,i in enumerate(indices):
        pp=P[d['physical_id'][i]];sx,sy=offset[d['offset_id'][i]];h=d['h'][i];tt=int(d['transform'][i]);c=float(d['courants'][i]);oracle=np.empty((7,7))
        if pp[0]<2:
            gg=EllipseGeometry(center=np.array([-pp[1]*np.cos(pp[3])/h-sx,-pp[2]*np.sin(pp[3])/h-sy]),semi_axes=pp[1:3]/h,angle=0.,fluid_inside=True,local_normal=np.zeros(2),local_curvature=0.)
            area=lambda xa,xb,ya,yb:geometry_rectangle_area(gg,xa,xb,ya,yb)
        else:area=lambda xa,xb,ya,yb:I.cpu_oracle(pp,h,sx,sy,xa,xb,ya,yb)
        for iy,Y in enumerate(range(3,-4,-1)):
            for ix,X in enumerate(range(-3,4)):oracle[iy,ix]=area(X-.5,X+.5,Y-.5,Y+.5)
        vol=area(-.5,.5,.5-c,.5) if tt&1 else area(.5-c,.5,-.5,.5)
        if tt&1:oracle=oracle.T[::-1,::-1]
        if tt&2:oracle=1-oracle;vol=c-vol
        ff=oracle[3,3];ratio=0. if c==0 else (1. if c==1 else vol/ff);qq=V.cpu_nearest(pp,h,sx,sy)
        if tt&2:qq=-qq
        max_cell=max(max_cell,float(abs(oracle.astype('f4')-d['stencils'][i]).max()));max_ratio=max(max_ratio,abs(ratio-d['flux_ratio'][i]));max_q=max(max_q,abs(qq-d['q'][i]))
        if checked%16==0:status('physical_independent_truth_check',checked=checked+1,total=len(indices))
    assert max_cell<2e-7 and max_ratio<5e-6 and max_q<1e-9,(max_cell,max_ratio,max_q)
    circle_mask=d['family']==0;pp=P[d['physical_id'][circle_mask]];off=offset[d['offset_id'][circle_mask]];hh=d['h'][circle_mask]
    cx=-pp[:,1]*np.cos(pp[:,3])/hh-off[:,0];cy=-pp[:,2]*np.sin(pp[:,3])/hh-off[:,1];tt=d['transform'][circle_mask];swap=(tt&1)>0;oldcx=cx.copy();cx[swap]=cy[swap];cy[swap]=oldcx[swap]
    circle_keys=G.geometry_keys(dict(cx=cx,cy=cy,radius_cells=pp[:,1]/hh));circle_overlap={}
    for split in ['train','validation']:
        with np.load(R/'data'/f'{split}.npz') as zz:other={key:zz[key] for key in ['cx','cy','radius_cells']}
        circle_overlap[split]=int(np.intersect1d(circle_keys,G.geometry_keys(other)).size);assert circle_overlap[split]==0
    atomic_npz(dest,**d)
    rec=dict(physical_circle_train_validation_overlap=circle_overlap,dataset='physical',path=str(dest.relative_to(W)),sha256=sha(dest),seed=PHYSICAL_SEED,samples=len(d['q']),samples_per_N=int(keep.sum()//len(NS)),resolutions=NS,physical_shapes=48,offsets=9,phase_transpose_variants=4,parent_geometry_overlap_with_prior_physical_tests=overlap,
        radius_filter=resolved_policy,radius_filter_removed_queries=removed,truth='Complete closed circle/ellipse cell and right/top strip area; sine graph integral;nearest in-cell signed curvature. Transpose uses old top strip, phase complement changes area and q sign.',
        oracle=dict(geometries=len(indices),cell_integrals=49*len(indices),fraction_max_abs=max_cell,flux_ratio_max_abs=max_ratio,q_max_abs=max_q),created_utc=now())
    save(report,rec);return rec

def native_runners():
    NV=module('_blind_native_validation_interface',R/'native_validation.py')
    library=R/'native_validation/libaphros_production.so'
    if library.exists():
        lib=ctypes.CDLL(str(library));fn=lib.aphros_native_eval;fn.argtypes=[ctypes.c_int,ctypes.c_void_p,ctypes.c_int,ctypes.c_int64,ctypes.c_void_p,ctypes.c_void_p,ctypes.c_void_p,ctypes.c_void_p]
        def prod(x,c,out,res=None,it=None):return fn(0,x.ctypes.data,x.shape[1],len(x),c.ctypes.data,out.ctypes.data,None,None)
        build=read(R/'native_validation/production_build.json');assert sha(library)==build['library_sha256']
    else:
        NV.OUT=B/'native_production_build';NV.OUT.mkdir(parents=True,exist_ok=True);prod,build=NV.production_driver();library=NV.OUT/'libaphros_production.so'
    native=NV.NATIVE
    validation=read(R/'native_validation/summary.json');assert validation['state']=='complete'
    assert sha(library)==validation['identity']['production_library_sha256']
    assert all(sha(native/f'build/lib{backend}_local.so')==validation['identity']['libraries'][backend] for backend in ['aphros','simliquid'])
    return [dict(model_id='native_aphros_production_flux',operator='flux',backend='aphros',runtime='native_cpu_f64',library_sha256=sha(library),run=prod,
                 operator_boundary='raw f32 VOF promoted to f64;original normal+PLIC reconstruction+production flux q=CFL,dt=1;defined at zero CFL without output replacement'),
            dict(model_id='native_simliquid_flux',operator='flux',backend='simliquid',runtime='native_cpu_f64',library_sha256=sha(native/'build/libsimliquid_local.so'),run=NV.native_function('simliquid','flux'),operator_boundary='raw VOF;original native normal+PLIC flux'),
            dict(model_id='native_aphros_particles',operator='curvature',backend='aphros',runtime='native_cpu_f64',library_sha256=sha(native/'build/libaphros_local.so'),run=NV.native_function('aphros','curvature'),operator_boundary='raw VOF;original normal+neighbor PLIC fragments+particle string;not HF and not relabelled PLIC curvature'),
            dict(model_id='native_simliquid_particles',operator='curvature',backend='simliquid',runtime='native_cpu_f64',library_sha256=sha(native/'build/libsimliquid_local.so'),run=NV.native_function('simliquid','curvature'),operator_boundary='raw VOF;original native particle curvature;not HF')]

def infer_neural(model,kind,d):
    import torch
    k=int(model['stencil_size']);a=(7-k)//2;x=np.ascontiguousarray(d['stencils'][:,a:a+k,a:a+k]);c=np.ascontiguousarray(d['courants'],dtype='f4');chunks=[]
    with gpu_lease('resolved-blind-'+model['model_id'],formal_timing=False):
        status('neural_blind_accuracy',model_id=model['model_id'],dataset_samples=len(x),runtime=kind)
        if kind=='onnx_cuda':
            import onnxruntime as ort
            opts=ort.SessionOptions();opts.intra_op_num_threads=1;opts.inter_op_num_threads=1;opts.graph_optimization_level=ort.GraphOptimizationLevel.ORT_ENABLE_ALL;opts.execution_mode=ort.ExecutionMode.ORT_SEQUENTIAL
            limit=int(torch.cuda.get_device_properties(0).total_memory*.05)
            engine=ort.InferenceSession(str(pathof(model['onnx_path'])),sess_options=opts,providers=[('CUDAExecutionProvider',{'use_tf32':'0','gpu_mem_limit':str(limit),'arena_extend_strategy':'kSameAsRequested'})])
            assert engine.get_providers()[0]=='CUDAExecutionProvider';assert len(engine.get_outputs())==1
            for first in range(0,len(x),2048):
                xx=x[first:first+2048];cc=c[first:first+2048];feed={it.name:(cc if len(it.shape)==1 else xx) for it in engine.get_inputs()};chunks.append(engine.run(None,feed)[0].reshape(-1))
        elif kind=='torchscript_cuda':
            engine=torch.jit.load(str(pathof(model['torchscript_path'])),map_location='cuda').eval()
            with torch.inference_mode():
                for first in range(0,len(x),2048):
                    xx=torch.as_tensor(x[first:first+2048],device='cuda');cc=torch.as_tensor(c[first:first+2048],device='cuda');y=engine(xx,cc) if model['operator']=='flux' else engine(xx);chunks.append(y.cpu().numpy().reshape(-1))
        else:raise ValueError('Unregistered locked runtime '+kind)
        result=np.concatenate(chunks);torch.cuda.synchronize();del engine;gc.collect();torch.cuda.empty_cache();torch.cuda.synchronize()
    assert result.shape==(len(x),)
    return result

def run_case(dataset,d,rec,identity,runner=None):
    mid=rec['model_id'];dest=B/'predictions'/dataset/(mid+'.npz');meta=dest.with_suffix('.json')
    binding=dict(dataset_sha256=sha(B/'datasets'/(dataset+'.npz')),deployment_lock_sha256=identity['deployment_lock_sha256'],model_id=mid,runtime=rec['runtime'],artifact_sha256=rec.get('artifact_sha256',rec.get('library_sha256')))
    if dest.exists() and meta.exists():
        old=read(meta)
        if old.get('state')=='complete' and old.get('identity')==binding and old.get('prediction_sha256')==sha(dest):return old
    if runner is None:prediction=infer_neural(rec,rec['runtime'],d);residual=np.full(len(prediction),np.nan);iterations=np.full(len(prediction),-1,dtype='i4')
    else:
        status('native_blind_accuracy',model_id=mid,dataset=dataset)
        x=np.ascontiguousarray(d['stencils'],dtype='f8');c=np.ascontiguousarray(d['courants'],dtype='f8');prediction=np.full(len(x),np.nan);residual=np.full(len(x),np.nan);iterations=np.full(len(x),-1,dtype='i4')
        code=runner(x,c,prediction,residual,iterations)
        if code:raise RuntimeError((mid,'native return code',code))
    # Failure outputs remain raw NaN/Inf; no endpoint patching, truth fallback or sample dropping.
    valid=np.isfinite(prediction);atomic_npz(dest,prediction=prediction,valid=valid,query_id=d['query_id'],residual=residual,iterations=iterations)
    saved=dict(state='complete',identity=binding,record={k:v for k,v in rec.items() if k!='run'},dataset=dataset,prediction_path=str(dest.relative_to(W)),prediction_sha256=sha(dest),samples=len(prediction),valid_samples=int(valid.sum()),failures=int((~valid).sum()),formal_timing=False,speed_evidence=False,finished_utc=now())
    save(meta,saved);return saved

def error_stats(error,complete):
    if not len(error) or not complete:return None
    return float(np.mean(np.abs(error)))

def scope_record(d,p,valid,mask,rec,dataset,N,shape,rbin,subset,full_mask=None):
    required=int(mask.sum());available=mask&valid;nv=int(available.sum());complete=required==nv;op=rec['operator'];prediction=p[available];f32=d['stencils'][available,3,3].astype('f8');h=d['h'][available]
    row=dict(dataset=dataset,model_id=rec['model_id'],group=rec.get('group',rec['model_id']),model_type=rec['model_type'],operator=op,seed=rec.get('seed'),backend=rec.get('backend'),runtime=rec['runtime'],N=int(N),h=float(h[0]) if len(h) else (1. if N==-1 else 1/N),shape=shape,radius_bin=rbin,subset=subset,samples=required,valid_samples=nv,failures=required-nv,coverage=nv/required if required else None,
             original_samples=int(full_mask.sum()) if full_mask is not None else required,removed_from_original=int(full_mask.sum())-required if full_mask is not None else 0,
             original_failures=int((full_mask&~valid).sum()) if full_mask is not None else required-nv,formal_timing=False,speed_evidence=False,
             q_mae=None,available_only_q_mae=None,canonical_q_mae=None,kappa_mae=None,available_only_kappa_mae=None,ratio_mae=None,available_only_ratio_mae=None,flux_volume_mae=None,available_only_flux_volume_mae=None,exact_cell_transport_mae=None,available_only_exact_cell_transport_mae=None,physical_transport_area_mae=None,available_only_physical_transport_area_mae=None,max_abs_error=None)
    if op=='curvature':
        error=prediction-d['q'][available];physical=error/h
        row.update(q_mae=error_stats(error,complete),available_only_q_mae=error_stats(error,True),canonical_q_mae=error_stats(prediction-d['curvature'][available],complete),kappa_mae=error_stats(physical,complete),available_only_kappa_mae=error_stats(physical,True),max_abs_error=float(abs(physical).max()) if nv else None)
    else:
        ratio=prediction-d['flux_ratio'][available];canonical=prediction-d['flux'][available].astype('f8');weighted=f32*canonical;transport=f32*prediction-d['swept_volume'][available]
        row.update(ratio_mae=error_stats(ratio,complete),available_only_ratio_mae=error_stats(ratio,True),flux_volume_mae=error_stats(weighted,complete),available_only_flux_volume_mae=error_stats(weighted,True),exact_cell_transport_mae=error_stats(transport,complete),available_only_exact_cell_transport_mae=error_stats(transport,True),physical_transport_area_mae=error_stats(transport*h*h,complete),available_only_physical_transport_area_mae=error_stats(transport*h*h,True),max_abs_error=float(abs(ratio).max()) if nv else None)
    return row

def aggregate(datasets,cases,identity):
    rows=[];common_rows=[];endpoints=[];common_info=[];mask_arrays={}
    for name,d in datasets.items():
        relevant=[c for c in cases if c['dataset']==name];loaded={c['record']['model_id']:dict(np.load(pathof(c['prediction_path']))) for c in relevant}
        for op in ['flux','curvature']:
            cc=[c for c in relevant if c['record']['operator']==op];shared=np.ones(len(d['q']),bool)
            for c in cc:shared&=loaded[c['record']['model_id']]['valid']
            if name=='physical':
                ns=len(NS);per=len(shared)//ns;assert len(shared)==ns*per;shared=np.tile(shared.reshape(ns,per).all(0),ns)
            mask_arrays[name+'_'+op+'_common']=shared
            common_info.append(dict(dataset=name,operator=op,methods=len(cc),samples=len(shared),common_samples=int(shared.sum()),removed=int((~shared).sum()),crossN_common=(name=='physical')))
            for case in cc:
                rec=case['record'];arrays=loaded[rec['model_id']];p=arrays['prediction'];valid=arrays['valid'];assert np.array_equal(arrays['query_id'],d['query_id'])
                for N in np.unique(d['N']):
                    for family in [-1,*np.unique(d['family']).tolist()]:
                        base=(d['N']==N)&((d['family']==family) if family>=0 else True);shape='all' if family<0 else ['circle','ellipse','weak_sine'][family]
                        for label,rb in radius_bins(d['radius_cells']):
                            mask=base&rb
                            rows.append(scope_record(d,p,valid,mask,rec,name,int(N),shape,label,'full'))
                            common_rows.append(scope_record(d,p,valid,mask&shared,rec,name,int(N),shape,label,'all_methods_common_finite',mask))
                for scope,mask in [('all',np.ones(len(p),bool)),('courant_zero',d['courants']==0),('courant_one',d['courants']==1),('courant_low_nonzero',(d['courants']>0)&(d['courants']<.001)),('courant_interior',(d['courants']>0)&(d['courants']<1))]:
                    z=scope_record(d,p,valid,mask,rec,name,-1,'all','all','full');z['scope']=scope;endpoints.append(z)
    csvout(B/'metrics.csv',rows);csvout(B/'common_finite_metrics.csv',common_rows);save(B/'endpoint_scopes.json',dict(records=endpoints,no_endpoint_prediction_replacement=True))
    atomic_npz(B/'common_finite_masks.npz',**mask_arrays)
    save(B/'summary.json',dict(state='complete',scope='New blinded accuracy only;no formal timing or Pareto success claim from this file',formal_timing=False,speed_evidence=False,identity=identity,
        datasets=[dict(dataset=k,samples=len(v['q']),shapes=np.unique(v['shape']).tolist(),N=np.unique(v['N']).tolist()) for k,v in datasets.items()],models=[c['record'] for c in cases],case_records=len(cases),neural_models=24,native_models=4,
        metrics_schema=dict(identity=['dataset','model_id','group','model_type','operator','seed','backend','runtime','N','h','shape','radius_bin','subset'],full='Full predefined support,MAE undefined on any failed output;available_only statistics separate',common='All planned methods finite intersection;physical repeats matched across N;counts+removed reported',curvature='q analyticF64 and physical_kappa error=error_q/h;canonical_q_mae auxiliary',flux='flux_volume_mae=mean(rawF32center*abs(pred-canonicalF32ratio)),matchingvalidation;ratio_mae vs exact full strip ratio;exact_cell_transport_mae=abs(rawF32center*pred-exactnormalizedstrip),physical_transport_area_mae adds explicit h² only as secondary',circle='N=-1,h=1 cell-normalized radius sweep;do not plot as physical refinement'),common_support=common_info,
        all_seeds_retained=True,no_post_test_model_selection=True,deployment_choices=read(R/'deployment_selection_locked.json')['chosen_groups'],finished_utc=now()))
    return len(rows),len(common_rows)

def main(args):
    B.mkdir(parents=True,exist_ok=True)
    if args.describe:
        print(json.dumps(dict(gate='deployment_selection_locked.json',models=24,circle_count=CIRCLE_COUNT,circle_seed=CIRCLE_SEED,physical_seed=PHYSICAL_SEED,resolutions=NS,outputs='blind_test/datasets/manifest.json,metrics.csv,common_finite_metrics.csv,summary.json,completion.json',starts_gpu=False)));return
    if args.check_gate_only:
        exists=(R/'deployment_selection_locked.json').exists();status('gate_check_only',deployment_lock_exists=exists,test_generated=False,gpu_started=False);print(json.dumps(dict(deployment_lock_exists=exists,generated=False,gpu_started=False)));return
    lock=(B/'evaluation.lock').open('a+')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:print('Existing blind accuracy controller retained');return
    (B/'controller.pid').write_text(str(os.getpid()))
    while not (R/'deployment_selection_locked.json').exists():
        status('waiting_for_deployment_selection_lock',test_generated=False,gpu_started=False)
        if not args.wait:return
        time.sleep(15)
    manifest_path=pathof(read(R/'deployment_selection_locked.json').get('model_manifest_path',str((R/'inference/model_manifest.json').relative_to(W))))
    while not manifest_path.exists() or not read(manifest_path).get('complete'):
        status('waiting_for_all24_exports',test_generated=False,gpu_started=False)
        if not args.wait:raise RuntimeError('Deployment lock exists but complete export manifest missing')
        time.sleep(15)
    selection,models,identity=locked_plan()
    if (B/'completion.json').exists():
        complete=read(B/'completion.json')
        if complete.get('state')=='complete' and complete.get('identity')==identity:status('complete',cached=True,models=24);return
    import torch
    torch.set_num_threads(1);torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False;torch.cuda.set_per_process_memory_fraction(.10)
    status('generating_new_blind_data',deployment_lock_sha256=identity['deployment_lock_sha256'],models=24)
    manifests=[circle_data(identity),physical_data(identity)];assert_gate(identity)
    save(B/'datasets/manifest.json',dict(complete=True,identity=identity,files=manifests,all_new_after_deployment_lock=True,raw_input_keys=dict(flux=['stencils','courants'],curvature=['stencils']),geometry_id='Parent physical geometry,kept across N/offset/phase transforms;circle unique base_id',query_id='Unique row identity within each dataset',radius_bins=['4_to_8','8_to_32','32_to_128','128_to_512','above_512'],created_utc=now()))
    datasets={m['dataset']:dict(np.load(pathof(m['path']))) for m in manifests};cases=[]
    for dname,d in datasets.items():
        for native in native_runners():
            rec={k:v for k,v in native.items() if k!='run'};rec.update(model_type='native',group=native['model_id'],seed=None)
            cases.append(run_case(dname,d,rec,identity,native['run']));assert_gate(identity)
    for index,m in enumerate(models):
        kind=selection['model_runtimes'][m['model_id']];key='onnx' if kind=='onnx_cuda' else 'torchscript';assert sha(pathof(m[key+'_path']))==m[key+'_sha256']
        roles=[role for role,group in selection['chosen_groups'][m['operator']].items() if isinstance(group,str) and group==m['group']]
        rec=dict(model_id=m['model_id'],group=m['group'],operator=m['operator'],seed=m['seed'],stencil_size=m['stencil_size'],runtime=kind,artifact_sha256=m[key+'_sha256'],model_type='neural',backend='GPU',deployment_roles=roles,onnx_path=m.get('onnx_path'),torchscript_path=m['torchscript_path'])
        for dname,d in datasets.items():cases.append(run_case(dname,d,rec,identity));assert_gate(identity)
        status('evaluating_locked_models',models_completed=index+1,total_models=24,prediction_cases=len(cases))
    save(B/'prediction_manifest.json',dict(complete=True,identity=identity,cases=cases,formal_timing=False,speed_evidence=False))
    full,common=aggregate(datasets,cases,identity);assert_gate(identity)
    save(B/'completion.json',dict(state='complete',identity=identity,neural_models=24,neural_seeds=SEEDS,native_methods=4,prediction_cases=len(cases),datasets=2,metric_rows=full,common_metric_rows=common,formal_timing=False,speed_evidence=False,finished_utc=now(),elapsed_seconds=time.monotonic()-START))
    status('complete',models=24,native_methods=4,prediction_cases=len(cases),completion='blind_test/completion.json',formal_timing=False)
    print(json.dumps(dict(state='complete',models=24,native_methods=4,metric_rows=full,common_metric_rows=common,formal_timing=False)),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--wait',action='store_true');parser.add_argument('--check-gate-only',action='store_true');parser.add_argument('--describe',action='store_true');args=parser.parse_args()
    try:main(args)
    except BaseException as error:
        B.mkdir(parents=True,exist_ok=True);failure=dict(state='failed',error=repr(error),traceback=traceback.format_exc(),pid=os.getpid(),updated_utc=now())
        with (B/'failures.jsonl').open('a') as stream:stream.write(json.dumps(failure)+'\n')
        save(B/'status.json',failure);raise
