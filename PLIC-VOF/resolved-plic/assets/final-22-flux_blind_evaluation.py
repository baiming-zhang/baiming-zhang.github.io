"""Fixed three-seed flux-only independent accuracy; hard pre-generation gate."""
from pathlib import Path
import os,sys
sys.dont_write_bytecode=True
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
import argparse,fcntl,importlib.util,json,time,traceback
import numpy as np
R=Path(__file__).resolve().parent;W=R.parent;B=R/'flux_blind_test'
sys.path.insert(0,str(W));sys.path.insert(0,str(R))
spec=importlib.util.spec_from_file_location('_flux_blind_shared_math',R/'blind_evaluation.py');base=importlib.util.module_from_spec(spec);sys.modules[spec.name]=base;spec.loader.exec_module(base)
base.B=B
save=base.save;read=base.read;sha=base.sha;status=base.status
def pathof(value):
    text=str(value).replace('\\','/')
    if len(text)>2 and text[1:3]==':/':text='/mnt/'+text[0].lower()+text[2:]
    p=Path(text);return p if p.is_absolute() else W/p
base.pathof=pathof
GROUP='flux_compact_s5_w64_b16__sixteen_courants';SEEDS=[4101,4102,4103]
LOCK=R/'deployment_selection_locked.json'

def assert_gate(identity=None):
    if not LOCK.exists():raise RuntimeError('New flux blind geometry forbidden before root deployment lock')
    if identity and sha(LOCK)!=identity['deployment_lock_sha256']:raise RuntimeError('Deployment changed after blind registration;do not reuse/resample')
base.assert_gate=assert_gate

def cpu_checks():
    # Index transform check, not a generated blind geometry or model inference.
    y,x=np.meshgrid(np.arange(3,-4,-1),np.arange(-3,4),indexing='ij');coded=10*x+y
    transformed=coded.T[::-1,::-1]
    assert np.array_equal(transformed,10*y+x)
    f=np.array([.05,.25,.5,.75,.95]);c=np.array([0.,.01,.5,.9,1.]);v=np.minimum(f,c)*.7+np.maximum(0,f+c-1)*.3
    comp=c-v;assert np.max(abs(c-comp-v))<1e-15
    assert np.all(comp>=np.maximum(0,(1-f)+c-1)-1e-15) and np.all(comp<=np.minimum(1-f,c)+1e-15)
    save(B/'cpu_checks.json',dict(passed=True,generated_test_data=False,gpu_started=False,checks=['top-row-first transpose+reverse axes maps new right strip to old top strip','phase complement swept area c-v and capacity identity','all fixed seed/config/runtime assertions registered'],created_utc=base.now()))

def locked_plan():
    assert_gate();lock=read(LOCK);manifest_path=pathof(lock['model_manifest_path']);manifest=read(manifest_path)
    assert manifest['complete'];models=[m for m in manifest['models'] if m['model_id'] in lock['model_ids']]
    assert len(models)==3 and len(lock['model_ids'])==3 and sorted(m['seed'] for m in models)==SEEDS
    assert all(m['group']==GROUP and m['operator']=='flux' and m['stencil_size']==5 and m['config']['width']==64 and m['config']['num_bins']==16 and m['status']=='ready' for m in models)
    assert set(lock['model_runtimes'])==set(lock['model_ids']) and set(lock['model_runtimes'].values())=={'cuda_graph'}
    adapter=pathof(lock['runtime_adapter_path']);assert sha(adapter)==lock['runtime_adapter_sha256']
    for m in models:
        assert m['torchscript_status']=='ready' and sha(pathof(m['torchscript_path']))==m['torchscript_sha256']
        assert sha(pathof(m['checkpoint_path']))==m['checkpoint_sha256']
    identity=dict(deployment_lock_sha256=sha(LOCK),protocol_sha256=sha(R/'flux_blind_protocol.json'),script_sha256=sha(__file__),shared_math_sha256=sha(R/'blind_evaluation.py'),export_manifest_path=str(manifest_path.relative_to(W)),export_manifest_sha256=sha(manifest_path),runtime_adapter_path=str(adapter.relative_to(W)),runtime_adapter_sha256=sha(adapter),circle_seed=base.CIRCLE_SEED,circle_count=base.CIRCLE_COUNT,physical_seed=base.PHYSICAL_SEED,resolutions=base.NS,
                  models=[dict(model_id=m['model_id'],seed=m['seed'],runtime='cuda_graph',artifact_sha256=m['torchscript_sha256'],checkpoint_sha256=m['checkpoint_sha256']) for m in sorted(models,key=lambda m:m['seed'])])
    p=B/'registration.json'
    if p.exists():assert read(p)['identity']==identity,'Blind binding changed; preserve test,do not silently regenerate'
    else:save(p,dict(identity=identity,registered_utc=base.now(),new_after_root_deployment_lock=True,curvature_excluded=True,selection_source='Frozen validation/development scout only;no test selection'))
    return lock,sorted(models,key=lambda m:m['seed']),identity

def add_support(datasets,manifests):
    reports=[]
    for name,d in datasets.items():
        local=d['radius_cells']>=4-1e-12
        if name=='circle':structure=d['radius_cells'].copy()
        else:
            p=d['physical_parameters'][d['physical_id']];a=p[:,1];b=p[:,2]
            rr=np.where(p[:,0]==0,a,np.minimum(a,b)**2/np.maximum(a,b))
            sine=p[:,0]==2;rr[sine]=1/(np.abs(p[sine,5])*p[sine,6]**2)
            structure=rr/d['h']
        resolved=local&(structure>=4-1e-12)
        if name=='physical':
            per=len(resolved)//len(base.NS);assert len(resolved)==per*len(base.NS)
            common=np.tile(resolved.reshape(len(base.NS),per).all(0),len(base.NS))
        else:common=resolved.copy()
        old_support=d.get('resolved_allN_common')
        if old_support is not None:assert np.array_equal(old_support,common) and np.array_equal(d['resolved'],resolved)
        d.update(local_resolved=local,structure_min_radius_cells=structure,resolved=resolved,resolved_allN_common=common)
        path=B/'datasets'/f'{name}.npz'
        if old_support is None:base.atomic_npz(path,**d)
        manifest=next(m for m in manifests if m['dataset']==name);manifest['sha256']=sha(path)
        save(path.with_suffix('.json'),manifest)
        rows=[]
        for N in np.unique(d['N']):
            for fam in [-1,*np.unique(d['family']).tolist()]:
                m=(d['N']==N)&((d['family']==fam) if fam>=0 else True);n=int(m.sum())
                rows.append(dict(N=int(N),shape='all' if fam<0 else ['circle','ellipse','weak_sine'][fam],planned=n,local_resolved=int((m&local).sum()),structure_resolved=int((m&(structure>=4-1e-12)).sum()),primary_resolved=int((m&resolved).sum()),primary_coverage=float((m&resolved).sum()/n),allN_common_resolved=int((m&common).sum())))
        reports.append(dict(dataset=name,samples=len(resolved),physical_queries_removed=0,perN=rows))
    save(B/'resolution_coverage.json',dict(records=reports,raw_queries_kept=True,primary_requires='local and conservative structure R/h>=4',crossN='Fixed original sample_id+geometry. allN common support separately;perN aggregate supports may differ and cannot alone prove convergence.'))
    return reports

def infer_graph(model,kind,d):
    assert kind=='cuda_graph'
    identity=ACTIVE_IDENTITY;assert_gate(identity)
    adapter=base.module('_locked_flux_graph_adapter',pathof(identity['runtime_adapter_path']))
    # Adapter API is fixed before deployment lock,including any tail padding.
    with base.gpu_lease('flux-independent-accuracy-'+model['model_id'],formal_timing=False):
        import torch
        try:
            torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
            x=np.ascontiguousarray(d['stencils'][:,1:6,1:6],dtype='f4');c=np.ascontiguousarray(d['courants'],dtype='f4')
            runner=adapter.FluxCUDAGraphAdapter(pathof(model['torchscript_path']),x,c,mode='cpu_gpu_cpu',expected_artifact_sha256=model['artifact_sha256'],warmup=12)
            result=runner.predict_numpy()
            result=np.asarray(result).reshape(-1)
            assert result.shape==(len(d['courants']),)
            if hasattr(runner,'close'):runner.close()
            del runner
            return result
        finally:
            if torch.cuda.is_initialized():torch.cuda.synchronize()
base.infer_neural=infer_graph
ACTIVE_IDENTITY=None

def aggregate(datasets,cases,identity):
    rows=[];coverage=[];mask_arrays={}
    for name,d in datasets.items():
        cc=[c for c in cases if c['dataset']==name];loaded={c['record']['model_id']:dict(np.load(pathof(c['prediction_path']))) for c in cc}
        common_finite=np.logical_and.reduce([loaded[c['record']['model_id']]['valid'] for c in cc])
        if name=='physical':
            per=len(common_finite)//len(base.NS);common_finite=np.tile(common_finite.reshape(len(base.NS),per).all(0),len(base.NS))
        mask_arrays[name+'_all_methods_finite']=common_finite
        mask_arrays[name+'_resolved']=d['resolved'];mask_arrays[name+'_resolved_allN_common']=d['resolved_allN_common']
        for case in cc:
            rec=case['record'];a=loaded[rec['model_id']];assert np.array_equal(a['query_id'],d['query_id']);pred=a['prediction'];valid=a['valid']
            for N in np.unique(d['N']):
                for fam in [-1,*np.unique(d['family']).tolist()]:
                    scope=(d['N']==N)&((d['family']==fam) if fam>=0 else True);shape='all' if fam<0 else ['circle','ellipse','weak_sine'][fam]
                    for rbin,rmask in base.radius_bins(d['radius_cells']):
                        original=scope&rmask
                        for subset,condition in [('planned_full',np.ones(len(pred),bool)),('resolved_primary',d['resolved']),('resolved_allN_common',d['resolved_allN_common']),('resolved_common_finite_conditional',d['resolved']&common_finite),('resolved_allN_common_finite_conditional',d['resolved_allN_common']&common_finite)]:
                            row=base.scope_record(d,pred,valid,original&condition,rec,name,int(N),shape,rbin,subset,original)
                            row['resolution_coverage']=row['samples']/row['original_samples'] if row['original_samples'] else None;rows.append(row)
            coverage.append(dict(dataset=name,model_id=rec['model_id'],planned=len(pred),failures=int((~valid).sum()),coverage=float(valid.mean()),resolved_planned=int(d['resolved'].sum()),resolved_failures=int((d['resolved']&~valid).sum())))
    base.csvout(B/'metrics.csv',rows);base.atomic_npz(B/'support_masks.npz',**mask_arrays)
    save(B/'summary.json',dict(state='complete',identity=identity,group=GROUP,seeds=SEEDS,neural_models=3,native_models=2,curvature_tested=False,records=len(rows),coverage=coverage,formal_timing=False,speed_evidence=False,validation_selection_frozen=True,no_post_test_tuning=True,
        primary='resolved_primary exact_cell_transport_mae;ratio_mae and canonical flux_volume_mae retained',physical_support='Original allN fixed geometries retained. PerN resolved coverage varies;fixed allN resolved subset separate. Full and common finite conditional metrics never conflated.',metric_schema=['dataset','model_id','group','seed','N','h','shape','radius_bin','subset','samples','failures','coverage','original_samples','resolution_coverage','exact_cell_transport_mae','ratio_mae','flux_volume_mae','physical_transport_area_mae'],finished_utc=base.now()))
    return len(rows)

def main(args):
    global ACTIVE_IDENTITY
    B.mkdir(parents=True,exist_ok=True);cpu_checks()
    if args.check_gate_only:
        status('gate_check_only',deployment_lock_exists=LOCK.exists(),generated_test_data=False,gpu_started=False);print(json.dumps(dict(deployment_lock_exists=LOCK.exists(),generated=False,gpu_started=False)));return
    handle=(B/'evaluation.lock').open('a+')
    try:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:print('Existing flux blind evaluator retained');return
    (B/'controller.pid').write_text(str(os.getpid()))
    while not LOCK.exists():
        status('waiting_for_root_deployment_lock',generated_test_data=False,gpu_started=False)
        if not args.wait:return
        time.sleep(15)
    lock,models,identity=locked_plan();ACTIVE_IDENTITY=identity
    complete=B/'completion.json'
    if complete.exists() and read(complete).get('identity')==identity and read(complete).get('complete'):status('complete',cached=True);return
    import torch
    torch.set_num_threads(1);torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False;torch.cuda.set_per_process_memory_fraction(.15)
    manifests=[base.circle_data(identity),base.physical_data(identity,resolved_policy='preserve_all')]
    datasets={m['dataset']:dict(np.load(pathof(m['path']))) for m in manifests};support=add_support(datasets,manifests)
    save(B/'datasets/manifest.json',dict(complete=True,identity=identity,files=manifests,resolution_coverage=support,new_after_root_lock=True,curvature_tested=False))
    cases=[]
    for name,d in datasets.items():
        for native in [v for v in base.native_runners() if v['operator']=='flux']:
            rec={k:v for k,v in native.items() if k!='run'};rec.update(model_type='native',group=native['model_id'],seed=None)
            cases.append(base.run_case(name,d,rec,identity,native['run']));assert_gate(identity)
        for m in models:
            rec=dict(model_id=m['model_id'],group=m['group'],operator='flux',seed=m['seed'],stencil_size=5,runtime='cuda_graph',artifact_sha256=m['torchscript_sha256'],torchscript_path=m['torchscript_path'],model_type='neural',backend='GPU')
            cases.append(base.run_case(name,d,rec,identity));assert_gate(identity)
    save(B/'prediction_manifest.json',dict(complete=True,identity=identity,cases=cases,formal_timing=False))
    count=aggregate(datasets,cases,identity);assert_gate(identity)
    save(complete,dict(complete=True,state='complete',identity=identity,neural_models=3,seeds=SEEDS,native_methods=2,prediction_cases=10,metrics_rows=count,curvature_tested=False,formal_timing=False,speed_evidence=False,finished_utc=base.now()));status('complete',prediction_cases=10,formal_timing=False)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--wait',action='store_true');ap.add_argument('--check-gate-only',action='store_true');args=ap.parse_args()
    try:main(args)
    except BaseException as e:status('failed',error=repr(e),traceback=traceback.format_exc());raise
