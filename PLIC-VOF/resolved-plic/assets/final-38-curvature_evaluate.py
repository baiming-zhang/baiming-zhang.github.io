"""One registered four-method/three-seed validation matrix; no training or timing."""
from pathlib import Path
import os,sys
sys.dont_write_bytecode=True
for name in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[name]='1'
import hashlib,importlib.util,json,statistics,time,traceback
import numpy as np
import torch
R=Path(__file__).resolve().parent;P=R.parent;W=P.parent
sys.path.insert(0,str(R));sys.path.insert(0,str(W/'operator-ablation-20261122/shared'))
from curvature_models import CurvatureProjection,METHODS,d4_views,view_counts
from gpu_lease import gpu_lease
SEEDS=[4101,4102,4103]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def save(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False));tmp.replace(p)
def now():return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
def status(state,**kw):save(R/'curvature_status.json',dict(state=state,pid=os.getpid(),updated_utc=now(),**kw))
def build(config):
    name='_curvature_fixed_original'
    if name not in sys.modules:
        spec=importlib.util.spec_from_file_location(name,W/'operator-ablation-20261122/03-data-stencil/models.py');m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m)
    return sys.modules[name].OperatorMLP(config)
def source(seed):return P/'stage2_courant_objective/groups/curvature__unweighted_q_MSE'/f'seed{seed}/best.pt'
def metrics(pred,y,weak,plane):
    valid=np.isfinite(pred);e=np.abs(pred.astype(np.float64)-y.astype(np.float64));n=len(y);out=dict(samples=n,valid_samples=int(valid.sum()),failures=int((~valid).sum()),coverage=float(valid.mean()),q_mae=float(e.mean()) if valid.all() else None,available_only_q_mae=float(e[valid].mean()) if valid.any() else None)
    for name,mask in [('weak',weak),('plane',plane)]:
        out[name+'_samples']=int(mask.sum());out[name+'_mae']=float(e[mask].mean()) if mask.any() and valid[mask].all() else None
    return out
def cpu_check(checkpoints):
    grid=torch.arange(49,dtype=torch.float64).reshape(1,7,7)/64;views=d4_views(grid)
    assert len({tuple(v.reshape(-1).tolist()) for v in views})==8
    assert all(torch.equal(v[:,3,3],grid[:,3,3]) for v in views)
    assert torch.equal(1-(1-grid),grid)
    for a in views:
        assert {tuple(v.reshape(-1).tolist()) for v in d4_views(a)}=={tuple(v.reshape(-1).tolist()) for v in views}
    gen=torch.Generator().manual_seed(202611270601);x=torch.rand(11,7,7,generator=gen,dtype=torch.float64);rows=[]
    with torch.no_grad():
        for seed,ck in zip(SEEDS,checkpoints):
            raw=build(ck['config']).double().eval();raw.load_state_dict(ck['model']);params=sum(p.numel() for p in raw.parameters())
            for method in METHODS:
                m=CurvatureProjection(raw,method);assert sum(p.numel() for p in m.parameters())==params
                baseline=m(x);assert torch.isfinite(baseline).all()
                manual=torch.stack([raw(v) for v in (d4_views(x) if 'd4' in method else [x])]).mean(0)
                if 'phase' in method:manual=.5*(manual-torch.stack([raw(1-v) for v in (d4_views(x) if 'd4' in method else [x])]).mean(0))
                diff=float((m(x)-manual).abs().max());assert diff<1e-12
                spatial=max(float((m(v)-baseline).abs().max()) for v in d4_views(x)) if 'd4' in method else None
                odd=float((m(1-x)+baseline).abs().max()) if 'phase' in method else None
                if spatial is not None:assert spatial<1e-12
                if odd is not None:assert odd<1e-12
                assert torch.count_nonzero(m(torch.cat((torch.zeros_like(x),torch.ones_like(x)),0)))==0
                rows.append(dict(seed=seed,method=method,params=params,manual_average_max_abs=diff,d4_invariance_max_abs=spatial,complement_odd_max_abs=odd,views=view_counts(method)))
    save(R/'curvature_cpu_checks.json',dict(passed=True,device='CPU_float64',signed_q_convention='Spatial D4 leaves labels unchanged; phase complement negates labels. Verified group closure/center mapping and projected parity, no geometric q computation.',checks=rows,created_utc=now()))
def main():
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    protocol=read(R/'curvature_protocol.json');assert protocol['candidates']==list(METHODS)
    sources=[source(s) for s in SEEDS];identity=dict(protocol_sha256=sha(R/'curvature_protocol.json'),sources=[dict(seed=s,path=str(p.relative_to(W)),sha256=sha(p)) for s,p in zip(SEEDS,sources)],script_sha256=sha(__file__),models_sha256=sha(R/'curvature_models.py'))
    oldpath=W/'operator-ablation-20261122/03-data-stencil/data/enhanced_validation.npz';newpath=P/'data/validation.npz'
    identity['datasets']={n:dict(path=str(p.relative_to(W)),sha256=sha(p),samples=count) for n,p,count in [('old',oldpath,16355),('circle',newpath,16384)]}
    if (R/'curvature_completion.json').exists() and read(R/'curvature_completion.json').get('identity')==identity:print('cached complete');return
    checkpoints=[torch.load(p,map_location='cpu',weights_only=False) for p in sources]
    status('cpu_preflight',gpu_started=False);cpu_check(checkpoints)
    datasets={}
    for name,path in [('old',oldpath),('circle',newpath)]:
        with np.load(path) as z:
            keep=(abs(z['curvature'])<=.25) if name=='old' else np.ones(len(z['curvature']),bool)
            datasets[name]=dict(stencils=z['stencils'][keep],truth=z['curvature'][keep],radius=z['radius_cells'][keep] if 'radius_cells' in z else None)
    assert all(len(datasets[n]['truth'])==identity['datasets'][n]['samples'] for n in datasets)
    refs=read(P/'groups/curvature_continue_original/reference_validation.json')['old']
    native=read(P/'native_validation/summary.json');native_rows=[r for r in native['records'] if r['operator']=='curvature'] if 'records' in native else [r for r in native['methods'] if r['operator']=='curvature']
    assert native['identity']['dataset_sha256']==identity['datasets']['circle']['sha256']
    # Native summary shape is checked against its saved primary metric, never reconstructed from truth geometry.
    native_values=[]
    for r in native_rows:
        scope=r.get('metrics',r.get('scopes',r));scope=scope.get('all',scope)
        value=scope.get('q_mae',scope.get('mae'));native_values.append(dict(model_id=r.get('model_id',r.get('backend')+'_particles'),q_mae=value,coverage=scope.get('coverage'),failures=scope.get('failures'),prediction_path=r.get('prediction_path')))
    assert native_values and all(r['q_mae'] is not None and r['coverage']==1 for r in native_values),native_values
    target=min(r['q_mae'] for r in native_values)
    rows=[];status('waiting_gpu_lease',gpu_started=False)
    with gpu_lease('stage3-curvature-fixed-symmetry-validation',formal_timing=False):
        try:
            torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False;torch.cuda.set_per_process_memory_fraction(.15)
            status('validation_running',gpu_started=True)
            for i,(seed,ck) in enumerate(zip(SEEDS,checkpoints)):
                raw=build(ck['config']).eval().cuda();raw.load_state_dict(ck['model']);params=sum(p.numel() for p in raw.parameters())
                for method in METHODS:
                    wrapper=CurvatureProjection(raw,method).eval();scopes={};radius_metrics={}
                    for name,d in datasets.items():
                        pred=[]
                        with torch.inference_mode():
                            for start in range(0,len(d['truth']),512):pred.append(wrapper(torch.from_numpy(d['stencils'][start:start+512]).cuda()).cpu().numpy())
                        pred=np.concatenate(pred);np.savez_compressed(R/f'curvature_predictions_{seed}_{method}_{name}.npz',prediction=pred)
                        y=d['truth'];scopes[name]=metrics(pred,y,abs(y)<.001,y==0)
                        if name=='circle':
                            for lo,hi in [(4,8),(8,32),(32,128),(128,512.000001)]:
                                mask=(d['radius']>=lo)&(d['radius']<hi);radius_metrics[f'{lo}_to_{int(hi)}']=metrics(pred[mask],y[mask],abs(y[mask])<.001,y[mask]==0)
                    a=scopes['old'];guard=refs[i];eligible=all(a[f] is not None and a[f]<=1.1*guard[g]+1e-7 for f,g in [('q_mae','primary_mae'),('weak_mae','weak_mae'),('plane_mae','plane_mae')]) and a['coverage']==1 and scopes['circle']['coverage']==1
                    if method=='raw':assert abs(scopes['circle']['q_mae']-ck['selection']['circle_mae'])<2e-6
                    row=dict(seed=seed,method=method,source_step=ck['step'],params=params,views=view_counts(method),guard_pass=bool(eligible),reaches_native=bool(eligible and scopes['circle']['q_mae']<=target),metrics=scopes,circle_radius_bins=radius_metrics)
                    rows.append(row);save(R/'curvature_partial.json',dict(identity=identity,records=rows,complete=False));status('validation_running',seed=seed,method=method,records_completed=len(rows));print(json.dumps(dict(seed=seed,method=method,mae=scopes['circle']['q_mae'],guard=eligible)),flush=True)
                del wrapper,raw;torch.cuda.empty_cache()
        finally:
            if torch.cuda.is_initialized():torch.cuda.synchronize()
    assert all(sha(source(s))==identity['sources'][i]['sha256'] for i,s in enumerate(SEEDS))
    groups=[]
    for method in METHODS:
        rr=[r for r in rows if r['method']==method];vals=[r['metrics']['circle']['q_mae'] for r in rr];groups.append(dict(method=method,seeds=SEEDS,q_mae_mean=statistics.mean(vals),q_mae_std=statistics.stdev(vals),all_guard_pass=all(r['guard_pass'] for r in rr),all_reach_native=all(r['reaches_native'] for r in rr),views=view_counts(method)))
    summary=dict(state='complete',identity=identity,records=rows,groups=groups,native_reference=native_values,native_q_mae_target=target,formal_timing=False,validation_only=True,new_training=False,blind_test_generated=False,parameters_unchanged=True,guard='Original enhanced7,10%+1e-7 overall/weak/plane,full coverage',next='Only all-seed passing non-raw candidates may receive real runtime measurement; no estimated timing from views.',finished_utc=now())
    save(R/'curvature_summary.json',summary);save(R/'curvature_completion.json',dict(complete=True,state='complete',identity=identity,records=12,groups=4,training=False,blind_test_generated=False,formal_timing=False,finished_utc=now()));status('complete',records=12,gpu_released=True)
if __name__=='__main__':
    try:main()
    except BaseException as e:status('failed',error=repr(e),traceback=traceback.format_exc());raise
