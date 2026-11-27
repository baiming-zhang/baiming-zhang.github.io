"""One bounded frozen-weight CUDA-Graph deployment scout;development only."""
from pathlib import Path
import os,sys
sys.dont_write_bytecode=True
for key in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[key]='1'
import fcntl,gc,hashlib,importlib.util,json,random,time,traceback
import numpy as np
import torch
R=Path(__file__).resolve().parent;P=R.parent;W=P.parent;OUT=R/'runtime_results'
spec=importlib.util.spec_from_file_location('_runtime_reused_validation_helpers',P/'validation_pareto.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
sys.path.insert(0,str(W/'operator-ablation-20261122/shared'))
from gpu_lease import gpu_lease
WARMUP=10;REPEATS=31;STREAM=None;MAX_CAPTURE=30;MAX_RUN=300
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,obj):h.save(path,obj)
def status(state,**extra):save(R/'runtime_status.json',dict(state=state,pid=os.getpid(),updated_utc=h.now(),**extra))
def metric(p,y,f):return h.metric(p,y,f,'flux')
class Runner:
    def __init__(self,model,x,c,mode,stream):
        self.model=model;self.mode=mode;self.stream=stream;self.graph=torch.cuda.CUDAGraph()
        self.host=torch.from_numpy(np.ascontiguousarray(x,dtype=np.float32)).pin_memory();self.hc=torch.from_numpy(np.ascontiguousarray(c,dtype=np.float32)).pin_memory()
        self.ho=torch.empty(len(x),dtype=torch.float32,pin_memory=True)
        stream.wait_stream(torch.cuda.current_stream())
        with torch.cuda.stream(stream):self.x=self.host.to('cuda',non_blocking=True);self.c=self.hc.to('cuda',non_blocking=True)
        stream.synchronize();start=time.perf_counter()
        for _ in range(12):self.invoke('torchscript_eager')
        with torch.inference_mode(),torch.cuda.graph(self.graph,stream=stream):self.graph_output=self.forward()
        self.capture_seconds=time.perf_counter()-start
        assert self.capture_seconds<=MAX_CAPTURE,('capture_budget_exceeded',self.capture_seconds)
        self.invoke('cuda_graph')
    def forward(self):
        if self.mode=='cpu_gpu_cpu':self.x.copy_(self.host,non_blocking=True);self.c.copy_(self.hc,non_blocking=True)
        out=self.model(self.x,self.c)
        if self.mode=='cpu_gpu_cpu':self.ho.copy_(out,non_blocking=True)
        return out
    def invoke(self,runtime):
        with torch.inference_mode(),torch.cuda.stream(self.stream):
            if runtime=='cuda_graph':self.graph.replay()
            else:self.eager_output=self.forward()
        self.stream.synchronize()
    def prediction(self,runtime):
        if self.mode=='cpu_gpu_cpu':return self.ho.numpy().astype(np.float64,copy=True)
        out=self.graph_output if runtime=='cuda_graph' else self.eager_output
        return out.detach().cpu().numpy().astype(np.float64,copy=True)
    def replace_inputs(self,x,c):
        self.host.copy_(torch.from_numpy(np.ascontiguousarray(x,dtype=np.float32)));self.hc.copy_(torch.from_numpy(np.ascontiguousarray(c,dtype=np.float32)))
        if self.mode=='resident':
            with torch.cuda.stream(self.stream):self.x.copy_(self.host,non_blocking=True);self.c.copy_(self.hc,non_blocking=True)
            self.stream.synchronize()
def constraints(p,x,grid):
    p=p.reshape(len(x),len(grid));f=x[:,x.shape[1]//2,x.shape[2]//2].astype(np.float64);v=f[:,None]*p
    lower=np.maximum(0,f[:,None]+grid[None,:]-1);upper=np.minimum(f[:,None],grid[None,:])
    checks=dict(all_finite=bool(np.isfinite(p).all()),zero_endpoint=float(abs(p[:,0]).max()),one_endpoint=float(abs(p[:,-1]-(f>0)).max()),
        monotonicity=float(np.maximum(0,p[:,:-1]-p[:,1:]).max()),capacity=float(np.maximum(np.maximum(0,lower-v),np.maximum(0,v-upper)).max()))
    assert checks['all_finite'] and max(checks[k] for k in ['zero_endpoint','one_endpoint','monotonicity','capacity'])<=2e-6,checks
    return checks
def validate(model,m,val):
    k=m['stencil_size'];a=(7-k)//2;x=np.ascontiguousarray(val['stencils'][:,a:a+k,a:a+k]);c=val['courants'];checks=[];origin=None
    for mode in ['resident','cpu_gpu_cpu']:
        runner=Runner(model,x,c,mode,STREAM);runner.invoke('torchscript_eager');expected=runner.prediction('torchscript_eager');runner.invoke('cuda_graph');actual=runner.prediction('cuda_graph')
        assert np.isfinite(expected).all() and np.isfinite(actual).all();diff=float(abs(actual-expected).max());assert diff<=2e-5
        # Changing data in the static buffers must change outputs correctly.
        permutation=np.roll(np.arange(len(x)),17);runner.replace_inputs(x[permutation],c[permutation]);runner.invoke('cuda_graph');probe=runner.prediction('cuda_graph')
        input_diff=float(abs(probe-expected[permutation]).max());assert input_diff<=2e-5
        checks.append(dict(mode=mode,full_validation_samples=len(x),graph_vs_same_weights_eager_max_abs=diff,input_update_probe_max_abs=input_diff,finite=True))
        if origin is None:origin=expected
        else:assert float(abs(origin-expected).max())<=2e-5
        del runner;gc.collect();torch.cuda.empty_cache()
    oldfile=P/'stage2_courant_objective/validation_pareto/accuracy'/f"{m['model_id']}.npz"
    old=np.load(oldfile)['prediction'];stage2diff=float(abs(origin-old).max());assert stage2diff<=2e-5
    sx=np.concatenate([x[:32],np.zeros((1,k,k),dtype=np.float32),np.ones((1,k,k),dtype=np.float32)])
    grid=np.linspace(0,1,33,dtype=np.float32);xx=np.repeat(sx,len(grid),axis=0);cc=np.tile(grid,len(sx));capacity=[]
    for mode in ['resident','cpu_gpu_cpu']:
        runner=Runner(model,xx,cc,mode,STREAM)
        for runtime in ['torchscript_eager','cuda_graph']:
            runner.invoke(runtime);capacity.append(dict(mode=mode,runtime=runtime,checks=constraints(runner.prediction(runtime),sx,grid)))
        del runner;gc.collect();torch.cuda.empty_cache()
    return dict(model_id=m['model_id'],all_checks_passed=True,validation_metrics=metric(origin,val['flux'].astype(np.float64),val['stencils'][:,3,3].astype(np.float64)),
        runtime_parity=checks,capacity=capacity,source_stage2_max_abs=stage2diff,validation_samples=16384,pure_frozen_model=True)
def native_measure(workload,backend,workers,x,c,truth,identity):
    os.sched_setaffinity(0,{0});indices=[np.arange(len(x))];obj=h.CPUCase(backend,'flux',workers,x,c,indices)
    try:
        for _ in range(WARMUP):obj(0)
        times=[];ratios=[]
        with h.Monitor() as monitor:
            for _ in range(REPEATS):
                t=time.perf_counter_ns();r=obj(0);times.append((time.perf_counter_ns()-t)/1e6);ratios.extend(r)
        pred=obj.banks[0][2].copy();metrics=metric(pred,truth,x[:,3,3].astype(np.float64))
        result=dict(identity=identity,workload_id=workload,backend=backend,operator='flux',native_method='actual_compiled_PLIC_with_normal_plane_reconstruction',workers=workers,
            cpu_pool=list(range(8)),worker_affinity_verified=True,batch=len(x),unique_real_grid_queries=len(x),duplicated_for_workload=False,precision='rawf32_promoted_float64',
            wall_ms=h.stats(times),raw_wall_ms=times,warmup=WARMUP,repeats=REPEATS,metrics=metrics,formal=False,speed_evidence_eligible=False,
            resource=h.machine(),worker_cpu_wall_ratio_median=float(np.median(ratios)),external_load=monitor.evidence(),library_sha256=obj.library_sha,
            boundary='CPU-ready promoted VOF -> original local normal/PLIC/flux -> completedCPU output;persistent pinned dispatch included')
        save(OUT/'native'/f'{workload}_{backend}_W{workers}.json',result);return result
    finally:obj.close()
def gpu_measure(model,m,workload,x,c,truth,identity):
    k=m['stencil_size'];a=(7-k)//2;xx=np.ascontiguousarray(x[:,a:a+k,a:a+k]);results=[]
    cases=['resident','cpu_gpu_cpu'];random.Random(202611270303+m['seed']+len(x)).shuffle(cases)
    for mode in cases:
        runner=Runner(model,xx,c,mode,STREAM);runner.invoke('torchscript_eager');expected=runner.prediction('torchscript_eager');runner.invoke('cuda_graph');graphpred=runner.prediction('cuda_graph')
        diff=float(abs(expected-graphpred).max());assert diff<=2e-5 and np.isfinite(graphpred).all()
        runtimes=['torchscript_eager','cuda_graph'];random.Random(m['seed']+len(x)+(mode=='resident')).shuffle(runtimes)
        for runtime in runtimes:
            for _ in range(WARMUP):runner.invoke(runtime)
            times=[]
            with h.Monitor() as monitor:
                for _ in range(REPEATS):
                    start=time.perf_counter_ns();runner.invoke(runtime);times.append((time.perf_counter_ns()-start)/1e6)
            pred=runner.prediction(runtime);assert float(abs(pred-expected).max())<=2e-5
            result=dict(identity=identity,model_id=m['model_id'],group=m['group'],seed=m['seed'],stencil_size=k,source_checkpoint_sha256=m['checkpoint_sha256'],
                runtime=runtime,mode=mode,workload_id=workload,batch=len(x),unique_real_grid_queries=len(x),duplicated_for_workload=False,
                graph_capture_includes_H2D_and_D2H=mode=='cpu_gpu_cpu',boundary='pinnedCPU input ->actual H2D->unchanged frozen model->actual D2H->completedCPU output' if mode=='cpu_gpu_cpu' else 'GPU-ready input->unchanged frozen model->completedGPU output',
                wall_ms=h.stats(times),raw_wall_ms=times,warmup=WARMUP,repeats=REPEATS,metrics=metric(pred,truth,x[:,3,3].astype(np.float64)),
                same_batch_eager_max_abs=diff,capture_seconds=runner.capture_seconds,stream_pointer=int(STREAM.cuda_stream),precision='float32;TF32off;noAMP',
                external_load=monitor.evidence(),formal=False,speed_evidence_eligible=False,development_grid_only=True)
            save(OUT/'gpu'/f"{m['model_id']}_{workload}_{mode}_{runtime}.json",result);results.append(result)
        del runner;gc.collect();torch.cuda.empty_cache()
    return results
def aggregate(models,accuracy,native,gpu,identity):
    groups=[];comparisons=[]
    for group in sorted({m['group'] for m in models}):
        ids=[m['model_id'] for m in models if m['group']==group];assert len(ids)==3
        for workload in sorted({t['workload_id'] for t in native}):
            for mode in ['resident','cpu_gpu_cpu']:
                for runtime in ['torchscript_eager','cuda_graph']:
                    points=[t for t in gpu if t['model_id'] in ids and t['workload_id']==workload and t['mode']==mode and t['runtime']==runtime];assert len(points)==3
                    ts=[p['wall_ms']['median'] for p in points];es=[p['metrics']['error'] for p in points]
                    values=dict(group=group,model_ids=ids,workload_id=workload,mode=mode,runtime=runtime,mean_ms=float(np.mean(ts)),seed_sd_ms=float(np.std(ts,ddof=1)),max_seed_ms=max(ts),
                        mean_volume_mae=float(np.mean(es)),seed_sd_volume_mae=float(np.std(es,ddof=1)),max_seed_volume_mae=max(es),
                        validation_volume_mae_mean=float(np.mean([accuracy[i]['validation_metrics']['error'] for i in ids])),formal=False)
                    groups.append(values)
                    refs=[t for t in native if t['workload_id']==workload];fast=min(refs,key=lambda t:t['wall_ms']['median'])
                    comparisons.append(dict(group=group,workload_id=workload,mode=mode,runtime=runtime,native_backend=fast['backend'],native_workers=fast['workers'],
                        native_ms=fast['wall_ms']['median'],native_volume_mae=fast['metrics']['error'],neural=values,
                        all_seed_accuracy_reaches_native=max(es)<=fast['metrics']['error'],all_seed_time_fits_native=max(ts)<=fast['wall_ms']['median'],
                        scout_candidate_dominance=bool(max(es)<=fast['metrics']['error'] and max(ts)<=fast['wall_ms']['median']),formal_speed_evidence=False))
    improvements=[]
    for g in [g for g in groups if g['runtime']=='cuda_graph']:
        eager=next(e for e in groups if e['group']==g['group'] and e['workload_id']==g['workload_id'] and e['mode']==g['mode'] and e['runtime']=='torchscript_eager')
        improvements.append(dict(group=g['group'],workload_id=g['workload_id'],mode=g['mode'],eager_mean_ms=eager['mean_ms'],graph_mean_ms=g['mean_ms'],observed_ratio=eager['mean_ms']/g['mean_ms'],scout_only=True))
    result=dict(state='complete',identity=identity,attempt='one_CUDA_graph_launch_optimization',validation_only=True,development_grid_only=True,
        blind_test_generated=False,retraining=False,formal_speed_evidence=False,all_six_frozen_models_retained=True,accuracy=accuracy,native=native,groups=groups,comparisons=comparisons,
        same_weights_runtime_improvements=improvements,stop='Bounded attempt complete;no automatic further runtime search. Promising scout points still need new-seed independent real grid and uncontaminated formal measurements.')
    save(R/'runtime_summary.json',result);return result
def main():
    global STREAM
    OUT.mkdir(exist_ok=True);lock=(R/'runtime.lock').open('a+')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:print('Existing runtime scout retained');return
    (R/'runtime.pid').write_text(str(os.getpid()));protocol=json.loads((R/'runtime_protocol.json').read_text())
    manifest=P/'stage2_courant_objective/inference/model_manifest.json';allmodels=json.loads(manifest.read_text())['models']
    models=[m for m in allmodels if m['group'] in protocol['frozen_groups']];assert len(models)==6 and all(m['torchscript_status']=='ready' for m in models)
    for group in protocol['frozen_groups']:assert sorted(m['seed'] for m in models if m['group']==group)==[4101,4102,4103]
    source=W/protocol['development_workload']['fixture'];identity=dict(protocol_sha256=sha(R/'runtime_protocol.json'),script_sha256=sha(__file__),source_manifest_sha256=sha(manifest),
        validation_sha256=sha(P/'data/validation.npz'),development_fixture_sha256=sha(source),helper_sha256=sha(P/'validation_pareto.py'),frozen_models=[dict(model_id=m['model_id'],torchscript_sha256=m['torchscript_sha256'],checkpoint_sha256=m['checkpoint_sha256']) for m in models])
    if (R/'runtime_completion.json').exists() and json.loads((R/'runtime_completion.json').read_text()).get('identity')==identity:print('Matching completed runtime scout retained');return
    with np.load(P/'data/validation.npz',allow_pickle=False) as z:val={k:np.array(z[k],copy=True) for k in ['stencils','courants','flux']}
    with np.load(source,allow_pickle=False) as z:
        x=np.ascontiguousarray(z['stencils'][:,2:9,2:9]);truth=np.array(z['flux_ratio_truth'],copy=True)
    assert x.shape==(127090,7,7) and x.dtype==np.float32
    workloads=[('real_grid_N4096_unique65536',65536),('real_grid_N4096_all127090',127090)]
    for m in models:
        assert sha(W/m['torchscript_path'])==m['torchscript_sha256'] and sha(W/m['checkpoint_path'])==m['checkpoint_sha256']
    native=[];gpu=[];accuracy={};jobs=[(wid,b,backend,w) for wid,b in workloads for backend in ['aphros','simliquid'] for w in [1,4,8]]
    random.Random(202611270304).shuffle(jobs);begin=time.perf_counter()
    for wid,b,backend,workers in jobs:
        status('native_CPU_scout',workload_id=wid,backend=backend,workers=workers,formal=False)
        with gpu_lease('stage3-runtime-native-scout',formal_timing=False):native.append(native_measure(wid,backend,workers,x[:b],np.full(b,.5,dtype=np.float32),truth[:b],identity))
    torch.set_num_threads(1);torch.set_num_interop_threads(1);torch.cuda.set_per_process_memory_fraction(.10);torch.cuda.reset_peak_memory_stats()
    torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False;torch.set_float32_matmul_precision('highest');os.sched_setaffinity(0,{8})
    random.Random(202611270305).shuffle(models)
    for m in models:
        status('waiting_shared_GPU_lease',model_id=m['model_id'],formal=False)
        with gpu_lease('stage3-frozen-flux-CUDAgraph-scout',formal_timing=False):
            STREAM=torch.cuda.Stream();assert STREAM.cuda_stream!=0
            model=torch.jit.load(str(W/m['torchscript_path']),map_location='cuda').eval()
            status('validating_frozen_runtime',model_id=m['model_id'],formal=False);checks=validate(model,m,val);accuracy[m['model_id']]=checks
            save(OUT/'accuracy'/f"{m['model_id']}.json",checks)
            for wid,b in workloads:
                status('runtime_scout',model_id=m['model_id'],workload_id=wid,formal=False)
                gpu.extend(gpu_measure(model,m,wid,x[:b],np.full(b,.5,dtype=np.float32),truth[:b],identity))
            torch.cuda.synchronize();del model;gc.collect();torch.cuda.empty_cache()
        assert time.perf_counter()-begin<=MAX_RUN,('bounded_runtime_budget_exceeded',time.perf_counter()-begin)
    result=aggregate(models,accuracy,native,gpu,identity)
    save(R/'runtime_completion.json',dict(state='complete',identity=identity,models=6,native_cases=len(native),GPU_cases=len(gpu),formal=False,
        peak_allocated_bytes=torch.cuda.max_memory_allocated(),peak_reserved_bytes=torch.cuda.max_memory_reserved(),allocation_fraction_cap=.10,elapsed_seconds=time.perf_counter()-begin))
    status('complete',models=6,formal=False,summary='runtime_summary.json')
    print(json.dumps(dict(state='complete',runtime_improvements=result['same_weights_runtime_improvements'],comparisons=result['comparisons']),indent=2))
if __name__=='__main__':
    try:main()
    except Exception:status('stopped_bounded_runtime_attempt',traceback=traceback.format_exc(),formal=False,no_further_runtime_search=True);raise
