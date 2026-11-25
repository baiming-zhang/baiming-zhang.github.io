"""Resumable GPU neural vs native CPU deployment measurements on real VOF grids."""
from pathlib import Path
import os,sys
sys.dont_write_bytecode=True
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
import json,time,subprocess,hashlib,gc,traceback,fcntl,random,threading,concurrent.futures,platform
import numpy as np
from large_grid_operator_fixtures import R,ROOT,save

NATIVE=ROOT/'operator-width-native-baselines/02-native-bench'
sys.path.insert(0,str(NATIVE));sys.path.insert(0,str(ROOT/'operator-ablation-20261122/shared'))
from native_predict import native_function
from nn_registry import registry
from gpu_lease import gpu_lease

PROTOCOL='actual_grid_gpu32_cpu64_operator_v1'
WARMUP=10;REPEATS=93;CHUNK=65536
START=time.time()

def status(state,**extra):save(R/'status.json',dict(state=state,pid=os.getpid(),elapsed_seconds=time.time()-START,**extra))
def stats(a):
    a=np.array(a,float);return dict(count=len(a),median=float(np.median(a)),p05=float(np.quantile(a,.05)),p95=float(np.quantile(a,.95)))
def smi():
    def query(args):
        r=subprocess.run(['/usr/lib/wsl/lib/nvidia-smi',*args],capture_output=True,text=True)
        return dict(code=r.returncode,output=r.stdout.strip(),error=r.stderr.strip())
    return dict(gpu=query(['--query-gpu=name,memory.total,memory.used,utilization.gpu,temperature.gpu,clocks.sm','--format=csv,noheader,nounits']),processes=query(['--query-compute-apps=pid,process_name,used_memory','--format=csv,noheader,nounits']))
def foreign_processes(s):
    out=[]
    for line in s['processes']['output'].splitlines():
        parts=line.split(',')
        try:pid=int(parts[0].strip())
        except (ValueError,IndexError):continue
        if pid!=os.getpid():out.append(line)
    return out
def wait_idle():
    first=None
    while True:
        sample=smi();foreign=foreign_processes(sample)
        try:util=float(sample['gpu']['output'].split(',')[3].strip());idle=sample['gpu']['code']==0 and util<=5 and not foreign
        except Exception:idle=False
        if idle:
            if first is None:first=time.monotonic()
            if time.monotonic()-first>=90:return sample
        else:first=None
        status('waiting_for_uncontended_gpu',waiting_seconds=time.time()-START,external_gpu_processes=foreign,gpu=sample,reason='Formal deployment timing waits for90s idle and no foreign CUDA process;CPU fixture/accuracy preparation is complete.')
        time.sleep(10)
def load(case):return dict(np.load(ROOT/case['source']))
def metrics(pred,d,op,n):
    good=np.isfinite(pred);answer=dict(samples=len(pred),valid_samples=int(good.sum()),coverage=float(good.mean()))
    if op=='curvature':values={'q':(pred,d['q_truth']),'physical_kappa':(pred*n,d['kappa_truth'])}
    else:values={'ratio':(pred,d['flux_ratio_truth']),'flux_density':(d['center_fraction']*pred/.5,d['flux_density_truth'])}
    for name,(p,t) in values.items():
        error=abs(p[good]-t[good]);mean=float(error.mean()) if len(error) else None
        answer[name]=dict(mae=mean if good.all() else None,available_mae=mean,max_abs=float(error.max()) if len(error) else None)
    return answer

class Native:
    def __init__(self,backend,op,x,workers):
        self.run=native_function(backend,op);self.x=x;self.c=np.full(len(x),.5);self.out=np.empty(len(x));self.workers=workers
        self.edges=np.linspace(0,len(x),workers+1,dtype=int);self.pool=None
        if workers>1:
            allowed=list(range(os.cpu_count()));counter=iter(allowed[:workers]);lock=threading.Lock()
            def initialize():
                with lock:core=next(counter)
                os.sched_setaffinity(0,{core})
            self.pool=concurrent.futures.ThreadPoolExecutor(max_workers=workers,initializer=initialize)
    def part(self,i):
        a,b=self.edges[i:i+2];begin=time.thread_time_ns();t=time.perf_counter_ns()
        code=self.run(self.x[a:b],self.c[a:b],self.out[a:b])
        if code:raise RuntimeError(('native return',code))
        return dict(cpu_ns=time.thread_time_ns()-begin,wall_ns=time.perf_counter_ns()-t)
    def __call__(self):
        return list(self.pool.map(self.part,range(self.workers))) if self.pool else [self.part(0)]
    def close(self):
        if self.pool:self.pool.shutdown(wait=True)

def cpu_point(case,backend,op,workers,d,formal):
    key=f"{case['scene']}_N{case['N']}_{backend}_{op}_cpu{workers}";dest=R/('cpu_timings' if formal else 'cpu_preflight')/(key+'.json')
    if dest.exists():return json.loads(dest.read_text())
    x=np.ascontiguousarray(d['stencils'][:,2:9,2:9],dtype=np.float64);obj=Native(backend,op,x,workers)
    try:
        obj();pred=obj.out.copy();truthfile=R/'predictions'/f"{case['scene']}_N{case['N']}_{backend}_{op}.npz"
        truthfile.parent.mkdir(exist_ok=True)
        if truthfile.exists():
            baseline=np.load(truthfile)['prediction'];assert np.array_equal(np.isfinite(baseline),np.isfinite(pred))
            assert np.allclose(baseline,pred,atol=1e-12,rtol=1e-12,equal_nan=True)
        else:np.savez_compressed(truthfile,prediction=pred)
        result=dict(protocol=PROTOCOL,case=case,backend=backend,operator=op,device='CPU',workers=workers,precision='float64',metrics=metrics(pred,d,op,case['N']),prediction=str(truthfile.relative_to(ROOT)),formal=formal)
        if formal:
            for _ in range(WARMUP):obj()
            samples=[];cpu_ratios=[]
            for _ in range(REPEATS):
                begin=time.perf_counter_ns();slices=obj();samples.append((time.perf_counter_ns()-begin)/1e6)
                cpu_ratios.extend(v['cpu_ns']/max(v['wall_ns'],1) for v in slices if v['wall_ns']>10_000_000)
            result.update(wall_ms=stats(samples),raw_wall_ms=samples,worker_cpu_wall_ratio_median=float(np.median(cpu_ratios)) if cpu_ratios else None,
                measurement_boundary='Prepared CPU VOF to completed CPU output;native geometry and persistent parallel dispatch included',
                scheduling='one pinned thread per worker;one/8/24 workers independently measured',loadavg=list(os.getloadavg()))
        save(dest,result);return result
    finally:obj.close()

def prepare():
    from large_grid_operator_fixtures import main as fixtures
    if not (R/'fixture_manifest.json').exists():fixtures()
    cases=json.loads((R/'fixture_manifest.json').read_text())['cases'];complete=[]
    for case in cases:
        d=load(case)
        for backend in ['aphros','simliquid']:
            for op in ['flux','curvature']:
                status('cpu_accuracy_preflight',case=case['scene'],N=case['N'],backend=backend,operator=op,completed=len(complete))
                complete.append(cpu_point(case,backend,op,1,d,False))
        del d;gc.collect()
    save(R/'preparation_complete.json',dict(state='complete',cases=len(cases),native_accuracy_records=len(complete),preparation_seconds=time.time()-START))
    return cases

def gpu_measure(m,case,d,session,stream,torch):
    prefix=f"{m['model_id']}_{case['scene']}_N{case['N']}";paths=[R/'gpu_timings'/(prefix+'_'+mode+'.json') for mode in ['resident','cpu_gpu_cpu']]
    if all(p.exists() for p in paths):return
    k=int(m['stencil_size']);a=(11-k)//2;x=np.ascontiguousarray(d['stencils'][:,a:a+k,a:a+k],dtype=np.float32);b=len(x)
    host=torch.from_numpy(x).pin_memory();hc=torch.full((b,),.5,dtype=torch.float32,pin_memory=True);ho=torch.empty(b,dtype=torch.float32,pin_memory=True)
    chunks=[]
    with torch.cuda.stream(stream):
        for start in range(0,b,CHUNK):
            stop=min(start+CHUNK,b);xx=host[start:stop].to('cuda',non_blocking=True);cc=hc[start:stop].to('cuda',non_blocking=True);out=torch.empty(stop-start,device='cuda')
            binding=session.io_binding()
            for inp in session.get_inputs():
                val=cc if len(inp.shape)==1 else xx
                binding.bind_input(inp.name,'cuda',0,np.float32,tuple(val.shape),val.data_ptr())
            binding.bind_output(session.get_outputs()[0].name,'cuda',0,np.float32,(stop-start,),out.data_ptr())
            chunks.append((start,stop,xx,cc,out,binding))
    stream.synchronize()
    def invoke(mode):
        with torch.cuda.stream(stream):
            for start,stop,xx,cc,out,binding in chunks:
                if mode=='cpu_gpu_cpu':
                    xx.copy_(host[start:stop],non_blocking=True)
                    if m['operator']=='flux':cc.copy_(hc[start:stop],non_blocking=True)
                session.run_with_iobinding(binding)
                if mode=='cpu_gpu_cpu':ho[start:stop].copy_(out,non_blocking=True)
        stream.synchronize()
    invoke('cpu_gpu_cpu');pred=ho.numpy().copy();accuracy=metrics(pred,d,m['operator'],case['N'])
    predfile=R/'predictions'/(prefix+'.npz');predfile.parent.mkdir(exist_ok=True);np.savez_compressed(predfile,prediction=pred)
    for mode,dest in zip(['resident','cpu_gpu_cpu'],paths):
        if dest.exists():continue
        before=smi()
        if foreign_processes(before):raise InterruptedError('Foreign CUDA process appeared before formal timing')
        for _ in range(WARMUP):invoke(mode)
        samples=[]
        for _ in range(REPEATS):
            t=time.perf_counter_ns();invoke(mode);samples.append((time.perf_counter_ns()-t)/1e6)
        after=smi();qualified=not foreign_processes(after)
        result=dict(protocol=PROTOCOL,case=case,model=m,operator=m['operator'],device='GPU',mode=mode,precision='float32 input/native artifact arithmetic',formal=qualified,
            wall_ms=stats(samples),raw_wall_ms=samples,metrics=accuracy,prediction=str(predfile.relative_to(ROOT)),chunk_cap=CHUNK,chunks=len(chunks),warmup=WARMUP,repeats=REPEATS,
            input_location='GPU resident' if mode=='resident' else 'preallocated pinned CPU buffer',output_location='GPU completed' if mode=='resident' else 'CPU returned',
            includes_transfers=mode=='cpu_gpu_cpu',input_gathering_excluded=True,device_before=before,device_after=after)
        if qualified:save(dest,result)
        else:
            save(R/'attempts'/(dest.stem+'_'+str(int(time.time()))+'.json'),result)
            raise InterruptedError('Foreign CUDA process appeared;case retained as nonformal attempt')
    del chunks,host,hc,ho;gc.collect();torch.cuda.empty_cache()

def formal(cases):
    import torch
    import importlib.util
    hpath=ROOT/'operator-ablation-20261122/02-inference-pareto/harness.py'
    spec=importlib.util.spec_from_file_location('validated_large_grid_session',hpath);h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h);h.HERE=R
    torch.set_num_threads(1);torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    os.sched_setaffinity(0,{0})
    stream=torch.cuda.Stream();models=registry();random.Random(202611250402).shuffle(models)
    cpu_jobs=[(c,b,o,w) for c in cases for b in ['aphros','simliquid'] for o in ['flux','curvature'] for w in [1,8,min(24,os.cpu_count())]]
    random.Random(202611250403).shuffle(cpu_jobs)
    for i,(case,backend,op,workers) in enumerate(cpu_jobs):
        status('formal_cpu_deployment',case=case['scene'],N=case['N'],backend=backend,operator=op,workers=workers,completed=i,total=len(cpu_jobs))
        d=load(case);cpu_point(case,backend,op,workers,d,True);del d;gc.collect()
    for mi,m in enumerate(models):
        if foreign_processes(smi()):raise InterruptedError('Foreign GPU process appeared between models')
        session=h.session(ROOT/m['onnx_path'],stream);assert len(session.get_outputs())==1
        for case in cases:
            status('formal_gpu_deployment',model=m['model_id'],model_index=mi,total_models=len(models),case=case['scene'],N=case['N'])
            d=load(case);gpu_measure(m,case,d,session,stream,torch);del d;gc.collect()
        del session;gc.collect();torch.cuda.empty_cache()
    save(R/'completion.json',dict(state='complete',cpu_cases=len(cpu_jobs),gpu_points=len(models)*len(cases)*2,models=len(models),N=[512,1024,2048,4096],formal_compute_complete=True))
    status('compute_complete',cpu_cases=len(cpu_jobs),gpu_points=len(models)*len(cases)*2,report_pending=True)

def main():
    R.mkdir(exist_ok=True);lock=(R/'controller.lock').open('a+')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:print('Existing controller owns this study');return
    (R/'controller.pid').write_text(str(os.getpid()));cases=prepare()
    if (R/'completion.json').exists():status('compute_complete',report_pending=True);return
    while True:
        idle=wait_idle()
        with gpu_lease('large-grid-gpu-cpu',formal_timing=True):
            if foreign_processes(smi()):continue
            try:formal(cases);break
            except InterruptedError as error:status('formal_timing_interrupted',reason=str(error),completed_cases_retained=True)

if __name__=='__main__':
    R.mkdir(exist_ok=True)
    sys.stdout=(R/'controller.log').open('a',buffering=1);sys.stderr=sys.stdout
    try:main()
    except Exception as error:status('failed',error=repr(error),traceback=traceback.format_exc());raise
