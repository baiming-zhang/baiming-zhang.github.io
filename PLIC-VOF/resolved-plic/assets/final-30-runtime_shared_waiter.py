"""One full preregistered timing-order correction;frozen model/data/runtime."""
from pathlib import Path
import os,sys
sys.dont_write_bytecode=True
for key in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
import argparse,fcntl,gc,hashlib,importlib.util,json,random,threading,time,traceback
import numpy as np
R=Path(__file__).resolve().parent;P=R.parent;W=P.parent;OUT=R/'runtime_shared_results'
spec=importlib.util.spec_from_file_location('_formal_reused_helpers',P/'validation_pareto.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
sys.path.insert(0,str(W/'operator-ablation-20261122/shared'));from gpu_lease import gpu_lease
WARMUP=10;REPEATS=31;ORDER_SEED=202611270452
class Contamination(RuntimeError):pass
class StableMonitor:
    """Before snapshot precedes warmup;only logical9 polls asynchronously."""
    def __init__(self):
        self.samples=[];self.phase='before_warmup';self.stop_event=threading.Event();self.thread=threading.Thread(target=self.loop,daemon=True)
    def loop(self):
        os.sched_setaffinity(0,{9})
        while not self.stop_event.wait(.25):
            phase=self.phase;begin=time.perf_counter_ns();sample=h.smi()
            sample.update(phase_at_query_start=phase,phase_at_query_end=self.phase,query_start_ns=begin,query_end_ns=time.perf_counter_ns(),telemetry_thread_affinity=sorted(os.sched_getaffinity(0)))
            self.samples.append(sample)
    def __enter__(self):
        self.before=h.smi();self.thread.start();self.phase='warmup';return self
    def __exit__(self,*args):
        self.phase='after_measurement';self.stop_event.set();self.thread.join();self.after=h.smi()
    def evidence(self):return dict(before=self.before,during=self.samples,after=self.after,telemetry_logical_core=9,poll_interval_seconds=.25)
    def clean(self):return all(clean(s) for s in [self.before,*self.samples,self.after])
def continuous_warmup(call):
    start=time.perf_counter_ns();count=0
    while count<10 or time.perf_counter_ns()-start<500_000_000:
        call();count+=1
    end=time.perf_counter_ns();return count,(end-start)/1e9,end

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def save(path,obj):h.save(path,obj)
def pathof(value):
    value=str(value).replace('\\','/')
    if len(value)>2 and value[1]==':':value='/mnt/'+value[0].lower()+value[2:]
    p=Path(value);return p if p.is_absolute() else W/p
def status(stage,**extra):save(R/'runtime_shared_status.json',dict(state=stage,pid=os.getpid(),updated_utc=h.now(),**extra))
def metric(pred,truth,f,swept):
    base=h.metric(pred,truth,f,'flux');good=np.isfinite(pred)
    exact=float(np.mean(np.abs(f[good]*pred[good]-swept[good]))) if good.any() else None
    base.update(ratio_weighted_volume_mae=base['error'],metric='exact_cell_transport_mae',error=exact if good.all() else None,
        exact_cell_transport_mae=exact if good.all() else None,transported_volume_mae=exact if good.all() else None,available_only_exact_cell_transport_mae=exact,
        definition='mean(abs(rawf32_center_fraction*predicted_ratio-exact_swept_cell_area));unitcell area1,no h^2 factor')
    return base
def clean(s):return s['telemetry_ok'] and not s['foreign_processes']
def idle(s):
    try:return clean(s) and bool(s['gpu']['stdout']) and all(float(row.split(',')[2])<=5 for row in s['gpu']['stdout'].splitlines())
    except (ValueError,IndexError):return False
def wait_idle(seconds,once):
    begin=time.monotonic();first=None;history=[]
    while True:
        s=h.smi();history.append(s);first=(first if first is not None else time.monotonic()) if idle(s) else None
        elapsed=time.monotonic()-first if first is not None else 0
        status('waiting_independent_GPU_idle',idle_seconds=elapsed,required_seconds=90,waited_seconds=time.monotonic()-begin,latest=s,formal_started=False)
        if elapsed>=90:return history
        if once or time.monotonic()-begin>=seconds:
            save(OUT/'last_wait_attempt.json',dict(state='waiting_external_GPU_or_idle_window',samples=history,formal_started=False));return None
        time.sleep(min(5,max(.1,seconds-(time.monotonic()-begin))))
def frozen_inputs():
    lockpath=P/'deployment_selection_locked.json';manifestpath=P/'flux_blind_test/large_grid/manifest.json'
    if not lockpath.exists() or not manifestpath.exists():
        status('waiting_locked_new_grid_manifest',lock_exists=lockpath.exists(),fixture_manifest_exists=manifestpath.exists(),GPU_started=False);return None
    lock=read(lockpath);manifest=read(manifestpath)
    assert lock['operators']==['flux'] and lock['seeds']==[4101,4102,4103] and lock['all_seeds_retained']
    assert manifest['state']=='complete' and manifest['identity']['seed']==202611270400
    assert manifest['identity']['deployment_lock_sha256']==sha(lockpath)
    assert sha(pathof(lock['runtime_adapter_path']))==lock['runtime_adapter_sha256']
    assert sha(pathof(lock['runtime_core_path']))==lock['runtime_core_sha256']
    assert sha(pathof(lock['test_protocol_path']))==lock['test_protocol_sha256']
    models=sorted(lock['models'],key=lambda m:m['seed']);assert [m['seed'] for m in models]==[4101,4102,4103]
    for m in models:
        assert sha(pathof(m['torchscript_path']))==m['torchscript_sha256'] and sha(pathof(m['checkpoint_path']))==m['checkpoint_sha256']
    cases=manifest['cases'];assert len(cases)==10
    assert {(c['scene'],c['N']) for c in cases}=={(scene,n) for scene in ['single_drop','many_drops'] for n in [512,1024,2048,4096,8192]}
    for case in cases:
        source=pathof(case['source']);assert 'flux_blind_test/large_grid/' in str(source).replace('\\','/')
        assert case['genuine_field'] and case['minimum_radius_cells']>=4 and sha(source)==case['path_sha256']
    identity=dict(deployment_lock_sha256=sha(lockpath),new_grid_manifest_sha256=sha(manifestpath),test_protocol_sha256=lock['test_protocol_sha256'],
        adapter_sha256=lock['runtime_adapter_sha256'],core_sha256=lock['runtime_core_sha256'],formal_script_sha256=sha(__file__),formal_plan_sha256=sha(R/'runtime_formal_plan.json'),helper_sha256=sha(P/'validation_pareto.py'),timing_recovery_protocol_sha256=sha(R/'timing_recovery_protocol.json'),preliminary_completion_sha256=sha(R/'runtime_formal_results/completion.json'),shared_execution_policy_sha256=sha(R/'shared_execution_policy.json'))
    return lock,models,cases,identity
def load_fields(cases):
    result={}
    for case in cases:
        with np.load(pathof(case['source']),allow_pickle=False) as z:
            x=np.ascontiguousarray(z['stencils'][:,2:9,2:9],dtype=np.float32);truth=np.asarray(z['flux_ratio_truth'],dtype=np.float64)
            swept=np.asarray(z['flux_density_truth'],dtype=np.float64)*.5
            ix=z['ix'];iy=z['iy'];assert len(np.unique(np.stack([ix,iy],1),axis=0))==len(x)
        assert len(x)==case['active_interface_cells'] and x.shape[1:]==(7,7)
        result[f"{case['scene']}_N{case['N']}"]=dict(case=case,x=x,c=np.full(len(x),.5,dtype=np.float32),truth=truth,swept=swept)
    return result
def cpu_case(fid,d,backend,workers,identity):
    target=OUT/'cpu_timings'/f'{fid}_{backend}_W{workers}.json'
    if target.exists():
        old=read(target)
        if old['identity']==identity and old.get('measurement_valid',False):return old
    os.sched_setaffinity(0,{0});before=h.machine();obj=h.CPUCase(backend,'flux',workers,d['x'],d['c'],[np.arange(len(d['x']))])
    try:
        samples=[];ratios=[]
        with StableMonitor() as monitor:
            warm_count,warm_seconds,warm_end=continuous_warmup(lambda:obj(0))
            monitor.phase='measurement';first_start=None
            for _ in range(REPEATS):
                t=time.perf_counter_ns()
                if first_start is None:first_start=t
                r=obj(0);samples.append((time.perf_counter_ns()-t)/1e6);ratios.extend(r)
            gap_ns=first_start-warm_end
        pred=obj.banks[0][2].copy();metrics=metric(pred,d['truth'],d['x'][:,3,3].astype(np.float64),d['swept'])
        record=dict(identity=identity,fixture_id=fid,scene=d['case']['scene'],N=d['case']['N'],actual_queries=len(pred),total_grid_cells=d['case']['total_grid_cells'],device='CPU',backend=backend,workers=workers,
            native_method='actual_original_PLIC_normal_plane_flux',precision='same_rawf32_promoted_float64',cpu_pool=list(range(8)),worker_affinity_verified=True,
            warmup=WARMUP,repeats=REPEATS,wall_ms=h.stats(samples),raw_wall_ms=samples,metrics=metrics,formal=False,measurement_valid=True,execution_policy='user_authorized_shared_GPU',exclusive_hardware_claim=False,
            external_load=monitor.evidence(),resource_before=before,resource_after=h.machine(),worker_cpu_wall_ratio_median=float(np.median(ratios)),library_sha256=obj.library_sha,
            duplicated_queries=False,boundary='CPU-ready promotedVOF->requirednative normal+PLIC+flux->completedCPU output;persistentdispatch included')
        record.update(batch=len(pred),error=metrics['error'],median_ms=record['wall_ms']['median'],p05_ms=record['wall_ms']['p05'],p95_ms=record['wall_ms']['p95'])
        record.update(warmup=warm_count,warmup_actual_calls=warm_count,warmup_seconds=warm_seconds,warmup_minimum_calls=10,warmup_minimum_seconds=.5,warmup_to_measure_gap_ns=gap_ns,telemetry_thread_affinity=[9],timing_round='single_complete_method_recovery_primary')
        save(target,record)
        assert record['repeats']==31 and metrics['coverage']==1,'Incomplete CPU measurement'
        return record
    finally:obj.close()
def gpu_case(fid,d,m,mode,adapter_module,identity):
    target=OUT/'gpu_timings'/f"{fid}_{m['model_id']}_{mode}.json"
    if target.exists():
        old=read(target)
        if old['identity']==identity and old.get('measurement_valid',False):return old
    probe=h.smi()
    if not probe['telemetry_ok']:raise RuntimeError('GPU telemetry unavailable before shared capture')
    os.sched_setaffinity(0,{8});x=np.ascontiguousarray(d['x'][:,1:6,1:6]);before=h.machine()
    cold=time.perf_counter();adapter=adapter_module.FluxCUDAGraphAdapter(pathof(m['torchscript_path']),x,d['c'],mode=mode,expected_artifact_sha256=m['torchscript_sha256'],warmup=12)
    cold_seconds=time.perf_counter()-cold
    try:
        samples=[]
        with StableMonitor() as monitor:
            warm_count,warm_seconds,warm_end=continuous_warmup(adapter.run)
            monitor.phase='measurement';first_start=None
            for _ in range(REPEATS):
                start=time.perf_counter_ns()
                if first_start is None:first_start=start
                last=adapter.run();samples.append((time.perf_counter_ns()-start)/1e6)
            gap_ns=first_start-warm_end
        pred=last.detach().cpu().numpy().copy();metrics=metric(pred.astype(np.float64),d['truth'],d['x'][:,3,3].astype(np.float64),d['swept'])
        record=dict(identity=identity,fixture_id=fid,scene=d['case']['scene'],N=d['case']['N'],actual_queries=len(pred),total_grid_cells=d['case']['total_grid_cells'],device='GPU',model_id=m['model_id'],seed=m['seed'],mode=mode,
            runtime='locked_CUDAgraph_FP32',artifact_sha256=m['torchscript_sha256'],checkpoint_sha256=m['checkpoint_sha256'],adapter_metadata=adapter.metadata,
            capture_seconds=adapter.metadata['capture_seconds'],capture_excluded_from_hot_latency=True,warmup=WARMUP,repeats=REPEATS,wall_ms=h.stats(samples),raw_wall_ms=samples,
            metrics=metrics,formal=False,measurement_valid=True,execution_policy='user_authorized_shared_GPU',exclusive_hardware_claim=False,external_load=monitor.evidence(),resource_before=before,resource_after=h.machine(),duplicated_queries=False,
            boundary=adapter.metadata['boundary'],fixed_actual_shape=True,padding=False)
        record.update(batch=len(pred),error=metrics['error'],median_ms=record['wall_ms']['median'],p05_ms=record['wall_ms']['p05'],p95_ms=record['wall_ms']['p95'],
            cold_load_and_capture_seconds=cold_seconds,capture_cost_definition='capture_seconds includes12 preparation calls+graphcapture;full constructor cold_load_and_capture_seconds also includes artifactload and buffers',
            resident_output_CPU_copy_after_timer=mode=='resident')
        record.update(warmup=warm_count,warmup_actual_calls=warm_count,warmup_seconds=warm_seconds,warmup_minimum_calls=10,warmup_minimum_seconds=.5,warmup_to_measure_gap_ns=gap_ns,telemetry_thread_affinity=[9],timing_round='single_complete_method_recovery_primary')
        save(target,record)
        assert record['repeats']==31 and metrics['coverage']==1,'Incomplete GPU measurement'
        return record
    finally:adapter.close();del adapter;gc.collect();adapter_module.torch.cuda.empty_cache()
def aggregate(cases,models,cpu,gpu,identity):
    groups=[];iso=[]
    for case in cases:
        fid=f"{case['scene']}_N{case['N']}";refs=[r for r in cpu if r['fixture_id']==fid];assert len(refs)==6
        for mode in ['resident','cpu_gpu_cpu']:
            points=[r for r in gpu if r['fixture_id']==fid and r['mode']==mode];assert len(points)==3
            times=[p['wall_ms']['median'] for p in points];errors=[p['metrics']['error'] for p in points]
            g=dict(fixture_id=fid,scene=case['scene'],N=case['N'],actual_queries=case['active_interface_cells'],mode=mode,model_ids=[p['model_id'] for p in points],
                mean_ms=float(np.mean(times)),seed_sd_ms=float(np.std(times,ddof=1)),max_seed_ms=max(times),min_seed_ms=min(times),
                volume_mae_mean=float(np.mean(errors)) if all(e is not None for e in errors) else None,volume_mae_seed_sd=float(np.std(errors,ddof=1)) if all(e is not None for e in errors) else None,
                max_seed_volume_mae=max(errors) if all(e is not None for e in errors) else None,all_seed_full_coverage=all(p['metrics']['coverage']==1 for p in points),formal=False)
            g.update(batch=case['active_interface_cells'],metric='exact_cell_transport_mae',error_mean=g['volume_mae_mean'],error_seed_sd=g['volume_mae_seed_sd'],time_mean_ms=g['mean_ms'],time_seed_sd_ms=g['seed_sd_ms'],member_ids=g['model_ids'])
            groups.append(g)
            for backend in ['aphros','simliquid','fastest_tested_native']:
                choices=refs if backend=='fastest_tested_native' else [r for r in refs if r['backend']==backend]
                reference=min(choices,key=lambda r:r['wall_ms']['median']);error=reference['metrics']['error'];complete=reference['metrics']['coverage']==1 and g['all_seed_full_coverage']
                iso.append(dict(fixture_id=fid,scene=case['scene'],N=case['N'],actual_queries=case['active_interface_cells'],mode=mode,reference_scope=backend,native_backend=reference['backend'],native_workers=reference['workers'],native_ms=reference['wall_ms']['median'],native_volume_mae=error,
                    group=g,equal_time_reachable=bool(complete and g['mean_ms']<=reference['wall_ms']['median']),equal_accuracy_reachable=bool(complete and g['volume_mae_mean']<=error),
                    all_seeds_fit_native_time=bool(complete and g['max_seed_ms']<=reference['wall_ms']['median']),all_seeds_reach_native_error=bool(complete and g['max_seed_volume_mae']<=error),
                    full_task_coverage=complete,formal=False,no_interpolation=True))
    result=dict(state='complete',identity=identity,operator='flux',formal_speed_evidence=False,measurement_valid=True,execution_policy='user_authorized_shared_GPU',exclusive_hardware_claim=False,only_frozen_three_seed_candidate=True,models=models,cases=cases,actual_unique_queries=True,
        native_points=cpu,per_seed_error_time=gpu,groups=groups,equal_time_equal_error=iso,warmup_minimum_calls=WARMUP,warmup_minimum_seconds=.5,repeats=REPEATS,timing_role='sole_primary_corrected_round_regardless_outcome',old_round_role='runtime_formal_results preserved as preliminary method diagnosis;never take minima across rounds',
        primary_metric='exact_cell_transport_mae=mean(abs(rawf32_center*pred_ratio-exact_geometric_swept_cell_area));no h^2',auxiliary_metric='ratio_weighted_volume_mae=f32center*abs(pred_ratio-analytic_ratio_truth),retained separately',
        limitations=['GPUFP32 versus nativeCPU64 is the prespecified deployment comparison,not precision-normalized algorithmic superiority.','CPU compute workers use0..7 logical-vCPU pool on shared host,not exclusive physical cores.','Every raw query/output failure is retained;incomplete tasks cannot claim an error target.','Cold capture/load excluded from hot single-call latency and capture_seconds separately reported.','Independent fixtures created only after immutable deployment lock;no test-based model,runtime,seed or batch changes.'])
    save(OUT/'summary.json',result);return result
def main(args):
    OUT.mkdir(exist_ok=True);guard=(R/'runtime_shared.lock').open('a+')
    try:fcntl.flock(guard,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:print('Existing bounded formal waiter retained');return
    (R/'runtime_shared.pid').write_text(str(os.getpid()));inputs=frozen_inputs()
    if inputs is None:return
    lock,models,cases,identity=inputs;completion=OUT/'completion.json'
    assert set(range(10)).issubset(h.ALLOWED_AFFINITY),'Require compute0..8 andtelemetry9'
    recovery=read(R/'timing_recovery_protocol.json');assert recovery['registered_before_execution']
    if completion.exists() and read(completion)['identity']==identity:
        status('complete',cached=True,formal_speed_evidence=False,measurement_valid=True,execution_policy='user_authorized_shared_GPU',exclusive_hardware_claim=False);return
    registration=OUT/'registration.json'
    if registration.exists():assert read(registration)['identity']==identity,'Formal binding changed;preserve artifacts and resolve explicitly'
    else:save(registration,dict(identity=identity,cases=cases,model_ids=[m['model_id'] for m in models],created_utc=h.now(),independent_new_seed=202611270400))
    history=[h.smi()]
    assert read(R/'shared_execution_policy.json')['foreign_GPU_processes_allowed']
    status('starting_user_authorized_shared_measurement',idle_gate_required=False,execution_policy='user_authorized_shared_GPU')
    with gpu_lease('frozen-flux-user-authorized-shared-measurement',formal_timing=False):
        probe=h.smi()
        if not probe['telemetry_ok']:raise RuntimeError('GPU telemetry unavailable before shared execution')
        save(OUT/'shared_execution.json',dict(user_authorized=True,idle_gate_required=False,exclusive_load=False,samples=history))
        fields=load_fields(cases);cpu=[];gpu=[]
        adapterpath=pathof(lock['runtime_adapter_path']);spec=importlib.util.spec_from_file_location('_final_locked_graph_adapter',adapterpath);a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
        a.torch.set_num_threads(1);a.torch.set_num_interop_threads(1)
        jobs=[]
        for fid in fields:
            jobs.extend(('cpu',fid,b,w) for b in ['aphros','simliquid'] for w in [1,4,8])
            jobs.extend(('gpu',fid,m,mode) for m in models for mode in ['resident','cpu_gpu_cpu'])
        random.Random(ORDER_SEED).shuffle(jobs);save(OUT/'case_order.json',dict(seed=ORDER_SEED,jobs=jobs))
        for i,(kind,fid,first,second) in enumerate(jobs):
            status('shared_measurement_running',completed=i,total=len(jobs),fixture_id=fid,device=kind)
            if kind=='cpu':cpu.append(cpu_case(fid,fields[fid],first,second,identity))
            else:gpu.append(gpu_case(fid,fields[fid],first,second,a,identity))
        aggregate(cases,models,cpu,gpu,identity)
        save(completion,dict(state='complete',identity=identity,CPU_configurations=len(cpu),GPU_configurations=len(gpu),models=3,fixtures=10,formal_speed_evidence=False,measurement_valid=True,execution_policy='user_authorized_shared_GPU',exclusive_hardware_claim=False,
            measured_invocations=len(jobs)*REPEATS,warmups=sum(t['warmup_actual_calls'] for t in cpu+gpu),warmup_minimum_calls=10,warmup_minimum_seconds=.5,summary='runtime_shared_results/summary.json'))
    status('complete',formal_speed_evidence=False,measurement_valid=True,execution_policy='user_authorized_shared_GPU',exclusive_hardware_claim=False,summary='runtime_shared_results/summary.json')
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--max-wait-seconds',type=float,default=300);parser.add_argument('--once',action='store_true');args=parser.parse_args();assert args.max_wait_seconds>0
    try:main(args)
    except Contamination as error:status('waiting_external_load_retry',reason=str(error),formal_completion=False);sys.exit(0)
    except Exception:status('failed',traceback=traceback.format_exc(),formal_completion=False);raise
