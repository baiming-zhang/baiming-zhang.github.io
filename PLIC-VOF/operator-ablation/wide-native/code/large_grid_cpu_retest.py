"""CPU-only retest with recorded host activity and an optional fixed-core policy."""
from pathlib import Path
import os,sys,json,time,threading,subprocess,collections,traceback
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'gpu-cpu-large-grid-20261125'
PINNED='--pinned' in sys.argv
RET=BASE/('cpu-retest-pinned' if PINNED else 'cpu-retest')
DEPENDENCIES=BASE/'cpu-retest/dependencies'
POOL=list(range(8))

def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(data,indent=2),encoding='utf-8')
    for i in range(30):
        try:temp.replace(path);return
        except PermissionError:
            if i==29:raise
            time.sleep(.02)

def read(path):
    try:return json.loads(path.read_text(encoding='utf-8'))
    except (OSError,ValueError):return {}

def host_monitor():
    sys.path.insert(0,str(DEPENDENCIES))
    import psutil
    own=os.getpid();previous={};samples=collections.deque(maxlen=120);events=RET/'host-load.jsonl'
    while True:
        per_core=psutil.cpu_percent(interval=1,percpu=True);now=time.time();current={};external=[]
        for proc in psutil.process_iter(['pid','name','cpu_times','create_time']):
            try:
                data=proc.info
                if data['pid']==0 or not data['cpu_times']:continue
                key=(data['pid'],data['create_time']);cpu=data['cpu_times'].user+data['cpu_times'].system
                current[key]=(cpu,now)
                if key not in previous:continue
                before,then=previous[key];pct=max(0,100*(cpu-before)/max(now-then,.1))/len(per_core)
                if data['pid']!=own and not (data['name'] or '').lower().startswith('vmmem'):
                    external.append(dict(pid=data['pid'],name=data['name'],percent=pct))
            except (psutil.NoSuchProcess,psutil.AccessDenied,TypeError):pass
        previous=current
        sample=dict(epoch=now,total_percent=sum(per_core)/len(per_core),per_core=per_core,
            external_percent=sum(p['percent'] for p in external),largest_external_percent=max((p['percent'] for p in external),default=0),
            top_external=sorted(external,key=lambda x:x['percent'],reverse=True)[:5])
        samples.append(sample)
        window=[s for s in samples if now-s['epoch']<=92]
        enough=len(window)>=85 and now-window[0]['epoch']>=89
        quiet=False
        if enough:
            totals=sorted(s['total_percent'] for s in window)
            per_mean=[sum(s['per_core'][i] for s in window)/len(window) for i in range(len(per_core))]
            quiet=sum(totals)/len(totals)<=5 and totals[int(.95*(len(totals)-1))]<=10 and max(per_mean)<=25
        sample.update(idle_90_seconds=quiet,idle_window_seconds=now-window[0]['epoch'],
            timing_load_ok=sample['external_percent']<=3 and sample['largest_external_percent']<=2)
        write(RET/'host_cpu.json',sample)
        with events.open('a',encoding='utf-8') as f:f.write(json.dumps(sample)+'\n')

def windows_main():
    RET.mkdir(exist_ok=True)
    sys.path.insert(0,str(DEPENDENCIES));import psutil
    pidfile=RET/'host-controller.pid'
    if pidfile.exists():
        try:
            other=psutil.Process(int(pidfile.read_text()))
            if other.is_running() and any('large_grid_cpu_retest.py' in x for x in other.cmdline()):
                print('Existing host controller retained');return
        except (psutil.Error,ValueError):pass
    pidfile.write_text(str(os.getpid()))
    threading.Thread(target=host_monitor,daemon=True).start()
    linux_script='/mnt/c/Users/ASUS/Desktop/Neural PLIC-VOF/agent_experiments/large_grid_cpu_retest.py'
    with (RET/'worker.log').open('ab',buffering=0) as log:
        child=subprocess.Popen(['wsl.exe','-d','cold-surge-ubuntu-24.04','--','/root/.venvs/neural-plic-vof/bin/python','-B',linux_script,'--worker',*(['--pinned'] if PINNED else [])],stdout=log,stderr=log)
        write(RET/'host-controller.json',dict(pid=os.getpid(),wsl_process_pid=child.pid,state='running'))
        code=child.wait()
        write(RET/'host-controller.json',dict(pid=os.getpid(),wsl_process_pid=child.pid,state='complete' if code==0 else 'failed',exit_code=code))
        if code:raise RuntimeError('CPU worker failed; see retained worker.log')

def proc_times():
    result={}
    for folder in Path('/proc').iterdir():
        if not folder.name.isdigit() or int(folder.name)==os.getpid():continue
        try:
            fields=(folder/'stat').read_text().rsplit(')',1)[1].split()
            result[int(folder.name)]=(int(fields[11])+int(fields[12]))/os.sysconf('SC_CLK_TCK')
        except (OSError,ValueError,IndexError):pass
    return result

def vcpu_times():
    result={}
    for line in Path('/proc/stat').read_text().splitlines():
        parts=line.split()
        if not parts or not parts[0].startswith('cpu') or not parts[0][3:].isdigit():continue
        core=int(parts[0][3:]);values=list(map(int,parts[1:9]));total=sum(values)
        result[core]=(total-values[3]-values[4],total)
    return result

def worker_main():
    import fcntl,random,gc,numpy as np
    RET.mkdir(exist_ok=True);lock=(RET/'worker.lock').open('a+')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:print('Existing CPU retest worker retained');return
    (RET/'worker.pid').write_text(str(os.getpid()))
    import large_grid_gpu_cpu as h
    h.R=RET
    if PINNED:
        original_native=h.Native
        class PinnedNative(original_native):
            def part(self,i):
                affinity=os.sched_getaffinity(0)
                assert len(affinity)==1 and affinity.issubset(POOL),('worker affinity',affinity)
                return super().part(i)
        h.Native=PinnedNative
    started=time.time();linux=dict(epoch=0,external_core_fraction=100);samples=collections.deque(maxlen=100)
    def sample_linux():
        try:os.sched_setaffinity(0,{max(os.sched_getaffinity(0))})
        except OSError:pass
        previous=proc_times();previous_vcpu=vcpu_times();then=time.monotonic()
        while True:
            time.sleep(1);current=proc_times();now=time.monotonic()
            fraction=sum(max(0,v-previous.get(pid,v)) for pid,v in current.items())/(now-then)
            vcpu=vcpu_times();usage={}
            for core in POOL:
                a,b=vcpu[core];pa,pb=previous_vcpu[core];usage[str(core)]=100*(a-pa)/max(b-pb,1)
            linux.update(epoch=time.time(),external_core_fraction=fraction,allocated_vcpu_percent=usage)
            samples.append(dict(linux));previous=current;previous_vcpu=vcpu;then=now
    threading.Thread(target=sample_linux,daemon=True).start();os.sched_setaffinity(0,{0})
    cases=read(BASE/'fixture_manifest.json')['cases']
    worker_counts=[1,4,8] if PINNED else [1,8,min(24,os.cpu_count())]
    jobs=[(c,b,o,w) for c in cases for b in ['aphros','simliquid'] for o in ['flux','curvature'] for w in worker_counts]
    random.Random(202611250601).shuffle(jobs)
    def state(name,**kwargs):write(RET/'status.json',dict(state=name,pid=os.getpid(),elapsed_seconds=time.time()-started,**kwargs))
    def ready(initial=False):
        host=read(RET/'host_cpu.json');fresh=time.time()-host.get('epoch',0)<5 and time.time()-linux['epoch']<5
        if not fresh:return False
        if PINNED:return True
        if initial:return host.get('idle_90_seconds',False) and len(samples)>=85 and max(s['external_core_fraction'] for s in list(samples)[-85:])<=.25
        return host.get('timing_load_ok',False) and linux['external_core_fraction']<=.25
    if (RET/'completion.json').exists():state('complete');return
    while not ready(initial=True):
        state('waiting_for_quiet_cpu',completed=len(list((RET/'cpu_timings').glob('*.json'))),total=96,
            host=read(RET/'host_cpu.json'),linux=dict(linux),remaining_seconds=None)
        time.sleep(5)
    write(RET/('affinity_policy_start.json' if PINNED else 'initial_idle_gate.json'),dict(state='passed',host=read(RET/'host_cpu.json'),linux=list(samples),whole_host_idle_required=not PINNED,cpu_pool=POOL if PINNED else None))
    for case,backend,op,workers in jobs:
        key=f"{case['scene']}_N{case['N']}_{backend}_{op}_cpu{workers}"
        dest=RET/'cpu_timings'/(key+'.json')
        if dest.exists():
            existing=read(dest)
            if existing.get('retest') and existing.get('formal'):continue
            attempts=RET/'attempts';attempts.mkdir(exist_ok=True)
            dest.rename(attempts/(key+'_unqualified_'+str(time.time_ns())+'.json'))
        while True:
            while not ready():
                state('waiting_for_quiet_cpu',case=key,completed=len(list((RET/'cpu_timings').glob('*.json'))),total=96,host=read(RET/'host_cpu.json'),linux=dict(linux),remaining_seconds=None)
                time.sleep(2)
            state('formal_cpu_retest',case=key,completed=len(list((RET/'cpu_timings').glob('*.json'))),total=96)
            before=read(RET/'host_cpu.json');linux_before=dict(linux);begin=time.time()
            d=h.load(case);record=h.cpu_point(case,backend,op,workers,d,True);del d;gc.collect()
            after=read(RET/'host_cpu.json');linux_after=dict(linux)
            contaminated=not ready() or not before.get('timing_load_ok',False)
            record.update(retest=True,host_before=before,host_after=after,linux_before=linux_before,linux_after=linux_after,measurement_start_epoch=begin,measurement_end_epoch=time.time(),formal=not contaminated)
            if PINNED:
                contaminated=not ready();record.update(formal=not contaminated,cpu_resource_policy='fixed_pool_shared_host',cpu_pool=POOL,worker_count_set=worker_counts,benchmark_affinity_verified=True,scheduling='one pinned logical-vCPU per compute worker,within0..7;1/4/8 workers measured;shared host,not exclusive physical cores')
            baseline=np.load(BASE/'predictions'/f"{case['scene']}_N{case['N']}_{backend}_{op}.npz")['prediction']
            pred=np.load(ROOT/record['prediction'])['prediction']
            assert np.allclose(baseline,pred,rtol=1e-12,atol=1e-12,equal_nan=True)
            write(dest,record)
            if not contaminated:break
            attempts=RET/'attempts';attempts.mkdir(exist_ok=True)
            dest.rename(attempts/(key+'_'+str(time.time_ns())+'.json'))
            state('load_contaminated_attempt_retained',case=key)
            time.sleep(2)
    write(RET/'completion.json',dict(state='complete',cpu_configurations=96,gpu_reused=True,initial_results_preserved=True,elapsed_seconds=time.time()-started,cpu_pool=POOL if PINNED else None,worker_counts=worker_counts,whole_host_idle_required=not PINNED))
    state('complete',completed=96,total=96)

if __name__=='__main__':
    sys.dont_write_bytecode=True
    try:worker_main() if '--worker' in sys.argv else windows_main()
    except Exception as error:
        write(RET/('worker_failure.json' if '--worker' in sys.argv else 'host_failure.json'),dict(error=repr(error),traceback=traceback.format_exc()))
        raise
