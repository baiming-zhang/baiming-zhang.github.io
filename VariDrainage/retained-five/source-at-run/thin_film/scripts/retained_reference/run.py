"""Frozen same-candidate references, held native handles, no automatic replay."""
from pathlib import Path
import argparse,csv,datetime,hashlib,json,os,subprocess,time
import numpy as np
HERE=Path(__file__).resolve().parent;R=HERE.parents[2];O=R/'results/retained_reference_20261113';RAW=R/'results/local-runs/retained_reference_20261113';LOG=R/'results/logs/retained_reference_20261113'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,v):p.write_text(json.dumps(v,indent=2,ensure_ascii=True,allow_nan=False)+'\n',encoding='utf-8')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def atomic(p,v):
    tmp=p.with_suffix('.tmp');save(tmp,v);tmp.replace(p)
def rows(p):
    if not p.exists():return []
    result=[]
    with p.open(encoding='utf-8') as stream:
        for r in csv.DictReader(stream):
            try:
                if None in r or any(v is None or v=='' for v in r.values()):continue
                data={k:float(v) for k,v in r.items()};assert all(np.isfinite(v) for v in data.values());result.append(data)
            except (ValueError,AssertionError):continue
    return result
def progress(root):
    h=rows(root/'history.csv');m=rows(root/'mechanics.csv');return dict(last_complete_history_step=int(h[-1]['step']) if h else 0,last_complete_mechanics_step=int(m[-1]['step']) if m else 0,last_force_ratio=m[-1]['g1_N']/m[-1]['target_N'] if m else None)
ENV={k:v for k,v in os.environ.items() if not k.startswith('VF_')};ENV['VF_THREADS']='2'
DLLS=sorted((Path(ENV['LOCALAPPDATA'])/'.xmake/packages').rglob('*.dll'));ENV['PATH']=';'.join(sorted({str(p.parent) for p in DLLS}))+';'+ENV['PATH']
def prepare():
    assert not O.exists() and not RAW.exists() and not LOG.exists()
    parent=R/'results/retained_motion_20261113';control=read(parent/'status.json');coverage=read(parent/'audit_coverage.json');p=read(parent/'protocol.json');assert coverage['actual_single_step_controls_qualified'] and len(control['cases'])==3 and all(c['passed'] for c in control['cases'])
    for rec in control['cases']:
        assert all(sha(R/f['path'])==f['sha256'] for f in rec['files']);d=parent/'audits/cases'/rec['name'];seal=read(d/'verified.json');assert all(sha(d/n)==h for n,h in seal['outputs'].items())
    assert all(sha(R/f)==h for f,h in p['sources'].items()) and all(sha(R/f)==h for f,h in p['original_sources'].items())
    build=read(parent/'build.json');exe=R/build['binary'];assert sha(exe)==build['sha256']
    O.mkdir();RAW.mkdir();LOG.mkdir();cases=[]
    for dx,N,cap1,cap2 in [(20,144,300,600),(10,544,1200,3600)]:
        for tag,dt,stride,cap in [('00625',.0625,20,cap1),('003125',.03125,40,cap2)]:cases.append(dict(name=f'u{dx}_dt{tag}',mode='uniform',dx_um=dx,N=N,dt_us=dt,T_us=10,radius_um=500,amplitude=0,c_bulk_mM=6,preparation='variational',save_stride=stride,refine_start_um=200,observer=True,timeout_s=cap))
    files=[Path(__file__),HERE/'PLAN.md',HERE/'monitor.py']
    save(O/'protocol.json',dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),prepared_utc=now(),source_hashes={**p['original_sources'],**p['sources'],**{str(f.relative_to(R)):sha(f) for f in files}},binary_sha256=sha(exe),executable=str(exe),control_coverage_sha256=sha(parent/'audit_coverage.json'),cases=cases,whole_native_cap_s=6600,monitor_interval_s=1800,runtime_dlls={str(f):sha(f) for f in DLLS},environment={'VF_THREADS':'2','caller_VF':'excluded'},gates=dict(force_ratio=1.000001,balance=1e-10,energy=1e-7,preparation=1e-8,time=.02,space=.1),scope='New same-candidate20/10um references at10us; old original trajectories not replayed; fullC1 incomplete'))
    save(O/'status.json',dict(state='prepared',cases=[],skipped=[],full_C1_complete=False));print('Prepared four same-candidate references; no native process launched.')
def run():
    p=read(O/'protocol.json');s=read(O/'status.json');assert s['state']=='prepared';assert all(sha(R/f)==h for f,h in p['source_hashes'].items());assert sha(Path(p['executable']))==p['binary_sha256'];assert all(sha(Path(f))==h for f,h in p['runtime_dlls'].items())
    s.update(state='running_matrix',started_utc=now(),parent_pid=os.getpid());save(O/'status.json',s);whole=time.monotonic();ready=True
    for case in p['cases']:
        name=case['name'];left=p['whole_native_cap_s']-(time.monotonic()-whole)
        if not ready or left<5:s['skipped'].append(dict(name=name,reason='Earlier native failure or whole native cap'));continue
        root=RAW/name;assert not root.exists();root.mkdir();log=LOG/(name+'.log');monitor=LOG/('monitor_'+name+'.jsonl');assert not log.exists() and not monitor.exists()
        cmd=[p['executable'],'uniform',str(case['dx_um']),str(case['dt_us']),'10','500','0',str(root),'6','variational',str(case['save_stride']),'200','1'];start=time.monotonic();cap=min(case['timeout_s'],left);timed_out=False
        with log.open('w',encoding='utf-8') as stream,monitor.open('w',encoding='utf-8') as watch:
            proc=subprocess.Popen(cmd,cwd=R,env=ENV,stdout=stream,stderr=subprocess.STDOUT)
            atomic(O/'live.json',dict(state='running',case=name,native_pid=proc.pid,parent_pid=os.getpid(),command=cmd,output_directory=str(root),started_utc=now(),expected_steps=round(10/case['dt_us']),timeout_s=cap))
            watch.write(json.dumps(dict(utc=now(),event='started',native_pid=proc.pid))+'\n');watch.flush();next_check=p['monitor_interval_s']
            while proc.poll() is None:
                elapsed=time.monotonic()-start
                if elapsed>=cap:timed_out=True;proc.kill();proc.wait();break
                if elapsed>=next_check:watch.write(json.dumps(dict(utc=now(),event='half_hour_check',native_pid=proc.pid,handle_alive=proc.poll() is None,elapsed_s=elapsed,**progress(root)))+'\n');watch.flush();next_check+=p['monitor_interval_s']
                time.sleep(1)
            rc=proc.wait();watch.write(json.dumps(dict(utc=now(),event='terminal',native_pid=proc.pid,returncode=rc,timed_out=timed_out,**progress(root)))+'\n')
        rec=dict(name=name,parameters=case,command=cmd,returncode=rc,timed_out=timed_out,wall_s=time.monotonic()-start,terminal_utc=now(),completed=False,native_gate=False,log=str(log.relative_to(R)),log_sha256=sha(log),monitor_log=str(monitor.relative_to(R)),monitor_log_sha256=sha(monitor))
        if (root/'summary.json').exists():
            m=read(root/'summary.json');rec.update(summary=m,preparation=read(root/'preparation.json'),mesh=read(root/'mesh.json'));rec['completed']=not timed_out and rc==0 and m['steps']==round(10/case['dt_us']) and m['failures']==0 and m['N']==case['N']
            if rec['completed']:
                force=rows(root/'mechanics.csv');ratio=max(r['g1_N']/r['target_N'] for r in force);rec['maximum_force_target_ratio']=ratio;rec['native_gate']=len(force)==m['steps'] and ratio<=1.000001 and max(m[k] for k in ['max_V_res','max_S_res','max_regionV_res','max_regionS_res'])<=1e-10 and m['max_budget_excess']<=1e-7 and rec['preparation']['pressure_residual']<=1e-8
        if (root/'failure.json').exists():rec['failure']=read(root/'failure.json')
        rec['files']=[dict(path=str(f.relative_to(R)),bytes=f.stat().st_size,sha256=sha(f)) for f in sorted(root.iterdir()) if f.is_file()];s['cases'].append(rec);ready=rec['native_gate'];save(O/'status.json',s);atomic(O/'live.json',dict(state='terminal_case',case=name,terminal_utc=now()));print(name,rc,'native',rec['native_gate'],round(rec['wall_s'],3),flush=True)
    assert all(sha(R/f)==h for f,h in p['source_hashes'].items()) and sha(Path(p['executable']))==p['binary_sha256'];s.update(state='completed',terminal_utc=now(),compute_wall_s=time.monotonic()-whole);save(O/'status.json',s);atomic(O/'live.json',dict(state='terminal_batch',terminal_utc=now(),native_cases=len(s['cases'])))
if __name__=='__main__':
    a=argparse.ArgumentParser();g=a.add_mutually_exclusive_group(required=True);g.add_argument('--prepare',action='store_true');g.add_argument('--run',action='store_true');args=a.parse_args();prepare() if args.prepare else run()
