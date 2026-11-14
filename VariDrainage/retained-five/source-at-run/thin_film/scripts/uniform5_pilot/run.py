"""Monitored bounded actual uniform5 pilot. Never restart a native attempt."""
from pathlib import Path
import importlib.util,json,subprocess,datetime,time,os,csv,sys
spec=importlib.util.spec_from_file_location('base',Path(__file__).parents[1]/'moving_voronoi/run.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
R=base.R;O=R/'results/uniform5_pilot_20261113';RAW=R/'results/local-runs/uniform5_pilot_20261113';LOG=R/'results/logs/uniform5_pilot_20261113';save,sha=base.save,base.sha
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def atomic(p,data):
    tmp=p.with_suffix('.tmp');save(tmp,data);tmp.replace(p)
def complete_rows(p):
    if not p.exists():return []
    with p.open(encoding='utf-8') as stream:
        rows=[]
        for r in csv.DictReader(stream):
            try:
                if None in r or any(v is None or v=='' for v in r.values()):continue
                values={k:float(v) for k,v in r.items()};assert all(__import__('math').isfinite(v) for v in values.values());rows.append(values)
            except (ValueError,AssertionError):continue
        return rows
def progress(dest):
    h=complete_rows(dest/'history.csv');m=complete_rows(dest/'mechanics.csv');return dict(last_complete_history_step=int(h[-1]['step']) if h else 0,last_complete_mechanics_step=int(m[-1]['step']) if m else 0,last_force_ratio=(m[-1]['g1_N']/m[-1]['target_N']) if m else None)
def prepare():
    for p in [O,RAW,LOG]:p.mkdir(parents=True,exist_ok=True)
    assert not (O/'protocol.json').exists();control=json.loads((R/'results/local_assembly_20261113/status.json').read_text());audit=json.loads((R/'results/local_assembly_20261113/audit_attempts.json').read_text());cp=json.loads((R/'results/local_assembly_20261113/protocol.json').read_text())
    assert control['state']==audit['state']=='completed' and len(control['operators'])==4 and all(c['passed'] for c in control['operators']) and len(control['cases'])==3 and all(c['passed'] for c in control['cases']) and audit['all_eligible_cases_audited'];assert sum(int(t['counts'][-1][0]) for t in control['tests'])==56 and all(t['passed'] for t in control['tests'])
    exe=base.BIN/'tf_local_assembly_bench.exe';assert sha(exe)==cp['binaries']['tf_local_assembly_bench']
    for name,digest in cp['sources'].items():assert sha(R/name)==digest,name
    for c in control['operators']+control['cases']:
        for f in c['files']:assert sha(R/f['path'])==f['sha256']
    for t in control['tests']:assert sha(R/t['log'])==t['log_sha256']
    refs=json.loads((R/'results/global_space_20261113/status.json').read_text());refs=[c for c in refs['cases'] if c['name'] in ['u10_dt00625','u10_dt003125']];assert len(refs)==2 and all(c['native_gate'] for c in refs)
    for c in refs:
        for f in c['files']:assert sha(R/f['path'])==f['sha256']
    cases=[]
    for name,dt,stride,cap in [('u5_dt00625',.0625,20,18000),('u5_dt003125',.03125,40,28000)]:cases.append(dict(name=name,mode='uniform',dx_um=5,dt_us=dt,T_us=10,radius_um=500,amplitude=0,c_bulk_mM=6,preparation='variational',save_stride=stride,refine_start_um=200,observer=True,timeout_s=cap))
    sources=[Path(__file__),Path(__file__).with_name('README.md'),Path(__file__).with_name('monitor.py')]
    p=dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),source_hashes={**cp['sources'],**{str(p.relative_to(R)):sha(p) for p in sources}},binary_sha256=sha(exe),executable=str(exe),qualified_control_protocol_sha256=sha(R/'results/local_assembly_20261113/protocol.json'),qualified_control_status_sha256=sha(R/'results/local_assembly_20261113/status.json'),qualified_audit_sha256=sha(R/'results/local_assembly_20261113/audit_attempts.json'),reused_candidate_unit_checks=56,new_unit_executions=0,reference_cases=refs,cases=cases,whole_native_cap_s=46000,monitor_interval_s=1800,gates=dict(balance=1e-10,energy=1e-7,force_ratio=1.000001,preparation=1e-8,time=.02,space=.1),environment={'VF_THREADS':'2','caller_VF':'excluded'},scope='Actual uniform5 at10us with qualified mathematically equivalent local assembly; original50/200us/full C1 remain required')
    save(O/'protocol.json',p);save(O/'status.json',dict(state='prepared',cases=[],skipped=[],full_C1_complete=False));print('Prepared monitored uniform5 pilot, no native process launched.')
def run():
    p=json.loads((O/'protocol.json').read_text());state=json.loads((O/'status.json').read_text());assert state['state']=='prepared'
    for name,digest in p['source_hashes'].items():assert sha(R/name)==digest,name
    assert sha(Path(p['executable']))==p['binary_sha256'];authorization=json.loads((O/'monitor_authorization.json').read_text());assert authorization['status']=='ACTIVE' and authorization['interval_minutes']==30 and authorization['native_tool_verified'];age=(datetime.datetime.now(datetime.timezone.utc)-datetime.datetime.fromisoformat(authorization['verified_utc'])).total_seconds();assert 0<=age<=900
    state.update(state='running_matrix',started_utc=now(),monitor_id=authorization['id']);atomic(O/'status.json',state);used=0
    for index,c in enumerate(p['cases']):
        if used>=p['whole_native_cap_s']:state['skipped'].extend(dict(name=g['name'],reason='whole native budget') for g in p['cases'][index:]);break
        dest=RAW/c['name'];assert not dest.exists();dest.mkdir();log=LOG/(c['name']+'.log');assert not log.exists();cmd=[p['executable'],c['mode']]+[str(c[k]) for k in ['dx_um','dt_us','T_us','radius_um','amplitude']]+[str(dest),str(c['c_bulk_mM']),c['preparation'],str(c['save_stride']),str(c['refine_start_um']),'1'];start=time.monotonic();start_utc=now();rc=None;checks=LOG/('monitor_'+c['name']+'.jsonl');assert not checks.exists()
        with log.open('w',encoding='utf-8') as stream,checks.open('w',encoding='utf-8') as monitoring:
            child=subprocess.Popen(cmd,cwd=R,env=base.ENV,stdout=stream,stderr=subprocess.STDOUT);live=dict(state='running',case=c['name'],native_pid=child.pid,parent_pid=os.getpid(),command=cmd,started_utc=start_utc,expected_steps=round(c['T_us']/c['dt_us']),output_directory=str(dest),monitor_id=authorization['id']);atomic(O/'live.json',live);monitoring.write(json.dumps(dict(utc=now(),event='native_started',native_pid=child.pid,handle_alive=child.poll() is None,**progress(dest)))+'\n');monitoring.flush();next_check=start+1800;limit=min(c['timeout_s'],p['whole_native_cap_s']-used)
            while True:
                left=limit-(time.monotonic()-start)
                if left<=0:child.kill();child.wait();rc='timeout';break
                try:rc=child.wait(timeout=min(60,left));break
                except subprocess.TimeoutExpired:
                    if time.monotonic()>=next_check:
                        monitoring.write(json.dumps(dict(utc=now(),event='half_hour_check',native_pid=child.pid,handle_alive=child.poll() is None,elapsed_s=time.monotonic()-start,**progress(dest)))+'\n');monitoring.flush();next_check+=1800
            monitoring.write(json.dumps(dict(utc=now(),event='native_terminal',native_pid=child.pid,returncode=rc,elapsed_s=time.monotonic()-start,**progress(dest)))+'\n');monitoring.flush()
        elapsed=time.monotonic()-start;used+=elapsed;rec=dict(name=c['name'],parameters=c,command=cmd,returncode=rc,native_pid=child.pid,started_utc=start_utc,terminal_utc=now(),wall_s=elapsed,log=str(log.relative_to(R)),log_sha256=sha(log),monitor_log=str(checks.relative_to(R)),monitor_log_sha256=sha(checks),completed=False,native_gate=False)
        for filename,key in [('summary.json','summary'),('preparation.json','preparation'),('mesh.json','mesh'),('failure.json','failure')]:
            if (dest/filename).exists():rec[key]=json.loads((dest/filename).read_text())
        m=rec.get('summary',{});rec['completed']=rc==0 and m.get('steps')==round(c['T_us']/c['dt_us']) and m.get('failures')==0 and m.get('N')==2112
        if rec['completed']:
            stats=complete_rows(dest/'mechanics.csv');ratio=max(r['g1_N']/r['target_N'] for r in stats);rec['maximum_force_target_ratio']=ratio;rec['native_gate']=len(stats)==m['steps'] and ratio<=1.000001 and max(m[k] for k in ['max_V_res','max_S_res','max_regionV_res','max_regionS_res'])<=1e-10 and m['max_budget_excess']<=1e-7 and rec['preparation']['pressure_residual']<=1e-8
        rec['files']=[dict(path=str(f.relative_to(R)),bytes=f.stat().st_size,sha256=sha(f)) for f in sorted(dest.iterdir()) if f.is_file()];state['cases'].append(rec);state['compute_wall_s']=used;atomic(O/'status.json',state);atomic(O/'live.json',dict(**{k:v for k,v in live.items() if k!='state'},state='terminal_case',returncode=rc,terminal_utc=now()));print(c['name'],rc,rec['native_gate'],round(elapsed,2),flush=True)
        if not rec['native_gate']:state['skipped'].extend(dict(name=g['name'],reason='prior native case failed; no retry') for g in p['cases'][index+1:]);break
    state.update(state='completed',terminal_utc=now(),compute_wall_s=used);atomic(O/'status.json',state);atomic(O/'live.json',dict(state='terminal_batch',terminal_utc=now(),native_cases=len(state['cases']),monitor_id=authorization['id']));print('Native batch terminal; independent audits/time-space comparisons still required.',flush=True)
if __name__=='__main__':prepare() if '--prepare' in sys.argv else run()
