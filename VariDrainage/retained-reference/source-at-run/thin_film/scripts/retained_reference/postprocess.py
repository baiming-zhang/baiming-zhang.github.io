"""Independent immutable terminal subsets; no current native frame reads."""
from pathlib import Path
import argparse,datetime,json,subprocess,sys,time
from run import R,O,RAW,LOG,ENV,read,save,sha
HERE=Path(__file__).resolve().parent;PLAN=O/'postprocess_protocol.json'
def sources():
    return {str(f.relative_to(R)):sha(f) for f in [Path(__file__),HERE/'POSTPROCESS.md',HERE.parent/'global_space/audit_ready.py',HERE.parent/'mechanical_work/analyze.py',HERE.parent/'curved_edge/analyze.py',HERE.parent/'prepared_motion/analyze.py']}
def eligible(name):
    state=read(O/'status.json');records=[r for r in state['cases'] if r['name']==name];assert len(records)==1,'Case not terminal; leave live frames untouched'
    rec=records[0];assert rec['completed'] and rec['native_gate'] and rec['returncode']==0 and not rec['timed_out']
    live=read(O/'live.json');assert not (live.get('state')=='running' and live.get('case')==name)
    native=read(O/'protocol.json');case=next(c for c in native['cases'] if c['name']==name);m=rec['summary'];assert m['N']==case['N'] and m['steps']==m['requested_steps']==round(10/case['dt_us']) and abs(m['T_s']-10e-6)<1e-15
    assert {p.name for p in (RAW/name).iterdir() if p.is_file()}=={Path(f['path']).name for f in rec['files']}
    assert all(sha(R/f['path'])==f['sha256'] for f in rec['files'])
    for key in ['log','monitor_log']:assert sha(R/rec[key])==rec[key+'_sha256']
    return rec
def adapter():
    source=HERE.parent/'global_space/audit_ready.py';code=source.read_text(encoding='utf-8')
    changes={'global_space_audit_inputs_20261113':'retained_reference_audit_inputs_20261113','24previously executed unit checks verified via frozen hashes; not rerun here':'24candidate fixture checks verified via frozen hashes; not rerun by this auditor'}
    for old,new in changes.items():assert code.count(old)==1;code=code.replace(old,new)
    ns=dict(__file__=str(Path(__file__)),__name__='retained_reference_terminal_auditor');exec(compile(code,str(source),'exec'),ns);return ns
def prepare():
    assert not PLAN.exists();p=read(O/'protocol.json');assert all(sha(R/f)==h for f,h in p['source_hashes'].items()) and sha(Path(p['executable']))==p['binary_sha256'];adapter()
    state=read(O/'status.json');done={r['name'] for r in state['cases']};pending=[r['name'] for r in p['cases'] if r['name'] not in done]
    for name in pending:
        try:eligible(name)
        except AssertionError:pass
        else:raise AssertionError('Unrecorded/live case admitted')
    save(PLAN,dict(prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources=sources(),native_protocol_sha256=sha(O/'protocol.json'),per_case_cap_s=180,pending_cases_refused=pending,new_independent_audits_executed=0,scope='Processing-order extension only; original independent formulas and native/comparison gates unchanged',full_C1_complete=False));print('Prepared terminal-only audit; pending cases refused',pending)
def one(name):
    p=read(PLAN);assert p['sources']==sources() and p['native_protocol_sha256']==sha(O/'protocol.json');rec=eligible(name);d=adapter()['one'](rec)
    seal=read(d/'verified.json');assert all(sha(d/n)==h for n,h in seal['outputs'].items())
    fa=read(d/'saved_formula_audit.json')['audits'][name];ma=read(d/'mechanical_audit.json')[name];assert fa['passed'] and ma and all(r['passed'] for r in ma)
    save(d/'wrapper_verified.json',dict(protocol_sha256=sha(PLAN),sources=sources(),inherited_seal_sha256=sha(d/'verified.json'),scalar_saved_steps=len(fa['audited_saved_steps']),mechanical_saved_steps=len(ma),passed=True,scope='Independent saved discrete formulas, not continuous physical flux/traction qualification'))
def ready():
    p=read(PLAN);assert p['sources']==sources();state=read(O/'status.json');results=[];new=0
    for rec in state['cases']:
        if not (rec['completed'] and rec['native_gate']):continue
        name=rec['name'];eligible(name);d=O/'audits/cases'/name;log=LOG/('audit_'+name+'.log')
        if not (d/'wrapper_verified.json').exists():
            assert not d.exists() and not log.exists(),'Preserve unsealed earlier attempt; no retry';cmd=[sys.executable,str(Path(__file__)),'--case',name];start=time.monotonic()
            with log.open('w',encoding='utf-8') as stream:
                try:rc=subprocess.run(cmd,cwd=R,env=ENV,stdout=stream,stderr=subprocess.STDOUT,timeout=p['per_case_cap_s']).returncode
                except subprocess.TimeoutExpired:rc='timeout'
            save(O/('audit_attempt_'+name+'.json'),dict(command=cmd,returncode=rc,wall_s=time.monotonic()-start,log=str(log.relative_to(R)),log_sha256=sha(log)));assert rc==0 and (d/'wrapper_verified.json').exists(),'Independent audit failed; retain evidence';new+=1
        w=read(d/'wrapper_verified.json');seal=read(d/'verified.json');assert w['passed'] and w['sources']==sources() and w['protocol_sha256']==sha(PLAN) and w['inherited_seal_sha256']==sha(d/'verified.json') and all(sha(d/n)==h for n,h in seal['outputs'].items());results.append(dict(name=name,**w))
    report=dict(native_state=state['state'],audited_cases=results,new_independent_case_audits_executed=new,verified_terminal_case_count=len(results),scalar_saved_steps=sum(r['scalar_saved_steps'] for r in results),mechanical_saved_steps=sum(r['mechanical_saved_steps'] for r in results),all4_native_cases_independently_audited=state['state']=='completed' and not state['skipped'] and len(results)==len(state['cases'])==4,full_C1_complete=False,scope='Terminal passing subsets only; no comparison/time/space qualification')
    save(O/'audit_coverage.json',report);print(json.dumps({k:v for k,v in report.items() if k!='audited_cases'},indent=2))
if __name__=='__main__':
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True);g.add_argument('--prepare',action='store_true');g.add_argument('--ready',action='store_true');g.add_argument('--case');a=parser.parse_args();prepare() if a.prepare else one(a.case) if a.case else ready()
