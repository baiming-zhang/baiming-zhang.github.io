"""Sealed independent audits of passing completed controls, no native replay."""
from pathlib import Path
import subprocess,sys,time
from run import R,O,LOG,read,save,sha,ENV
def main():
    state=read(O/'status.json');assert state['state']=='completed'
    assert not (O/'audit_coverage.json').exists()
    p=read(O/'protocol.json');assert all(sha(R/f)==h for f,h in p['sources'].items())
    attempts=[];scalar=0;mechanical=0
    for rec in state['cases']:
        if not rec['passed']:continue
        for f in rec['files']:assert sha(R/f['path'])==f['sha256']
        name=rec['name'];log=LOG/('audit_'+name+'.log');d=O/'audits/cases'/name;assert not log.exists() and not d.exists()
        cmd=[sys.executable,str(Path(__file__).with_name('audit_case.py')),name];start=time.monotonic()
        with log.open('w',encoding='utf-8') as stream:
            try:rc=subprocess.run(cmd,cwd=R,env=ENV,stdout=stream,stderr=subprocess.STDOUT,timeout=180).returncode
            except subprocess.TimeoutExpired:rc='timeout'
        attempt=dict(name=name,command=cmd,returncode=rc,wall_s=time.monotonic()-start,log=str(log.relative_to(R)),log_sha256=sha(log));attempts.append(attempt);save(O/'audit_attempts.json',attempts)
        assert rc==0,'Preserve failed independent audit; no retry'
        wrapper=read(d/'wrapper_verified.json');seal=read(d/'verified.json');assert wrapper['passed'] and wrapper['sealed_outputs_sha256']==sha(d/'verified.json') and all(sha(d/n)==h for n,h in seal['outputs'].items())
        fa=read(d/'saved_formula_audit.json')['audits'][name];ma=read(d/'mechanical_audit.json')[name];assert fa['passed'] and ma and all(r['passed'] for r in ma);scalar+=len(fa['audited_saved_steps']);mechanical+=len(ma)
    passed=len(attempts)==len(state['cases'])==len(p['cases'])==3 and not state['skipped'] and all(c['passed'] for c in state['cases']) and scalar==mechanical==3
    save(O/'audit_coverage.json',dict(completed_control_count=len(state['cases']),audited_controls=len(attempts),scalar_saved_steps=scalar,mechanical_saved_steps=mechanical,actual_single_step_controls_qualified=passed,full_C1_complete=False,scope='Three actual single-step controls and independent discrete audits only; not actual78 replay, time/space or physical flux validation'))
    print((O/'audit_coverage.json').read_text(encoding='utf-8'))
if __name__=='__main__':main()
