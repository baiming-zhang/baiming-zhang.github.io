"""Audit terminal cases once; seal input/source/output hashes for final reuse."""
from pathlib import Path
import json,hashlib,shutil,importlib.util
from run import R,O,RAW,save,sha
SCOPE='24previously executed unit checks verified via frozen hashes; not rerun here'
def signature(rec):
    sources=[Path(__file__).parents[1]/'mechanical_work/analyze.py',Path(__file__).parents[1]/'curved_edge/analyze.py',Path(__file__).parents[1]/'prepared_motion/analyze.py']
    return dict(record_sha256=hashlib.sha256(json.dumps(rec,sort_keys=True,allow_nan=False).encode()).hexdigest(),raw={f['path']:f['sha256'] for f in rec['files']},auditor_sources={str(p.relative_to(R)):sha(p) for p in sources})
def one(rec):
    assert rec['completed'] and rec['native_gate']
    for f in rec['files']:assert sha(R/f['path'])==f['sha256']
    d=O/'audits/cases'/rec['name'];seal=d/'verified.json';key=signature(rec)
    if seal.exists():
        old=json.loads(seal.read_text());assert old['inputs']==key
        for name,digest in old['outputs'].items():assert sha(d/name)==digest
        print('Reused sealed audit',rec['name'],flush=True);return d
    assert not d.exists(),'Unsealed audit attempt exists; preserve and inspect before retry'
    d.mkdir(parents=True)
    # The inherited auditor enumerates RAW when writing its manifest. Keep
    # this immutable copy separate from the still-running fine native case.
    snapshot=R/'results/local-runs/global_space_audit_inputs_20261113'/rec['name']/'native';snapshot.mkdir(parents=True);shutil.copytree(RAW/rec['name'],snapshot/rec['name'])
    for f in rec['files']:assert sha(snapshot/rec['name']/Path(f['path']).name)==f['sha256']
    save(d/'status.json',dict(state='completed',tests=[],cases=[rec],reused_native_checks=24,audited_subset_only=True,native_checks_scope=SCOPE))
    spec=importlib.util.spec_from_file_location('terminal_case_mechanical',Path(__file__).parents[1]/'mechanical_work/analyze.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.O=d;m.RAW=snapshot;m.old.RAW=snapshot;m.main()
    for fn in ['status.json','saved_formula_audit.json','summary.json']:
        p=d/fn;t=json.loads(p.read_text());t['native_checks_scope']=SCOPE;save(p,t)
    fa=json.loads((d/'saved_formula_audit.json').read_text());ma=json.loads((d/'mechanical_audit.json').read_text());assert fa['audits'][rec['name']]['passed'];assert ma[rec['name']] and all(r['passed'] for r in ma[rec['name']])
    save(seal,dict(inputs=key,outputs={p.name:sha(p) for p in d.iterdir() if p.is_file()},scope='Independent saved-step scalar/mechanical formulas only; immutable copy, no native rerun.'))
    return d
def ready(state):
    cases=[r for r in state['cases'] if r.get('completed') and r.get('native_gate')];paths=[one(r) for r in cases];a=O/'audits';a.mkdir(exist_ok=True);formula={};mechanical={};regions={};diags={}
    for d in paths:
        formula.update(json.loads((d/'saved_formula_audit.json').read_text())['audits']);mechanical.update(json.loads((d/'mechanical_audit.json').read_text()));regions.update(json.loads((d/'mechanical_regions.json').read_text()));diags.update(json.loads((d/'summary.json').read_text())['observer_diagnostics'])
    save(a/'saved_formula_audit.json',dict(attempts=len(cases),completed=len(cases),audits=formula,comparisons=[],native_checks=24,native_checks_scope=SCOPE,scope='Terminal-case independent saved-step audits, source/raw/output hashes sealed per case.'))
    save(a/'mechanical_audit.json',mechanical);save(a/'mechanical_regions.json',regions);save(a/'summary.json',dict(observer_diagnostics=diags,native_checks_scope=SCOPE,full_C1_complete=False));protocol=json.loads((O/'protocol.json').read_text());planned=[r['name'] for r in protocol['guards']+protocol['conditional_matrix']];terminal=[r['name'] for r in state['cases']];skipped=[r['name'] for r in state.get('skipped',[])];pending=[n for n in planned if n not in terminal and n not in skipped];save(a/'coverage.json',dict(native_state_at_audit=state['state'],audited_cases=[r['name'] for r in cases],not_audited=[r['name'] for r in state['cases'] if r not in cases]+pending,planned_but_not_terminal=pending,skipped=skipped,native_matrix_may_still_be_running=state['state']!='completed',saved_scalar_steps=sum(len(r['audited_saved_steps']) for r in formula.values()),saved_mechanical_steps=sum(len(v) for v in mechanical.values()),scope='Only terminal native cases; no qualification of an unfinished time pair.'))
    return a
if __name__=='__main__':ready(json.loads((O/'status.json').read_text()))
