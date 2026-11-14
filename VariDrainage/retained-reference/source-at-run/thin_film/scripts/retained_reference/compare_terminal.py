"""Original moving-overlap metrics, all four same-candidate native cases only."""
from pathlib import Path
import argparse,datetime,importlib.util,json
import numpy as np
import shapely
from run import R,O,read,save,sha
import postprocess as audits
HERE=Path(__file__).resolve().parent;PLAN=O/'comparison_protocol.json'
def sources():return {str(f.relative_to(R)):sha(f) for f in [Path(__file__),HERE/'COMPARE.md',HERE/'postprocess.py',HERE.parent/'prepared_motion/analyze.py',HERE.parent/'curved_edge/analyze.py']}
def load():
    s=importlib.util.spec_from_file_location('retained_original_overlap',HERE.parent/'prepared_motion/analyze.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def prepare():
    assert not PLAN.exists() and not (O/'comparisons.json').exists()
    old=read(R/'results/uniform5_pilot_20261113/comparison_protocol.json')
    for f in [HERE.parent/'prepared_motion/analyze.py',HERE.parent/'curved_edge/analyze.py']:assert old['sources'][str(f.relative_to(R))]==sha(f)
    ap=read(audits.PLAN);assert ap['sources']==audits.sources();load()
    save(PLAN,dict(prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),native_state_at_prepare=read(O/'status.json')['state'],sources=sources(),native_protocol_sha256=sha(O/'protocol.json'),postprocess_protocol_sha256=sha(audits.PLAN),original_comparison_protocol_sha256=sha(R/'results/uniform5_pilot_20261113/comparison_protocol.json'),original_baseline_reproduction=old['baseline_reproduction'],numerical_environment=dict(numpy=np.__version__,shapely=shapely.__version__),new_native_executions=0,full_C1_complete=False))
    print('Prepared original metric comparison; no new reference comparisons computed.')
def eligible():
    s=read(O/'status.json');assert s['state']=='completed','Native matrix still live; no frame reads/output'
    assert not s['skipped'] and len(s['cases'])==4
    p=read(PLAN);assert p['sources']==sources() and p['native_protocol_sha256']==sha(O/'protocol.json') and p['postprocess_protocol_sha256']==sha(audits.PLAN)
    coverage=read(O/'audit_coverage.json');assert coverage['all4_native_cases_independently_audited'] and coverage['scalar_saved_steps']==coverage['mechanical_saved_steps']==36
    seals={}
    for rec in s['cases']:
        assert audits.eligible(rec['name'])==rec;d=O/'audits/cases'/rec['name'];w=read(d/'wrapper_verified.json');a=read(d/'verified.json');assert w['passed'] and w['sources']==audits.sources() and w['protocol_sha256']==sha(audits.PLAN) and w['inherited_seal_sha256']==sha(d/'verified.json') and all(sha(d/f)==h for f,h in a['outputs'].items());seals[rec['name']]=sha(d/'verified.json')
    return s,p,seals
def run():
    s,p,seals=eligible();assert not (O/'comparisons.json').exists(),'Preserve existing comparison'
    prep=load();cache={}
    for rec in s['cases']:
        root=next(R/f['path'] for f in rec['files'] if Path(f['path']).name=='states.csv').parent;m=rec['summary'];initial=prep.old.read(root/'initial_states.csv');states=prep.old.read(root/'states.csv');end=states[(states['step']==m['steps'])&(states['stage']==2)];assert len(end)==len(initial)==m['N'];cache[rec['name']]=(initial,states,m)
    rows=[]
    for a,b,kind in [('u20_dt00625','u20_dt003125','time'),('u10_dt00625','u10_dt003125','time'),('u20_dt00625','u10_dt00625','space'),('u20_dt003125','u10_dt003125','space')]:
        pair=prep.compare(a,b,kind,cache);assert len(pair)==1 and abs(pair[0]['t_s']-10e-6)<1e-15;rows.extend(pair)
    time_pass=all(r['passed'] for r in rows if r['kind']=='time');space_pass=all(r['passed'] for r in rows if r['kind']=='space')
    save(O/'status_at_native_completion.json',s);save(O/'comparison_inputs.json',dict(protocol_sha256=sha(PLAN),sources=sources(),raw={f['path']:f['sha256'] for rec in s['cases'] for f in rec['files']},audit_seals=seals,numerical_environment=p['numerical_environment']))
    save(O/'comparisons.json',rows);save(O/'convergence_summary.json',dict(first_checkpoint_time_pairs_passed=time_pass,first_checkpoint_space_rows_passed=space_pass,first_checkpoint_space_qualified=time_pass and space_pass,missing_original_checkpoints_us=[50,200],new_candidate5um_pair_completed=False,physical_normal_flux_traction_validated=False,full_C1_complete=False,scope='Same-candidate20/10um at10us only; original overlap norms and failed rows retained'))
    print(json.dumps(rows,indent=2))
if __name__=='__main__':
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True);g.add_argument('--prepare',action='store_true');g.add_argument('--run',action='store_true');g.add_argument('--readiness',action='store_true');a=parser.parse_args()
    if a.prepare:prepare()
    elif a.run:run()
    else:
        try:eligible()
        except (AssertionError,FileNotFoundError) as e:print(json.dumps(dict(ready=False,reason=str(e),comparisons_exist=(O/'comparisons.json').exists(),new_native_executions=0)))
        else:print(json.dumps(dict(ready=True,comparisons_exist=(O/'comparisons.json').exists(),new_native_executions=0)))
