"""Apply only preregistered geometry gates to the fixed original/candidate matrix."""
import csv,json
from run import R,O,LOCAL,read,save,sha,source_check

def main():
    assert not (O/'summary.json').exists()
    p=read(O/'protocol.json');assert read(O/'status.json')['state']=='completed';assert p['change']==source_check()
    assert all(sha(R/f)==h for f,h in p['original_sources'].items())
    raw=read(O/'raw_manifest.json');assert all(sha(R/f)==v['sha256'] for f,v in raw.items())
    results=[]
    for case in p['cases']:
        d=LOCAL/'native'/case['name'];rows=list(csv.DictReader((d/'observations.csv').open(encoding='utf-8')));assert [r['variant'] for r in rows]==['original','candidate']
        for r in rows:
            for key in r:
                if key not in ['variant','message']:r[key]=float(r[key])
        original,candidate=rows;eq=read(d/'comparison.json')
        accepted=lambda r:r['completed']==1 and r['normal_closure']<=p['normal_closure_gate'] and r['area_closure']<=p['area_closure_gate'] and r['orthogonality_certificate']<=p['orthogonality_certificate_gate'] and r['area_import_error']<=p['area_import_gate']
        expected_rejection=case['name']=='short_face';original_gate=(original['completed']==0 and 'normal_closure=' in original['message']) if expected_rejection else accepted(original)
        equality=True if expected_rejection else eq.get('both_admitted',False) and eq['cells_exact'] and eq['areas_exact']
        if case['name'] in ['unperturbed','resolved_face']:equality=equality and eq['pairs_exact'] and eq['conductances_exact']
        results.append(dict(case=case,observations=rows,equality=eq,expected_original_outcome=bool(original_gate),candidate_geometry_admitted=bool(accepted(candidate)),comparison_contract_satisfied=bool(equality),passed=bool(original_gate and accepted(candidate) and equality)))
    summary=dict(fixed_geometry_matrix_admitted=all(r['passed'] for r in results),cases=results,candidate_admitted_count=sum(r['candidate_geometry_admitted'] for r in results),paired_native_geometry_observations=14,new_native_trajectory_steps=0,short_face_candidate_closure=results[1]['observations'][1]['normal_closure'],actual_snapshots_all_graphs_exact=all(r['equality'].get('conductances_exact',False) and r['equality'].get('pairs_exact',False) for r in results if r['case']['kind']=='saved_actual_geometry'),original_default_changed=False,actual_failed78_reproduced=False,derivative_checks_executed=False,dynamic_convergence_qualified=False,full_C1_complete=False,scope='Geometry-only admission, not operator/dynamic/physical qualification')
    save(O/'summary.json',summary);print(json.dumps({k:v for k,v in summary.items() if k!='cases'},indent=2))
if __name__=='__main__':main()
