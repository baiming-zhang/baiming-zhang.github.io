"""Preregistered physical-output equivalence, separate from raw byte flags."""
from pathlib import Path
import json,numpy as np
def check(reference,candidate,dx_m):
    results={};good=True
    for filename in ['initial_states.csv','states.csv','initial_faces.csv','faces.csv']:
        a=np.atleast_1d(np.genfromtxt(reference/filename,delimiter=',',names=True));b=np.atleast_1d(np.genfromtxt(candidate/filename,delimiter=',',names=True));entry=dict(same_structure=a.dtype.names==b.dtype.names and len(a)==len(b),errors={})
        if not entry['same_structure']:entry['passed']=False;results[filename]=entry;good=False;continue
        velocity_scale=max(1e-30,float(np.max(np.hypot(a['ux_m_s'],a['uy_m_s'])))) if 'ux_m_s' in a.dtype.names else None
        for key in a.dtype.names:
            if key in ['step','stage','id','i','j','bath']:
                err=0. if np.array_equal(a[key],b[key]) else 1.;tol=0.
            else:
                scale=max(1e-30,float(np.max(np.abs(a[key]))));tol=1e-10
                if key in ['x_m','y_m']:scale=dx_m
                elif key in ['ux_m_s','uy_m_s']:scale=velocity_scale;tol=1e-8
                elif key=='F_m3':tol=1e-8
                err=float(np.max(np.abs(a[key]-b[key])))/scale
            entry['errors'][key]=dict(normalized_max=err,threshold=tol,passed=bool(np.isfinite(err) and err<=tol))
        entry['passed']=all(r['passed'] for r in entry['errors'].values());good &= entry['passed'];results[filename]=entry
    for filename in ['summary.json','mesh.json','preparation.json']:
        a=json.loads((reference/filename).read_text());b=json.loads((candidate/filename).read_text());entry=dict(same_keys=set(a)==set(b),fields={});ignored={'wall_s','max_V_res','max_S_res','max_regionV_res','max_regionS_res','max_budget_excess'} if filename=='summary.json' else {'pressure_residual','iterations'} if filename=='preparation.json' else set()
        for key,value in a.items():
            if key in ignored:entry['fields'][key]=dict(reference=value,candidate=b.get(key),diagnostic_only=True);continue
            if isinstance(value,(str,bool)) or filename=='mesh.json':passed=value==b.get(key);entry['fields'][key]=dict(equal=passed,passed=passed)
            else:
                tol=1e-8 if key=='max_displacement_m' else 1e-10;err=abs(float(b.get(key,float('inf')))-value)/max(1e-30,abs(value));entry['fields'][key]=dict(normalized_error=err,threshold=tol,passed=bool(np.isfinite(err) and err<=tol))
        entry['passed']=entry['same_keys'] and all(v.get('passed',True) for v in entry['fields'].values());good &= entry['passed'];results[filename]=entry
    return dict(passed=bool(good),checks=results,scope='Preregistered numerical output equivalence, clock/residual diagnostics separately checked by native gates; native work/force rows independently audited. Not a time/space/physical-flux qualification.')
