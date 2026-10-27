"""Retain089 and target only its two remaining CAD-distance failures."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3';A=B/'copv_c04_0000/attempt089_curved_torus56_refined'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,v):
 assert not p.exists(),'Preserve '+str(p)
 p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
g=read(A/'geometry_check/review.json');p=read(A/'preflight/review.json');failed=[f for f in g['faces'] if not f['passed']];fixed=next(f for f in g['faces'] if f['surface']==56)
assert fixed['passed'] and {f['surface'] for f in failed}=={92,96} and all(f['type']=='Torus' and f['part_ids']==[10] for f in failed)
assert not p['candidate_FE_validation_allowed'] and p['resource_guard_passed']
now=datetime.now(timezone.utc).isoformat()
put(A/'local_geometry_hold_v1.json',dict(created=now,failed_faces=failed,repaired_face56=fixed,geometry_review_sha256=sha(A/'geometry_check/review.json'),preflight_sha256=sha(A/'preflight/review.json'),GPU_allowed=False,reason='Upper torus56 CAD screen now passes. Lower titanium torus92/96 still off CAD after guarded curving. Restore those boundary maps before any fixed-boundary interior quality optimization. All original attempts retained.'))
c=read(B/'worker88_config.json');c.update(created=now,worker_state='worker90_progress.json',mesh_attempt='attempt090_lower_tori_refined',curved_attempt='attempt091_curved_lower_tori_refined',scope='Retain passing56 target.15mm; change only remaining lower-Ti torus92/96 surface and shared curves .25→.18mm, bounded by unchanged size/cost/CAD/quality gates. Same0000 local-accuracy research; no stress or production acceptance implied.')
c['environment']['COPV_C04_TORUS_FACE_H_JSON']='{"56":0.15,"92":0.18,"96":0.18}'
put(B/'worker90_config.json',c)
put(B/'same_geometry_mesh_sensitivity_v1/lower_tori_repair_design_v1.json',dict(created=now,failed_attempt='attempt089_curved_torus56_refined',hold_sha256=sha(A/'local_geometry_hold_v1.json'),evidence='56 distance .0249546→.00626266mm passed after local refinement;089 residual failures only92 and96 at .021388/.021467mm against .015mm.',changes_mm={'56':[.15,.15],'92':[.25,.18],'96':[.25,.18]},hypothesis='Finer owning-face and shared-curve boundary resolution may reduce guarded projection retreat. .18mm is a sizing trial, not an acceptance tolerance; independently check actual CAD distance and full native Jacobians/SICN. No guarantee or inferred stress convergence.',config_sha256=sha(B/'worker90_config.json'),production_eligible=False))
print(json.dumps(dict(config='worker90_config.json',started=False)))
