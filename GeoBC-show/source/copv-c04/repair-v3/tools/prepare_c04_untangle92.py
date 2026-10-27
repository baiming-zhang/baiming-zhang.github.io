"""Retain091 and freeze a fixed-CAD, volume-only untangling repair attempt."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,ast
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3';C=B/'copv_c04_0000';A=C/'attempt091_curved_lower_tori_refined';P=C/'attempt090_lower_tori_refined'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,v):
 assert not p.exists(),'Preserve '+str(p)
 p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
g=read(A/'geometry_check/review.json');p=read(A/'preflight/review.json');bad=[f for f in g['faces'] if not f['passed']];assert len(bad)==1 and bad[0]['surface']==52 and not p['candidate_FE_validation_allowed'];assert read(B/'sicn_negative_J_validation_v1/review.json')['passed']
now=datetime.now(timezone.utc).isoformat()
put(A/'local_geometry_hold_v1.json',dict(created=now,failed_faces=bad,geometry_review_sha256=sha(A/'geometry_check/review.json'),preflight_sha256=sha(A/'preflight/review.json'),GPU_allowed=False,reason='091 retains retreated CAD nodes at upperTi torus52. Do not fix this incorrect boundary in a shape patch. Next attempt starts freshly from full CAD projection and permits only volume-owned nodes to move.'))
c=read(B/'worker78_config.json');c.update(created=now,worker_state='worker92_progress.json',mesh_attempt=P.name,parent_mesh_sha256=sha(P/'mesh/mesh.npz'),curved_attempt='attempt092_full_CAD_untangle',case_id='copv_c04_0000',purpose='local_stress_accuracy_repair_only',scope='No remeshing: reuse090 original affine parent and freshly recompute full CAD P2 projection. Only complete-star volume-owned nodes may move; every CAD node is fixed at its true full projection. Explicit intermediate negative-J repair permitted only by validated signed objective; full native positiveJ/SICN and actualCAD/volume/Tri6/interface/resource checks still decide eligibility. Original091 retained; never fix its retreated boundary.')
c['environment'].pop('COPV_C04_PROJECTED_CHECKPOINT');c['environment'].update(COPV_C04_UNTANGLE_FULL_CAD='1',COPV_C04_SICN_NATIVE_CHECK_EVERY='10');c['curver_timeout_seconds']=1500
ast.parse((W/'notes/curve_c04_repair_candidate.py').read_text(encoding='utf-8'));put(B/'worker92_config.json',c)
put(B/'same_geometry_mesh_sensitivity_v1/full_CAD_untangle_design_v1.json',dict(created=now,failed_attempt=A.name,hold_sha256=sha(A/'local_geometry_hold_v1.json'),parent_affine_mesh_sha256=c['parent_mesh_sha256'],config_sha256=sha(B/'worker92_config.json'),curver_sha256=sha(W/'notes/curve_c04_repair_candidate.py'),negative_J_derivative_proof_sha256=sha(B/'sicn_negative_J_validation_v1/review.json'),method='Fresh full CAD projection; signed-SICN and signed determinant penalty, volume-owned complete-star coordinate variables only. Boundary coordinates cannot retreat. Invalid intermediate is never a solver mesh.',guards='50k patch cells/20k free nodes,1ring/400iter/1500s; full native checks and original independent screens. No acceptance gate relaxed.',production_eligible=False))
print(json.dumps(dict(config='worker92_config.json',started=False)))
