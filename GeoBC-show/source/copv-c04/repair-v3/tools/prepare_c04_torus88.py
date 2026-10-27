"""Retain087 CAD failure and prepare a bounded same-CAD torus refinement trial."""
from pathlib import Path
from datetime import datetime,timezone
import ast,hashlib,json
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3';C=B/'copv_c04_0000';A=C/'attempt087_curved_local_accuracy'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,v):
 assert not p.exists(),str(p)+' already exists; retain it'
 p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
g=read(A/'geometry_check/review.json');r=read(A/'preflight/review.json');bad=[f for f in g['faces'] if not f['passed']]
assert len(bad)==1 and bad[0]['surface']==56 and bad[0]['type']=='Torus' and bad[0]['part_ids']==[9]
assert not r['candidate_FE_validation_allowed'] and r['resource_guard_passed'] and r['material_volume_passed']
now=datetime.now(timezone.utc).isoformat()
put(A/'local_geometry_hold_v1.json',dict(created=now,case_id='copv_c04_0000',failed_faces=bad,geometry_review_sha256=sha(A/'geometry_check/review.json'),preflight_sha256=sha(A/'preflight/review.json'),GPU_allowed=False,reason='Guarded curving retained off-CAD nodes on upper titanium torus56. Do not fix these erroneous boundary coordinates in an interior SICN optimization. Preserve087 and restore geometry in a new attempt.'))
cfg=read(B/'worker86_config.json');cfg.update(created=now,worker_state='worker88_progress.json',mesh_attempt='attempt088_torus56_refined',curved_attempt='attempt089_curved_torus56_refined',scope='Same0000 original CAD/material/load and086 local-accuracy recipe; only failed torus56 surface and its shared boundary curves .25→.15mm. Other torus faces retain .25mm. This tests whether finer boundary triangles permit valid CAD projection without retreat; no guarantee. Original quality/CAD/volume/load/resource gates unchanged; no fixed-wrong-boundary optimization.')
cfg['environment']['COPV_C04_TORUS_FACE_H_JSON']='{"56":0.15}'
for source in ['publication/GeoBC-Neural-Operator/tools/mesh_copv_c04_guarded_repair.py','notes/run_c04_repair_candidate.py']:
 ast.parse((W/source).read_text(encoding='utf-8'))
put(B/'worker88_config.json',cfg)
put(B/'same_geometry_mesh_sensitivity_v1/torus56_repair_design_v1.json',dict(created=now,baseline_failed_attempt='attempt087_curved_local_accuracy',geometry_hold_sha256=sha(A/'local_geometry_hold_v1.json'),original_CAD_sha256=g['CAD_sha256'],changes={'torus56_surface_and_boundary_target_mm':[.25,.15]},unchanged='All086 local accuracy bands/root sizing, CAD/materials/physics/quality and resource criteria',diagnosis='Only CAD face56 failed its distance screen;302 volume elements remain below SICN. Geometry failure precedes any fixed-boundary SICN repair. Boundary refinement is an unproven mesh repair hypothesis, not a stress-accuracy result.',config_sha256=sha(B/'worker88_config.json'),source_sha256={p:sha(W/p) for p in ['publication/GeoBC-Neural-Operator/tools/mesh_copv_c04_guarded_repair.py','notes/run_c04_repair_candidate.py']},production_eligible=False))
print(json.dumps(dict(config='worker88_config.json',original_attempt_retained=True,started=False)))
