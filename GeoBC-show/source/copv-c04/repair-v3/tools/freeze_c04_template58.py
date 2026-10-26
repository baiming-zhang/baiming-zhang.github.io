"""Freeze the first all-screened SAME-0000 mesh; preserve the idle queue history."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import hashlib,json
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,v):
 assert not p.exists(),f'Preserve existing record: {p}'
 p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
stamp=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()
candidate=B/'copv_c04_0000/attempt058_two_cell_patch56';r=read(candidate/'preflight/review.json')
assert r['original_mesh_screen_passed'] and r['candidate_FE_validation_allowed'] and not r['production_eligible']
assert sha(candidate/'mesh/mesh.npz')==r['mesh_sha256']
paths=[candidate/'mesh/mesh.npz',candidate/'preflight/review.json',candidate/'geometry_check/review.json',candidate/'result.json',B/'worker52_config.json',B/'worker56_config.json',B/'worker58_config.json',B/'sicn_objective_validation_v1/review.json',W/'notes/c04_shape_objective.py']
record=dict(schema='c04_mesh_template_v1',created=stamp,case_id=r['case_id'],candidate=str(candidate.relative_to(W)).replace('\\','/'),original_mesh_screen_passed=True,production_eligible=False,recipe=['052 independent linear mesh','053 CAD projection with safe-local pullback','056 fixed-CAD interior analytic SICN patch','058 recentered two-ring polish'],original_gates=read(B/'protocol.json')['original_quality_gate'],evidence_sha256={str(p.relative_to(W)).replace('\\','/'):sha(p) for p in paths},scope='Only this exact SAME-0000 mesh passed all screening gates. GPU-P2 fields and same-geometry sensitivity remain pending. Applying this recipe to another frozen ID requires independent checks.')
write(B/'mesh_template_058_v1.json',record)
q=read(B/'gpu_validation_queue.json');assert q['status']=='waiting_resources' and not q['runs'] and q['new_FE_labels']==0
write(B/'gpu_selection_before058_v1.json',q)
q.setdefault('candidate_history',[]).append(dict(candidate=q['candidate'],selected_until=stamp,runs=0,reason='Superseded before any run by SAME-ID mesh058 passing all original mesh screens. Original candidate retained.'))
q.update(candidate=record['candidate'],updated=stamp,selection_reason='First original-mesh-screen-passing SAME0000 mesh; GPU fields/convergence pending',mesh_template_file='data/copv/c04_repair_v3/mesh_template_058_v1.json')
tmp=B/'gpu_validation_queue.tmp';tmp.write_text(json.dumps(q,indent=2)+'\n');tmp.replace(B/'gpu_validation_queue.json')
cfg=read(B/'worker52_config.json');cfg.update(created=stamp,case_id='copv_c04_0011',worker_state='worker60_progress.json',mesh_attempt='attempt060_template58_linear',curved_attempt='attempt061_template58_curved',scope='Apply the frozen 0000 sizing recipe to original failed ID0011 with its own unchanged CAD/material dimensions. Independent actual checks required; this does not replace its original failure or guarantee acceptance.')
write(B/'worker60_config.json',cfg)
print(json.dumps(dict(template=record['candidate'],next_case=cfg['case_id'])))
