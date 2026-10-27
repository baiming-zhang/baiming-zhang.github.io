"""Declare an additional screened same-CAD q8 accuracy study, without launching it."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib
from c04_repair_contracts import build,physical_identity
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3';O=B/'same_geometry_mesh_sensitivity_v1';candidate=B/'copv_c04_0000/attempt092_full_CAD_untangle';baseline=B/'copv_c04_0000/attempt058_two_cell_patch56'
def put(p,v):
 assert not p.exists(),'Preserve '+str(p)
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
record=read(B/'screened_mesh_copv_c04_0000_attempt092_full_CAD_untangle_v1.json');r=read(candidate/'preflight/review.json');assert r['original_mesh_screen_passed'] and r['candidate_FE_validation_allowed'] and record['mesh_sha256']==sha(candidate/'mesh/mesh.npz')
assert not list(candidate.glob('local_geometry_hold*.json'))
contracts=[build(c,8) for c in [baseline,candidate]];normalized=[]
for c in contracts:
 c=dict(c);c['mesh_sha256']='SAME_CAD_DIFFERENT_SCREENED_MESH';normalized.append(physical_identity(c))
assert normalized[0]==normalized[1]
reviews=[read(c/'mesh/mesh_review_curved_checked.json') for c in [baseline,candidate]];assert all(reviews[0][k]==reviews[1][k] for k in ['source_CAD_sha256','source_model_sha256','parameters_mm'])
now=datetime.now(timezone.utc).isoformat();protocol=O/'accuracy092_protocol_v1.json'
put(protocol,dict(created=now,baseline='attempt058_two_cell_patch56',additional_reference='attempt085_two_cell_patch84',variant=candidate.name,purpose='Local stress accuracy research after observed235-probe sensitivity; not an independent geometry or replacement of old fields.',lineage=['attempt086_local_accuracy','attempt087_curved_local_accuracy','attempt088_torus56_refined','attempt089_curved_torus56_refined','attempt090_lower_tori_refined','attempt091_curved_lower_tori_refined','attempt092_full_CAD_untangle'],changes='Root min .25→.18mm/plateau.25/growth.6,local Ti9/10 hotspot bands.65mm,torus56.15/92and96.18mm,then fixed-full-CAD volume-only untangling. Exact same CAD/material/load/BC convention.',assembly_degrees=[8],solver='FrozenV4 FP64 with original resource/actualCSR/CUDA guards',probe_plan='Reuse exact baseline058235 probe IDs/physical locations and material-specific trueP2 inverse/gradients; compare092 against058 and retained085 values at those same IDs. Preserve every failed/ambiguous probe ID; no substitution or zero-fill. Report raw and6DOF rigid-aligned displacement,unweighted tensor point-sampling L2(notvolume/trutherror),energy and nativeGP extrema separately.',acceptance='Original mesh screens must all pass before GPU. Numerical residual/work/native ID checks after GPU; scientific local stress sensitivity remains an explicit separate hold. Three nonuniform meshes do not by themselves establish h-limit or peak convergence.',production_eligible=False))
put(O/'accuracy092_planned_q8_contract.json',contracts[1]);put(O/'accuracy092_declaration_v1.json',dict(created=now,candidate=record['candidate'],mesh_sha256=record['mesh_sha256'],frozen_record_sha256=sha(B/'screened_mesh_copv_c04_0000_attempt092_full_CAD_untangle_v1.json'),protocol_sha256=sha(protocol),independent_geometry_increment=0,production_eligible=False,original_CAD_and_numerical_physics_equal=True,actual_field_status='not_started'))
queue=B/('gpu_validation_queues/copv_c04_0000_accuracy_'+record['mesh_sha256'][:12]+'.json')
put(queue,dict(schema='c04_curved_GPU_validation_queue_v1',case_id='copv_c04_0000',candidate=record['candidate'],status='waiting_resources',purpose='same_geometry_mesh_sensitivity_only',assembly_degrees=[8],production_eligible=False,new_FE_labels=0,runs=[],created=now,contract_schema='c04_curved_P2_validation_contract_v2',priority_after='Remaining own-case0011,0002,0018 paired validations; retain one C04 GPU lane.',study_protocol=str(protocol.relative_to(W)).replace('\\','/')))
print(json.dumps(dict(queue=str(queue.relative_to(W)),new_GPU_jobs=0,independent_geometry_increment=0)))
