"""Declare the screened sensitivity mesh without replacing the canonical baseline."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import json
from c04_repair_contracts import build,physical_identity
from c04_curved_point_fields import read,sha
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3';case=B/'copv_c04_0000';out=B/'same_geometry_mesh_sensitivity_v1';now=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat();baseline=case/'attempt058_two_cell_patch56';candidate=case/'attempt085_two_cell_patch84'
def put(p,v):
 assert not p.exists(),'Preserve frozen versions';p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
frozen=B/'screened_mesh_copv_c04_0000_attempt085_two_cell_patch84_v1.json';record=read(frozen);assert record['original_mesh_screen_passed'];contracts=[build(c,8) for c in [baseline,candidate]]
ids=[]
for c in contracts:c=dict(c);c['mesh_sha256']='SAME_CAD_DIFFERENT_SCREENED_MESH';ids.append(physical_identity(c))
assert ids[0]==ids[1]
reviews=[read(c/'mesh/mesh_review_curved_checked.json') for c in [baseline,candidate]];assert all(reviews[0][k]==reviews[1][k] for k in ['source_CAD_sha256','source_model_sha256','parameters_mm'])
put(out/'variant085_planned_q8_contract.json',contracts[1]);put(out/'variant085_declaration_v1.json',dict(schema='c04_mesh_sensitivity_variant_declaration_v1',created=now,baseline=str(baseline.relative_to(W)).replace('\\','/'),variant=record['candidate'],variant_mesh_sha256=record['mesh_sha256'],frozen_mesh_record_sha256=sha(frozen),protocol_sha256=sha(out/'protocol.json'),lineage=['attempt082_same_geometry_root_coarsening','attempt083_curved_root_coarsening','attempt084_sicn_patch83','attempt085_two_cell_patch84'],original_CAD_and_numerical_physics_equal=True,independent_geometry_increment=0,production_eligible=False,actual_field_status='not_started',scope='Additional nonuniform numerical mesh of SAME0000; baseline058 and all historical records retained. Full original mesh screening passed; actual q8 field and common-point stress/mesh sensitivity still pending.'))
queue=B/'gpu_validation_queues/copv_c04_0000_sensitivity_77a904cd977e.json';put(queue,dict(schema='c04_curved_GPU_validation_queue_v1',case_id='copv_c04_0000',candidate=record['candidate'],status='waiting_resources',purpose='same_geometry_mesh_sensitivity_only',assembly_degrees=[8],production_eligible=False,new_FE_labels=0,runs=[],created=now,contract_schema='c04_curved_P2_validation_contract_v2',priority_after='Complete primary0000 pair first, then this q8-only sensitivity variant, then0011/0002/0018; single C04 GPU lane.',study_protocol='data/copv/c04_repair_v3/same_geometry_mesh_sensitivity_v1/protocol.json'))
put(B/'point_evaluator_validation_v1/diagnosis_v1.json',dict(created=now,failure_preserved=True,reason='Initial inverse physical residual tolerance1e-8mm could stop too early for reference containment at thin curved element edge GP locations. Actual14-point output rule contains edge points. No FEA error or field data modification.',repair='Tighten inverse residual tolerance to1e-11mm without relaxing reference containment; immutable validation_v2 independently matched504 actual q6 GP values.',passed_validation='data/copv/c04_repair_v3/point_evaluator_validation_v2/review.json'))
print(json.dumps(dict(declared=record['candidate'],q8_only=True,independent_geometry_increment=0,new_GPU_jobs=0)))
