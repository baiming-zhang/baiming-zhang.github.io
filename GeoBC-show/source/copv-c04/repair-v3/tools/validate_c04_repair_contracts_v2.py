"""Check per-ID contracts and record a transparent legacy descriptive erratum."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import ast,hashlib,json,shutil
from c04_repair_contracts import build,physical_identity
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,v):
 assert not p.exists(),'Never overwrite validation/contract evidence';p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')
q=read(B/'gpu_validation_queue.json');assert q['status']!='running';out=B/'contract_metadata_validation_v2';assert not out.exists();out.mkdir();stamp=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat();candidate=W/q['candidate'];oldrun=candidate/'GPU_q6_attempt001';old=read(oldrun/'executed_contract.json');summary=read(oldrun/'summary.json');protocol=read(oldrun/'quadrature_protocol.json');command=(oldrun/'run.sh').read_text();assert summary['status']=='completed' and summary['displacement_degree']==2 and 'curved quadratic' in summary['geometry_mapping'] and protocol['assembly_degree']==6 and protocol['output_points']==14 and '--degree 2' in command and '--assembly-quadrature 6' in command;assert sha(B/'gpu_solver_sources_v3/solver.py')==sha(W/'publication/GeoBC-Neural-Operator/tools/copv_fenicsx_c04_repair.py')==summary['solver_sha256'];assert sha(oldrun/'executed_contract.json')==summary['contract_sha256'];rows=[]
for item in read(B/'screened_mesh_registry.json')['meshes']:
 c=W/item['candidate'];contracts=[build(c,d) for d in [6,8]];assert all(v['case_id']==item['case_id'] and v['parent_case_id']==item['case_id'] and v['mesh_sha256']==item['mesh_sha256'] for v in contracts);assert physical_identity(contracts[0])==physical_identity(contracts[1])
 for d,v in zip([6,8],contracts):put(out/(item['case_id']+'_planned_q'+str(d)+'.json'),v)
 rows.append(dict(case_id=item['case_id'],candidate=item['candidate'],physical_input_identity_sha256=physical_identity(contracts[0])))
new=build(candidate,8);assert physical_identity(old)==physical_identity(new),'Numerical problem differs from completed q6'
errata=dict(schema='c04_legacy_contract_metadata_errata_v1',created=stamp,legacy_contract_sha256=summary['contract_sha256'],affected_runs=[str((candidate/n).relative_to(W)).replace('\\','/') for n in ['GPU_q6_attempt001','GPU_q8_attempt001']],correction='Legacy C02 discretization/geometry/equivalent-cap descriptive text remained in copied metadata. Actual run CLI, frozen executed solver, summary and raw native schema prove isoparametric P2 with actual quadratic cap-edge integration; q6 used assembly degree6 and output14GP. Original executed contracts and outputs remain immutable.',actual_completed_run=dict(displacement_degree=2,assembly_degree=6,output_degree=4,output_points=14,geometry_mapping=summary['geometry_mapping'],cap_convention=summary['cap_convention'],solver_sha256=summary['solver_sha256'],run_script_sha256=sha(oldrun/'run.sh')),physical_input_identity_sha256=physical_identity(old),next_contract_schema='c04_curved_P2_validation_contract_v2',numerical_solver_source_unchanged=True,old_outputs_rewritten=False,production_eligible=False)
put(B/'legacy_contract_metadata_errata_v1.json',errata);shutil.copyfile(W/'notes/c04_repair_contracts.py',out/'builder_validated.py');put(out/'review.json',dict(passed=True,scope='Metadata and exact physical-input identity regression only; not a new FEA run or convergence result.',created=stamp,legacy_q6_to_v2_q8_same_physics=True,legacy_q6_contract_immutable_sha256=summary['contract_sha256'],solver_source_unchanged_sha256=summary['solver_sha256'],per_case_contracts=rows,builder_sha256=sha(W/'notes/c04_repair_contracts.py')))
folder=B/'gpu_validation_queues'
for item in read(B/'screened_mesh_registry.json')['meshes']:
 if item['case_id']=='copv_c04_0000':continue
 path=folder/(item['case_id']+'_'+item['mesh_sha256'][:12]+'.json');put(path,dict(schema='c04_curved_GPU_validation_queue_v1',case_id=item['case_id'],candidate=item['candidate'],status='waiting_resources',production_eligible=False,new_FE_labels=0,runs=[],created=stamp,contract_schema='c04_curved_P2_validation_contract_v2',priority_after='Primary0000 q6/q8 validation must complete first; CLI enforces this and a single C04 GPU lane'))
print(json.dumps(dict(passed=True,cases=[r['case_id'] for r in rows],old_q6_physics_unchanged=True,legacy_contract_preserved=True,secondary_GPU_jobs_started=0)))
