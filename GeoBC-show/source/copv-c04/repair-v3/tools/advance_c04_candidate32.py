"""New lower-cost curved-geometry candidate; preserve every earlier attempt."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import json,os,subprocess,sys,time
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3';C=B/'copv_c04_0000';start=time.perf_counter();state=dict(pid=os.getpid(),started=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),status='active',new_FE_labels=0,new_production_labels=0)
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(phase,**kw):
 state.update(kw,phase=phase,seconds=time.perf_counter()-start,updated=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat());p=B/'worker32_progress.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(state,indent=2)+'\n');t.replace(p)
env=os.environ.copy();env['PYTHONPATH']=os.pathsep.join([str(W/'publication/GeoBC-Neural-Operator'),str(W/'.rod_deps'),str(W/'.local_deps')]);env.update(COPV_C04_ROOT_BAND_MM='3',COPV_C04_MATCH_ROOT_CIRCLES='1',COPV_C04_MAX_SURFACE_NODES='700000',COPV_C04_COLLAR_FRONTAL='1',COPV_C04_FRONTAL_PARTS='7,8',COPV_C04_LINEAR_OPTIMIZER='Netgen',COPV_C04_ADJACENT_ROOT_BAND='1',COPV_C04_ADJACENT_ROOT_H='.5',COPV_C04_MAIN_CHORD_MM='.045',COPV_C04_CURVE_DET_ONLY='1')
for key in ['COPV_C04_SURFACE_CHECKPOINT','COPV_C04_PROJECTED_CHECKPOINT']:env.pop(key,None)
def run(phase,args,limit):
 save(phase)
 with (B/(phase+'_execution.log')).open('x',encoding='utf-8') as f:return subprocess.run([sys.executable,'-X','utf8','-B',*args],cwd=W,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=limit).returncode
try:
 rc=run('attempt032_lower_cost_root',['notes/run_c04_repair_candidate.py','--case','copv_c04_0000','--attempt','attempt032_lower_cost_root','--partition','data/copv/c04_pilot_v1/cases/copv_c04_0000/partition_v2','--h','6','--shell-h','5.5','--bond-h','3','--root-h','.25','--max-cells','3300000','--backend','linux'],1150)
 parent=C/'attempt032_lower_cost_root'
 if rc or not (parent/'mesh/mesh.npz').exists():save('mesh_failed_retained',status='failed_retained',reason='Mesh stage did not export a positive-J native candidate. Preserve all artifacts.',returncode=rc);raise SystemExit(1)
 rc=run('attempt033_curved_lower_cost',['notes/curve_c04_repair_candidate.py','--parent',str(parent),'--output',str(C/'attempt033_curved_lower_cost'),'--safe-local'],1500)
 candidate=C/'attempt033_curved_lower_cost'
 if rc or not (candidate/'result.json').exists() or read(candidate/'result.json')['status']!='curved_positive_J_mesh_preflight_only':save('curved_failed_retained',status='failed_retained',returncode=rc);raise SystemExit(1)
 rc=run('attempt033_geometry',['notes/check_c04_curved_geometry.py',str(candidate)],600)
 if rc:save('geometry_check_failed_retained',status='failed_retained',returncode=rc);raise SystemExit(1)
 rc=run('attempt033_preflight',['notes/preflight_c04_curved_candidate.py',str(candidate)],600)
 if rc:save('preflight_failed_retained',status='failed_retained',returncode=rc);raise SystemExit(1)
 r=read(candidate/'preflight/review.json');save('CPU_candidate_review_finished',status='CPU_candidate_ready' if r['candidate_FE_validation_allowed'] else 'CPU_candidate_rejected_retained',candidate=str(candidate.relative_to(W)),candidate_FE_validation_allowed=r['candidate_FE_validation_allowed'],original_volume_SICN_target_passed=r['original_volume_SICN_target_passed'],scope='No automatic GPU selection, no production promotion; inspect actual resource, geometry and quality result.')
except SystemExit:raise
except BaseException as e:save(state.get('phase','worker_failed'),status='failed_retained',error=repr(e));raise
