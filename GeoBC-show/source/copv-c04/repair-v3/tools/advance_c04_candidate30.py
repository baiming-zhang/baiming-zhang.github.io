"""CPU-only optimizer diagnosis beyond the unchanged GPU candidate size ceiling."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import json,os,subprocess,sys,time
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3';C=B/'copv_c04_0000';start=time.perf_counter();state=dict(pid=os.getpid(),started=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),status='active',phase='attempt030_CPU_Netgen',new_FE_labels=0,new_production_labels=0)
def save(**kw):
 state.update(kw,seconds=time.perf_counter()-start,updated=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat());p=B/'worker30_progress.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(state,indent=2)+'\n');t.replace(p)
env=os.environ.copy();env['PYTHONPATH']=os.pathsep.join([str(W/'publication/GeoBC-Neural-Operator'),str(W/'.rod_deps'),str(W/'.local_deps')]);env.update(COPV_C04_ROOT_BAND_MM='3',COPV_C04_MATCH_ROOT_CIRCLES='1',COPV_C04_MAX_SURFACE_NODES='700000',COPV_C04_COLLAR_FRONTAL='1',COPV_C04_LINEAR_OPTIMIZER='Netgen');env.pop('COPV_C04_SURFACE_CHECKPOINT',None)
with (B/'resource_amendment_v4_CPU_only.json').open('x',encoding='utf-8') as f:json.dump(dict(created=state['started'],max_CPU_optimizer_cells=3300000,max_GPU_candidate_cells=2200000,max_GPU_P2_DOFs=10000000,reason='029 matching root circles reduced low-SICN count to 1102 but generated 3062535 positive-J cells, above GPU candidate ceiling. Preserve 029. CPU-only Netgen diagnosis on the exact retained linear mesh in NEW 030, to distinguish topology/size effects before designing a cheaper mesh. No GPU eligibility or original quality threshold changes.'),f,indent=2)
save()
try:
 with (B/'attempt030_execution.log').open('x',encoding='utf-8') as log:
  run=subprocess.run([sys.executable,'-X','utf8','-B','notes/run_c04_repair_candidate.py','--case','copv_c04_0000','--attempt','attempt030_CPU_Netgen_root','--partition','data/copv/c04_pilot_v1/cases/copv_c04_0000/partition_v2','--shell-h','4.5','--bond-h','3','--root-h','.25','--max-cells','3300000','--backend','linux','--linear-checkpoint',str(C/'attempt029_matched_root_volume/mesh/failed_mesh.msh')],cwd=W,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1100)
 result=json.loads((C/'attempt030_CPU_Netgen_root/state.json').read_text());save(status='CPU_stage_finished',phase=result['phase'],returncode=run.returncode,reason='CPU-only diagnosis; GPU size ceiling remains unchanged. Inspect quality and material volumes before designing next candidate.')
except BaseException as e:save(status='failed_retained',error=repr(e));raise
