"""Execute one frozen CPU config with exclusive attempts and bounded children."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse,json,os,subprocess,sys,time
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3';C=B/'copv_c04_0000'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True);args=ap.parse_args();cfg=read(args.config);start=time.perf_counter();state=dict(pid=os.getpid(),started=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),status='active',config=str(args.config),new_FE_labels=0,new_production_labels=0)
 progress=B/cfg['worker_state'];assert not progress.exists(),'Do not repeat an existing worker'
 for key in ['mesh_attempt','curved_attempt']:assert not (C/cfg[key]).exists(),'Never overwrite an attempt'
 env=os.environ.copy();env['PYTHONPATH']=os.pathsep.join([str(W/'publication/GeoBC-Neural-Operator'),str(W/'.rod_deps'),str(W/'.local_deps')])
 for key in list(env):
  if key.startswith('COPV_C04_'):del env[key]
 env.update(cfg['environment']);env['COPV_C04_CURVE_DET_ONLY']='1'
 def save(phase,**kw):
  state.update(kw,phase=phase,seconds=time.perf_counter()-start,updated=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat());temp=progress.with_suffix('.tmp');temp.write_text(json.dumps(state,indent=2)+'\n');temp.replace(progress)
 def run(phase,arguments,limit):
  save(phase)
  with (B/(phase+'_execution.log')).open('x',encoding='utf-8') as f:return subprocess.run([sys.executable,'-X','utf8','-B',*arguments],cwd=W,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=limit).returncode
 try:
  parameters=cfg['parameters'];command=['notes/run_c04_repair_candidate.py','--case','copv_c04_0000','--attempt',cfg['mesh_attempt'],'--partition','data/copv/c04_pilot_v1/cases/copv_c04_0000/partition_v2','--backend','linux']
  for key,value in parameters.items():command.extend(['--'+key,str(value)])
  if cfg.get('linear_checkpoint'):command.extend(['--linear-checkpoint',str(W/cfg['linear_checkpoint'])])
  rc=run(cfg['mesh_attempt'],command,1150);parent=C/cfg['mesh_attempt'];candidate=C/cfg['curved_attempt']
  if rc or not (parent/'mesh/mesh.npz').exists():save('mesh_failed_retained',status='failed_retained',returncode=rc);return
  rc=run(cfg['curved_attempt'],['notes/curve_c04_repair_candidate.py','--parent',str(parent),'--output',str(candidate),'--safe-local'],1500)
  if rc or not (candidate/'result.json').exists() or read(candidate/'result.json')['status']!='curved_positive_J_mesh_preflight_only':save('curved_failed_retained',status='failed_retained',returncode=rc);return
  rc=run(cfg['curved_attempt']+'_geometry',['notes/check_c04_curved_geometry.py',str(candidate)],600)
  if rc:save('geometry_failed_retained',status='failed_retained',returncode=rc);return
  rc=run(cfg['curved_attempt']+'_preflight',['notes/preflight_c04_curved_candidate.py',str(candidate)],600)
  if rc:save('preflight_failed_retained',status='failed_retained',returncode=rc);return
  r=read(candidate/'preflight/review.json');save('CPU_review_finished',status='CPU_candidate_ready' if r['candidate_FE_validation_allowed'] else 'CPU_candidate_rejected_retained',candidate=str(candidate.relative_to(W)),candidate_FE_validation_allowed=r['candidate_FE_validation_allowed'],original_volume_SICN_target_passed=r['original_volume_SICN_target_passed'])
 except BaseException as e:save(state.get('phase','failed'),status='failed_retained',error=repr(e));raise
if __name__=='__main__':main()
