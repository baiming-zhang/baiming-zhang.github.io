"""Separate CPU elastic extension of CAD boundary displacement, never pull CAD back."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import json,os,subprocess,sys,time
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3';C=B/'copv_c04_0000';out=C/'attempt038_elastic36';start=time.perf_counter();state=dict(pid=os.getpid(),started=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),status='active',new_FE_labels=0,new_production_labels=0);assert not out.exists();assert not (B/'worker38_progress.json').exists();env=os.environ.copy();env['PYTHONPATH']=os.pathsep.join([str(W/'publication/GeoBC-Neural-Operator'),str(W/'.rod_deps'),str(W/'.local_deps')])
for key in list(env):
 if key.startswith('COPV_C04_'):del env[key]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(phase,**kw):
 state.update(kw,phase=phase,seconds=time.perf_counter()-start,updated=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat());p=B/'worker38_progress.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(state,indent=2)+'\n');t.replace(p)
def run(phase,args,limit):
 save(phase)
 with (B/(phase+'_execution.log')).open('x',encoding='utf-8') as f:return subprocess.run([sys.executable,'-X','utf8','-B',*args],cwd=W,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=limit).returncode
try:
 rc=run('attempt038_elastic36',['notes/curve_c04_repair_candidate.py','--parent',str(C/'attempt036_graded_root'),'--output',str(out),'--high-order','--high-order-method','HighOrderElastic'],1500)
 if rc or not (out/'result.json').exists() or read(out/'result.json')['status']!='curved_positive_J_mesh_preflight_only':save('curved_failed_retained',status='failed_retained',returncode=rc);raise SystemExit(0)
 rc=run('attempt038_geometry',['notes/check_c04_curved_geometry.py',str(out)],600)
 if rc:save('geometry_failed_retained',status='failed_retained',returncode=rc);raise SystemExit(0)
 rc=run('attempt038_preflight',['notes/preflight_c04_curved_candidate.py',str(out)],600)
 if rc:save('preflight_failed_retained',status='failed_retained',returncode=rc);raise SystemExit(0)
 r=read(out/'preflight/review.json');save('CPU_review_finished',status='CPU_candidate_ready' if r['candidate_FE_validation_allowed'] else 'CPU_candidate_rejected_retained',candidate=str(out.relative_to(W)),candidate_FE_validation_allowed=r['candidate_FE_validation_allowed'])
except SystemExit:raise
except BaseException as e:save(state.get('phase','failed'),status='failed_retained',error=repr(e));raise
