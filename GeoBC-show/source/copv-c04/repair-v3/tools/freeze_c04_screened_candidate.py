"""Freeze a screened original-ID mesh and register it without changing active runs."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse,hashlib,json
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
 return h.hexdigest()
def write(p,v):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');t.replace(p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('candidate',type=Path);ap.add_argument('--config',action='append',type=Path,required=True);args=ap.parse_args();c=args.candidate.resolve();r=read(c/'preflight/review.json');case=r['case_id']
 assert case in {i['id'] for i in read(B/'protocol.json')['inputs']}
 assert c.parent==B/case and r['original_mesh_screen_passed'] and r['candidate_FE_validation_allowed'] and not r['production_eligible']
 assert not (c/'local_geometry_hold.json').exists() and sha(c/'mesh/mesh.npz')==r['mesh_sha256']
 configs=[p.resolve() for p in args.config]
 assert all(p.parent==B and read(p).get('case_id','copv_c04_0000')==case for p in configs)
 files=[c/'mesh/mesh.npz',c/'preflight/review.json',c/'geometry_check/review.json',c/'result.json',*configs,W/'notes/c04_shape_objective.py',B/'sicn_objective_validation_v1/review.json']
 name='screened_mesh_'+case+'_'+c.name+'_v1.json';target=B/name;assert not target.exists(),'Never overwrite a frozen mesh record'
 record=dict(schema='c04_screened_mesh_v1',created=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),case_id=case,candidate=str(c.relative_to(W)).replace('\\','/'),mesh_sha256=r['mesh_sha256'],original_mesh_screen_passed=True,production_eligible=False,GPU_status='not_started',evidence_sha256={str(p.relative_to(W)).replace('\\','/'):sha(p) for p in files},original_quality_targets=read(B/'protocol.json')['original_quality_gate'],scope='Exact original-ID mesh screening only; actual GPU fields and same-geometry sensitivity pending. Active GPU queue and original history untouched.')
 write(target,record)
 registry_path=B/'screened_mesh_registry.json';registry=read(registry_path) if registry_path.exists() else dict(schema='c04_screened_mesh_registry_v1',production_labels=0,meshes=[])
 assert not any(m['candidate']==record['candidate'] for m in registry['meshes'])
 registry['meshes'].append(dict(case_id=case,candidate=record['candidate'],mesh_sha256=r['mesh_sha256'],record_file=str(target.relative_to(W)).replace('\\','/'),GPU_status='not_started',production_eligible=False));registry['updated']=record['created'];write(registry_path,registry);print(json.dumps(dict(frozen=str(target.relative_to(W)),GPU_queue_unchanged=True)))
if __name__=='__main__':main()
