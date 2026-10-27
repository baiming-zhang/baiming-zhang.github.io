"""Retain a diagnosed resource interruption and allow a NEW pending-degree run."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import json
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3';Q=B/'gpu_validation_queue.json'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def exclusive(p,value):
 assert not p.exists(),'Do not overwrite diagnostic evidence';p.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
q=read(Q);assert q['status']=='failed_retained';r=q['runs'][-1];assert r['degree']==8 and r['status']=='failed_retained';out=W/r['output'];resource=read(out/'resources_before_assembly.json');assert resource['host_available_bytes']>=resource['host_required_bytes'] and resource['GPU_free_MiB']>=resource['GPU_required_MiB'] and resource['GPU_utilization_percent']>10
stamp=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat();evidence=dict(schema='c04_GPU_resource_interruption_v1',created=stamp,degree=8,failed_run=r['output'],resource_stage='before_assembly',actual_resources=resource,reason='GPU utilization exceeded unchanged idle guard after the CPU/domain preparation. Host RAM and GPU free-memory gates passed. No completed q8 field.',numerical_divergence_observed=False,completed_q6_retained=True,retry_policy='Wait for original resource gates; NEW q8 attempt002 only, never overwrite attempt001 or rerun completed q6. No GPU/host guard or source changes.')
exclusive(B/'gpu_queue_failure_q8a1_v1.json',q);exclusive(out/'resource_interruption_v1.json',evidence);r['error']=q.get('error');r['failure_evidence_file']=str((out/'resource_interruption_v1.json').relative_to(W)).replace('\\','/');q.setdefault('retry_decisions',[]).append(dict(created=stamp,failed_run=r['output'],action='Retain failure and wait before NEW q8 run',evidence=r['failure_evidence_file']));q.update(status='waiting_resources',updated=stamp,resource_reason='q8 attempt001 retained: GPU utilization became busy; pending NEW q8 attempt002');q.pop('error',None);t=Q.with_suffix('.tmp');t.write_text(json.dumps(q,indent=2)+'\n');t.replace(Q);print(json.dumps(dict(status=q['status'],completed_fields=q['new_FE_labels'],next_degree=8,next_attempt=2)))
