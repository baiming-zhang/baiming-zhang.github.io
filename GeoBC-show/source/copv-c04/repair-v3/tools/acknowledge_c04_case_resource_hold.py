"""Retain an own-case GPU resource stop and acknowledge one versioned retry."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def put(p,v):
 with p.open('x',encoding='utf-8') as f:f.write(json.dumps(v,indent=2)+'\n')
ap=argparse.ArgumentParser();ap.add_argument('--queue-file',type=Path,required=True);ap.add_argument('--expected-run',required=True);ap.add_argument('--stage',choices=['before_domain','before_assembly','before_CUDA'],required=True);a=ap.parse_args();Q=a.queue_file.resolve();assert Q.parent==B/'gpu_validation_queues'
q=read(Q);assert q['status']=='failed_retained';r=q['runs'][-1];assert r['status']=='failed_retained' and Path(r['output']).name==a.expected_run and r['degree']==8
out=W/r['output'];assert out.parent.resolve()==(W/q['candidate']).resolve();assert not (out/'summary.json').exists(),'Completed field must never be treated as resource failure'
assert hashlib.sha256((out/'solver_executed_snapshot.py').read_bytes()).hexdigest()=='c6ec8d3c92af624c13f32ed2d32a776a3dda489c95e2c71fdb48c04780335a57'
resource=read(out/('resources_'+a.stage+'.json'));assert resource['host_available_bytes']>=resource['host_required_bytes'] and resource['GPU_free_MiB']>=resource['GPU_required_MiB'] and resource['GPU_utilization_percent']>30
assert any(v['degree']==6 and v['status']=='completed' for v in q['runs']);assert 'Resources occupied or insufficient' in (out/'execution.log').read_text(encoding='utf-8',errors='replace')
stamp=datetime.now(timezone.utc).isoformat();case=q['case_id'];attempt=int(out.name.rsplit('attempt',1)[1]);snapshot=B/f'gpu_queue_failure_{case}_q8a{attempt}_v1.json'
put(snapshot,q);evidence=out/'resource_interruption_v1.json';put(evidence,dict(schema='c04_GPU_resource_interruption_v1',created=stamp,case_id=case,degree=8,failed_run=r['output'],resource_stage=a.stage,actual_resources=resource,utilization_limit_percent=30,reason='GPU utilization exceeded authorized30% resource gate; hostRAM and freeVRAM passed. No completed q8 field.',numerical_divergence_observed=False,completed_q6_retained=True,retry_policy=f'Retain this run, no guard/source change. Stable5x2s sample and all original resource gates before NEW q8 attempt{attempt+1:03d}; never repeat completed q6.'))
r['error']=q.get('error');r['failure_evidence_file']=str(evidence.relative_to(W)).replace('\\','/');q.setdefault('retry_decisions',[]).append(dict(created=stamp,failed_run=r['output'],action='Preserve failure; wait for resources before NEW pending-degree attempt',evidence=r['failure_evidence_file']));q.update(status='waiting_resources',updated=stamp,resource_reason=f'q8 attempt{attempt:03d} resource stop retained; NEW attempt{attempt+1:03d} only after stable resources');q.pop('error',None)
t=Q.with_suffix('.tmp');t.write_text(json.dumps(q,indent=2)+'\n',encoding='utf-8');t.replace(Q);print(json.dumps(dict(case_id=case,status=q['status'],completed_fields=q['new_FE_labels'],next_degree=8,next_attempt=attempt+1)))
