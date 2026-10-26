"""Single C04 CUDA validation lane, invoked on heartbeat only when resources idle."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse,hashlib,json,os,shlex,shutil,subprocess,time
W=Path(__file__).resolve().parents[1];P=W/'publication/GeoBC-Neural-Operator';B=W/'data/copv/c04_repair_v3';Q=B/'gpu_validation_queue.json'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def stamp():return datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()
def write(p,v):p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix('.tmp');t.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8');t.replace(p)
def lin(p):s=str(p.resolve()).replace('\\','/');return '/mnt/'+s[0].lower()+s[2:]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--launch-if-idle',action='store_true');args=ap.parse_args()
 state=read(Q) if Q.exists() else dict(schema='c04_curved_GPU_validation_queue_v1',candidate='data/copv/c04_repair_v3/copv_c04_0000/attempt027_local_J_from025',status='waiting_resources',production_eligible=False,new_FE_labels=0,runs=[])
 # A live owned controller, and all retained failures, require inspection first.
 if state['status'] in ['running','failed_retained','completed_candidate_validation']:print(json.dumps(state));return
 query=subprocess.run(['wsl','-d','cold-surge-ubuntu-24.04','--','nvidia-smi','--query-gpu=memory.free,utilization.gpu','--format=csv,noheader,nounits'],capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=30)
 if query.returncode:state.update(status='waiting_resources',resource_reason='GPU availability query failed',updated=stamp());write(Q,state);return
 lines=query.stdout.strip().splitlines();assert len(lines)==1;free,util=[int(x.strip()) for x in lines[0].split(',')];state.update(GPU_free_MiB=free,GPU_utilization_percent=util,updated=stamp())
 if free<24576 or util>10 or not args.launch_if_idle:state.update(status='waiting_resources');write(Q,state);print(json.dumps(state));return
 import numpy as np
 candidate=W/state['candidate'];preflight=read(candidate/'preflight/review.json');assert preflight['candidate_FE_validation_allowed'] and not preflight['production_eligible'];assert not (candidate/'local_geometry_hold.json').exists();assert sha(candidate/'mesh/mesh.npz')==preflight['mesh_sha256'];assert read(B/'quadrature_validation_v2/review.json')['passed'];assert shutil.disk_usage(W).free>=40*1024**3
 frozen=B/'gpu_solver_sources_v2'
 if not frozen.exists():
  frozen.mkdir();shutil.copyfile(P/'tools/copv_fenicsx_c04_repair.py',frozen/'solver.py')
  for name in ['__init__.py','copv_material.py','copv_structured_mesh.py','copv_geometry.py']:
   dest=frozen/'geobc'/name;dest.parent.mkdir(exist_ok=True);shutil.copyfile(P/'geobc'/name,dest)
  write(frozen/'source_sha256.json',{str(p.relative_to(frozen)):sha(p) for p in frozen.rglob('*.py')})
 state.update(status='running',pid=os.getpid(),started=stamp());write(Q,state)
 try:
  for degree in [6,8]:
   if any(r['degree']==degree and r['status']=='completed' for r in state['runs']):continue
   attempt=sum(r['degree']==degree for r in state['runs'])+1;out=candidate/f'GPU_q{degree}_attempt{attempt:03d}';out.mkdir(exist_ok=False);start=time.perf_counter();row=dict(degree=degree,status='running',output=str(out.relative_to(W)),started=stamp());state['runs'].append(row);write(Q,state)
   contract=read(W/'data/copv/c04_pilot_v1/validation_contract.json');params=read(candidate/'mesh/mesh_review_curved_checked.json')['parameters_mm'];contract.update(case_id='copv_c04_0000',parent_case_id='copv_c04_0000',geometry_parameters_mm=params,mesh_sha256=preflight['mesh_sha256'],dataset_release='quarantined_repair_validation_only');contract['load'].update(upper_bore_radius_mm=.68*params['boss_bore_radius_mm'],lower_bore_radius_mm=.55*params['boss_bore_radius_mm']);write(out/'executed_contract.json',contract);shutil.copyfile(frozen/'solver.py',out/'solver_executed_snapshot.py')
   runtime=W/'.runtime/fenicsx_cuda';command=[lin(runtime/'env/bin/python'),'-B',lin(out/'solver_executed_snapshot.py'),'--assembly-quadrature',str(degree),'--repair-preflight',lin(candidate/'preflight/review.json'),'--mesh',lin(candidate/'mesh/mesh.npz'),'--review',lin(candidate/'mesh/mesh_review_curved_checked.json'),'--surface',lin(candidate/'preflight/surface.npz'),'--contract',lin(out/'executed_contract.json'),'--output',lin(out),'--degree','2','--preconditioner','p1-gamg','--save-gp-local','--quarantined-research']
   command.append('--save-raw-gp')
   script='#!/bin/bash\nset -euo pipefail\nexport PETSC_OPTIONS="-use_gpu_aware_mpi 0 -no_signal_handler"\nexport OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4\nexport TMPDIR='+shlex.quote(lin(runtime/'tmp'))+'\nexport XDG_CACHE_HOME='+shlex.quote(lin(runtime/'cache'))+'\nexport PATH='+shlex.quote(lin(runtime/'env/bin'))+':"$PATH"\nexport PYTHONPATH='+shlex.quote(lin(frozen))+'\nexec timeout --signal=TERM --kill-after=20s 2400s '+shlex.join(command)+'\n';(out/'run.sh').write_bytes(script.encode())
   with (out/'execution.log').open('x',encoding='utf-8') as log:run=subprocess.run(['wsl','-d','cold-surge-ubuntu-24.04','--','bash',lin(out/'run.sh')],stdout=log,stderr=subprocess.STDOUT,timeout=2470)
   if run.returncode:raise RuntimeError(f'GPU q{degree} returned {run.returncode}; retain run and diagnose before new attempt')
   report=read(out/'summary.json');m=np.load(candidate/'mesh/mesh.npz');u=np.load(out/'nodal_fields.npz');loads=np.load(out/'native_load_and_bc.npz');gp=np.load(out/'GP_stress_local_MPa.npy',mmap_mode='r');assert np.array_equal(m['node_ids'],u['node_ids']) and np.array_equal(m['node_ids'],loads['node_ids']);assert gp.shape==(len(m['cells']),14,6);assert report['matrix_type']=='seqaijcusparse' and report['vector_type']=='seqcuda' and report['displacement_degree']==2 and report['candidate_verification_only'] and not report['production_eligible'];assert report['solver_sha256']==sha(out/'solver_executed_snapshot.py') and report['mesh_sha256']==preflight['mesh_sha256'];assert int(loads['locked_components'].sum())==6
   cell=np.load(out/'cell_display_fields.npz');assert np.array_equal(cell['cell_ids'],m['cell_ids']) and np.array_equal(cell['cell_part_ids'],m['cell_part_ids']);assert np.isfinite(u['U_mm']).all() and np.max(abs(u['U_mm'][loads['locked_components']]))<1e-8
   for first in range(0,len(gp),50000):assert np.isfinite(gp[first:first+50000]).all()
   row.update(status='completed',finished=stamp(),seconds=time.perf_counter()-start,summary=report);state['new_FE_labels']+=1;write(Q,state)
  state.update(status='completed_candidate_validation',finished=stamp());write(Q,state)
 except BaseException as error:
  state.update(status='failed_retained',updated=stamp(),error=repr(error))
  if state['runs'] and state['runs'][-1]['status']=='running':state['runs'][-1].update(status='failed_retained',seconds=time.perf_counter()-start)
  write(Q,state);raise
 print(json.dumps(state))
if __name__=='__main__':main()
