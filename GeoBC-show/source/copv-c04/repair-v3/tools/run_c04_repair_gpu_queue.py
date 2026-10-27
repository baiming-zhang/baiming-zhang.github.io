"""Single C04 CUDA validation lane, invoked on heartbeat only when resources idle."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse,ast,hashlib,json,os,shlex,shutil,subprocess,time
from c04_repair_contracts import build as build_case_contract,physical_identity
W=Path(__file__).resolve().parents[1];P=W/'publication/GeoBC-Neural-Operator';B=W/'data/copv/c04_repair_v3';Q=B/'gpu_validation_queue.json'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def stamp():return datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()
def write(p,v):p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix('.tmp');t.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8');t.replace(p)
def lin(p):s=str(p.resolve()).replace('\\','/');return '/mnt/'+s[0].lower()+s[2:]
def main():
 global Q
 ap=argparse.ArgumentParser();ap.add_argument('--launch-if-idle',action='store_true');ap.add_argument('--queue-file',type=Path);args=ap.parse_args()
 if args.queue_file:
  selected=args.queue_file.resolve();assert selected.parent==B/'gpu_validation_queues' and selected.exists(),'Only declared case-specific C04 queues are supported'
  if read(Q)['status']!='completed_candidate_validation':print(json.dumps(dict(status='waiting_primary_0000_validation')));return
  for p in selected.parent.glob('*.json'):
   if p.resolve()!=selected and read(p).get('status')=='running':print(json.dumps(dict(status='waiting_other_C04_GPU_lane',queue=str(p))));return
  Q=selected
 state=read(Q) if Q.exists() else dict(schema='c04_curved_GPU_validation_queue_v1',candidate='data/copv/c04_repair_v3/copv_c04_0000/attempt027_local_J_from025',status='waiting_resources',production_eligible=False,new_FE_labels=0,runs=[])
 # A live owned controller, and all retained failures, require inspection first.
 if state['status'] in ['running','failed_retained','completed_candidate_validation']:print(json.dumps(state));return
 query=subprocess.run(['wsl','-d','cold-surge-ubuntu-24.04','--','nvidia-smi','--query-gpu=memory.free,utilization.gpu','--format=csv,noheader,nounits'],capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=30)
 if query.returncode:state.update(status='waiting_resources',resource_reason='GPU availability query failed',updated=stamp());write(Q,state);return
 lines=query.stdout.strip().splitlines();assert len(lines)==1;free,util=[int(x.strip()) for x in lines[0].split(',')];state.update(GPU_free_MiB=free,GPU_utilization_percent=util,updated=stamp())
 if free<24576 or util>10 or not args.launch_if_idle:state.update(status='waiting_resources');write(Q,state);print(json.dumps(state));return
 import numpy as np
 candidate=W/state['candidate'];preflight=read(candidate/'preflight/review.json');assert preflight['candidate_FE_validation_allowed'] and not preflight['production_eligible'];assert not (candidate/'local_geometry_hold.json').exists();assert sha(candidate/'mesh/mesh.npz')==preflight['mesh_sha256'];assert read(B/'quadrature_validation_v2/review.json')['passed'];assert shutil.disk_usage(W).free>=40*1024**3
 assert read(B/'compact_frame_validation_v1/review.json')['passed']
 identity=read(B/'compact_frame_validation_v1/helper_identity.json');helper=next(n for n in ast.parse((P/'tools/copv_fenicsx_c04_repair.py').read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='compact_material_frame');assert hashlib.sha256(ast.dump(helper,include_attributes=False).encode()).hexdigest()==identity['helper_AST_sha256'];assert sha(B/'compact_frame_validation_v1/review.json')==identity['validation_review_sha256']
 # Apply the solver's existing before-domain RAM requirement before creating
 # a GPU run, so an active C04 CPU mesher causes waiting rather than a known
 # resource failure. The solver still checks real resources at every stage.
 degrees=state.get('assembly_degrees',[6,8]);assert degrees in [[6,8],[8]],'Only canonical 6/8 or declared mesh-sensitivity 8-only runs are supported'
 if degrees==[8]:assert state.get('purpose')=='same_geometry_mesh_sensitivity_only','8-only mode is reserved for the declared same-geometry sensitivity study'
 pending=[d for d in degrees if not any(r['degree']==d and r['status']=='completed' for r in state['runs'])]
 host_query=subprocess.run(['wsl','-d','cold-surge-ubuntu-24.04','--','cat','/proc/meminfo'],capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=30)
 available=[int(line.split()[1])*1024 for line in host_query.stdout.splitlines() if line.startswith('MemAvailable:')]
 if host_query.returncode or len(available)!=1:state.update(status='waiting_resources',resource_reason='WSL host availability query failed',updated=stamp());write(Q,state);return
 nc=int(read(candidate/'mesh/mesh_review_curved_checked.json')['cells']);required=8*1024**3+nc*max(pending,default=8)*1024
 state.update(host_available_bytes=available[0],host_before_domain_required_bytes=required)
 if available[0]<required:state.update(status='waiting_resources',resource_reason='WSL RAM below unchanged solver before-domain guard; wait for active CPU work',updated=stamp());write(Q,state);print(json.dumps(state));return
 state.pop('resource_reason',None)
 frozen=B/'gpu_solver_sources_v3'
 contract_review=read(B/'contract_metadata_validation_v2/review.json');assert contract_review['passed'] and sha(W/'notes/c04_repair_contracts.py')==contract_review['builder_sha256'],'Contract builder changed; preserve validation and declare a new version'
 if frozen.exists():assert sha(frozen/'solver.py')==sha(P/'tools/copv_fenicsx_c04_repair.py'),'Frozen solver differs; preserve it and declare a NEW source version'
 if not frozen.exists():
  frozen.mkdir();shutil.copyfile(P/'tools/copv_fenicsx_c04_repair.py',frozen/'solver.py')
  for name in ['__init__.py','copv_material.py','copv_structured_mesh.py','copv_geometry.py']:
   dest=frozen/'geobc'/name;dest.parent.mkdir(exist_ok=True);shutil.copyfile(P/'geobc'/name,dest)
  write(frozen/'source_sha256.json',{str(p.relative_to(frozen)):sha(p) for p in frozen.rglob('*.py')})
 state.update(status='running',pid=os.getpid(),started=stamp());write(Q,state)
 try:
  for degree in degrees:
   if any(r['degree']==degree and r['status']=='completed' for r in state['runs']):continue
   attempt=sum(r['degree']==degree for r in state['runs'])+1;out=candidate/f'GPU_q{degree}_attempt{attempt:03d}';out.mkdir(exist_ok=False);start=time.perf_counter();row=dict(degree=degree,status='running',output=str(out.relative_to(W)),started=stamp());state['runs'].append(row);write(Q,state)
   contract=build_case_contract(candidate,degree);row['case_id']=contract['case_id'];row['physical_input_identity_sha256']=physical_identity(contract);write(out/'executed_contract.json',contract);shutil.copyfile(W/'notes/c04_repair_contracts.py',out/'contract_builder_executed.py');shutil.copyfile(frozen/'solver.py',out/'solver_executed_snapshot.py');write(Q,state)
   runtime=W/'.runtime/fenicsx_cuda';command=[lin(runtime/'env/bin/python'),'-B',lin(out/'solver_executed_snapshot.py'),'--assembly-quadrature',str(degree),'--repair-preflight',lin(candidate/'preflight/review.json'),'--mesh',lin(candidate/'mesh/mesh.npz'),'--review',lin(candidate/'mesh/mesh_review_curved_checked.json'),'--surface',lin(candidate/'preflight/surface.npz'),'--contract',lin(out/'executed_contract.json'),'--output',lin(out),'--degree','2','--preconditioner','p1-gamg','--save-gp-local','--quarantined-research']
   command.append('--save-raw-gp')
   script='#!/bin/bash\nset -euo pipefail\nexport PETSC_OPTIONS="-use_gpu_aware_mpi 0 -no_signal_handler"\nexport OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4\nexport TMPDIR='+shlex.quote(lin(runtime/'tmp'))+'\nexport XDG_CACHE_HOME='+shlex.quote(lin(runtime/'cache'))+'\nexport PATH='+shlex.quote(lin(runtime/'env/bin'))+':"$PATH"\nexport PYTHONPATH='+shlex.quote(lin(frozen))+'\nexec timeout --signal=TERM --kill-after=20s 2400s '+shlex.join(command)+'\n';(out/'run.sh').write_bytes(script.encode())
   with (out/'execution.log').open('x',encoding='utf-8') as log:run=subprocess.run(['wsl','-d','cold-surge-ubuntu-24.04','--','bash',lin(out/'run.sh')],stdout=log,stderr=subprocess.STDOUT,timeout=2470)
   if run.returncode:raise RuntimeError(f'GPU q{degree} returned {run.returncode}; retain run and diagnose before new attempt')
   report=read(out/'summary.json');m=np.load(candidate/'mesh/mesh.npz');u=np.load(out/'nodal_fields.npz');loads=np.load(out/'native_load_and_bc.npz');gp=np.load(out/'GP_stress_local_MPa.npy',mmap_mode='r');assert np.array_equal(m['node_ids'],u['node_ids']) and np.array_equal(m['node_ids'],loads['node_ids']);assert gp.shape==(len(m['cells']),14,6);assert report['matrix_type']=='seqaijcusparse' and report['vector_type']=='seqcuda' and report['displacement_degree']==2 and report['candidate_verification_only'] and not report['production_eligible'];assert report['solver_sha256']==sha(out/'solver_executed_snapshot.py') and report['mesh_sha256']==preflight['mesh_sha256'];assert int(loads['locked_components'].sum())==6
   assert report['contract_sha256']==sha(out/'executed_contract.json') and report['case_id']==contract['case_id'] and report['geometry_parameters_mm']==contract['geometry_parameters_mm'] and report['assembly_quadrature_degree']==degree
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
