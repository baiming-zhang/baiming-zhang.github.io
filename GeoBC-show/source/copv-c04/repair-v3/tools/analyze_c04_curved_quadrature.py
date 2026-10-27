"""Compare completed same-mesh q6/q8 native fields; report differences, not release."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse,hashlib,json,os,time
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import numpy as np
from c04_repair_contracts import physical_identity
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('candidate',type=Path);ap.add_argument('--q6',type=Path,required=True);ap.add_argument('--q8',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();c=a.candidate.resolve();runs=[a.q6.resolve(),a.q8.resolve()];out=a.output.resolve();assert not out.exists(),'Preserve prior comparison attempts';s=[read(p/'summary.json') for p in runs];pre=read(c/'preflight/review.json')
 for p,r,d in zip(runs,s,[6,8]):
  assert r['status']=='completed' and r['mesh_sha256']==pre['mesh_sha256'] and r['matrix_type']=='seqaijcusparse' and r['vector_type']=='seqcuda' and r['displacement_degree']==2
  assert not r['production_eligible'] and r['candidate_verification_only'] and r['raw_GP_saved']
  q=read(p/'quadrature_protocol.json');assert q['assembly_degree']==d and q['output_degree']==4 and q['output_points']==14
 contracts=[read(p/'executed_contract.json') for p in runs]
 for contract,p,r,degree in zip(contracts,runs,s,[6,8]):
  if contract['discretization'].get('displacement_degree')==2:assert contract['discretization']['volume_quadrature_degree']==degree
  else:
   errata=read(c.parents[1]/'legacy_contract_metadata_errata_v1.json');assert r['contract_sha256']==errata['legacy_contract_sha256'] and errata['numerical_solver_source_unchanged']
 assert all(r['contract_sha256']==sha(p/'executed_contract.json') for p,r in zip(runs,s))
 assert s[0]['case_id']==s[1]['case_id']==pre['case_id'] and physical_identity(contracts[0])==physical_identity(contracts[1]) and s[0]['solver_sha256']==s[1]['solver_sha256'] and s[0]['preconditioner']==s[1]['preconditioner'] and s[0]['cap_convention']==s[1]['cap_convention']=='quadratic-opening-edge'
 m=np.load(c/'mesh/mesh.npz');nodal=[np.load(p/'nodal_fields.npz') for p in runs];loads=[np.load(p/'native_load_and_bc.npz') for p in runs];cell=[np.load(p/'cell_display_fields.npz') for p in runs]
 for n,l,cl in zip(nodal,loads,cell):
  assert np.array_equal(n['node_ids'],m['node_ids']) and np.array_equal(l['node_ids'],m['node_ids']) and np.array_equal(cl['cell_ids'],m['cell_ids']) and np.array_equal(cl['cell_part_ids'],m['cell_part_ids'])
  assert l['locked_components'].sum()==6 and np.isfinite(n['U_mm']).all()
 assert np.array_equal(loads[0]['locked_components'],loads[1]['locked_components'])
 load_diffs={}
 for key in ['consistent_load_N','pressure_load_N','upper_cap_load_N','lower_cap_load_N']:
  x,y=loads[0][key],loads[1][key];assert np.allclose(x,y,rtol=1e-10,atol=1e-9);load_diffs[key]=float(np.max(abs(x-y)))
 u0,u1=nodal[0]['U_mm'],nodal[1]['U_mm'];du=u0-u1;u_ref=np.linalg.norm(u1,axis=1);delta=np.linalg.norm(du,axis=1)
 displacement=dict(relative_nodal_L2=float(np.linalg.norm(du)/np.linalg.norm(u1)),max_vector_difference_mm=float(delta.max()),max_vector_difference_over_q8_peak=float(delta.max()/u_ref.max()),q6_peak_U_mm=float(np.linalg.norm(u0,axis=1).max()),q8_peak_U_mm=float(u_ref.max()),scope='Same native nodes. Nodal L2 is unweighted, not an FE volume norm.')
 arrays=[{k:np.load(p/k,mmap_mode='r') for k in ['GP_stress_local_MPa.npy','GP_stress_xyz_MPa.npy','GP_xyz_mm.npy','GP_volume_weights_mm3.npy']} for p in runs];nc=len(m['cell_ids']);pids=m['cell_part_ids'];parts={int(pid):dict(cells=0,volume_mm3=0.,local_tensor_diff_squared_volume=0.,local_tensor_q8_squared_volume=0.,global_tensor_diff_squared_volume=0.,global_tensor_q8_squared_volume=0.,max_abs_local_component_difference_MPa=0.,max_abs_global_component_difference_MPa=0.) for pid in np.unique(pids)};max_xyz=max_weight=0.;tensor=np.array([1,1,1,2,2,2])
 for ar in arrays:
  assert ar['GP_stress_local_MPa.npy'].shape==(nc,14,6) and ar['GP_stress_xyz_MPa.npy'].shape==(nc,14,6) and ar['GP_xyz_mm.npy'].shape==(nc,14,3) and ar['GP_volume_weights_mm3.npy'].shape==(nc,14)
 for first in range(0,nc,10000):
  sl=slice(first,min(first+10000,nc));w0,w1=[ar['GP_volume_weights_mm3.npy'][sl] for ar in arrays];xyz0,xyz1=[ar['GP_xyz_mm.npy'][sl] for ar in arrays];assert np.isfinite(w0).all() and np.isfinite(w1).all() and (w0>0).all() and (w1>0).all();assert np.allclose(w0,w1,rtol=1e-10,atol=1e-12) and np.allclose(xyz0,xyz1,rtol=0,atol=1e-9);max_xyz=max(max_xyz,float(np.max(abs(xyz0-xyz1))));max_weight=max(max_weight,float(np.max(abs(w0-w1))))
  labels=pids[sl]
  for pid in np.unique(labels):
   mask=labels==pid;rec=parts[int(pid)];rec['cells']+=int(mask.sum());rec['volume_mm3']+=float(w1[mask].sum())
  for name,frame in [('GP_stress_local_MPa.npy','local'),('GP_stress_xyz_MPa.npy','global')]:
   x,y=[ar[name][sl] for ar in arrays];assert np.isfinite(x).all() and np.isfinite(y).all();diff=x-y
   for pid in np.unique(labels):
    mask=labels==pid;rec=parts[int(pid)];rec[frame+'_tensor_diff_squared_volume']+=float(np.sum(w1[mask,...,None]*diff[mask]**2*tensor));rec[frame+'_tensor_q8_squared_volume']+=float(np.sum(w1[mask,...,None]*y[mask]**2*tensor));key='max_abs_'+frame+'_component_difference_MPa';rec[key]=max(rec[key],float(np.max(abs(diff[mask]))))
 for rec in parts.values():
  for frame in ['local','global']:rec[frame+'_tensor_volume_L2_relative']=float(np.sqrt(rec[frame+'_tensor_diff_squared_volume']/rec[frame+'_tensor_q8_squared_volume']))
 energy=[r['strain_energy_Nmm'] for r in s];report=dict(schema='c04_same_mesh_quadrature_field_comparison_v1',created=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),case_id=pre['case_id'],mesh_sha256=pre['mesh_sha256'],runs=[str(p) for p in runs],summary_sha256=[sha(p/'summary.json') for p in runs],assembly_degrees=[6,8],output_GP_points=14,displacement=displacement,strain_energy_Nmm=energy,strain_energy_relative_difference_to_q8=abs(energy[0]-energy[1])/abs(energy[1]),load_max_abs_component_differences_N=load_diffs,GP_coordinate_max_difference_mm=max_xyz,GP_volume_weight_max_difference_mm3=max_weight,parts=parts,production_eligible=False,scope='Actual native-node/14GP same-mesh integration-degree differences. q8 is a comparison reference, not known truth. No automatic production acceptance, mesh convergence, or point-peak stress accuracy claim. Tensor L2 includes factor2 for the three symmetric shear components; no averaging across materials.',seconds=time.perf_counter()-start)
 report.update(executed_contract_sha256=[r['contract_sha256'] for r in s],physical_input_identity_sha256=physical_identity(contracts[0]),contract_identity_scope='Per-run discretization degree and corrected legacy descriptive metadata may differ; each complete executed contract hash is bound to its summary, all physical inputs and actual native load/lock/geometry identities must match.')
 out.mkdir(parents=True,exist_ok=False);(out/'review.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['case_id','displacement','strain_energy_relative_difference_to_q8','production_eligible','seconds']}))
if __name__=='__main__':main()
