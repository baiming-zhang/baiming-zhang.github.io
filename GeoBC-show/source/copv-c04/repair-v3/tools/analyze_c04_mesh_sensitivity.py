"""Actual same-CAD q8 mesh comparison at common physical, material-specific probes."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse,json,time
from c04_curved_point_fields import Evaluator,read,sha,np
from c04_repair_contracts import physical_identity
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3'
def ordinary(value):
 if isinstance(value,np.ndarray):return value.tolist()
 if isinstance(value,np.generic):return value.item()
 raise TypeError(type(value).__name__)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--baseline',type=Path,required=True);ap.add_argument('--reference-run',type=Path,required=True);ap.add_argument('--variant',type=Path,required=True);ap.add_argument('--variant-run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert not a.output.exists();start=time.perf_counter()
 validation=read(B/'point_evaluator_validation_v2/review.json');assert validation['passed'] and validation['evaluator_sha256']==sha(W/'notes/c04_curved_point_fields.py')
 runs=[a.reference_run,a.variant_run];candidates=[a.baseline,a.variant];summaries=[read(p/'summary.json') for p in runs];reviews=[read(p/'mesh/mesh_review_curved_checked.json') for p in candidates];contracts=[read(p/'executed_contract.json') for p in runs]
 for c,p,r in zip(candidates,runs,summaries):
  assert r['assembly_quadrature_degree']==8 and r['output_quadrature_degree']==4 and r['matrix_type']=='seqaijcusparse' and r['vector_type']=='seqcuda' and r['case_id']=='copv_c04_0000' and not r['production_eligible']
  assert read(c/'preflight/review.json')['original_mesh_screen_passed'] and sha(p/'solver_executed_snapshot.py')==r['solver_sha256'] and sha(p/'executed_contract.json')==r['contract_sha256']
 assert summaries[0]['solver_sha256']==summaries[1]['solver_sha256'] and summaries[0]['mesh_sha256']!=summaries[1]['mesh_sha256']
 assert all(reviews[0][k]==reviews[1][k] for k in ['parameters_mm','source_CAD_sha256','source_model_sha256'])
 assert [{k:v for k,v in c.items() if k in ['part_id','material','angle','structural']} for c in reviews[0]['parts']]==[{k:v for k,v in c.items() if k in ['part_id','material','angle','structural']} for c in reviews[1]['parts']]
 normalized=[]
 for c in contracts:
  c=dict(c);c['mesh_sha256']='SAME_CAD_DIFFERENT_SCREENED_MESH';normalized.append(physical_identity(c))
 assert normalized[0]==normalized[1],'Different physical inputs beyond the mesh'
 for r in summaries:assert r['cap_convention']=='quadratic-opening-edge' and r['pressure_MPa']==1 and r['energy_work_relative_error']<1e-6
 e0=Evaluator(a.baseline,a.reference_run);gp=np.load(a.reference_run/'GP_xyz_mm.npy',mmap_mode='r');stress=np.load(a.reference_run/'GP_stress_xyz_MPa.npy',mmap_mode='r');probes=[];peak_rows=[]
 # Select actual baseline GP coordinates only after checking that they are
 # interior to their source cell. Edge/face gradients are discontinuous and
 # are not made into an arbitrary common-point stress sample.
 for pid in sorted(e0.parts):
  ids=np.flatnonzero(e0.pids==pid);sampled=ids[np.unique(np.linspace(0,len(ids)-1,4,dtype=int))];best=None
  for cell in sampled:
   for q in range(14):
    v,reason=e0.in_cell(gp[cell,q],int(cell));assert v is not None,reason;r=v['reference_xyz']
    if min(1-r.sum(),*r)<=1e-4:continue
    if q%2==0:probes.append(dict(probe_id=f'p{pid}-cell{int(e0.cell_ids[cell])}-GP{q}',part_id=pid,baseline_cell=int(cell),baseline_GP=q,xyz_mm=gp[cell,q].copy(),reference=v))
    norm=float(np.linalg.norm(stress[cell,q]))
    if best is None or norm>best[0]:best=(norm,int(cell),q,v)
  # Global raw-GP maxima remain separate mesh-dependent diagnostics.
  peak_rows.append(dict(part_id=pid,baseline=summaries[0]['part_statistics'][str(pid)],variant=summaries[1]['part_statistics'][str(pid)]))
  assert best is not None;_,cell,q,v=best
  key=f'p{pid}-cell{int(e0.cell_ids[cell])}-GP{q}'
  if not any(p['probe_id']==key for p in probes):probes.append(dict(probe_id=key,part_id=pid,baseline_cell=cell,baseline_GP=q,xyz_mm=gp[cell,q].copy(),reference=v))
 del e0
 e1=Evaluator(a.variant,a.variant_run);matched=[];failures=[]
 for p in probes:
  value,reason=e1.locate(p['xyz_mm'],p['part_id'])
  if value is not None:
   r=value['reference_xyz']
   if min(1-r.sum(),*r)<=1e-6:reason='variant_element_boundary_gradient_ambiguous';value=None
  if value is None:failures.append(dict(probe_id=p['probe_id'],part_id=p['part_id'],xyz_mm=p['xyz_mm'],reason=reason));continue
  matched.append(dict(**p,variant=value))
 a.output.mkdir(parents=True,exist_ok=False)
 (a.output/'point_location_review.json').write_text(json.dumps(dict(requested_probes=len(probes),matched_probes=len(matched),failed_or_ambiguous_probes=failures,scope='Every failed original probe is retained; no replacement or zero-filled field.'),default=ordinary,indent=2)+'\n')
 (a.output/'common_physical_probe_fields.json').write_text(json.dumps(matched,default=ordinary,indent=2)+'\n')
 assert len(matched)>=6,'Too few unambiguous common physical probes; point-location diagnostic retained, use a new version'
 points=np.array([p['xyz_mm'] for p in matched]);u0=np.array([p['reference']['U_mm'] for p in matched]);u1=np.array([p['variant']['U_mm'] for p in matched]);du=u1-u0;center=points.mean(0);x=points-center;design=np.zeros((len(x),3,6));design[:,:,:3]=np.eye(3)
 for i,(xx,yy,zz) in enumerate(x):design[i,:,3:]=[[0,zz,-yy],[-zz,0,xx],[yy,-xx,0]]
 coefficients,_,rank,_=np.linalg.lstsq(design.reshape(-1,6),du.ravel(),rcond=None);assert rank==6;rigid=(design@coefficients);aligned=du-rigid;parts={};tensor=np.array([1,1,1,2,2,2])
 for pid in sorted({p['part_id'] for p in matched}):
  rows=[p for p in matched if p['part_id']==pid];s0=np.array([p['reference']['stress_xyz_MPa'] for p in rows]);s1=np.array([p['variant']['stress_xyz_MPa'] for p in rows]);parts[pid]=dict(matched_probes=len(rows),stress_tensor_sample_L2_relative_to_baseline=float(np.sqrt(np.sum((s1-s0)**2*tensor)/np.sum(s0**2*tensor))),max_absolute_stress_component_difference_MPa=float(abs(s1-s0).max()),scope='Unweighted common-point sample norm; not a volume norm.')
 report=dict(schema='c04_same_geometry_nonuniform_mesh_sensitivity_v1',created=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),case_id='copv_c04_0000',baseline_mesh_sha256=summaries[0]['mesh_sha256'],variant_mesh_sha256=summaries[1]['mesh_sha256'],source_summary_sha256=[sha(p/'summary.json') for p in runs],same_CAD_and_numerical_physics=True,quadrature_degree=8,energy_relative_difference_to_baseline=abs(summaries[1]['strain_energy_Nmm']/summaries[0]['strain_energy_Nmm']-1),displacement=dict(raw_common_probe_L2_mm=float(np.linalg.norm(du)),rigid_aligned_common_probe_L2_mm=float(np.linalg.norm(aligned)),max_raw_vector_difference_mm=float(np.linalg.norm(du,axis=1).max()),max_rigid_aligned_vector_difference_mm=float(np.linalg.norm(aligned,axis=1).max()),least_squares_translation_mm=coefficients[:3],least_squares_infinitesimal_rotation=coefficients[3:],fitting_center_mm=center,scope='Equal-weight common physical probes. Raw and six-DOF rigid-aligned differences retained; no nearest node/GP substitution.'),parts=parts,requested_probes=len(probes),matched_probes=len(matched),failed_or_ambiguous_probes=failures,material_volume_and_raw_GP_peak_diagnostics=peak_rows,production_eligible=False,limitations=['Two nonuniform meshes, unchanged minimum root/fillet size. This does not establish global h-limit or point-peak convergence.','Every failed/ambiguous original probe retained, with no substituted point or zero-filled field.','Raw GP maxima may occur on material/element boundaries and are mesh-dependent diagnostics.','Same nominal pressure is checked; actual per-run boundary integration/load verification remains bound to each completed run.'],seconds=time.perf_counter()-start)
 (a.output/'review.json').write_text(json.dumps(report,default=ordinary,indent=2)+'\n');print(json.dumps(dict(matched=len(matched),retained_failures=len(failures),energy_relative_difference=report['energy_relative_difference_to_baseline'])))
if __name__=='__main__':main()
