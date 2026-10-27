"""Validate new P2 point/gradient evaluator against actual completed native GP fields."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse,json,time
from c04_curved_point_fields import Evaluator,read,sha,np
def main():
 ap=argparse.ArgumentParser();ap.add_argument('candidate',type=Path);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert not a.output.exists(),'Preserve validation attempts';a.output.mkdir(parents=True);start=time.perf_counter()
 try:
  e=Evaluator(a.candidate,a.run);coords=np.load(a.run/'GP_xyz_mm.npy',mmap_mode='r');local=np.load(a.run/'GP_stress_local_MPa.npy',mmap_mode='r');glob=np.load(a.run/'GP_stress_xyz_MPa.npy',mmap_mode='r');rows=[];diff2=ref2=0.;maxabs=maxres=0.
  for pid in sorted(e.parts):
   ids=np.flatnonzero(e.pids==pid)
   for cell in ids[np.unique(np.linspace(0,len(ids)-1,3,dtype=int))]:
    for gp in range(14):
     value,reason=e.in_cell(coords[cell,gp],int(cell));assert value is not None,(pid,int(cell),gp,reason)
     for key,original in [('stress_local_MPa',local[cell,gp]),('stress_xyz_MPa',glob[cell,gp])]:
      delta=value[key]-original;diff2+=float(delta@delta);ref2+=float(original@original);maxabs=max(maxabs,float(abs(delta).max()))
     maxres=max(maxres,value['physical_residual_mm']);rows.append(dict(native_cell_id=int(e.cell_ids[cell]),part_id=pid,GP_index=gp))
  relative=float(np.sqrt(diff2/ref2));assert relative<1e-6 and maxabs<1e-4 and maxres<=1e-8,(relative,maxabs,maxres)
  report=dict(passed=True,created=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),scope='New actual quadratic inverse geometry, P2 displacement gradient and material stress reconstruction compared with saved real q6 native 14GP fields at 3 deterministic cells per material. This validates the evaluator only; no new FEA, quadrature or mesh-convergence result.',points=len(rows),sample_native_IDs=rows,combined_stress_component_relative_L2=relative,max_absolute_stress_component_difference_MPa=maxabs,max_inverse_geometry_residual_mm=maxres,mesh_sha256=e.summary['mesh_sha256'],source_summary_sha256=sha(a.run/'summary.json'),evaluator_sha256=sha(Path(__file__).with_name('c04_curved_point_fields.py')),production_eligible=False,seconds=time.perf_counter()-start)
  (a.output/'review.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='sample_native_IDs'}))
 except BaseException as error:
  (a.output/'failure.json').write_text(json.dumps(dict(error=repr(error),seconds=time.perf_counter()-start),indent=2)+'\n');raise
if __name__=='__main__':main()
