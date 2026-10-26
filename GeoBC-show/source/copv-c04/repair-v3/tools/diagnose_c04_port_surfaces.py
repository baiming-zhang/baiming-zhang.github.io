"""Use existing measured bad-face quality and exact CAD dimensions once."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse,hashlib,json,time
import gmsh
W=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('candidate',type=Path);a=ap.parse_args();c=a.candidate.resolve();o=c/'port_surface_diagnosis_v1';o.mkdir(exist_ok=False);start=time.perf_counter();r=read(c/'geometry_check/review.json');review=read(c/'mesh/mesh_review.json');parent=W/review['parent_attempt'];state=read(parent/'state.json');cad=Path(state['partition_file']);assert hashlib.sha256(cad.read_bytes()).hexdigest()==state['partition_sha256'];rows=[]
 gmsh.initialize();gmsh.option.setNumber('General.NumThreads',1)
 try:
  gmsh.model.occ.importShapes(str(cad));gmsh.model.occ.synchronize()
  for f in r['faces']:
   if f['min_surface_SICN']>=.05:continue
   tag=f['surface'];area=gmsh.model.occ.getMass(2,tag);curves=[]
   for d,k in gmsh.model.getBoundary([(2,tag)],combined=False):
    k=abs(k);curves.append(dict(curve=k,type=gmsh.model.getType(1,k),length_mm=gmsh.model.occ.getMass(1,k),bbox_mm=list(gmsh.model.getBoundingBox(1,k))))
   perimeter=sum(k['length_mm'] for k in curves);rows.append(dict(**f,CAD_area_mm2=area,CAD_perimeter_mm=perimeter,narrow_strip_width_estimate_2A_over_P_mm=2*area/perimeter,CAD_bbox_mm=list(gmsh.model.getBoundingBox(2,tag)),CAD_center_mm=list(gmsh.model.occ.getCenterOfMass(2,tag)),curves=curves))
 finally:gmsh.finalize()
 report=dict(schema='c04_bad_port_CAD_dimensions_v1',created=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),case_id=r['case_id'],CAD_sha256=state['partition_sha256'],measured_geometry_review_sha256=hashlib.sha256((c/'geometry_check/review.json').read_bytes()).hexdigest(),bad_faces=rows,seconds=time.perf_counter()-start,scope='Existing actual curved surface SICN, exact CAD area/perimeter/edges. 2A/P is a strip-sizing heuristic, not a rigorous minimum width or quality guarantee. No CAD or mesh changed.')
 (o/'review.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps([dict(surface=f['surface'],type=f['type'],parts=f['part_ids'],width_est_mm=f['narrow_strip_width_estimate_2A_over_P_mm'],edge_lengths=[k['length_mm'] for k in f['curves']]) for f in rows]))
if __name__=='__main__':main()
