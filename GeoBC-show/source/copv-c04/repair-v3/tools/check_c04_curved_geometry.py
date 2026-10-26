"""Full boundary-node and sampled Tri6-interior distance to owning CAD surfaces."""
from pathlib import Path
import argparse,hashlib,json,shutil,time
import numpy as np
import gmsh
W=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def write(p,d):p.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);ap.add_argument('--samples-per-face',type=int,default=128);args=ap.parse_args();root=args.root.resolve();out=root/'geometry_check';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'checker_executed.py');review=read(root/'mesh/mesh_review.json');case=read(root/'result.json')['case_id'];parent=W/review['parent_attempt'];state=read(parent/'state.json');partition=Path(state['partition_file']) if 'partition_file' in state else W/'data/copv/c04_pilot_v1/cases'/case/'partition_v2/structural_partition.brep';assert hashlib.sha256(partition.read_bytes()).hexdigest()==state['partition_sha256'];start=time.perf_counter();rows=[]
 gmsh.initialize();gmsh.option.setNumber('General.NumThreads',1)
 try:
  gmsh.model.occ.importShapes(str(partition));gmsh.model.occ.synchronize();cad={tag:gmsh.model.getType(2,tag) for _,tag in gmsh.model.getEntities(2)};gmsh.merge(str(root/'mesh/mesh.msh'));assert all(gmsh.model.getType(2,tag)==typ for tag,typ in cad.items());parts={int(t):int(p) for _,p in gmsh.model.getPhysicalGroups(3) for t in gmsh.model.getEntitiesForPhysicalGroup(3,p)};owners={}
  for volume,part in parts.items():
   for dim,face in gmsh.model.getBoundary([(3,volume)],combined=False):owners.setdefault(abs(face),set()).add(part)
  tags,coords,_=gmsh.model.mesh.getNodes();tags=np.asarray(tags);xyz=np.asarray(coords).reshape(-1,3);order=np.argsort(tags);sortedtags=tags[order];lam=np.array([[1/3]*3,[.6,.2,.2],[.2,.6,.2],[.2,.2,.6]]);basis=np.c_[lam*(2*lam-1),4*lam[:,0]*lam[:,1],4*lam[:,1]*lam[:,2],4*lam[:,2]*lam[:,0]]
  def distance(surface,points):
   worst=0.
   for first in range(0,len(points),5000):
    p=points[first:first+5000];q,_=gmsh.model.getClosestPoint(2,surface,p.ravel());q=np.asarray(q).reshape(-1,3);assert np.isfinite(q).all();worst=max(worst,float(np.linalg.norm(q-p,axis=1).max()))
   return worst
  for _,surface in gmsh.model.getEntities(2):
   types,ids,cc=gmsh.model.mesh.getElements(2,surface)
   if not len(types):continue
   assert list(types)==[9];native=np.asarray(cc[0]).reshape(-1,6);nodeids,positions,_=gmsh.model.mesh.getNodes(2,surface,includeBoundary=True);positions=np.asarray(positions).reshape(-1,3);node_error=distance(surface,positions);indices=order[np.searchsorted(sortedtags,native)];sample=indices[np.unique(np.linspace(0,len(indices)-1,min(args.samples_per_face,len(indices)),dtype=int))];points=np.einsum('qn,fni->fqi',basis,xyz[sample]).reshape(-1,3);sample_error=distance(surface,points);quality=np.asarray(gmsh.model.mesh.getElementQualities(ids[0],'minSICN'));partids=sorted(owners[surface]);tol=min([.015]+[.003 if p in [7,8] else .0075 if p==2 else .015 for p in partids]);rows.append(dict(surface=surface,type=cad[surface],part_ids=partids,nodes=len(nodeids),triangles=len(native),max_node_CAD_distance_mm=node_error,sampled_interior_points=len(points),max_sampled_interior_CAD_distance_mm=sample_error,project_screen_tolerance_mm=tol,min_surface_SICN=float(quality.min()),below_005=int((quality<.05).sum()),passed=bool(max(node_error,sample_error)<=tol)));write(out/'progress.json',dict(completed_faces=len(rows),seconds=time.perf_counter()-start,latest=rows[-1]))
  record=dict(schema='c04_local_CAD_screen_v1',case_id=case,mesh_sha256=hashlib.sha256((root/'mesh/mesh.npz').read_bytes()).hexdigest(),CAD_sha256=hashlib.sha256(partition.read_bytes()).hexdigest(),method='All nodes classified on each CAD surface, plus four Tri6 interior reference points on up to the declared number of deterministically selected triangles per CAD face; Euclidean distance to Gmsh getClosestPoint on owning surface.',tolerance_basis='Additional project candidate-screening tolerances, not industry standards: collar .003mm (10% of nominal .03mm onset excess), adhesive .0075mm (5% of nominal .15mm bondline), other surfaces .015mm. This sampled interior test is not a rigorous global geometry error bound.',samples_per_face=args.samples_per_face,faces=rows,max_node_CAD_distance_mm=max(r['max_node_CAD_distance_mm'] for r in rows),max_sampled_interior_CAD_distance_mm=max(r['max_sampled_interior_CAD_distance_mm'] for r in rows),min_surface_SICN=min(r['min_surface_SICN'] for r in rows),surface_below_005=sum(r['below_005'] for r in rows),passed=bool(all(r['passed'] for r in rows)),production_eligible=False,seconds=time.perf_counter()-start);write(out/'review.json',record);print(json.dumps({k:v for k,v in record.items() if k!='faces'}),flush=True)
 finally:gmsh.finalize()
if __name__=='__main__':main()
