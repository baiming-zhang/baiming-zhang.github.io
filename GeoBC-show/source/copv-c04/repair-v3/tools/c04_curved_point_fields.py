"""Actual native TET10 P2 evaluation; candidate search never substitutes a field."""
from pathlib import Path
import hashlib,json,os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import numpy as np
from scipy.spatial import cKDTree
from c04_shape_objective import reference_derivatives,EDGES
from geobc.copv_structured_mesh import analytic_segments
from geobc.copv_material import engineering_stiffness
PAIRS=[(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)]
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
 return h.hexdigest()
def basis(r):
 lam=np.r_[1-np.sum(r),r];return np.r_[lam*(2*lam-1),[4*lam[a]*lam[b] for a,b in EDGES]]
def inverse(x,point,tolerance_mm=1e-11):
 """Positive-J Newton inversion; return a reason instead of an extrapolated value."""
 try:
  r=np.linalg.solve((x[1:4]-x[0]).T,point-x[0])
  for iteration in range(25):
   d=reference_derivatives(r[None])[0];j=x.T@d
   if np.linalg.det(j)<=0:return None,'nonpositive_local_J'
   delta=basis(r)@x-point
   if np.linalg.norm(delta)<=tolerance_mm:
    if min(1-r.sum(),*r)<-1e-8:return None,'outside_reference_tetrahedron'
    return r,None
   step=np.linalg.solve(j,delta);r-=step
   if not np.isfinite(r).all() or np.linalg.norm(r)>20:return None,'Newton_diverged'
  return None,'Newton_iteration_limit'
 except np.linalg.LinAlgError:return None,'singular_J'
class Evaluator:
 def __init__(self,candidate,run):
  self.c=Path(candidate);self.run=Path(run);self.summary=read(self.run/'summary.json');self.contract=read(self.run/'executed_contract.json');pre=read(self.c/'preflight/review.json');review=read(self.c/'mesh/mesh_review_curved_checked.json')
  assert self.summary['status']=='completed' and self.summary['displacement_degree']==2 and self.summary['raw_GP_saved'] and not self.summary['production_eligible']
  assert pre['original_mesh_screen_passed'] and pre['mesh_sha256']==self.summary['mesh_sha256']==sha(self.c/'mesh/mesh.npz') and self.summary['contract_sha256']==sha(self.run/'executed_contract.json')
  m=np.load(self.c/'mesh/mesh.npz');n=np.load(self.run/'nodal_fields.npz');assert np.array_equal(n['node_ids'],m['node_ids']);self.xyz=m['xyz'];self.cells=m['cells'];self.cell_ids=m['cell_ids'];self.pids=m['cell_part_ids'];self.u=n['U_mm'];assert np.isfinite(self.u).all();self.parts={int(p['part_id']):p for p in review['parts'] if p['structural']};self.matrices={}
  for name,p in self.contract['materials'].items():
   if 'E_MPa' in p:
    e=p['E_MPa'];nu=p['nu'];g=e/(2*(1+nu));self.matrices[name]=engineering_stiffness(E1=e,E2=e,E3=e,G12=g,G13=g,G23=g,nu12=nu,nu13=nu,nu23=nu)
   else:self.matrices[name]=engineering_stiffness(**{k:p[k] for k in ['E1','E2','E3','G12','G13','G23','nu12','nu13','nu23']})
  rz,normal=analytic_segments(self.summary['geometry_parameters_mm'],(4,2500,16,2500,4));self.profile_tree=cKDTree(rz);self.tangent=np.c_[-normal[:,1],normal[:,0]];self.search={}
 def axes(self,point,pid):
  part=self.parts[int(pid)]
  if part['material']!='carbon':return np.eye(3)
  radius=np.linalg.norm(point[:2]);assert radius>0;_,near=self.profile_tree.query([radius,point[2]]);t=self.tangent[near];er=np.r_[point[:2]/radius,0.];et=np.array([-er[1],er[0],0.]);mer=-t[0]*er;mer[2]-=t[1];normal=np.cross(mer,et);a=np.deg2rad(part.get('angle',0));q=np.stack([np.cos(a)*mer+np.sin(a)*et,-np.sin(a)*mer+np.cos(a)*et,normal],axis=-1);assert abs(np.linalg.det(q)-1)<1e-10;return q
 def in_cell(self,point,cell):
  conn=self.cells[cell];x=self.xyz[conn];r,reason=inverse(x,np.asarray(point))
  if reason:return None,reason
  d=reference_derivatives(r[None])[0];grad=self.u[conn].T@d@np.linalg.inv(x.T@d);strain=.5*(grad+grad.T);pid=int(self.pids[cell]);q=self.axes(point,pid);local=q.T@strain@q;ev=np.array([local[i,j]*(1 if i==j else 2) for i,j in PAIRS]);sl=self.matrices[self.parts[pid]['material']]@ev;tensor=np.zeros((3,3))
  for v,(i,j) in zip(sl,PAIRS):tensor[i,j]=tensor[j,i]=v
  sg=q@tensor@q.T
  return dict(cell_index=int(cell),native_cell_id=int(self.cell_ids[cell]),part_id=pid,reference_xyz=r,U_mm=basis(r)@self.u[conn],stress_local_MPa=sl,stress_xyz_MPa=np.array([sg[i,j] for i,j in PAIRS]),physical_residual_mm=float(np.linalg.norm(basis(r)@x-point))),None
 def locate(self,point,pid,max_candidates=128):
  """Bounded centroid candidate search, followed by exact curved containment.

  A missed cell is retained as a location failure; no nearest-cell field is returned.
  """
  pid=int(pid)
  if pid not in self.search:
   ids=np.flatnonzero(self.pids==pid);centers=np.empty((len(ids),3))
   for first in range(0,len(ids),10000):centers[first:first+10000]=np.mean(self.xyz[self.cells[ids[first:first+10000],:4]],axis=1)
   self.search[pid]=(ids,cKDTree(centers))
  ids,tree=self.search[pid];_,near=tree.query(point,k=min(max_candidates,len(ids)))
  for item in np.atleast_1d(near):
   value,reason=self.in_cell(point,int(ids[item]))
   if value is not None:return value,None
  return None,'no_containing_same_material_TET10_among_'+str(min(max_candidates,len(ids)))+'_candidates'
