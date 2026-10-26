"""Check SICN against native Gmsh and derivatives against finite differences."""
from pathlib import Path
import hashlib,json,time
import gmsh,numpy as np
from c04_shape_objective import IDEAL,IDEAL_INV,sicn_and_gradient,reference_derivatives,objective_gradient,EDGES
W=Path(__file__).resolve().parents[1];O=W/'data/copv/c04_repair_v3/sicn_objective_validation_v1';O.mkdir(exist_ok=False);start=time.perf_counter();rng=np.random.default_rng(17);mat=np.stack([np.eye(3),IDEAL,np.array([[1,.2,.1],[0,.8,.1],[.1,.1,.6]])]);vertices=np.concatenate([np.vstack([np.zeros(3),m.T])+np.array([3*k,0,0]) for k,m in enumerate(mat)]);gmsh.initialize();gmsh.model.add('sicn_probe');gmsh.option.setNumber('Mesh.Renumber',0)
try:
 gmsh.model.addDiscreteEntity(3,1);gmsh.model.mesh.addNodes(3,1,np.arange(1,len(vertices)+1),vertices.ravel());gmsh.model.mesh.addElementsByType(1,4,[1,2,3],np.arange(1,len(vertices)+1));native=np.array(gmsh.model.mesh.getElementQualities([1,2,3],'minSICN'));calculated=sicn_and_gradient(mat@IDEAL_INV)[0];difference=float(np.max(abs(native-calculated)));assert difference<1e-12,(native,calculated)
finally:gmsh.finalize()
a=mat[2]@IDEAL_INV;q,g,_,_=sicn_and_gradient(a);fd=np.zeros((3,3));h=1e-6
for i in range(3):
 for j in range(3):
  plus=a.copy();minus=a.copy();plus[i,j]+=h;minus[i,j]-=h;fd[i,j]=(sicn_and_gradient(plus)[0]-sicn_and_gradient(minus)[0])/(2*h)
matrix_error=float(np.max(abs(g-fd)));assert matrix_error<1e-7
corners=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]]);base=np.vstack([corners,*[.5*(corners[i]+corners[j]) for i,j in EDGES]]);cells=np.arange(10).reshape(1,10);free=np.array([4,5,6,7,8,9]);dN=reference_derivatives([[.1,.2,.3],[.25,.25,.25],[0,0,0],[1,0,0]]);x=rng.normal(size=18)*.01;ref=np.array([1/np.linalg.det(IDEAL)]);scale=np.full(6,.5);loss,gradient=objective_gradient(x,base,cells,free,dN,.95,ref,scale);finite=np.zeros_like(x)
for k in range(len(x)):
 p=x.copy();m=x.copy();p[k]+=h;m[k]-=h;finite[k]=(objective_gradient(p,base,cells,free,dN,.95,ref,scale)[0]-objective_gradient(m,base,cells,free,dN,.95,ref,scale)[0])/(2*h)
geometry_error=float(np.max(abs(gradient-finite)));assert geometry_error<1e-7
report=dict(passed=True,linear_Gmsh_SICN=native.tolist(),linear_formula_SICN=calculated.tolist(),max_native_difference=difference,max_matrix_derivative_error=matrix_error,max_coordinate_objective_derivative_error=geometry_error,helper_sha256=hashlib.sha256((W/'notes/c04_shape_objective.py').read_bytes()).hexdigest(),seconds=time.perf_counter()-start,scope='Three linear tetrahedra native metric comparison and finite-difference derivatives of sampled curved TET10 objective. Not a COPV mesh acceptance or physical solve.');(O/'review.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
