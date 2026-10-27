"""New-regime derivative check for fixed-CAD repair starting with inverted TET10."""
from pathlib import Path
import hashlib,json
import numpy as np
import c04_shape_objective as o
W=Path(__file__).resolve().parents[1];out=W/'data/copv/c04_repair_v3/sicn_negative_J_validation_v1';assert not out.exists()
v=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]]);base=np.vstack([v,*[.5*(v[a]+v[b]) for a,b in o.EDGES]]);base[4]=[.5,0,.8]
points=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[.25,.25,.25],[.1,.2,.3]]);dn=o.reference_derivatives(points);cells=np.arange(10).reshape(1,10);free=np.array([4,5,6,7,8,9]);x=np.zeros(18);ref=np.array([1/np.linalg.det(o.IDEAL)]);scale=np.ones(6)
j=np.einsum('ni,qnj->qij',base,dn);det=np.linalg.det(j);assert det.min()<0 and det.max()>0
loss,g=o.objective_gradient(x,base,cells,free,dn,.07,ref,scale,movement_weight=1e-7);fd=np.empty_like(x);h=1e-6
for k in range(len(x)):
 p=x.copy();m=x.copy();p[k]+=h;m[k]-=h;fd[k]=(o.objective_gradient(p,base,cells,free,dn,.07,ref,scale,movement_weight=1e-7)[0]-o.objective_gradient(m,base,cells,free,dn,.07,ref,scale,movement_weight=1e-7)[0])/(2*h)
error=float(np.max(abs(fd-g)));relative=float(np.linalg.norm(fd-g)/np.linalg.norm(g));assert np.isfinite(loss) and error<1e-6 and relative<1e-6
step=1e-4/max(1,float(np.linalg.norm(g)));new=o.objective_gradient(x-step*g,base,cells,free,dn,.07,ref,scale,movement_weight=1e-7)[0];assert new<loss
out.mkdir();r=dict(passed=True,min_test_J=float(det.min()),max_test_J=float(det.max()),max_coordinate_gradient_error=error,relative_gradient_error=relative,loss_before=loss,loss_after_small_descent=new,helper_sha256=hashlib.sha256(Path(o.__file__).read_bytes()).hexdigest(),scope='Additional negative-J regime finite-difference derivative/descent check for unchanged validated helper. Not a COPV untangling result or mesh acceptance. Final full native Jacobian/SICN/CAD checks remain mandatory.');(out/'review.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
