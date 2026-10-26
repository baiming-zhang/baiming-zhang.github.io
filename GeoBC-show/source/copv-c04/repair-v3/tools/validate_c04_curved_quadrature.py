"""Small CPU-only FEniCSx check of generated assembly/output frame separation."""
from pathlib import Path
import ast,json,time,sysconfig
_cfg=[sysconfig.get_config_vars()]
try:
 from setuptools._distutils import sysconfig as _dc
 _cfg.append(_dc.get_config_vars())
except ImportError:pass
for cfg in _cfg:cfg.update(CFLAGS='-O2 -fPIC',LDFLAGS='',CC='/usr/bin/gcc',CXX='/usr/bin/g++',LDSHARED='/usr/bin/gcc -shared',CPPFLAGS='')
import numpy as np
from scipy.spatial import cKDTree
from mpi4py import MPI
import basix,basix.ufl,ufl
from dolfinx import mesh,fem
from dolfinx.fem import petsc as fp
W=Path(__file__).resolve().parents[1];O=W/'data/copv/c04_repair_v3/quadrature_validation_v2';O.mkdir(exist_ok=False);start=time.perf_counter();source=W/'publication/GeoBC-Neural-Operator/tools/copv_fenicsx_c04_repair.py';s=source.read_text();main=next(n for n in ast.parse(s).body if isinstance(n,ast.FunctionDef) and n.name=='main');frame=next(n for n in main.body if isinstance(n,ast.FunctionDef) and n.name=='material_axes');scope=dict(np=np,part_local=np.array([3]),carbon={3},tree=cKDTree([[10.,0.]]),tangent=np.array([[1.,0.]]),angle_by_id={3:25.});exec(compile(ast.fix_missing_locations(ast.Module(body=[frame],type_ignores=[])),str(source),'exec'),scope);axes_fn=scope['material_axes']
# One genuinely curved isoparametric TET10 with positive Jacobians.
corners=np.array([[10.,0.,0.],[12.,0.,0.],[10.,1.,0.],[10.,0.,1.]]);edges=[(0,1),(1,2),(2,0),(0,3),(2,3),(1,3)];xyz=np.vstack([corners,*[.5*(corners[a]+corners[b]) for a,b in edges]]);xyz[4,2]+=.03;xyz[7,1]+=.02
perm=[0,1,2,3]+[4+next(i for i,e in enumerate(edges) if set(e)==set(edge)) for edge in basix.topology(basix.CellType.tetrahedron)[1]];element=basix.ufl.element('Lagrange','tetrahedron',2,shape=(3,),lagrange_variant=basix.LagrangeVariant.equispaced);domain=mesh.create_mesh(MPI.COMM_SELF,np.array([perm],np.int64),ufl.Mesh(element),xyz);V=fem.functionspace(domain,element);field=fem.Function(V);field.interpolate(lambda x:np.vstack([x[0],2*x[1],3*x[2]]));scalar=basix.create_element(basix.ElementFamily.P,basix.CellType.tetrahedron,2,basix.LagrangeVariant.equispaced);nodes=domain.geometry.x[domain.geometry.dofmap[0]];C=np.diag([11.,13.,17.,3.,5.,7.]);ids=np.array([0],np.int32);rows=[]
def eng(e):return ufl.as_vector([e[0,0],e[1,1],e[2,2],2*e[1,2],2*e[0,2],2*e[0,1]])
for degree in [4,6,8]:
 qp,qw=basix.make_quadrature(basix.CellType.tetrahedron,degree);tab=scalar.tabulate(1,qp);positions=np.einsum('qi,ij->qj',tab[0,:,:,0],nodes);frames=axes_fn(positions[None],ids)[0];Qspace=fem.functionspace(domain,basix.ufl.quadrature_element('tetrahedron',value_shape=(3,3),degree=degree,scheme='default'));Q=fem.Function(Qspace);Q.x.array.reshape(-1,9)[Qspace.dofmap.cell_dofs(0)]=frames.reshape(-1,9);strain=eng(Q.T*ufl.sym(ufl.grad(field))*Q);dx=ufl.Measure('dx',domain=domain,metadata={'quadrature_degree':degree});energy=float(fem.assemble_scalar(fem.form(.5*ufl.dot(strain,ufl.as_matrix(C)*strain)*dx)));u=ufl.TrialFunction(V);v=ufl.TestFunction(V);a=fem.form(ufl.dot(eng(Q.T*ufl.sym(ufl.grad(v))*Q),ufl.as_matrix(C)*eng(Q.T*ufl.sym(ufl.grad(u))*Q))*dx);A=fp.assemble_matrix(a);A.assemble();action=A.createVecRight();A.mult(field.x.petsc_vec,action);matrix_energy=.5*field.x.petsc_vec.dot(action);assert abs(matrix_energy/energy-1)<1e-10;action.destroy();A.destroy()
 # Output uses its own 14-point frame field, independent of assembly degree.
 outqp,outqw=basix.make_quadrature(basix.CellType.tetrahedron,4);outtab=scalar.tabulate(1,outqp);outpos=np.einsum('qi,ij->qj',outtab[0,:,:,0],nodes);outaxes=axes_fn(outpos[None],ids)[0];Qoutspace=fem.functionspace(domain,basix.ufl.quadrature_element('tetrahedron',value_shape=(3,3),degree=4,scheme='default'));Qout=fem.Function(Qoutspace);Qout.x.array.reshape(-1,9)[Qoutspace.dofmap.cell_dofs(0)]=outaxes.reshape(-1,9);stress=fem.Expression(ufl.as_matrix(C)*eng(Qout.T*ufl.sym(ufl.grad(field))*Qout),outqp).eval(domain,ids).reshape(14,6);local=np.swapaxes(outaxes,-1,-2)@np.diag([1.,2.,3.])@outaxes;strain_np=np.stack([local[:,0,0],local[:,1,1],local[:,2,2],2*local[:,1,2],2*local[:,0,2],2*local[:,0,1]],axis=-1);error=float(np.max(abs(stress-strain_np@C.T)));assert error<1e-9
 jac=np.stack([np.einsum('qi,ij->qj',outtab[j,:,:,0],nodes) for j in [1,2,3]],axis=-1);assert np.linalg.det(jac).min()>0
 rows.append(dict(assembly_degree=degree,assembly_points=len(qp),output_points=14,energy=energy,matrix_energy=float(matrix_energy),max_output_stress_error=error,curved_volume=float(np.linalg.det(jac)@outqw)))
report=dict(passed=True,scope='One curved TET10 CPU FEniCSx unit probe, analytic affine physical displacement, anisotropic spatially varying material frame. Not a COPV GPU solve or convergence certification.',rows=rows,seconds=time.perf_counter()-start);(O/'review.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
