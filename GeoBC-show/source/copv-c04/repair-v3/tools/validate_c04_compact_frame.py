"""Independent full-frame vs compact-frame FEniCSx matrix/stress equivalence."""
from pathlib import Path
import ast,json,time,sysconfig
for cfg in [sysconfig.get_config_vars()]:cfg.update(CFLAGS='-O2 -fPIC',LDFLAGS='',CC='/usr/bin/gcc',CXX='/usr/bin/g++',LDSHARED='/usr/bin/gcc -shared',CPPFLAGS='')
try:
 from setuptools._distutils import sysconfig as _dc
 _dc.get_config_vars().update(CFLAGS='-O2 -fPIC',LDFLAGS='',CC='/usr/bin/gcc',CXX='/usr/bin/g++',LDSHARED='/usr/bin/gcc -shared',CPPFLAGS='')
except ImportError:pass
import numpy as np
from mpi4py import MPI
import basix,basix.ufl,ufl
from dolfinx import mesh,fem
from dolfinx.fem import petsc as fp
W=Path(__file__).resolve().parents[1];O=W/'data/copv/c04_repair_v3/compact_frame_validation_v1';O.mkdir(exist_ok=False);start=time.perf_counter();source=W/'publication/GeoBC-Neural-Operator/tools/copv_fenicsx_c04_repair.py';tree=ast.parse(source.read_text());helper=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='compact_material_frame');scope=dict(ufl=ufl);exec(compile(ast.fix_missing_locations(ast.Module(body=[helper],type_ignores=[])),str(source),'exec'),scope)
corners=np.array([[10.,0.,0.],[12.,0.,0.],[10.,1.,0.],[10.,0.,1.]]);edges=[(0,1),(1,2),(2,0),(0,3),(2,3),(1,3)];xyz=np.vstack([corners,*[.5*(corners[a]+corners[b]) for a,b in edges]]);xyz[4,2]+=.03;xyz[7,1]+=.02;perm=[0,1,2,3]+[4+next(i for i,e in enumerate(edges) if set(e)==set(edge)) for edge in basix.topology(basix.CellType.tetrahedron)[1]];element=basix.ufl.element('Lagrange','tetrahedron',2,shape=(3,),lagrange_variant=basix.LagrangeVariant.equispaced);domain=mesh.create_mesh(MPI.COMM_SELF,np.array([perm],np.int64),ufl.Mesh(element),xyz);V=fem.functionspace(domain,element);field=fem.Function(V);field.interpolate(lambda x:np.vstack([x[0],2*x[1],3*x[2]]));u=ufl.TrialFunction(V);v=ufl.TestFunction(V);scalar=basix.create_element(basix.ElementFamily.P,basix.CellType.tetrahedron,2,basix.LagrangeVariant.equispaced);nodes=domain.geometry.x[domain.geometry.dofmap[0]];C=ufl.as_matrix(np.diag([11.,13.,17.,3.,5.,7.]));rows=[]
def eng(e):return ufl.as_vector([e[0,0],e[1,1],e[2,2],2*e[1,2],2*e[0,2],2*e[0,1]])
for degree in [4,6,8]:
 qp,qw=basix.make_quadrature(basix.CellType.tetrahedron,degree);tab=scalar.tabulate(0,qp);pos=np.einsum('qi,ij->qj',tab[0,:,:,0],nodes);phi=.2+pos[:,2]*.3;tv=np.c_[np.cos(phi),np.sin(phi)];r=np.linalg.norm(pos[:,:2],axis=-1);er=np.c_[pos[:,:2]/r[:,None],np.zeros(len(pos))];et=np.c_[-er[:,1],er[:,0],np.zeros(len(pos))];m=-tv[:,0,None]*er;m[:,2]-=tv[:,1];normal=np.cross(m,et)
 for mask in [0.,1.]:
  angle=np.deg2rad(25.);ca=np.cos(angle);sa=np.sin(angle);axes=np.stack([ca*m+sa*et,-sa*m+ca*et,normal],axis=-1);axes=mask*axes+(1-mask)*np.eye(3);Qspace=fem.functionspace(domain,basix.ufl.quadrature_element('tetrahedron',value_shape=(3,3),degree=degree,scheme='default'));Q=fem.Function(Qspace);Q.x.array.reshape(-1,9)[Qspace.dofmap.cell_dofs(0)]=axes.reshape(-1,9);Tspace=fem.functionspace(domain,basix.ufl.quadrature_element('tetrahedron',value_shape=(2,),degree=degree,scheme='default'));T=fem.Function(Tspace);T.x.array.reshape(-1,2)[Tspace.dofmap.cell_dofs(0)]=tv;Fspace=fem.functionspace(domain,basix.ufl.element('DG','tetrahedron',0,shape=(3,)));F=fem.Function(Fspace);F.x.array.reshape(-1,3)[Fspace.dofmap.cell_dofs(0)]=np.array([[mask,ca,sa]]);Qc=scope['compact_material_frame'](domain,T,F);dx=ufl.Measure('dx',domain=domain,metadata={'quadrature_degree':degree});matrices=[];stresses=[];packed=[]
  for frame in [Q,Qc]:
   def strain(w):return eng(frame.T*ufl.sym(ufl.grad(w))*frame)
   a=fem.form(ufl.dot(strain(v),C*strain(u))*dx);A=fp.assemble_matrix(a);A.assemble();i,j,values=A.getValuesCSR();dense=np.zeros((30,30))
   for row in range(30):dense[row,j[i[row]:i[row+1]]]=values[i[row]:i[row+1]]
   matrices.append(dense);stresses.append(fem.Expression(C*strain(field),qp).eval(domain,np.array([0],np.int32)));packed.append(sum(x.nbytes for x in fem.pack_coefficients(a).values()));A.destroy()
  error=float(np.linalg.norm(matrices[0]-matrices[1])/np.linalg.norm(matrices[0]));sigma=float(np.max(abs(stresses[0]-stresses[1])));assert error<1e-12 and sigma<1e-10
  rows.append(dict(assembly_degree=degree,points=len(qp),carbon_mask=int(mask),relative_matrix_difference=error,max_stress_difference=sigma,full_frame_packed_bytes=packed[0],compact_frame_packed_bytes=packed[1]))
report=dict(passed=True,scope='One curved TET10, independent numeric full 9-component frames vs generated compact UFL frame, spatially varying meridian direction, anisotropic stiffness and metal identity branch. All matrix entries and local stresses compared; not a COPV solve or convergence certification.',rows=rows,seconds=time.perf_counter()-start);(O/'review.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
