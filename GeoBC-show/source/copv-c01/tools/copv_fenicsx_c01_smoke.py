"""FEniCSx CUDA geometry-import smoke: axisymmetric TWO-layer isotropic elasticity.

This deliberately limited C01 integration check is not the C02 orthotropic
layup solver, a closed-end vessel validation, or a mesh-convergence claim.
Units mm/N/MPa. Bottom axial displacement fixed; top neck is traction-free.
"""
import argparse,os,time,json,sysconfig
from pathlib import Path
os.environ.setdefault('PETSC_OPTIONS','-use_gpu_aware_mpi 0')
for cfg in [sysconfig.get_config_vars()]:cfg.update(CFLAGS='-O2 -fPIC',LDFLAGS='',CC='/usr/bin/gcc',CXX='/usr/bin/g++',LDSHARED='/usr/bin/gcc -shared',CPPFLAGS='')
try:
    from setuptools._distutils import sysconfig as dcfg
    dcfg.get_config_vars().update(CFLAGS='-O2 -fPIC',LDFLAGS='',CC='/usr/bin/gcc',CXX='/usr/bin/g++',LDSHARED='/usr/bin/gcc -shared',CPPFLAGS='')
except ImportError:pass
import numpy as np
from scipy.spatial import cKDTree
from mpi4py import MPI
from petsc4py import PETSc
import basix,basix.ufl,ufl,dolfinx
from dolfinx import mesh,fem
from dolfinx.fem import petsc as fp

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();start=time.perf_counter()
    if not PETSc.Sys.hasExternalPackage('cuda'):raise RuntimeError('PETSc CUDA required; CPU fallback refused')
    src=np.load(args.input);xy=src['xy'];cells=src['cells'];ij=src['node_ij'];e=src['cell_ij']
    element=basix.ufl.element('Lagrange','quadrilateral',2,shape=(2,),lagrange_variant=basix.LagrangeVariant.equispaced)
    points=basix.create_element(basix.ElementFamily.P,basix.CellType.quadrilateral,2,basix.LagrangeVariant.equispaced).points
    _,perm=cKDTree((src['gmsh_reference']+1)/2).query(points)
    domain=mesh.create_mesh(MPI.COMM_SELF,cells[:,perm].astype(np.int64),ufl.Mesh(element),xy)
    V=fem.functionspace(domain,element);dist,node_order=cKDTree(xy).query(V.tabulate_dof_coordinates()[:,:2])
    if dist.max()>1e-8 or len(np.unique(node_order))!=len(xy):raise RuntimeError('Original-node bijection failed')
    DG=fem.functionspace(domain,('DG',0));E=fem.Function(DG);nu=fem.Function(DG)
    orig=domain.topology.original_cell_index
    for i,source in enumerate(orig):
        dof=DG.dofmap.cell_dofs(i);E.x.array[dof]=69000. if e[source,0]==0 else 120000.;nu.x.array[dof]=.33 if e[source,0]==0 else .28
    domain.topology.create_connectivity(1,2);boundary=mesh.exterior_facet_indices(domain.topology)
    geom=mesh.entities_to_geometry(domain,1,boundary,False)
    _,nodes=cKDTree(xy).query(domain.geometry.x[geom,:2]);param=ij[nodes]
    ids=[];values=[]
    for facet,pr in zip(boundary,param):
        value=1 if np.all(pr[:,0]==0) else 2 if np.all(pr[:,1]==0) else 3
        ids.append(facet);values.append(value)
    order=np.argsort(ids);tags=mesh.meshtags(domain,1,np.array(ids,np.int32)[order],np.array(values,np.int32)[order])
    ds=ufl.Measure('ds',domain=domain,subdomain_data=tags);dx=ufl.Measure('dx',domain=domain,metadata={'quadrature_degree':6})
    r=ufl.SpatialCoordinate(domain)[0];u=ufl.TrialFunction(V);v=ufl.TestFunction(V)
    def eps(w):return ufl.as_tensor([[w[0].dx(0),0,(w[0].dx(1)+w[1].dx(0))/2],[0,w[0]/r,0],[(w[0].dx(1)+w[1].dx(0))/2,0,w[1].dx(1)]])
    mu=E/(2*(1+nu));lam=E*nu/((1+nu)*(1-2*nu))
    a=fem.form(2*np.pi*r*ufl.inner(lam*ufl.tr(eps(u))*ufl.Identity(3)+2*mu*eps(u),eps(v))*dx)
    normal=ufl.FacetNormal(domain);L=fem.form(2*np.pi*r*ufl.dot(-normal,v)*ds(1))
    bottom=tags.find(2);dofs=fem.locate_dofs_topological(V.sub(1),1,bottom);bc=fem.dirichletbc(PETSc.ScalarType(0),dofs,V.sub(1))
    A=fp.assemble_matrix(a,bcs=[bc]);A.assemble();b=fp.assemble_vector(L);fp.apply_lifting(b,[a],bcs=[[bc]]);b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES,mode=PETSc.ScatterMode.REVERSE);fp.set_bc(b,[bc])
    cpu=A.createVecRight();direct=PETSc.KSP().create(MPI.COMM_SELF);direct.setOperators(A);direct.setType('preonly');direct.getPC().setType('lu');direct.solve(b,cpu)
    gpu=A.convert('aijcusparse');rhs=PETSc.Vec().createSeq(len(b.array),comm=MPI.COMM_SELF);rhs.setType('cuda');rhs.array[:]=b.array;answer=rhs.duplicate()
    ksp=PETSc.KSP().create(MPI.COMM_SELF);ksp.setOperators(gpu);ksp.setType('cg');ksp.getPC().setType('gamg');ksp.setTolerances(rtol=1e-11,max_it=5000);tick=time.perf_counter();ksp.solve(rhs,answer);solve_seconds=time.perf_counter()-tick
    if ksp.getConvergedReason()<=0:raise RuntimeError(f'CUDA solve failed {ksp.getConvergedReason()}')
    residual=rhs.duplicate();gpu.mult(answer,residual);energy=float(answer.dot(residual));work=float(answer.dot(rhs));residual.axpy(-1,rhs)
    relative_residual=residual.norm()/rhs.norm();diff=float(np.linalg.norm(answer.array-cpu.array)/np.linalg.norm(cpu.array))
    if gpu.getType()!='seqaijcusparse' or answer.getType()!='seqcuda' or relative_residual>1e-8 or diff>1e-7:raise RuntimeError((gpu.getType(),answer.getType(),relative_residual,diff))
    args.output.mkdir(parents=True,exist_ok=True);original=np.empty((len(xy),2));original[node_order]=answer.array.reshape(-1,2)
    np.savez_compressed(args.output/'fields.npz',xy=xy,U_rz_mm=original,cells=cells)
    report={'status':'passed','scope':'C01 geometry import and CUDA algebra check only; 2-layer isotropic, two-component axisym Q2, open necks, bottom uz=0, top traction-free, inner pressure 1 MPa; NOT C02 physical COPV/layup validation','dolfinx':dolfinx.__version__,'petsc':list(PETSc.Sys.getVersion()),'backend':'CPU FFCx assembly + PETSc CUDA CG/GAMG','matrix_type':gpu.getType(),'vector_type':answer.getType(),'nodes':len(xy),'cells':len(cells),'dofs':len(answer.array),'iterations':ksp.getIterationNumber(),'converged_reason':int(ksp.getConvergedReason()),'relative_residual':relative_residual,'CPU_GPU_U_relative_L2':diff,'energy_work_relative_error':abs(energy-work)/abs(work),'max_U_mm':float(np.linalg.norm(original,axis=1).max()),'inner_area_mm2':float(fem.assemble_scalar(fem.form(2*np.pi*r*ds(1)))),'original_node_mapping_max_error_mm':float(dist.max()),'cuda_solve_seconds':solve_seconds,'total_seconds':time.perf_counter()-start,'material_E_MPa':[69000,120000],'material_nu':[.33,.28]}
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':main()
