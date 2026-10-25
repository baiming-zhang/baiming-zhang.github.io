"""C02 verification solver: three-component axisymmetric orthotropic elasticity.

This geometry-matched code verification specimen excludes detailed bosses.
It retains twist, true Q2 mapping, GP material rotation and equivalent caps.
"""
import argparse,hashlib,json,os,sysconfig,time
from pathlib import Path
os.environ.setdefault('PETSC_OPTIONS','-use_gpu_aware_mpi 0')
_cfg=[sysconfig.get_config_vars()]
try:
    from setuptools._distutils import sysconfig as _dc
    _cfg.append(_dc.get_config_vars())
except ImportError:pass
for cfg in _cfg:cfg.update(CFLAGS='-O2 -fPIC',LDFLAGS='',CC='/usr/bin/gcc',CXX='/usr/bin/g++',LDSHARED='/usr/bin/gcc -shared',CPPFLAGS='')
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree
from mpi4py import MPI
from petsc4py import PETSc
import basix,basix.ufl,ufl,dolfinx
from dolfinx import mesh,fem
from dolfinx.fem import petsc as fp

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();start=time.perf_counter();out=args.output;out.mkdir(parents=True,exist_ok=True)
    if not PETSc.Sys.hasExternalPackage('cuda'):raise RuntimeError('CUDA required; CPU fallback refused')
    src=np.load(args.input);xy=src['xy'];cells=src['cells'];ij=src['node_ij'];ci=src['cell_ij'];matrices=src['C_material'];params=src['parameters'];bore=params[3]
    ge=basix.ufl.element('Lagrange','quadrilateral',2,shape=(2,),lagrange_variant=basix.LagrangeVariant.equispaced)
    scalar=basix.create_element(basix.ElementFamily.P,basix.CellType.quadrilateral,2,basix.LagrangeVariant.equispaced)
    _,perm=cKDTree((src['gmsh_reference']+1)/2).query(scalar.points)
    domain=mesh.create_mesh(MPI.COMM_SELF,cells[:,perm].astype(np.int64),ufl.Mesh(ge),xy)
    V=fem.functionspace(domain,basix.ufl.element('Lagrange','quadrilateral',2,shape=(3,),lagrange_variant=basix.LagrangeVariant.equispaced))
    dist,node_order=cKDTree(xy).query(V.tabulate_dof_coordinates()[:,:2]);assert dist.max()<1e-8 and len(np.unique(node_order))==len(xy)
    cell_order=domain.topology.original_cell_index;qp,qw=basix.make_quadrature(basix.CellType.quadrilateral,5);tab=scalar.tabulate(1,qp)
    coordinates=domain.geometry.x[domain.geometry.dofmap,:2];ds=np.einsum('qi,eij->eqj',tab[2,:,:,0],coordinates);dt=np.einsum('qi,eij->eqj',tab[1,:,:,0],coordinates)
    jac=dt[...,0]*ds[...,1]-ds[...,0]*dt[...,1]
    if jac.min()<=0:raise RuntimeError('Nonpositive geometry Jacobian')
    tangent=ds/np.linalg.norm(ds,axis=-1)[...,None];meridian=-np.stack([tangent[...,0],np.zeros(tangent.shape[:-1]),tangent[...,1]],axis=-1);circ=np.zeros_like(meridian);circ[...,1]=1;normal=np.cross(meridian,circ)
    ply=ci[cell_order,0];angle=np.zeros(len(ply));angle[ply>0]=np.deg2rad(np.array([55,-55,-55,55]))[ply[ply>0]-1]
    c=np.cos(angle)[:,None,None];s=np.sin(angle)[:,None,None];axes=np.stack([c*meridian+s*circ,-s*meridian+c*circ,normal],axis=-1)
    if abs(np.linalg.det(axes)-1).max()>1e-10:raise RuntimeError('Invalid material frame')
    Qspace=fem.functionspace(domain,basix.ufl.quadrature_element('quadrilateral',value_shape=(3,3),degree=5,scheme='default'));Q=fem.Function(Qspace)
    for e in range(len(cell_order)):Q.x.array.reshape(-1,9)[Qspace.dofmap.cell_dofs(e)]=axes[e].reshape(-1,9)
    Cspace=fem.functionspace(domain,basix.ufl.element('DG','quadrilateral',0,shape=(6,6)));C=fem.Function(Cspace)
    for e in range(len(cell_order)):C.x.array.reshape(-1,36)[Cspace.dofmap.cell_dofs(e)]=matrices[int(ply[e]>0)].reshape(1,36)
    domain.topology.create_connectivity(1,2);boundary=mesh.exterior_facet_indices(domain.topology);geom=mesh.entities_to_geometry(domain,1,boundary,False);_,nodes=cKDTree(xy).query(domain.geometry.x[geom,:2]);values=[]
    for a in ij[nodes]:values.append(1 if np.all(a[:,0]==0) else 2 if np.all(a[:,1]==0) else 3 if np.all(a[:,1]==ij[:,1].max()) else 4)
    order=np.argsort(boundary);tags=mesh.meshtags(domain,1,boundary[order].astype(np.int32),np.array(values,np.int32)[order]);dx=ufl.Measure('dx',domain=domain,metadata={'quadrature_degree':5});dsurf=ufl.Measure('ds',domain=domain,subdomain_data=tags,metadata={'quadrature_degree':5});r=ufl.SpatialCoordinate(domain)[0]
    def eps(u):
        g=ufl.as_tensor([[u[0].dx(0),-u[1]/r,u[0].dx(1)],[u[1].dx(0),u[0]/r,u[1].dx(1)],[u[2].dx(0),0,u[2].dx(1)]]);return ufl.sym(g)
    def eng(e):return ufl.as_vector([e[0,0],e[1,1],e[2,2],2*e[1,2],2*e[0,2],2*e[0,1]])
    def localstrain(u):return eng(Q.T*eps(u)*Q)
    u=ufl.TrialFunction(V);v=ufl.TestFunction(V);a=fem.form(2*np.pi*r*ufl.dot(localstrain(v),C*localstrain(u))*dx)
    cap=np.pi*bore**2;areas={name:float(fem.assemble_scalar(fem.form(2*np.pi*r*dsurf(i)))) for i,name in [(2,'bottom'),(3,'top')]};n=ufl.FacetNormal(domain)
    L=fem.form(2*np.pi*r*(-n[0]*v[0]-n[1]*v[2])*dsurf(1)-2*np.pi*r*(cap/areas['bottom'])*v[2]*dsurf(2)+2*np.pi*r*(cap/areas['top'])*v[2]*dsurf(3))
    A=fp.assemble_matrix(a);A.assemble();rhs=fp.assemble_vector(L);rhs.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES,mode=PETSc.ScatterMode.REVERSE)
    K=sp.csr_matrix(A.getValuesCSR()[::-1],shape=A.getSize());F=rhs.array.copy();x=V.tabulate_dof_coordinates()[:,:2];bottom=np.isclose(x[:,1],x[:,1].min());top=np.isclose(x[:,1],x[:,1].max());gauge=np.flatnonzero(bottom)[np.argmax(x[bottom,0])]
    fixed=np.r_[3*np.flatnonzero(bottom)+2,3*gauge+1];topdofs=3*np.flatnonzero(top)+2;mask=np.ones(len(F),bool);mask[fixed]=False;mask[topdofs]=False;free=np.flatnonzero(mask);rows=np.r_[free,topdofs];cols=np.r_[np.arange(len(free)),np.full(len(topdofs),len(free))];P=sp.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(F),len(free)+1));Kr=(P.T@K@P).tocsr();Fr=P.T@F
    host=PETSc.Mat().createAIJ(size=Kr.shape,csr=(Kr.indptr.astype(PETSc.IntType),Kr.indices.astype(PETSc.IntType),Kr.data),comm=MPI.COMM_SELF);gpu=host.convert('aijcusparse');b=PETSc.Vec().createSeq(len(Fr),comm=MPI.COMM_SELF);b.setType('cuda');answer=b.duplicate();answer.set(0)
    direct=PETSc.KSP().create(MPI.COMM_SELF);direct.setOperators(host);direct.setType('preonly');direct.getPC().setType('lu');bhost=host.createVecRight();bhost.array[:]=Fr;reference=bhost.duplicate();direct.solve(bhost,reference)
    ksp=PETSc.KSP().create(MPI.COMM_SELF);ksp.setOperators(gpu);ksp.setType('cg');ksp.getPC().setType('gamg');ksp.setTolerances(rtol=1e-11,max_it=20000);history=[];peak=None;previous=np.zeros(len(Fr));factors=np.r_[np.arange(1,11)/10,np.arange(9,-1,-1)/10]
    for i,factor in enumerate(factors):
        b.array[:]=factor*Fr;answer.array[:]=previous if factor else 0.;ksp.setInitialGuessNonzero(bool(factor and np.linalg.norm(previous)>0));tick=time.perf_counter();ksp.solve(b,answer)
        if ksp.getConvergedReason()<=0:raise RuntimeError(f'CUDA solve failed state {i}: {ksp.getConvergedReason()}')
        solved=answer.array.copy();previous=solved;U=P@solved;residual=Kr@solved-factor*Fr;den=max(np.linalg.norm(factor*Fr),1e-30);rel=np.linalg.norm(residual)/den if factor else np.linalg.norm(U)
        if factor and rel>1e-8:raise RuntimeError('True residual gate failed')
        history.append({'state':i+1,'pressure_MPa':float(factor),'iterations':ksp.getIterationNumber(),'reason':int(ksp.getConvergedReason()),'seconds':time.perf_counter()-tick,'relative_residual_or_zero_U_norm':float(rel),'max_U_mm':float(np.linalg.norm(U.reshape(-1,3),axis=1).max())})
        if factor==1.:peak=U.copy();peak_reduced=solved.copy()
    field=fem.Function(V);field.x.array[:]=peak;expr=C*localstrain(field)
    localstress=fem.Expression(expr,qp).eval(domain,np.arange(len(cell_order),dtype=np.int32)).reshape(len(cell_order),len(qp),6)
    stress_tensor=np.zeros((len(cell_order),len(qp),3,3));pairs=[(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)]
    for j,(a0,a1) in enumerate(pairs):stress_tensor[...,a0,a1]=localstress[...,j];stress_tensor[...,a1,a0]=localstress[...,j]
    globalstress=axes@stress_tensor@np.swapaxes(axes,-1,-2);stress=np.stack([globalstress[...,i,j] for i,j in pairs],axis=-1);qxy=np.einsum('qi,eij->eqj',tab[0,:,:,0],coordinates);weight=jac*qw[None,:]*2*np.pi*qxy[...,0]
    original=np.empty((len(xy),3));original[node_order]=peak.reshape(-1,3);ordered_stress=np.empty_like(stress);ordered_qxy=np.empty_like(qxy);ordered_weights=np.empty_like(weight);ordered_stress[cell_order]=stress;ordered_qxy[cell_order]=qxy;ordered_weights[cell_order]=weight
    R=K@peak-F;energy=.5*peak@(K@peak);work=peak@F;diff=np.linalg.norm(peak_reduced-reference.array)/np.linalg.norm(reference.array)
    np.savez_compressed(out/'fields.npz',xyz_rz=xy,displacement_rtz=original,stress_rtz_voigt=ordered_stress,quadrature_rz=ordered_qxy,quadrature_volume_weight=ordered_weights,cells=cells)
    report={'status':'completed_verification_specimen','scope':'Closed-end orthotropic 3-component axisymmetric geometry-matched code verification; detailed boss/seal assembly is separate','dolfinx':dolfinx.__version__,'petsc':list(PETSc.Sys.getVersion()),'backend':'CPU FFCx assembly + PETSc CUDA CG/GAMG','matrix_type':gpu.getType(),'vector_type':answer.getType(),'nodes':len(xy),'cells':len(cells),'states':history,'CPU_GPU_U_relative_L2':float(diff),'strain_energy_Nmm':float(energy),'external_work_Nmm':float(work),'energy_work_relative_error':float(abs(2*energy/work-1)),'top_master_extension_mm':float(peak_reduced[-1]),'max_abs_utheta_mm':float(abs(original[:,1]).max()),'min_mapping_Jacobian':float(jac.min()),'material_frame_det_max_error':float(abs(np.linalg.det(axes)-1).max()),'cap_force_per_end_N':float(cap),'cap_area_mm2':areas,'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),'total_seconds':time.perf_counter()-start}
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':main()
