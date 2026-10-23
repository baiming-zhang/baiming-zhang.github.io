"""Independent curved P2 FE solver for the elastic branch of CCX finite-strain J2.

This is deliberately fail-closed on yielding; it does not substitute a different
plasticity law. The case contract and all accepted states record that limitation.

The CUDA backend is hybrid in the precise FEniCSx/PETSc sense: FFCx element
integration and sparse assembly run on the host, while the reduced saddle-point
linear system, Krylov iterations, GAMG preconditioner and sparse products run on
the NVIDIA GPU through PETSc ``AIJCUSPARSE``/``CUDA`` objects.  CUDA is mandatory
when requested; there is no silent CPU fallback.
"""
import argparse,hashlib,json,os,time
from pathlib import Path
import sysconfig
# This solver is intentionally single-rank (COMM_SELF). The conda MPICH build
# is not CUDA-aware, so prevent PETSc from aborting during device discovery;
# no inter-rank device transfer exists in this execution mode.
os.environ.setdefault('PETSC_OPTIONS','-use_gpu_aware_mpi 0')
# Conda embeds an unquoted '-isystem PREFIX/include' in Python's compiler flags.
# A Windows-mounted prefix may contain spaces. Keep includes as proper argv items
# supplied by CFFI and use the system compiler; no relocation or root writes.
configs=[sysconfig.get_config_vars()]
try:
    from setuptools._distutils import sysconfig as dist_config
    configs.append(dist_config.get_config_vars())
except ImportError:pass
for config in configs:
    config.update(CFLAGS='-O2 -fPIC',LDFLAGS='',CC='/usr/bin/gcc',CXX='/usr/bin/g++',LDSHARED='/usr/bin/gcc -shared',CPPFLAGS='')
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree
from mpi4py import MPI
from petsc4py import PETSc
import basix,basix.ufl,ufl,dolfinx
from dolfinx import fem,mesh
from dolfinx.fem import petsc as fp

def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2),encoding='utf-8');tmp.replace(path)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def constitutive(u):
    F=ufl.variable(ufl.Identity(3)+ufl.grad(u));J=ufl.det(F);C=F.T*F
    E=115100.;nu=.28;mu=E/(2*(1+nu));K=E/(3*(1-2*nu))
    W=mu/2*(J**(-2/3)*ufl.tr(C)-3)+K/4*(J*J-1-2*ufl.ln(J))
    P=ufl.diff(W,F);sigma=P*F.T/J
    return W,sigma,J
def make_domain(xyz,cells):
    # Native Gmsh TET10 edges: 01,12,20,03,23,13.
    # Basix tetrahedron entities use 23,13,12,03,02,01.
    edges=[(0,1),(1,2),(2,0),(0,3),(2,3),(1,3)]
    topology=basix.topology(basix.CellType.tetrahedron)
    perm=[0,1,2,3]+[4+next(i for i,e in enumerate(edges) if set(e)==set(edge)) for edge in topology[1]]
    coordinate=basix.ufl.element('Lagrange','tetrahedron',2,shape=(3,),lagrange_variant=basix.LagrangeVariant.equispaced)
    domain=mesh.create_mesh(MPI.COMM_SELF,np.asarray(cells[:,perm],dtype=np.int64),ufl.Mesh(coordinate),np.asarray(xyz,dtype=np.float64))
    V=fem.functionspace(domain,coordinate)
    coords=V.tabulate_dof_coordinates();dist,idx=cKDTree(xyz).query(coords)
    assert len(coords)==len(xyz) and len(np.unique(idx))==len(xyz) and dist.max()<1e-8
    return domain,V,idx,float(dist.max())
def reduction(normal,node_order):
    """Sparse exact radial elimination; mean rows use two Lagrange multipliers."""
    n=normal[node_order];rows=[];cols=[];vals=[];k=0
    for i,v in enumerate(n):
        if np.linalg.norm(v)>.5:
            p=int(np.argmax(np.abs(v[:2])));q=1-p;ratio=0. if abs(v[q])<1e-9 else -v[q]/v[p]
            rows.extend([3*i+q,3*i+p,3*i+2]);cols.extend([k,k,k+1]);vals.extend([1.,ratio,1.]);k+=2
        else:
            rows.extend([3*i,3*i+1,3*i+2]);cols.extend([k,k+1,k+2]);vals.extend([1.,1.,1.]);k+=3
    P=sp.csr_matrix((vals,(rows,cols)),shape=(len(n)*3,k));P.eliminate_zeros()
    return P
def rigid_near_modes(P,xyz):
    """Project the six 3-D rigid modes into the radially reduced coordinates."""
    centered=np.asarray(xyz)-np.asarray(xyz).mean(axis=0);fields=[]
    for axis in np.eye(3):fields.append(np.broadcast_to(axis,centered.shape).copy())
    for axis in np.eye(3):fields.append(np.cross(np.broadcast_to(axis,centered.shape),centered))
    gram=np.asarray(P.power(2).sum(axis=0)).ravel();assert np.all(gram>0)
    candidates=np.column_stack([(P.T@field.ravel())/gram for field in fields])
    modes,_=np.linalg.qr(candidates,mode='reduced')
    assert modes.shape[1]==6 and np.max(np.abs(modes.T@modes-np.eye(6)))<1e-10
    return modes
def petsc_csr(A,backend='cpu'):
    result=PETSc.Mat().createAIJ(size=A.shape,csr=(A.indptr.astype(PETSc.IntType),A.indices.astype(PETSc.IntType),A.data),comm=MPI.COMM_SELF)
    if backend=='cuda':
        device=result.convert(PETSc.Mat.Type.AIJCUSPARSE);result.destroy();result=device
    return result
def require_cuda_backend():
    """Prove that this PETSc build can execute CUDA work, not merely name types."""
    if not PETSc.Sys.hasExternalPackage('cuda'):
        raise RuntimeError('CUDA backend requested but PETSc was built without CUDA; refusing CPU fallback')
    vector=PETSc.Vec().createSeq(2,comm=MPI.COMM_SELF);matrix=None;device_matrix=None;answer=None
    try:
        vector.setType(PETSc.Vec.Type.CUDA);vector.set(2.);norm=float(vector.norm())
        matrix=PETSc.Mat().createAIJ(size=(2,2),csr=(np.array([0,1,2],PETSc.IntType),np.array([0,1],PETSc.IntType),np.ones(2)),comm=MPI.COMM_SELF)
        device_matrix=matrix.convert(PETSc.Mat.Type.AIJCUSPARSE);answer=vector.duplicate();device_matrix.mult(vector,answer)
        if device_matrix.getType()!='seqaijcusparse' or vector.getType()!='seqcuda' or abs(float(answer.norm())-norm)>1e-12:
            raise RuntimeError('PETSc CUDA execution probe returned unexpected types or values')
        return dict(petsc_cuda_package=True,probe_vector_type=vector.getType(),probe_matrix_type=device_matrix.getType(),probe_norm=norm)
    finally:
        for item in [answer,device_matrix,matrix,vector]:
            if item is not None:item.destroy()
def solve_direct(reduced,Cmat,Ct,zero,residual,g):
    nest=PETSc.Mat().createNest([[reduced,Ct],[Cmat,zero]],comm=MPI.COMM_SELF);nest.assemble();aug=nest.convert('aij')
    ksp=PETSc.KSP().create(MPI.COMM_SELF);ksp.setOperators(aug);ksp.setType('preonly');pc=ksp.getPC();pc.setType('lu');pc.setFactorSolverType('mumps');ksp.setFromOptions()
    rhs=aug.createVecRight();ans=rhs.duplicate();rhs.array[:len(residual)]=-residual;rhs.array[len(residual):]=-g
    ksp.solve(rhs,ans);reason=int(ksp.getConvergedReason());dq=ans.array[:len(residual)].copy()
    report=dict(ksp_type=ksp.getType(),pc_type=pc.getType(),reason=reason,iterations=int(ksp.getIterationNumber()),matrix_type=aug.getType(),vector_type=rhs.getType())
    for item in [ksp,rhs,ans,aug,nest]:item.destroy()
    if reason<0:raise RuntimeError('PETSc MUMPS factorization/solve failed')
    return dq,report
def solve_iterative(reduced,Cmat,Ct,zero,residual,g,backend,near_modes,linear_atol):
    """Exact two-constraint solve using a low-rank augmented GPU operator.

    For ``K dq + C.T dl = b, C dq = h`` use ``A=K+gamma*C.T*C``.
    Three solves with nonsingular ``A`` give the exact bordered solution via a
    2x2 Schur complement. ``C.T*C`` is applied matrix-free, so no dense global
    matrix is formed and every large operation remains a CUDA sparse kernel.
    """
    matrix=reduced;constraint=Cmat;constraint_t=Ct
    if backend=='cuda':
        matrix=reduced.convert(PETSc.Mat.Type.AIJCUSPARSE)
        constraint=Cmat.convert(PETSc.Mat.Type.AIJCUSPARSE);constraint_t=Ct.convert(PETSc.Mat.Type.AIJCUSPARSE)
    k_scale=max(float(matrix.norm(PETSc.NormType.INFINITY)),1.)
    c_scale=max(float(constraint.norm(PETSc.NormType.FROBENIUS)),1.)
    gamma=k_scale/(c_scale*c_scale)
    class AugmentedOperator:
        def __init__(self,K,C,Ct,weight):
            self.K=K;self.C=C;self.Ct=Ct;self.weight=weight
            self.small=C.createVecLeft();self.large=K.createVecLeft()
        def mult(self,mat,x,y):
            self.K.mult(x,y);self.C.mult(x,self.small);self.small.scale(self.weight)
            self.Ct.mult(self.small,self.large);y.axpy(1.,self.large)
        def destroy(self,mat):self.small.destroy();self.large.destroy()
    context=AugmentedOperator(matrix,constraint,constraint_t,gamma)
    operator=PETSc.Mat().createPython(matrix.getSize(),context=context,comm=MPI.COMM_SELF)
    operator.setVecType(PETSc.Vec.Type.CUDA if backend=='cuda' else PETSc.Vec.Type.SEQ);operator.setUp()
    operator.setOption(PETSc.Mat.Option.SYMMETRIC,True);operator.setOption(PETSc.Mat.Option.SPD,True)
    preconditioner=matrix.copy();preconditioner_shift=max(k_scale*1e-10,1e-12);preconditioner.shift(preconditioner_shift)
    augmentation_diagonal=np.zeros(len(residual))
    for row in range(2):
        columns,values=Cmat.getRow(row);augmentation_diagonal[columns]+=gamma*np.square(values)
    diagonal=matrix.createVecLeft();nz=np.flatnonzero(augmentation_diagonal)
    diagonal.setValues(nz.astype(PETSc.IntType),augmentation_diagonal[nz]);diagonal.assemble()
    preconditioner.setDiagonal(diagonal,addv=PETSc.InsertMode.ADD_VALUES);preconditioner.assemble()
    preconditioner.setOption(PETSc.Mat.Option.SYMMETRIC,True);preconditioner.setOption(PETSc.Mat.Option.SPD,True)
    mode_vectors=[]
    for values in near_modes.T:
        vector=matrix.createVecRight();vector.setValues(np.arange(len(values),dtype=PETSc.IntType),values);vector.assemble();mode_vectors.append(vector)
    near_nullspace=PETSc.NullSpace().create(constant=False,vectors=mode_vectors,comm=MPI.COMM_SELF)
    preconditioner.setNearNullSpace(near_nullspace)
    prefix=f'rod_{backend}_';options=PETSc.Options()
    configured={
        'ksp_type':'cg','ksp_rtol':'1e-8','ksp_atol':'1e-10','ksp_max_it':'1200','ksp_norm_type':'unpreconditioned',
        'pc_type':'gamg','pc_gamg_type':'agg','pc_gamg_threshold':'0.08','pc_gamg_agg_nsmooths':'0',
        'mg_levels_ksp_type':'chebyshev','mg_levels_pc_type':'jacobi','mg_coarse_ksp_type':'preonly','mg_coarse_pc_type':'lu'}
    for key,value in configured.items():options[prefix+key]=value
    ksp=PETSc.KSP().create(MPI.COMM_SELF);ksp.setOptionsPrefix(prefix);ksp.setOperators(operator,preconditioner);ksp.setFromOptions()
    rhs=matrix.createVecRight();rhs.setValues(np.arange(len(residual),dtype=PETSc.IntType),-residual);rhs.assemble()
    small_rhs=constraint.createVecLeft();basis=small_rhs.duplicate();column_rhs=rhs.duplicate()
    solutions=[];reasons=[];iterations=[];initial_residuals=[];final_residuals=[]
    for index in [-1,0,1]:
        current_rhs=rhs
        if index>=0:
            basis.set(0.);basis.setValue(index,1.);basis.assemble();constraint_t.mult(basis,column_rhs);current_rhs=column_rhs
            # Gauge columns have O(E) norms and stagnate at an absolute O(1e-2)
            # algebraic residual on the large mesh. This is below 2e-7 relative
            # to their scale; the physical Newton RHS keeps the strict tolerance.
            ksp.setTolerances(rtol=1e-8,atol=max(2e-7*c_scale,1e-10),max_it=1200)
        else:ksp.setTolerances(rtol=1e-8,atol=linear_atol,max_it=1200)
        answer=current_rhs.duplicate();answer.set(0.);ksp.setConvergenceHistory(reset=True);ksp.solve(current_rhs,answer)
        convergence=np.asarray(ksp.getConvergenceHistory(),dtype=float);reason=int(ksp.getConvergedReason())
        reasons.append(reason);iterations.append(int(ksp.getIterationNumber()));initial_residuals.append(float(convergence[0]) if len(convergence) else None);final_residuals.append(float(convergence[-1]) if len(convergence) else None)
        if reason<0:
            final=float(convergence[-1]) if len(convergence) else None
            target=max(2e-7*c_scale,1e-10) if index>=0 else linear_atol
            answer.destroy();raise RuntimeError(f'PETSc {backend} augmented solve {index} failed with reason {reason}, iterations={ksp.getIterationNumber()}, final_residual={final}, absolute_target={target}')
        solutions.append(answer)
    constraint.mult(solutions[0],small_rhs);cb=small_rhs.getArray(readonly=True).copy();schur=np.empty((2,2))
    for j in range(2):constraint.mult(solutions[j+1],small_rhs);schur[:,j]=small_rhs.getArray(readonly=True)
    h=-g;multiplier=np.linalg.solve(schur,cb-h);answer=solutions[0].copy()
    for j in range(2):answer.axpy(-multiplier[j],solutions[j+1])
    constraint.mult(answer,small_rhs);constraint_error=float(np.max(np.abs(small_rhs.getArray(readonly=True)-h)))
    if constraint_error>1e-7:raise RuntimeError(f'Bordered CUDA solve constraint error {constraint_error}')
    dq=answer.getArray(readonly=True).copy()
    report=dict(ksp_type=ksp.getType(),pc_type=ksp.getPC().getType(),reason=min(reasons),iterations=sum(iterations),subsolve_iterations=iterations,
        matrix_type=matrix.getType(),vector_type=rhs.getType(),initial_residual=max(x for x in initial_residuals if x is not None),
        final_residual=max(x for x in final_residuals if x is not None),preconditioner_shift=preconditioner_shift,
        augmentation_gamma=gamma,schur_condition=float(np.linalg.cond(schur)),constraint_error=constraint_error,near_nullspace_modes=len(mode_vectors))
    for key in configured:
        try:del options[prefix+key]
        except KeyError:pass
    for item in [answer,*solutions,ksp,rhs,small_rhs,basis,column_rhs,near_nullspace,*mode_vectors,diagonal,preconditioner,operator]:item.destroy()
    if matrix is not reduced:matrix.destroy();constraint.destroy();constraint_t.destroy()
    if backend=='cuda' and (report['matrix_type']!='seqaijcusparse' or report['vector_type']!='seqcuda'):
        raise RuntimeError(f'CUDA solve unexpectedly used {report["matrix_type"]}/{report["vector_type"]}; refusing result')
    return dq,report
def build_forms(domain,V,u):
    W,sigma,J=constitutive(u);v=ufl.TestFunction(V);du=ufl.TrialFunction(V)
    # degree-2 tetrahedron quadrature is the same symmetric 4-point rule as C3D10.
    points,weights=basix.make_quadrature(basix.CellType.tetrahedron,2)
    assert len(points)==4 and np.max(abs(weights-1/24))<1e-14
    dx=ufl.Measure('dx',domain=domain,metadata={'quadrature_degree':2})
    R=ufl.derivative(W*dx,u,v);K=ufl.derivative(R,u,du)
    opts={'cache_dir':os.environ['XDG_CACHE_HOME']+'/ffcx_safe_flags_v2','timeout':180}
    return fem.form(R,jit_options=opts),fem.form(K,jit_options=opts),sigma,J,points
def evaluate_gp(domain,V,u,sigma,J,points,node_order,inputs):
    count=domain.topology.index_map(3).size_local;entities=np.arange(count,dtype=np.int32)
    sigma_value=fem.Expression(sigma,points).eval(domain,entities).reshape(count,4,3,3)
    jac=fem.Expression(J,points).eval(domain,entities).reshape(count,4)
    dev=sigma_value-np.trace(sigma_value,axis1=2,axis2=3)[:,:,None,None]*np.eye(3)/3
    vm=np.sqrt(1.5*np.sum(dev*dev,axis=(2,3)))
    # Recover native cell ordering from the exact four corner IDs, not centroid proximity.
    native={tuple(sorted(map(int,c[:4]))):i for i,c in enumerate(inputs['cells'])}
    domain.topology.create_connectivity(3,0)
    native_cells=[]
    for i in entities:
        dofs=V.dofmap.cell_dofs(i);key=tuple(sorted(map(int,node_order[dofs[:4]])));native_cells.append(native[key])
    native_cells=np.asarray(native_cells)
    physical_gp=fem.Expression(ufl.SpatialCoordinate(domain),points).eval(domain,entities).reshape(count,4,3)
    return sigma_value,vm,jac,native_cells,physical_gp
def solve(inputs,output,benchmark=False,backend='cpu',linear_solver='auto'):
    start=time.time();output.mkdir(parents=True,exist_ok=True)
    def status(phase,**kw):write(output/'status.json',dict(phase=phase,pid=os.getpid(),elapsed_seconds=time.time()-start,**kw))
    if backend not in {'cpu','cuda'}:raise ValueError(f'Unsupported backend {backend!r}')
    if linear_solver=='auto':linear_solver='iterative' if backend=='cuda' else 'direct'
    if backend=='cuda' and linear_solver!='iterative':raise ValueError('CUDA backend requires the scalable iterative solver')
    cuda_review=require_cuda_backend() if backend=='cuda' else dict(petsc_cuda_package=bool(PETSc.Sys.hasExternalPackage('cuda')))
    status('creating_curved_P2_domain');domain,V,order,err=make_domain(inputs['xyz'],inputs['cells'])
    assert V.dofmap.index_map_bs==3
    u=fem.Function(V);Psc=reduction(inputs['radial_normal'],order);near_modes=rigid_near_modes(Psc,inputs['xyz'][order]);P=petsc_csr(Psc)
    force=inputs['force_N'][order].ravel();fr=Psc.T@force
    Cphysical=np.vstack([(Psc.T@inputs[k][order].ravel()) for k in ['mean_tangent','mean_axial']])
    # Exact row scaling changes only multiplier units, not admissible displacement.
    # Keep the global integral rows comparable to the stiffness in the saddle LU.
    row_scale=115100./np.linalg.norm(Cphysical,axis=1);C=Cphysical*row_scale[:,None]
    assert np.linalg.matrix_rank(C@C.T)==2
    # Constraint blocks remain on the host only until the monolithic matrix is
    # formed; the complete matrix is then transferred to cuSPARSE in one copy.
    Cmat=petsc_csr(sp.csr_matrix(C));Ct=petsc_csr(sp.csr_matrix(C.T));zero=petsc_csr(sp.csr_matrix((2,2)))
    q=np.zeros(Psc.shape[1]);Rform,Kform,sigma,J,points=build_forms(domain,V,u)
    write(output/'mapping_review.json',dict(dolfinx_version=dolfinx.__version__,petsc_version='.'.join(map(str,PETSc.Sys.getVersion())),nodes=len(order),native_mapping_max_error_mm=err,geometry_degree=2,space_degree=2,quadrature_points=points.tolist(),radial_reduced_dofs=len(q),integral_rows=2,backend=backend,linear_solver=linear_solver,cuda_review=cuda_review,assembly_device='cpu',linear_algebra_device='cuda' if backend=='cuda' else 'cpu'))
    A=fp.create_matrix(Kform);r=fp.create_vector(V)
    history=[];linear_history=[];peak=None;peak_stress=None;peak_vm=None;peak_cells=None;peak_gp=None;peak_internal=None
    factors=[.1] if benchmark else list(np.arange(1,11)/10)+list(np.arange(9,-1,-1)/10)
    for step,factor in enumerate(factors):
        for iteration in range(35):
            u.x.array[:]=Psc@q;u.x.scatter_forward()
            with r.localForm() as local:local.set(0.)
            fp.assemble_vector(r,Rform);r.ghostUpdate(addv=PETSc.InsertMode.ADD,mode=PETSc.ScatterMode.REVERSE)
            residual=Psc.T@(r.array-force*factor);g=C@q
            lam=np.linalg.solve(C@C.T,C@residual)
            free=residual-C.T@lam;rel=float(np.linalg.norm(free)/max(np.linalg.norm(fr)*max(factor,.1),1))
            status('newton',step=step,load_factor=float(factor),iteration=iteration,free_residual_relative=rel)
            if rel<2e-7 and np.max(abs(g/row_scale))<1e-8:break
            A.zeroEntries();fp.assemble_matrix(A,Kform);A.assemble()
            reduced=A.PtAP(P)
            # Radial restraints alone retain two rigid modes. Always solve the
            # exact saddle system; CUDA uses Schur/GAMG rather than a huge LU.
            if linear_solver=='direct':dq,linear=solve_direct(reduced,Cmat,Ct,zero,residual,g)
            else:dq,linear=solve_iterative(reduced,Cmat,Ct,zero,residual,g,backend,near_modes,max(np.linalg.norm(fr)*max(factor,.1),1.)*2e-7)
            linear.update(step=step,newton_iteration=iteration);linear_history.append(linear);write(output/'linear_solver_history.json',linear_history)
            # Load is bounded and smooth for this accepted elastic case. Refuse an unbounded Newton step.
            if not np.isfinite(dq).all() or np.max(abs(dq))>100:raise RuntimeError('Newton step not admissible; no fallback physics')
            q+=dq;reduced.destroy()
        else:raise RuntimeError('Newton iteration limit; no unconverged label exported')
        u.x.array[:]=Psc@q;u.x.scatter_forward()
        S,vm,jac,native_cells,physical_gp=evaluate_gp(domain,V,u,sigma,J,points,order,inputs)
        if jac.min()<=0:raise RuntimeError('Deformation gradient inverted')
        if vm.max()>=781.6:raise RuntimeError('Yield detected: full J2 return map required; elastic-branch comparison invalid')
        native_u=np.empty_like(inputs['xyz']);native_u[order]=u.x.array.reshape(-1,3)
        record=dict(step=step,load_factor=float(factor),iterations=iteration,min_J=float(jac.min()),max_von_Mises_MPa=float(vm.max()),PEEQ=0.,free_residual_relative=rel,mean_constraint_mm=(g/row_scale).tolist(),elapsed_seconds=time.time()-start)
        history.append(record);write(output/'accepted_history.json',history)
        np.savez_compressed(output/f'state_{step:02d}.npz',node_ids=inputs['node_ids'],U_mm=native_u,load_factor=factor)
        if factor==1.:
            peak=native_u.copy();peak_stress=S.copy();peak_vm=vm.copy();peak_cells=native_cells.copy();peak_gp=physical_gp.copy()
            peak_internal=np.empty_like(inputs['xyz']);peak_internal[order]=r.array.reshape(-1,3)
    if benchmark:
        peak=native_u.copy();peak_stress=S;peak_vm=vm;peak_cells=native_cells;peak_gp=physical_gp
        peak_internal=np.empty_like(inputs['xyz']);peak_internal[order]=r.array.reshape(-1,3)
    stress=np.empty_like(peak_stress);stress[peak_cells]=peak_stress
    vm_native=np.empty_like(peak_vm);vm_native[peak_cells]=peak_vm
    gp=np.empty_like(peak_gp);gp[peak_cells]=peak_gp
    np.savez_compressed(output/'fields.npz',node_ids=inputs['node_ids'],cell_ids=inputs['cell_ids'],U_peak_mm=peak,U_unload_mm=native_u,force_peak_N=inputs['force_N'],internal_force_peak_N=peak_internal,GP_xyz_mm=gp,GP_Cauchy_stress_MPa=stress,GP_von_Mises_MPa=vm_native,GP_PEEQ=np.zeros_like(vm_native))
    write(output/'summary.json',dict(seconds=time.time()-start,nodes=len(order),elements=len(inputs['cells']),max_U_mm=float(np.linalg.norm(peak,axis=1).max()),max_GP_von_Mises_MPa=float(vm_native.max()),max_PEEQ=0.,all_accepted_states_below_initial_yield=True,constitutive_scope='source-matched finite-strain J2 virgin elastic branch; abort on yield',fields_sha256=sha(output/'fields.npz'),same_native_curved_mesh=True,no_dataset_replacement=True,backend=backend,linear_solver=linear_solver,assembly_device='cpu',linear_algebra_device='cuda' if backend=='cuda' else 'cpu',cuda_review=cuda_review,total_ksp_iterations=sum(x['iterations'] for x in linear_history)))
    failure=output/'failure.json'
    if failure.exists():failure.unlink()
    status('complete_pending_cross_solver_review');A.destroy();r.destroy();P.destroy();Cmat.destroy();Ct.destroy();zero.destroy()
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',required=True);parser.add_argument('--output');parser.add_argument('--benchmark',action='store_true');parser.add_argument('--backend',choices=['cpu','cuda'],default='cpu');parser.add_argument('--linear-solver',choices=['auto','direct','iterative'],default='auto');args=parser.parse_args();root=Path(args.root);output=Path(args.output) if args.output else root/'fenicsx'
    contract=json.loads((root/'protocol.json').read_text(encoding='utf-8'));assert sha(root/'inputs.npz')==contract['inputs_sha256']
    try:solve(dict(np.load(root/'inputs.npz')),output,benchmark=args.benchmark,backend=args.backend,linear_solver=args.linear_solver)
    except Exception as exc:
        write(output/'failure.json',dict(error=repr(exc),time_unix=time.time(),no_accepted_full_case=True));raise
