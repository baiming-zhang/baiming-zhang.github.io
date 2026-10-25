"""Independent three-displacement axisymmetric Quad9 elasticity prototype.

Coordinates are (r,z), components (u_r,u_theta,u_z).  The cylindrical gradient
has columns r,theta,z and entries [[ur,r,-ut/r,ur,z], [ut,r,ur/r,ut,z],
[uz,r,0,uz,z]]. Symmetrization gives the six engineering strains below.
The missing theta derivatives are zero, not the circumferential displacement.
Only positive-radius meridians are supported. No production labels are made.
"""
from __future__ import annotations
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from skfem import MeshQuad2,ElementQuad2,ElementVector,Basis,FacetBasis,BilinearForm,LinearForm,asm
from geobc.copv_material import require_spd,fiber_axes


def strain_axisymmetric(u,r):
    """[rr,theta-theta,zz,theta-z,r-z,r-theta], engineering shear."""
    if np.any(np.asarray(r)<=0):raise ValueError('Axis r=0 is outside this prototype')
    g=u.grad
    return np.array([g[0,0],u[0]/r,g[2,1],g[1,1],g[0,1]+g[2,0],g[1,0]-u[1]/r])


def material_axes_rtz(tangent_rz,alpha_radians):
    """Right-handed [fiber,transverse,normal] columns in cylindrical r/theta/z.

    m=-increasing-z meridian tangent, c=+e_theta, n=m cross c. This matches
    the 3D vessel convention but is constructed independently in the r/z plane.
    """
    t=np.asarray(tangent_rz,float);t=t/np.linalg.norm(t,axis=-1)[...,None]
    m=np.stack([-t[...,0],np.zeros(t.shape[:-1]),-t[...,1]],axis=-1)
    c=np.zeros_like(m);c[...,1]=1
    n=np.cross(m,c)
    return fiber_axes(m,c,n,alpha_radians)


def quad9_meridian(nr,nz,ri,ro,length,warp=0.):
    """Exact-reference tensor Quad9 mapping, optionally curved at nodal samples.

    Warp applies r += warp*sin(pi*z/length) to every P2 geometry node; this
    is a curved patch specimen, not an exact-sine geometry claim.
    """
    import gmsh
    own=not bool(gmsh.isInitialized())
    if own:gmsh.initialize()
    try:ref=np.asarray(gmsh.model.mesh.getElementProperties(10)[4]).reshape(9,2)
    finally:
        if own:gmsh.finalize()
    qref=ElementQuad2().doflocs
    keys={tuple(np.rint(v+1).astype(int)):i for i,v in enumerate(ref)}
    permutation=np.array([keys[tuple(np.rint(v*2).astype(int))] for v in qref])
    shape=(2*nr+1,2*nz+1);ij=np.indices(shape).reshape(2,-1).T
    z=ij[:,1]*length/(2*nz);r=ri+ij[:,0]*(ro-ri)/(2*nr)+warp*np.sin(np.pi*z/length)
    xyz=np.c_[r,z];grid=np.arange(len(xyz)).reshape(shape)
    e=np.indices((nr,nz)).reshape(2,-1).T;off=np.rint(ref+1).astype(int)
    ids=2*e[:,None,:]+off;cells=grid[ids[:,:,0],ids[:,:,1]]
    mesh=MeshQuad2(doflocs=xyz.T,t=cells[:,permutation].T,sort_t=False)
    orig=np.full(len(xyz),-1,int);proposed=cells[:,permutation].T;dofs=mesh.dofs.element_dofs
    lo=np.full(len(xyz),len(xyz),int);hi=np.full(len(xyz),-1,int)
    np.minimum.at(lo,dofs.ravel(),proposed.ravel());np.maximum.at(hi,dofs.ravel(),proposed.ravel())
    if not np.array_equal(lo,hi):raise ValueError('Quad9 shared geometry nodes disagree')
    orig[:]=lo
    if len(np.unique(orig))!=len(xyz) or not np.array_equal(mesh.doflocs.T,xyz[orig]):raise ValueError('Quad9 original-node bijection failed')
    inverse=np.empty(len(xyz),int);inverse[orig]=np.arange(len(xyz))
    # Corner topology identifies meridian edges; field/geometry midsides remain intact.
    lookup={tuple(sorted(orig[nodes])):i for i,nodes in enumerate(mesh.facets.T)}
    groups={'inner':[],'outer':[],'bottom':[],'top':[]}
    for j in range(nz):
        groups['inner'].append(lookup[tuple(sorted([grid[0,2*j],grid[0,2*j+2]]))])
        groups['outer'].append(lookup[tuple(sorted([grid[-1,2*j],grid[-1,2*j+2]]))])
    for i in range(nr):
        groups['bottom'].append(lookup[tuple(sorted([grid[2*i,0],grid[2*i+2,0]]))])
        groups['top'].append(lookup[tuple(sorted([grid[2*i,-1],grid[2*i+2,-1]]))])
    return dict(mesh=mesh,xyz=xyz,cells=cells,original_to_sk=inverse,sk_to_original=orig,
                cell_parametric_indices=e,node_parametric_indices=ij,facets={k:np.array(v,int) for k,v in groups.items()})


def vector_basis(adapted,intorder=6):
    return Basis(adapted['mesh'],ElementVector(ElementQuad2(),dim=3),intorder=intorder)


def apply_C(C,strain):
    C=np.asarray(C)
    if C.ndim==2:return np.einsum('ij,jeq->ieq',C,strain)
    if C.ndim==3:return np.einsum('eij,jeq->ieq',C,strain)
    return np.einsum('eqij,jeq->ieq',C,strain)


def assemble_axisymmetric(basis,C):
    C=np.asarray(C,float);require_spd(C)
    if C.shape not in [(6,6),(basis.nelems,6,6),(basis.nelems,basis.X.shape[-1],6,6)]:raise ValueError('C does not match element/GP order')
    @BilinearForm
    def elasticity(u,v,w):
        eu=strain_axisymmetric(u,w.x[0]);ev=strain_axisymmetric(v,w.x[0])
        return 2*np.pi*w.x[0]*np.sum(apply_C(C,eu)*ev,axis=0)
    return asm(elasticity,basis).tocsr()


def boundary_load(basis,facets,traction,intorder=6):
    fb=FacetBasis(basis.mesh,basis.elem,facets=facets,intorder=intorder)
    @LinearForm
    def form(v,w):
        n=np.array([w.n[0],np.zeros_like(w.n[0]),w.n[1]])
        return 2*np.pi*w.x[0]*np.sum(traction(w.x,n)*v,axis=0)
    return asm(form,fb),fb


def volume_load(basis,bodyforce):
    @LinearForm
    def form(v,w):return 2*np.pi*w.x[0]*np.sum(bodyforce(w.x)*v,axis=0)
    return asm(form,basis)


def extract(adapted,basis,U,F,R,C):
    components=basis.split_indices()
    for dofs in components:
        if not np.array_equal(basis.doflocs[:,dofs],basis.mesh.doflocs):
            raise ValueError('Unexpected vector/scalar Quad9 geometry node order')
    def nodes(v):return np.column_stack([v[d] for d in components])[adapted['original_to_sk']]
    strain=strain_axisymmetric(basis.interpolate(U),basis.global_coordinates()[0]);stress=apply_C(C,strain)
    qx=np.moveaxis(np.asarray(basis.global_coordinates()),0,-1);weights=basis.dx*2*np.pi*qx[...,0]
    return dict(xyz_rz=adapted['xyz'],displacement_rtz=nodes(U),nodal_load_rtz=nodes(F),nodal_reaction_rtz=nodes(R),
                stress_rtz_voigt=np.moveaxis(stress,0,-1),strain_rtz_engineering=np.moveaxis(strain,0,-1),
                quadrature_rz=qx,quadrature_volume_weight=weights)


def solve_axial_mpc(adapted,basis,K,F,C):
    """Bottom uz=0, top common unknown uz; radial motion free, one utheta gauge.

    utheta(r_outer,z_bottom)=0 removes axial rigid rotation. This is the only
    admissible rigid rotation in the axisymmetric displacement subspace.
    """
    components=basis.split_indices();xyz=basis.mesh.doflocs.T
    bottom=np.isclose(xyz[:,1],xyz[:,1].min());top=np.isclose(xyz[:,1],xyz[:,1].max())
    gauge=np.flatnonzero(bottom)[np.argmax(xyz[bottom,0])]
    fixed=np.r_[components[2][bottom],components[1][gauge]]
    topdofs=components[2][top];mask=np.ones(basis.N,bool);mask[fixed]=False;mask[topdofs]=False
    free=np.flatnonzero(mask);rows=np.r_[free,topdofs];cols=np.r_[np.arange(len(free)),np.full(len(topdofs),len(free))]
    T=sparse.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(basis.N,len(free)+1))
    reduced=(T.T@K@T).tocsc();rhs=T.T@F;v=spsolve(reduced,rhs);U=T@v;R=K@U-F
    residual=float(np.linalg.norm(T.T@R)/max(np.linalg.norm(rhs),1e-30))
    if residual>1e-8 or not np.isfinite(U).all():raise RuntimeError('Axisymmetric reduced solve failed')
    data=extract(adapted,basis,U,F,R,C)
    diag=dict(reduced_relative_residual=residual,end_extension_mm=float(v[-1]),n_reduced_dofs=len(v),
        strain_energy_Nmm=float(.5*U@(K@U)),external_work_Nmm=float(U@F),constraint_reaction_work_Nmm=float(U@R),
        top_axial_jump_mm=float(np.ptp(U[topdofs])),gauge_utheta_mm=float(U[components[1][gauge]]),
        axial_force_balance_N=float(np.sum(data['nodal_load_rtz'][:,2]+data['nodal_reaction_rtz'][:,2])),
        axial_torque_balance_Nmm=float(np.sum(adapted['xyz'][:,0]*(data['nodal_load_rtz'][:,1]+data['nodal_reaction_rtz'][:,1]))),
        end_radial_displacement_prescribed=False,analytic_U_used_as_BC=False,
        node_association='Exact original meridian node rows; vectors in cylindrical r/theta/z',
        GP_association='Original Quad9 cell then Basis GP; volume weights already include 2*pi*r')
    return data,diag
