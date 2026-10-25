"""Minimal 3-D orthotropic P2 elasticity assembly, not a COPV pressure solver.

Constant, piecewise-element, or quadrature-point global engineering-Voigt
stiffness is supported. Actual plies and COPV BCs need later validation.
"""
from __future__ import annotations
import numpy as np
from scipy.sparse.linalg import spsolve
from skfem import Basis, FacetBasis, ElementVector, ElementTetP2, BilinearForm, LinearForm, asm
from skfem.helpers import sym_grad
from geobc.copv_material import require_spd


def engineering_strain(field):
    """skfem vector field -> strain [6,element,quadrature]."""
    e = sym_grad(field)
    return np.array([e[0,0], e[1,1], e[2,2], 2*e[1,2], 2*e[0,2], 2*e[0,1]])


def apply_stiffness(stiffness, strain):
    """Constant, element, or element/quadrature C times strain[6,E,Q]."""
    C=np.asarray(stiffness,float)
    if C.ndim==2:
        return np.einsum('ij,j...->i...',C,strain)
    if C.ndim==3:
        return np.einsum('eij,jeq->ieq',C,strain)
    return np.einsum('eqij,jeq->ieq',C,strain)


def assemble_stiffness(basis,stiffness):
    """Assemble K; element/GP C follows basis element and quadrature order."""
    C=np.asarray(stiffness,float)
    if C.shape not in ((6,6),(basis.nelems,6,6),(basis.nelems,basis.X.shape[-1],6,6)):
        raise ValueError('Expected C[6,6], C[E,6,6], or C[E,Q,6,6] matching basis')
    require_spd(C)
    @BilinearForm
    def elasticity(u,v,w):
        return np.sum(apply_stiffness(C,engineering_strain(u))*engineering_strain(v),axis=0)
    return asm(elasticity,basis).tocsr()


def solve_orthotropic(mesh, stiffness_global, fixed_facets, traction_facets, prescribed_displacement, traction):
    """Solve zero-body-force elasticity with prescribed displacement/traction.

    stiffness_global [6,6], [mesh.nelements,6,6] or [mesh.nelements,Q,6,6]
    uses order [11,22,33,23,13,12], engineering
    shear. prescribed_displacement(x) gets x[3,node] and returns same shape;
    traction(x,n) gets [3,facet,q] coordinates/normals and returns same shape.
    Returned nodal fields follow mesh.doflocs, including ALL quadratic nodes.
    Integration-point stress/strain follow [element,q,6]. Units mm,N,MPa.
    """
    C = np.asarray(stiffness_global, float)
    basis = Basis(mesh, ElementVector(ElementTetP2()), intorder=4)
    K = assemble_stiffness(basis,C)
    fb = FacetBasis(mesh, basis.elem, facets=traction_facets, intorder=5)
    @LinearForm
    def loading(v, w):
        return np.sum(traction(w.x, w.n)*v, axis=0)
    f = asm(loading, fb)
    fixed = basis.get_dofs(facets=fixed_facets).all()
    free = np.setdiff1d(np.arange(basis.N), fixed)
    components = basis.split_indices()
    u = np.zeros(basis.N)
    prescribed = prescribed_displacement(basis.doflocs)
    if prescribed.shape != basis.doflocs.shape:
        raise ValueError('Prescribed displacement must return shape [3,vectorDOF]')
    for i, dofs in enumerate(components):
        u[dofs] = prescribed[i, dofs]
    u[free] = spsolve(K[free][:,free].tocsc(), f[free]-K[free][:,fixed]@u[fixed])
    residual = K@u-f
    relative_residual = float(np.linalg.norm(residual[free])/max(np.linalg.norm(f[free]),1e-30))
    if not np.isfinite(u).all() or relative_residual > 1e-8:
        raise RuntimeError(f'Linear solve failed: residual {relative_residual}')
    reaction = np.zeros_like(u); reaction[fixed] = residual[fixed]
    def nodal(vector):
        for dofs in components:
            if not np.allclose(basis.doflocs[:,dofs], mesh.doflocs, rtol=0, atol=1e-13):
                raise ValueError('Unexpected scalar/vector P2 node order')
        return np.column_stack([vector[dofs] for dofs in components])
    U, F, R = nodal(u), nodal(f), nodal(reaction)
    strain = engineering_strain(basis.interpolate(u))
    stress = apply_stiffness(C,strain)
    xyz = mesh.doflocs.T
    info = {
        'n_nodes':int(len(xyz)), 'n_elements':int(mesh.nelements), 'n_vector_dofs':int(basis.N),
        'n_fixed_dofs':int(len(fixed)), 'n_free_dofs':int(len(free)),
        'free_relative_residual':relative_residual,
        'assembled_force_N':F.sum(0).tolist(), 'reaction_force_N':R.sum(0).tolist(),
        'force_balance_N':(F+R).sum(0).tolist(),
        'load_moment_origin_Nmm':np.cross(xyz,F).sum(0).tolist(),
        'reaction_moment_origin_Nmm':np.cross(xyz,R).sum(0).tolist(),
        'moment_balance_Nmm':np.cross(xyz,F+R).sum(0).tolist(),
        'strain_energy_Nmm':float(.5*u@(K@u)),
        'traction_work_Nmm':float(u@f), 'prescribed_boundary_reaction_work_Nmm':float(u@reaction),
        'integration_volume_mm3':float(basis.dx.sum()),
        'stress_association':'element_quadrature_point; engineering Voigt stress [E,Q,6], MPa',
        'strain_association':'element_quadrature_point; engineering Voigt strain [E,Q,6]',
    }
    return dict(xyz=xyz, displacement=U, nodal_load=F, nodal_reaction=R,
                stress=stress.transpose(1,2,0), strain=strain.transpose(1,2,0),
                quadrature_xyz=basis.global_coordinates().transpose(1,2,0), quadrature_weight=basis.dx), info
