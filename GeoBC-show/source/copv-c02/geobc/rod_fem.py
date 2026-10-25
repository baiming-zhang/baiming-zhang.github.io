"""Quadratic, curved tetrahedral linear elasticity; units mm, N, MPa.

Boundary group 1 is an ideal clamp (all translations zero), not a claim
about the tutorial's cylindrical support. Group 2 has constant global
traction, integrated consistently on the curved quadratic surface.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.sparse.linalg import spsolve
from skfem import MeshTet, MeshTet2, Basis, FacetBasis, ElementTetP2, ElementVector, LinearForm, asm
from skfem.models.elasticity import linear_elasticity


def quadratic_mesh(xyz, cells):
    """Map Gmsh TET10 to skfem scalar P2 DOFs without dropping curved nodes."""
    import gmsh
    own = not gmsh.isInitialized()
    if own:
        gmsh.initialize()
    try:
        reference = np.asarray(gmsh.model.mesh.getElementProperties(11)[4]).reshape(10, 3)
    finally:
        if own:
            gmsh.finalize()
    pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    edges = []
    for q in reference[4:]:
        distances = [np.linalg.norm(q - .5 * (reference[i] + reference[j])) for i, j in pairs]
        k = int(np.argmin(distances))
        if distances[k] > 1e-12:
            raise ValueError('Unexpected Gmsh TET10 reference ordering')
        edges.append(pairs[k])
    vertices = np.unique(cells[:, :4])
    old_to_vertex = np.full(len(xyz), -1, dtype=np.int64)
    old_to_vertex[vertices] = np.arange(len(vertices))
    linear = MeshTet(xyz[vertices].T, old_to_vertex[cells[:, :4]].T)
    mesh = MeshTet2.from_mesh(linear)
    edge_node = {}
    for c in cells:
        for k, (a, b) in enumerate(edges):
            key = tuple(sorted((int(c[a]), int(c[b]))))
            value = int(c[4 + k])
            if key in edge_node and edge_node[key] != value:
                raise ValueError('Inconsistent midside node shared by tetrahedra')
            edge_node[key] = value
    scalar_to_original = np.empty(mesh.doflocs.shape[1], dtype=np.int64)
    scalar_to_original[:len(vertices)] = vertices
    for e, (i, j) in enumerate(linear.edges.T):
        scalar_to_original[mesh.dofs.edge_dofs[0, e]] = edge_node[tuple(sorted((int(vertices[i]), int(vertices[j]))))]
    if len(np.unique(scalar_to_original)) != len(xyz) or len(scalar_to_original) != len(xyz):
        raise ValueError('P2 nodes do not map bijectively onto all original nodes')
    mesh.doflocs[:] = xyz[scalar_to_original].T
    return mesh, scalar_to_original, old_to_vertex


def boundary_facets(mesh, faces, face_groups, original_to_vertex, group):
    lookup = {tuple(sorted(f)): k for k, f in enumerate(mesh.facets.T)}
    ids = []
    for face in faces[face_groups == group]:
        key = tuple(sorted(original_to_vertex[face[:3]]))
        ids.append(lookup[key])
    ids = np.unique(ids)
    if not len(ids):
        raise ValueError(f'Empty surface group {group}')
    return ids


def assemble_solve(mesh, fixed_facets, load_facets, E=110000., nu=.28, force=(-6100., 0., 0.), fixed_values=None):
    """Solve an ideal clamp / consistent total-traction problem on MeshTet2."""
    start = time.perf_counter()
    basis = Basis(mesh, ElementVector(ElementTetP2()), intorder=4)
    mu = E / (2 * (1 + nu))
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    K = asm(linear_elasticity(lam, mu), basis).tocsr()
    fb = FacetBasis(mesh, basis.elem, facets=load_facets, intorder=5)
    area = float(fb.dx.sum())
    traction = np.asarray(force, dtype=float) / area
    @LinearForm
    def loading(v, w):
        return traction[0] * v[0] + traction[1] * v[1] + traction[2] * v[2]
    f = asm(loading, fb)
    fixed = basis.get_dofs(facets=fixed_facets).all()
    free = np.setdiff1d(np.arange(basis.N), fixed)
    u = np.zeros(basis.N)
    if fixed_values is not None:
        u[fixed] = np.asarray(fixed_values)[fixed]
    u[free] = spsolve(K[free][:, free].tocsc(), f[free] - K[free][:, fixed] @ u[fixed])
    residual = K @ u - f
    relative_residual = float(np.linalg.norm(residual[free]) / max(np.linalg.norm(f[free]), 1e-30))
    if not np.isfinite(u).all() or relative_residual > 1e-7:
        raise RuntimeError(f'Linear solve failed residual check: {relative_residual}')
    # Components are mapped by split_indices, never inferred from array reshape.
    components = basis.split_indices()
    xyz = mesh.doflocs.T
    U = np.column_stack([u[d] for d in components])
    F = np.column_stack([f[d] for d in components])
    for d in components:
        if not np.allclose(basis.doflocs[:, d], mesh.doflocs, rtol=0, atol=1e-12):
            raise ValueError('Unexpected vector/scalar DOF ordering')
    grad = basis.interpolate(u).grad
    strain = .5 * (grad + np.swapaxes(grad, 0, 1))
    stress = 2 * mu * strain
    trace = sum(strain[i, i] for i in range(3))
    for i in range(3):
        stress[i, i] += lam * trace
    # stress and strain: [element, quadrature point, row, column].
    sigma = stress.transpose(2, 3, 0, 1)
    deviator = sigma - np.trace(sigma, axis1=-2, axis2=-1)[..., None, None] * np.eye(3) / 3
    vm = np.sqrt(1.5 * np.sum(deviator ** 2, axis=(-1, -2)))
    reaction = np.zeros_like(residual)
    reaction[fixed] = residual[fixed]
    reaction_xyz = np.column_stack([reaction[d] for d in components])
    jac = mesh.mapping().detDF(basis.X)
    diagnostics = {
        'E_MPa': E, 'nu': nu, 'load_area_mm2': area,
        'prescribed_total_force_N': list(force), 'assembled_force_N': F.sum(0).tolist(),
        'reaction_force_N': reaction_xyz.sum(0).tolist(),
        'force_balance_relative': float(np.linalg.norm(F.sum(0) + reaction_xyz.sum(0)) / np.linalg.norm(force)),
        'load_moment_origin_Nmm': np.cross(xyz, F).sum(0).tolist(),
        'reaction_moment_origin_Nmm': np.cross(xyz, reaction_xyz).sum(0).tolist(),
        'moment_balance_Nmm': np.cross(xyz, F + reaction_xyz).sum(0).tolist(),
        'free_relative_residual': relative_residual, 'strain_energy_Nmm': float(.5 * u @ (K @ u)),
        'external_work_Nmm': float(u @ f), 'max_displacement_mm': float(np.linalg.norm(U, axis=1).max()),
        'traction_weighted_displacement_along_force_mm': float(u @ f / np.linalg.norm(force)),
        'max_quadrature_von_mises_MPa': float(vm.max()),
        'min_abs_quadrature_detJ': float(np.abs(jac).min()),
        'jacobian_sign_change_elements': int(np.sum(np.min(jac, axis=1) * np.max(jac, axis=1) <= 0)),
        'n_nodes': int(len(xyz)), 'n_elements': int(mesh.nelements), 'n_vector_dofs': int(basis.N),
        'solver': 'scipy.sparse.linalg.spsolve', 'seconds': time.perf_counter() - start,
        'stress_association': 'element_quadrature_point; stress[e,q,i,j], MPa',
        'support': 'ideal clamp: all xyz translations fixed on big bore',
        'load': 'constant global traction on selected half small bore, consistent P2 integration',
    }
    if diagnostics['jacobian_sign_change_elements']:
        raise ValueError('Curved element has Jacobian sign change at integration points')
    arrays = dict(displacement=U, nodal_load=F, nodal_reaction=reaction_xyz,
                  stress=sigma, strain=strain.transpose(2, 3, 0, 1), von_mises=vm,
                  quadrature_xyz=basis.global_coordinates().transpose(1, 2, 0), quadrature_weight=basis.dx)
    return arrays, diagnostics


def solve(meshpath, outpath, E=110000., nu=.28, force=(-6100., 0., 0.)):
    meshpath, outpath = Path(meshpath), Path(outpath)
    with np.load(meshpath) as data:
        source = {k: data[k] for k in data.files}
    xyz, cells = source['xyz'], source['cells']
    mesh, to_original, to_vertex = quadratic_mesh(xyz, cells)
    fixed = boundary_facets(mesh, source['faces'], source['face_group_ids'], to_vertex, 1)
    loaded = boundary_facets(mesh, source['faces'], source['face_group_ids'], to_vertex, 2)
    arrays, diagnostics = assemble_solve(mesh, fixed, loaded, E, nu, force)
    for name in ('displacement', 'nodal_load', 'nodal_reaction'):
        reordered = np.empty_like(arrays[name])
        reordered[to_original] = arrays[name]
        arrays[name] = reordered
    diagnostics['node_mapping_bijective'] = True
    diagnostics['max_displacement_over_bbox_x_length'] = float(np.linalg.norm(arrays['displacement'], axis=1).max() / np.ptp(xyz[:, 0]))
    outpath.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(outpath, **source, **arrays)
    outpath.with_suffix('.json').write_text(json.dumps(diagnostics, indent=2), encoding='utf-8')
    return diagnostics


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('meshpath'); p.add_argument('outpath')
    p.add_argument('--E', type=float, default=110000.); p.add_argument('--nu', type=float, default=.28)
    p.add_argument('--force', nargs=3, type=float, default=[-6100., 0., 0.])
    a = p.parse_args()
    print(json.dumps(solve(a.meshpath, a.outpath, a.E, a.nu, a.force), indent=2))
