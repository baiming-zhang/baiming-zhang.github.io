"""Assemble the explicitly approved 1.5 MPa bending recipe on each native mesh."""
import numpy as np
from skfem import Basis,ElementVector,ElementTetP2
from geobc.rod_fem import quadratic_mesh,boundary_facets
from geobc.rod_v2_loading import snapshot,assemble_case
from geobc.tet_reference_faces import reference_facet_basis

def assemble(m,p):
    mesh,to_original,to_vertex=quadratic_mesh(m['xyz'],m['cells'])
    basis=Basis(mesh,ElementVector(ElementTetP2()),intorder=5)
    loaded=boundary_facets(mesh,m['faces'],m['face_group_ids'],to_vertex,2)
    engine=dict(crank_radius_mm=40.,rod_center_distance_mm=p['L'],rpm=4500,piston_bore_mm=80.,piston_assembly_mass_kg=.45)
    snap=snapshot(60,1.5,engine);old=snapshot(60,2.5,engine)
    for key in ['a_small_local_mm_s2','omega_rod_rad_s','alpha_rod_rad_s2']:
        assert np.array_equal(np.asarray(snap[key]),np.asarray(old[key]))
    case=dict(id='combined_bending',width_skew=.3,snapshot=snap)
    f,load=assemble_case(basis,loaded,case,p['L'],p['Ts']/2-p['re'],facet_basis_factory=reference_facet_basis)
    ordered=np.column_stack([f[d] for d in basis.split_indices()]);native=np.empty_like(ordered);native[to_original]=ordered
    assert np.isfinite(native).all() and native.shape==m['xyz'].shape
    return native,load
