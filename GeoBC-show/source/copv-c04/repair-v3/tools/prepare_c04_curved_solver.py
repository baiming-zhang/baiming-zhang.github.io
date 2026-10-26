"""Version a candidate-validation solver; preserve the original C04 executable."""
from pathlib import Path
W=Path(__file__).resolve().parents[1];P=W/'publication/GeoBC-Neural-Operator';source=P/'tools/copv_fenicsx_c04_research.py';target=P/'tools/copv_fenicsx_c04_repair.py';s=source.read_text(encoding='utf-8')
s=s.replace("ap.add_argument('--mesh',type=Path,required=True)","ap.add_argument('--repair-preflight',type=Path,required=True);ap.add_argument('--mesh',type=Path,required=True)",1)
anchor="    if review['connected_components']!=1 or source['min_det_jac'].min()<=0:raise RuntimeError('Unaccepted mesh')"
insert="""    repair_preflight=json.loads(args.repair_preflight.read_text())
    if not repair_preflight['candidate_FE_validation_allowed'] or repair_preflight['production_eligible']:raise RuntimeError('Repair candidate has no geometry/topology/load validation permission')
    if repair_preflight['mesh_sha256']!=sha(args.mesh):raise RuntimeError('Repair preflight mesh identity mismatch')
    native_path=args.repair_preflight.parent/'native_boundary_vectors.npz'
    closure_path=args.repair_preflight.parent/'quadratic_load_closure.json'
    if sha(native_path)!=repair_preflight['native_vectors_sha256'] or sha(closure_path)!=repair_preflight['closure_sha256']:raise RuntimeError('Quadratic load preflight identity mismatch')
    independent_native=np.load(native_path);independent_boundary=json.loads(closure_path.read_text())
    if not independent_boundary['passed'] or not np.array_equal(independent_native['node_ids'],source['node_ids']):raise RuntimeError('Quadratic load preflight/node identity failed')
"""
assert anchor in s;s=s.replace(anchor,insert+anchor,1)
anchor="    write(out/'cap_definition.json',dict(selected=args.cap_convention,circular_peak_N=circular_forces,discrete_polygon_peak_N=polygon_forces))"
insert="""    forces=independent_boundary['cap_forces_N'].copy()
    for cap in ['top','bottom']:
        if abs(cap_areas[cap]/independent_boundary['cap_surface_areas_mm2'][cap]-1)>1e-9:raise RuntimeError('Independent curved cap-area identity failed')
    write(out/'cap_definition.json',dict(selected='quadratic_opening_edge_integral',circular_peak_N=circular_forces,corner_polygon_peak_N=polygon_forces,quadratic_opening_peak_N=forces))"""
assert anchor in s;s=s.replace(anchor,insert,1)
start=s.index('        # On an affine Tri6 face,');end=s.index('        actual_load=np.zeros_like(original);',start)
s=s[:start]+"        # Independent native vectors integrate the true quadratic face map.\n        loads=[independent_native[key] for key in ['pressure_load_N','upper_cap_load_N','lower_cap_load_N']]\n"+s[end:]
s=s.replace("basis='Actual P2 affine Tri6 surface shape-function integrals; no corner-only load transfer'","basis='Independent native quadratic Tri6 coordinate/basis quadrature; actual quadratic bore-edge cap area; compared with true FEniCSx-assembled RHS'")
s=s.replace('cap_convention=args.cap_convention,geometry_parameters_mm=params',"cap_convention='quadratic-opening-edge',geometry_parameters_mm=params")
s=s.replace("f'Straight-sided geometry, P{args.degree} displacement; curved mesh attempt rejected; convergence belongs to C03'","review['geometry_mapping']+'; repair candidate only, new same-geometry sensitivity not yet certified'")
s=s.replace("    write(out/'summary.json',report);status", "    report['repair_preflight']=repair_preflight;report['candidate_verification_only']=True;report['original_shape_target_unresolved']=not repair_preflight['original_volume_SICN_target_passed'];report['production_eligible']=False\n    write(out/'summary.json',report);status",1)
compile(s,str(target),'exec');target.write_text(s,encoding='utf-8');print('Prepared separate curved-geometry GPU P2 candidate-verification executable')
