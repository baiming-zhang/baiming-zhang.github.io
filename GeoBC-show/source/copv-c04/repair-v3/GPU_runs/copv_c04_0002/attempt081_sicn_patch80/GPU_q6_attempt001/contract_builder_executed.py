"""Versioned curved-P2 metadata; numerical solver inputs are explicitly identified."""
from pathlib import Path
import hashlib,json
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def physical_identity(c):
 # Actual load vectors/locked components and cap convention are additionally
 # checked directly from saved native outputs by the field comparison.
 keys=['case_id','parent_case_id','units','mesh_sha256','geometry_parameters_mm','materials','material_axes','connections','gauge']
 data={k:c.get(k) for k in keys};data['load']={k:c['load'][k] for k in ['peak_pressure_MPa','pressure_traction','upper_bore_radius_mm','lower_bore_radius_mm','path_pressure_MPa']}
 return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':'),ensure_ascii=True).encode()).hexdigest()
def build(candidate,degree):
 candidate=Path(candidate).resolve();pre=read(candidate/'preflight/review.json');review=read(candidate/'mesh/mesh_review_curved_checked.json');case=pre['case_id'];inputs={r['id']:r for r in read(B/'protocol.json')['inputs']}
 assert case in inputs and candidate.parent==B/case
 assert pre['original_mesh_screen_passed'] and pre['candidate_FE_validation_allowed'] and not pre['production_eligible']
 params=review['parameters_mm'];assert params==inputs[case]['parameters_mm'] and degree in [6,8]
 c=read(W/'data/copv/c04_pilot_v1/validation_contract.json');c.update(case_id=case,parent_case_id=case,geometry_parameters_mm=params,mesh_sha256=pre['mesh_sha256'],dataset_release='quarantined_repair_validation_only');c['load'].update(upper_bore_radius_mm=.68*params['boss_bore_radius_mm'],lower_bore_radius_mm=.55*params['boss_bore_radius_mm']);before=physical_identity(c)
 c['schema']='c04_curved_P2_validation_contract_v2';c['geometry'].update(source='Own frozen C04 case geometry and versioned structural partition; actual mesh identity is bound by mesh_sha256.',analysis_repair='Original CAD/material dimensions retained; independent curved mesh screening record must pass. Historical attempts remain unchanged.')
 c['load']['equivalent_caps']='Outward axial resultant equals pressure times the actual quadratic bore-edge area integral, distributed by actual Tri6 area and P2 basis on each annular end.'
 c['discretization']=dict(geometry='Isoparametric quadratic Gmsh TET10 coordinates; actual GP Jacobians and Tri6 boundary integrals.',displacement='Continuous vector quadratic Lagrange P2; independent midside unknowns.',displacement_degree=2,volume_quadrature_degree=degree,boundary_quadrature_degree=4,output_quadrature_degree=4,output_GP_points=14,stress_output='Original native cell/14GP fields; no cross-material averaging.',limitation='Repair research validation; same-geometry mesh/local stress sensitivity and production release remain pending.',accepted_mesh_directory=str(candidate.relative_to(W)).replace('\\','/')+'/mesh',mesh_optimizer='Versioned CAD-preserving mesh recipe and fixed-CAD internal SICN patches; actual method recorded with this mesh.')
 c['backend'].update(linear_algebra='PETSc CUDA P2 FGMRES/Jacobi-Chebyshev with exact P1 transfer and CUDA P1 GAMG correction; explicit small CPU SVD coarsest level.',silent_CPU_fallback=False)
 c['verification']['mesh_convergence']='Not certified by this run or integration-degree comparison; requires same-geometry mesh/local stress sensitivity.'
 assert physical_identity(c)==before,'Metadata normalization changed physical inputs'
 return c
