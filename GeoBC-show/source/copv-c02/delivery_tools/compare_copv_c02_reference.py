"""Compare independently assembled fields at exact original nodes/GPs."""
from pathlib import Path
import json
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1]/'data/copv/c02_structural_v1/reference'
a=np.load(ROOT/'skfem/solution.npz');b=np.load(ROOT/'fenicsx/fields.npz')
assert np.array_equal(a['xyz_rz'],b['xyz_rz'])
u=a['displacement_rtz'];v=b['displacement_rtz'];delta=v-u
stress=[];errors=[];weights=[];qerror=0.
for i in range(len(a['quadrature_rz'])):
    dist,idx=cKDTree(b['quadrature_rz'][i]).query(a['quadrature_rz'][i]);assert len(np.unique(idx))==len(idx)
    qerror=max(qerror,float(dist.max()));stress.append(a['stress_rtz_voigt'][i]);errors.append(b['stress_rtz_voigt'][i,idx]-a['stress_rtz_voigt'][i]);weights.append(a['quadrature_volume_weight'][i])
stress=np.array(stress);errors=np.array(errors);weights=np.array(weights);factor=np.array([1,1,1,2,2,2]);stress_error=float(np.sqrt(np.sum(errors**2*factor*weights[...,None])/np.sum(stress**2*factor*weights[...,None])))
report={'specimen':'same Quad9 geometry, actual orthotropic 3-component closed-end pressure, same axial master/gauge; independent skfem versus FEniCSx','U_relative_L2':float(np.linalg.norm(delta)/np.linalg.norm(u)),'U_max_abs_mm':float(abs(delta).max()),'utheta_relative_L2':float(np.linalg.norm(delta[:,1])/np.linalg.norm(u[:,1])),'utheta_max_abs_error_mm':float(abs(delta[:,1]).max()),'raw_GP_stress_volume_weighted_relative_Frobenius':stress_error,'GP_xyz_max_match_error_mm':qerror,'voigt_physical_shear_multiplicity':factor.tolist(),'passed':False}
report['passed']=bool(report['U_relative_L2']<1e-6 and stress_error<1e-5 and qerror<1e-8 and report['utheta_relative_L2']<1e-3)
(ROOT/'comparison.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if not report['passed']:raise RuntimeError('Independent code verification gate failed')
