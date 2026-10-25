"""Independent skfem closed-end, three-component layered COPV reference."""
from pathlib import Path
import argparse,json,time
import numpy as np
import gmsh
from geobc.copv_geometry import DEFAULT_PARAMETERS
from geobc.copv_axisymmetric_vessel import meridian_mesh,solve_vessel_2d
from geobc.copv_orthotropic_reference import local_stiffness

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    counts=(2,24,16,24,2);p=dict(DEFAULT_PARAMETERS)
    adapted=meridian_mesh(p,counts,1,1)
    gmsh.initialize();ref=np.array(gmsh.model.mesh.getElementProperties(10)[4]).reshape(9,2);gmsh.finalize()
    np.savez_compressed(a.output/'geometry.npz',xy=adapted['xyz'],cells=adapted['cells'],node_ij=adapted['node_parametric_indices'],cell_ij=adapted['cell_parametric_indices'],gmsh_reference=ref,parameters=np.array(list(p.values())),C_material=local_stiffness())
    result=solve_vessel_2d(a.output/'skfem',p,counts,liner_elements=1,wrap_per_ply=1,pressure=1.)
    print(json.dumps({'nodes':result['n_nodes'],'cells':result['n_quad9'],'seconds':result['seconds'],'max_utheta':result['max_abs_utheta_mm'],'diagnostics':result['diagnostics']}),flush=True)
