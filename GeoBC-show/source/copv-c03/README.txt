C03 detailed COPV research review, 2026-10-26.
Download the original h6/h4 NPZ meshes separately from the C03 page; hashes match the solver reports. No GP or displacement samples are fabricated.
Run in a CUDA-enabled FEniCSx/PETSc environment with this directory on PYTHONPATH:
python -m tools.copv_fenicsx_c03 --mesh INPUT_H4.npz --review inputs/h4/mesh_review.json --contract configs/copv_c03_mesh_v1.json --output OUTPUT --degree 2 --preconditioner p1-gamg --cap-convention circular
Degree is selected by the CLI and recorded in each run; the inherited C02 default-displacement prose in the executed JSON is not the CLI degree. Never edit that executed JSON when validating its hash.
CPU builds the FE matrices; CUDA performs fine and P1-coarse iterative algebra. The tiny last SVD is explicitly on CPU. --save-raw-gp can additionally retain all GP arrays; original-cell GP statistics are reconstructable from the native mesh, saved nodal field and source constitutive frames.
All current meshes are RESEARCH ONLY and remain quarantined from production by the original SICN gate. Global P2 stability does not qualify local peak stress or pressure-vessel strength.
