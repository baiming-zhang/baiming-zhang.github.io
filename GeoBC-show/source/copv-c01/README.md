# COPV C01 geometry delivery

15 CAD cases: baseline + lower/upper value of each of seven independent dimensions. Length units mm, volume mm^3. All dimensions other than R are ratios to R in the range configuration. Eligibility must be checked for every combination; the whole seven-dimensional box is not certified.

The inner profile uses a C1 cubic Bezier dome with axial tangents, reflected at z=0, a straight cylinder and two open cylindrical necks. Liner and wrap boundaries are constant analytical normal offsets, approximated by OCC interpolating splines. Profiles, dense curvature checks, continuous distances between sampled polyline segments, shared CAD interfaces, positive layer volumes and STEP volume round trips are recorded. Profile samples are not a global analytic clearance proof. No boss, cap/insert/contact, winding path, isotensoid, material qualification or manufactured vessel claim.

Reproduce geometry with Python 3.12 and numpy, scipy, gmsh 4.15.2:
```
python -m tools.build_copv_c01 --output results/c01
```
The script resolves the bundled configs relative to its own repository root. Run from this folder. The 3D viewer is a tessellation of the analytical profile, not a finite-element mesh; CAD interpolation errors are in each geometry_manifest.json.

GPU FEniCSx smoke requires a Linux DOLFINx 0.11.0 / PETSc 3.25.6 CUDA environment, scipy, a working C compiler and NVIDIA driver:
```
PETSC_OPTIONS='-use_gpu_aware_mpi 0' python -m tools.copv_fenicsx_c01_smoke --input fenicsx_geometry.npz --output fenicsx_smoke
```
FFCx integration/assembly runs on CPU; CG/GAMG algebra uses CUDA sparse matrices and vectors. CPU LU is used only as an independent equality check. CUDA is required; there is no CPU fallback. The smoke uses a 735-node/156-Quad9 meridian, two-component axisymmetric isotropic elasticity; E=69000/120000 MPa, nu=.33/.28, 1 MPa inner pressure, bottom uz=0, top neck traction-free, no pressure caps. It verifies geometry import and CUDA execution only; it is not the orthotropic C02 layup solver, a realistic closed-end pressure case, an accuracy/convergence study or a performance benchmark. C02/C03 require material/ply axes, end resultants and independent reference/convergence work.

Geometry-only STEP/BREP/meridian/manifests are under cases/. The complete numerical fields and runtime logs remain local; the geometry input and smoke summary are bundled. Original task dates/dependencies remain planning metadata; C01 was executed early on 2026-10-25. C01 geometry delivery is complete, with human review pending separately.
