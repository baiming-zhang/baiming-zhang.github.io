# C02 detailed COPV elastic benchmark

Accepted C01 construction: 12 structural solids (Al liner, epoxy bondline, four orthotropic ply groups, two polar collars, two Ti bosses, two steel inserts). Four O-rings and two PEEK rings remain display-only. Shared structural interfaces are perfectly bonded. No contact opening, seal compression, plasticity, damage, burst or rated-pressure qualification is calculated.

Units: mm, N, MPa, N mm. Pressure is 1 MPa on the wetted inner wall and port shoulders. Separate closed-end resultants p*pi*r^2 are applied once to the two annular ends: 326.8513 N and 213.8246 N. A rank-six 3-2-1 nodal gauge removes free rigid motion; it is not a mounting clamp. The independent discrete-surface closure check records the small force remainder caused by analytical circular cap areas versus polygonized bores.

The accepted P1 solve has 194,235 corner nodes and 582,705 displacement DOF on 914,613 tetrahedra. The coordinate/output mesh stores 1,380,976 TET10 points. Midside displacements are exact affine interpolation, not additional P2 unknowns. Geometry is straight. Full detailed P2 AMG attempts exceeded the current resource budget. The two sampled polar collars required an inner-interface fit to the exact outer-ply surface; source hashes and volume changes are recorded in partition_review.json. Mesh, P1/P2 and local peak-stress convergence belong to C03.

Material constants and assumptions are explicit in configs/copv_c02_structural_v1.json. The carbon lamina uses all nine 3-D constants, engineering shear and GP-dependent [fiber, transverse, normal] axes. Ply angles are [55,-55,-55,55] degrees; polar collars use +25/-25 degrees. These are research direction fields, not qualified winding paths. Ti E is inherited from R02, while nu=0.28 is the retained modeling assumption. Steel E=200 GPa is sourced, nu=0.30 is assumed. Adhesive E=3 GPa and nu=0.35 are both numerical assumptions.

Backend: DOLFINx 0.11.0, PETSc 3.25.6. CPU FFCx assembly; PETSc CUDA FGMRES/GAMG on full displacement operator, unsmoothed aggregation and Jacobi/Chebyshev levels; small coarsest SVD on CPU, hierarchy recorded in KSP view. Matrix/vector types: seqaijcusparse / seqcuda. Exact hierarchy, coarse solve and options are retained in evidence/fenicsx/ksp_view.txt and petsc_options.json. Failed solver attempts are not accepted results.

Validation: independent skfem versus GPU FEniCSx axisymmetric vessel U relative L2=4.72380392e-10, volume-weighted raw GP tensor difference=2.91597718e-09. A rotated straight-sided TET10 patch matches analytical fields in both CalculiX and GPU FEniCSx; CalculiX stress output is ASCII precision. These P2 verification specimens are distinct from the full detailed P1 solve. Neither check is a second full detailed-assembly CalculiX calculation.

Detailed solve: true residual=1.01809678e-10, energy/work error=8.91509089e-12, half-load independent solve error=9.79914841e-10. Full peak solve 28.01 s; accepted-run total 118.59 s (excludes earlier rejected CAD/mesh/solver attempts and independent reference preparations). Twenty detailed display states are exact linear superpositions, independently checked at half and zero load. The smaller axisymmetric reference performed twenty actual solves. No nonlinear unloading or hysteresis is claimed.

Raw original-element GP local/global stress, coordinates and weights are saved locally as float64 .npy files. Public stress colors use per-cell volume-weighted GP means, retain material/ply discontinuities and are not raw peak values. P1 displacement is evaluated for the FE cut section. No temporal frame blending is used. Von Mises is displayed only for metal parts; fiber stress is displayed for composites. Composite strength allowables/failure indices are not supplied.

Reproduce in a DOLFINx/PETSc CUDA environment with CPU HYPRE, NumPy, SciPy, Gmsh 4.15.2, scikit-fem and optional CalculiX 2.23:
```
python -m tools.mesh_copv_c02_partition --partition evidence/mesh/structural_partition.brep --review evidence/mesh/partition_review.json --output mesh_h6_netgen --h 6 --curvature 64 --algorithm 1 --straight-geometry --optimize-netgen
python -m tools.copv_fenicsx_structural --mesh mesh_h6_netgen/mesh.npz --review mesh_h6_netgen/mesh_review.json --contract configs/copv_c02_structural_v1.json --output fenicsx
```
Set PETSC_OPTIONS='-use_gpu_aware_mpi 0', OPENBLAS_NUM_THREADS=4 and OMP_NUM_THREADS=4 for this single-rank run. The download includes the repaired partition, source and small verification evidence; the volume mesh and large raw GP arrays are generated locally, not embedded in the website download. The exact executed solver is retained in evidence/fenicsx/solver_executed_snapshot.py and matches solver_sha256. The delivered tools/copv_fenicsx_structural.py only corrects a stale displacement-description string in report metadata; numerical code and saved fields are unchanged. CLI help on the reference and patch scripts gives their independent reproduction steps.

The live viewer has independent boundary-condition and load visibility buttons. A/B/C mark the actual saved lower-port nodes: A locks Ux/Uy/Uz, B locks Uy/Uz, C locks Uz in global coordinates. This is a numerical 3-2-1 gauge, not a fixed end face, wall or hinge. Orange highlights all 57,544 selected wetted facets; direction arrows are a 28-facet sample. Magenta highlights the actual annular end faces and displays the two equivalent cap resultants. Markers are shown through the model and follow the displayed displacement; arrow lengths are schematic. Loads and gauge reactions scale with the selected elastic state, and zero load has no force arrows. Visibility buttons do not change the solved boundary conditions. No FEM calculation was rerun for this display addition.

Re-export these locations after the mesh/solution commands above:
```
python -m tools.prepare_copv_c02_surface --mesh mesh_h6_netgen/mesh.npz --output surface
python -m tools.export_copv_c02_bc --root . --output visualization
```
The export verifies the mesh hash, six saved constraint components, selected facet counts and the integrated pressure/cap resultant against the accepted C02 evidence.
