# Observer-only mechanical work, frozen protocol

Following film_side_20261113: improved thickness coverage did not qualify motion.
Observe the actual mechanics rather than changing tolerances/materials blindly.
Build test_curved_mechanical_observer and tf_curved_edge_bench. Run run.py, inspect
guards, then --matrix with unchanged source/binary hashes. Never overwrite.

New optional final observer flag changes outputs only. Physical state, complete
before/mechanics/drainage snapshots and every native history row must reproduce
the frozen source six files BYTE-EXACT for each200us case. Run24checks (6new,
4preparation,14curved), an observer-off dt.5us200us control and1us observer guard.
Then replay dt1/.5/.25us cases, all200us, local20 fine-start200um. Same fixed
x240um physical cut, geometry, SDS6mM preparation, boundary and all solver options.
Each<=600s; guard extrapolation<=600s; whole computation<=3300s. No >2h task.

All accepted steps: native initial/final force gradient, actual grad tolerance,
Newton/CG/negative-curvature counts, kinetic change, Rayleigh action, inertial
backward-Euler loss, incremental objective and last accepted displacement.
At step1, every10us and final: reconstruct actual frozen metric with q=dt*u_mid
(avoids subtracting nearly equal positions), common energy gradient and derivative
of a fixed weighted SUM of the SAME node energies. Global reference coverage is
unchanged. Moving ownership is accounted separately, not differentiated twice.

Regional action identity (old region weights in force work, new weights in stored
energies): dK+dE+inertial_BE_loss+viscous_D+endpoint_remainder
-ownershipK-ownershipE-capillary_port-viscous_port = force_residual dot q.
The capillary port is derived from weighted regional energy gradients, not a
remaining balance residual. Endpoint remainder may be negative and is NOT called
physical dissipation. Viscous_D independently integrates full-gradient strains on
actual clipped OLD Voronoi face pyramids. Its port compares local positive D to
weighted nodal viscous work. These are discrete variational model quantities.

Separately integrate recovered viscous stress times reconstructed velocity on
the actual fixed physical cut (left trace on exactly coincident triangle edges).
This nonconforming reconstruction is NOT an independently validated continuum
traction law. Compare it to the discrete port; do not equate an exact algebraic
identity to physical convergence. Capillary Cauchy traction, normal flux, full
space/time convergence and experimental prediction remain open C1 requirements.

Independent SI/GEOS analysis must reconstruct mass/inertia, the frozen viscosity
operator, clipped dissipation and cut traction from recorded states. Native
weighted energy derivatives are tested against direct-energy finite differences;
do not infer them from work residuals. All regional audits are sampled-step only.
No new C2 kinetic rate and no empirical drainage/exponent fit.
