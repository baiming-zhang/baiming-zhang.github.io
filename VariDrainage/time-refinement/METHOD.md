# Actual moving-film temporal refinement: frozen falsification

Same native curved-area film-side start200um geometry, physical cut240um,
radius500um, SDS6mM variational preparation, same materials and solver options.
Use the unchanged tf_curved_edge_bench binary from mechanical_work_fixed_20261113.
All new runs enable observer only; no solver/model implementation change.

Guards:24native checks, dt.25us200us replay with nine trajectory/observer CSV
files byte-exact to obs_dt025, then one dt.125us step. Conditional long matrix:
dt.125/.0625/.03125us, each200us, save first/every10us/final. Each<=600s,
whole native budget3300s; estimate from frozen59.531s dt.25us observed case,
linear step scaling times1.15. Startup and unchanged all-step conservation,
energy and Newton force criteria must pass. Failed attempts are preserved.
Never overwrite or rerun a frozen completed native batch.

Compare adjacent time steps at10/50/200us using the ORIGINAL moving-Voronoi-cell
overlap weighted relative dh and displacement metric, each<2%. Supplemental
matched-ID and absolute errors do not replace the original gate. Add all-step
cumulative BE inertial loss and viscosity, and sampled mechanical/formula audits.
Only dt changes; this is a time-refinement study, not causal isolation of splitting
vs BE. No spatial or experimental qualification is implied by a time-pair pass.
No >2h task planned; no physical viscosity or rate fitted to improve results.
