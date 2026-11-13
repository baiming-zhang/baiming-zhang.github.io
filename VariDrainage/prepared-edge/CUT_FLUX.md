# Independent cut reconstruction, added analysis after execution

Positive normal points toward increasing x, from film into edge across x240um.
Independent SI curved-area pressure is reconstructed from saved V/S, not from
mass balance. Each cell obtains a distance-weighted least-squares pressure
gradient using native neighbour coordinates. Relative volume flux density is
-h_lag^3/(12mu)*grad(p), using mechanical-stage h as in the frozen native mobility.
GEOS clips the ACTUAL current cells with the physical vertical line; its segment
lengths integrate normal flux. Total cut length must equal160um to1e-8relative.
Material volume and surface estimates use h_lag*u_x and (S/A_projected)*u_x.
Velocities are recorded mechanical-stage velocities; no rate is inferred from
regional inventory residuals. Compare relative, material and S parts separately.

This is a NONCONFORMING recovered field. Cell mobility vs native face-averaged
mobility and recovered vs pair pressure gradients differ at finite resolution.
It is not exact Stokes flow, a new evolution closure, or a proven physical flux.
Relative L2 differences over saved positive checkpoints are diagnostics, without
post-hoc acceptance-threshold selection. Native regional ledgers remain exact
bookkeeping; they do not by themselves establish independent normal flux.

At local20 dt1us, recovery-vs-ledger differences are8.46% relative /29.47%material
and29.46%surface. At local10 these become2.02% /8.61% /8.59%. Uniform20 differs
2.50% /38.87% /38.80%. All figures use the complete saved10us checkpoint series.
No claim of normal-flux convergence: the dynamic spatial comparisons still fail,
and uniform10was skipped by predeclared cost gate. Regional mechanical traction/
kinetic power remains unverified. Initial/motion surface thermodynamics is not
an experimental soluble-surfactant kinetics model.
