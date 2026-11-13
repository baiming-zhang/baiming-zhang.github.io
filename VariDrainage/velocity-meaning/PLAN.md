# Conditional velocity interpretation at the physical cut

Preregistered before execution,2026-11-13. Reuse original frozen uniform20 and
uniform10 fine-dt trajectories at10us only; no native replay or live files.
Verify every raw input SHA256 from the original run records. Current uniform5
batch and all captured native code/source/binary remain unchanged.

Source facts: Curved_Edge_Mechanics moves x by dt*v while keeping V and S
attached to nodes. Curved_Edge_Coupled_Step then transports V with mobility
h^3/(12mu), leaves S fixed in that drainage substep, and transports momentum
with volume flux. The metric/kinetic energy uses rho*V*v and rho*V*|v|^2/2.
These discrete facts alone do not establish a unique continuum velocity.

Conditional question: IF recorded node velocity is interpreted as the common
tangential velocity us of both interfaces, and IF quasi-steady parallel-slab
lubrication is applicable, what bulk profile follows from the independently
recovered pressure gradient? Here h denotes FULL thickness (some papers use
half-thickness). Set a=-h^2*grad(p)/(12mu); then mean velocity is us+a and
u(z)=us+1.5*a*(1-4*z^2/h^2). Integrate through full thickness. The conditional
kinetic density is rho*h*(|us|^2+2us.a+6|a|^2/5)/2; its excess over the energy
of the mean velocity is rho*h*|a|^2/10. These are independently integrated
classical profiles, not modifications to native momentum or energy ledgers.

Use the existing SI pressure and GEOS cut reconstruction: stage2 pressure,
stage1 lagged h and node velocity, exactly as prepared_motion physical_cut.
Recompute both recovered relative and material volume rates at x240um and
match the frozen global_space physical_cut.json values to relative1e-10.
Check cut length160um to relative1e-8 and Gauss profile integrals to1e-12.
These gates check input/reconstruction/arithmetic consistency ONLY.

Report the h*cutlength weighted RMS of a versus us, signed integral rates,
conditional energies per longitudinal length and positive shear variance.
No threshold labels these quantities physically acceptable/unacceptable.
Do not infer a cause of spatial error, physical inertia error, converged
normal flux/traction or a new validated closure. Previous startup diagnostic
already shows quasi-steady response may be delayed; this lift ignores that
delay, sloping surfaces, unsteady terms, lateral viscous terms and open-edge
matching. Prefer an independent coupled surface/bulk reference if recovered
relative flow is not negligible; do not fit a scalar exchange rate.

Primary reference: Lalli & Giusti(2024),JFM986,A7,
doi:10.1017/jfm.2024.335,eqs2.38-2.39 distinguish surface velocity and internal
flow; surface stress balances additionally determine surface velocity. Magnetic
forcing/NPs, surface-viscosity values and their boundary choices are not inputs
here. Classical shear/extension decomposition is prior art, not our novelty.
Save public PDF locally with retrieval URL/SHA256; no redistribution needed.

Budget300s, no native simulation/training. White Nature red/blue PNG/vectorPDF.
Original C1 time2%/space10%,10/50/200us, physical normal flux/traction, exchange
coefficient/new thin-element ratio and global-fringe/heldout experiments remain.
