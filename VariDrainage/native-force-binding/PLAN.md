# Native primitive binding: flat surfactant force and nodal mass

Preregistered2026-11-13 before new probe execution. Use existing unmodified
tf_curved_response_probe, no rebuild/source/default change or trajectory run.
It exports native Differentiate_Curved_Edge gradients and nodal modal masses.
Verify executable against curved_response protocol, selected headers/probe
source against both that protocol and the running uniform5 capture. Reuse
frozen flat16/flat32 native geometry exports; do not rerun their old fixtures.

Create NEW manufactured flat states, h6um, no bath, same original SDS6mM
insoluble material. Preserve native x,y,A,V,zero velocity; set S=S0*(1+a*cos(kx)).
W320um, Y160um; dx20/10um; wavelengths80/160um; amplitudes1e-3/5e-4.
Eight NEW read-only observations. Perturbations are numerical linearisation
controls, not parameters fitted to native dynamics or experiments. Avoid the
two-particle cosine mode that vanishes on cell-centred sites. Independent
flat SI Langmuir relative energy supplies the import check, not fitted offsets.
If any importer fails, preserve it and stop further cases; no tiny-oracle chase.

For q_i=sin(k*x_i), continuum force projection is
F_cont=-2*E0*a*k*sum(A_i*q_i^2), using native GRADIENT sign. Rectangular
Voronoi area derivative predicts F_native/F_cont=sin(k*dx)/(k*dx) to first order.
Finite-amplitude Langmuir terms remain in the native code; do not linearise
its input physics. Native modal mass should equal rho*sum(V_i*q_i^2).
The probe's three epsilon rows all share the SAME baseline force/mass. Use
these quantities only; do not claim its additional Hessian/FD outputs passed.

Gates: original importer area/energy relative<=1e-9, no relaxation; E0/rho
mapping verified from original parameters; modal mass relative1e-12; native
gradient field projection vs serialized force1e-12; ratio to analytic sinc
relative2e-3; paired amplitude-normalized force difference2e-3. No posthoc
gate change. Preserve all rows, including failures. Report resolution effect
as a primitive operator result, not actual dynamic spatial convergence.

Scope: partial binding of NODE inertia and initial monolayer force only. It
does not establish whether nodes are interface or bulk material particles,
the complete native continuum limit, periodic/unbounded viscous damping,
relative drainage, time splitting or actual curved-edge dynamics. Native box
viscous walls differ from infinite Fourier reference; do not compare damping
here. The independent free-film reference remains17/19 and its failures stay.

Per-case60s, whole480s; no >2h task. Existing native run/monitor untouched.
All commands/environment/binary/input/source/log/raw hashes retained. Native
raw files kept byte-exact in bounded archive with verified restoration; white
Nature red/blue PNG/vectorPDF. Original actual C1 liquid/S/work, physical flux/
traction,time/space10/50/200us,exchange coefficient/thin-element ratio and
global-fringe/heldout experiments remain open. No C2 or correction force.
