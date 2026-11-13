# Conditional continuum-analogue comparison, not a native trajectory test

Preregistered2026-11-13 before execution. Read saved linear_freefilm curves
and all nine modal hashes only; do not replay that reference or native runs.
Reference qualification remains17/19: two strict order-monotonicity failures.
Do not turn it into a validated truth by comparing another model to it.

Question: IF native node motion is interpreted as common interface motion,
does a standard flat long-wave analogue of its motion-plus-relative-drainage
ingredients give the same linear dynamics as the resolved surface/bulk model?
The analogue is explicitly derived here, not extracted as an exact native
Jacobian. Actual curved geometry, discrete mesh, finite time splitting and
bath are absent. This comparison cannot attribute the native spatial error.

For full thickness h, eta=h perturbation, relative surface coverage g, and
node/interface speed u, let Keta=gamma*k^2/2, p=Keta*eta and M=h^3/(12mu).
Define rho*h*u_dot=-4mu*h*k^2*u-ik*(h*Keta*eta+2E*g),
eta_dot=-ik*h*u-M*k^2*Keta*eta, g_dot=-ik*u.
Relative flux=-ik*M*p, conditional bulk mean=u+relative_flux/h.
Here inertia uses node/interface speed; h^3/(12mu) relative drainage uses
free-energy pressure, and in-plane extensional dissipation uses4mu*h.
These are continuum assumptions motivated by the exact source ingredients,
not proof that this limit is the current native discretization. Momentum
transport by relative drainage is second order about a resting flat base.

Quadratic complex-amplitude energy is rho*h*|u|^2/2+Keta*|eta|^2/2+E*|g|^2.
Its decay is 4mu*h*k^2*|u|^2+M*k^2*|p|^2. This makes the model dissipative
but does not imply physical equivalence to full linear NS. Preserve all terms;
do not fit coefficients or a relaxation exponent to reduce differences.

Use exactly the previous h,rho,mu,gamma,E,g0,wavelengths40/80/160um and ten
saved observation times through200us. Propagate ONLY the new3x3 analogue with
its exact matrix exponential. Compare eta,g,u_surface,u_mean against EACH
saved N8/16/32 reference using complex-amplitude relative L2 over the same
times. Record raw differences without a pass/fail model-agreement threshold.
Report dependence on reference order, retain its two failed conditions, and
never substitute this for original actual TIME2%/SPACE10% gates.

Analogue controls: Hermitian generator identity relative1e-12; sampled energy
nonincrease within1e-10 of initial energy. These are arithmetic/stability
checks only, not a physical acceptance of the analogue. Retain failures.

Freeze source hashes and verify selected native headers against the running
uniform5 protocol: ThinFilm.h,MovingVoronoiFilm.h,CurvedEdgeMotion.h,CurvedEdge.h.
No native modification/execution, no reference evolution replay. Isolated
single-BLAS-thread Python, budget300s. Nature red/blue PNG/vectorPDF.
Full actual curved-edge interface V/S/work,physical flux/traction,time/space
10/50/200us,exchange coefficient,thin-element ratio and experiments remain.
