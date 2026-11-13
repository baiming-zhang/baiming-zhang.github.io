# Independent pressure/traction-driven local Stokes reference

Preregistered2026-11-13, before execution. This is a reusable fixed-thickness
bulk reference with prescribed stress ports, not the completed C1 edge model.
No live native frames, native replay or modification of captured uniform5 code.

Solve rho*u_t=G+mu*u_zz in z[-h/2,h/2], initially at rest. G=-p_x is a
prescribed longitudinal pressure gradient. Both surfaces have the SAME applied
tangential traction tau: mu*u_z(top)=tau, -mu*u_z(bottom)=tau. Their velocities
are outputs, not fixed values. Pressure at the strip ends and matching exterior
normal pressure p_ext(x)=p0-G*x support the flat fixed geometry. No normal
surface motion, surfactant evolution or constitutive stress closure is solved.
This is an externally driven stress-boundary problem, not an unforced free film.

Implement a cell-centred transverse finite-volume state propagator. Mass per
cell is rho*dz; interior shear stiffness is mu/dz. Neumann tractions load the
boundary cells. Exact modal advancement handles arbitrary piecewise-constant
G,tau histories, preserving the evolving bulk state. Surface face velocities
use the prescribed traction across each half boundary cell. Pressure work is
G*q, surface work is tau*(u_top+u_bottom), per longitudinal area. Include the
two boundary half-cell shear losses. Check each interval's energy balance and
integrated momentum against applied force h*G+2*tau, not a fitted relaxation.

Controls: (1) tau=0,G=1Pa/m gives exact accelerating plug u=G*t/rho;
(2) tau=-G*h/2 balances net force, so mean flux remains zero while shear and
surface motion develop; (3) reverse both ports after50us to test state memory
and signed input work. These are natural zero/cancellation cases, not parameters
chosen to fit the native result. Use the same frozen rho,mu,h=6um and17.3585um
as the preceding transverse_startup protocol. No new material coefficients.

Independent analytic balanced response: mean=0; steady
u=G*(h^2/12-z^2)/(2mu), surface speed=-G*h^2/(12mu). Its transient surface
response is 1-6/pi^2 sum_n>=1 exp(-4*pi^2*mu*n^2*t/(rho*h^2))/n^2.
For dimensionless h=rho=mu=G=1, test32/64/128/256 cells at dimensionless
modal times .1,1,4. Require monotonically decreasing maximum surface-response
error, finest<1e-3. Verify accelerating plug against its exact solution and
both physical-height histories: momentum relative error<1e-9, interval energy
relative error<1e-9, nonnegative dissipation within roundoff. Preserve failures.

Compare the SAME pressure input with earlier stationary-surface Dirichlet
startup, shear-free traction and balanced traction. Distinct fluxes illustrate
why pressure alone cannot specify the local operator before surface stress
and normal-boundary assumptions are fixed. Do not interpret the new propagator
as a replacement for native inertia or as proof of a native error. Classical
Stokes flow, FV and stress/velocity ports are prior art, not claimed novelty.

Primary context: Lalli/Giusti2024JFM986A7,eq2.7-2.16 and2.38-2.39 (local
PDF/hash already saved), separate internal flow from interfacial velocity and
require surface momentum balance. Our unsteady flat externally stressed strip
is independently derived; no magnetic, surface-viscosity or fitted inputs used.

Budget300s, no native/training task. PNG/vectorPDF, white Nature red/blue.
Keep actual C1 ledgers, local normal flux/traction, time/space10/50/200us,
exchange coefficient/thin-element ratio and global-fringe/heldout experiments
open. The eventual edge reference must resolve geometry and determine stress
from its physical surface law; this flat strip cannot do that.
