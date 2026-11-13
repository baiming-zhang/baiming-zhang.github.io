# Conditional transverse startup diagnostic

Preregistered 2026-11-13, before diagnostic execution. This is an independent
canonical PDE benchmark, not a new native run, physical-film validation or a
replacement for any C1 gate. Do not modify the running uniform5 batch.

Question: under prescribed equal stationary tangential interface velocities,
how quickly does a pressure-driven liquid slab approach its steady relative
flux? Use rho and mu from the frozen original uniform10 fine-dt case, h=6um
(initial central-film modelling choice) and maximum non-bath initial thickness
from that case (state-derived representative thick location, not new truth).
Fix pressure gradient G=1Pa/m solely as a linear probe; report normalized flux
and traction, which do not depend on this amplitude. No fitted rate.

For z in [-h/2,h/2], solve rho*u_t=G+mu*u_zz, u(+-h/2)=0, u(z,0)=0.
Let tau=rho*h^2/(pi^2*mu). The independently derived analytic sine/cosine
series gives q/q_ss=1-96/pi^4 sum_odd exp(-n^2*t/tau)/n^4,
q_ss=G*h^3/(12*mu). Top traction divided by its signed steady value is
1-8/pi^2 sum_odd exp(-n^2*t/tau)/n^2. At t=0 both responses are exactly zero.
Keep the full analytic series; tau is a mode time, not an imposed single-rate law.

Verify with cell-centred finite volumes, ghost-cell zero Dirichlet boundaries
and exact eigen evolution, N=32,64,128,256 at t/tau=.1,1,4. Require decreasing
maximum normalized flux error and N256 error<1e-3. Independently Gauss-integrate
analytic u,u_z,u_t at the same times to check dK/dt=G*q-mu*integral(u_z^2)
relative to pressure power, gate1e-9. Also compare 128/256 odd-mode truncations,
gate1e-10 at positive samples. These gates verify this canonical calculation,
not native dynamics. Preserve any failed gate.

Report response at10/50/200us, analytical time to 99% steady flux (root of the
full series), white Nature red/blue PNG and vector PDF. Do not change the
native velocity definition or replace its traction with this one. In a free
film, tangential traction must be supplied by Marangoni/surface stress and
possibly extensional/inertial dynamics; imposed stationary surface velocities
are a different boundary-value problem. Actual slab deformation, lateral
gradients, surfactant transport and open-edge matching are absent here.

Primary context: van Nierop, Scheid & Stone, JFM602(2008),119-127,
doi:10.1017/S0022112008000955, Appendix A1-A11. It derives free-interface
stress conditions and extensional balance, not this transient constrained slab.
Do not take its fitted surface viscosity or adsorption estimates as inputs.
Download the public author copy locally with URL, retrieval date and SHA256;
retain bibliographic notes and provenance. No claim of a novel classical PDE.

Budget: analysis300s; no native execution or training. Original C1 time/space,
physical normal-flux/traction, exchange coefficient/thin-element ratio and
global-fringe/heldout experiments all remain required.
