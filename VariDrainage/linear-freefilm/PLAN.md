# Linear free-film surface/bulk reference

Preregistered2026-11-13 before execution. Extend beyond the prescribed-stress
strip: symmetric free surfaces, evolving full-thickness perturbation eta and
insoluble relative surface coverage perturbation g. No prescribed tangential
stress or fixed surface velocity. Exterior pressure is uniform; air inertia,
surface viscosity, diffusion, gravity and DLVO are absent. No real curved edge,
bath boundary, nonlinear deformation or C1 coupling is claimed.

Linearise about h0=6um, rho/mu from the frozen original uniform10 fine-dt case.
Choose nearest initial node to(80,80)um, derive gamma0 and per-interface
Langmuir elasticity E=RT*Gamma_inf*theta/(1-theta) from its frozen SI state.
No parameter fitted to a wave speed, native trajectory or experimental exponent.
Use wavelengths40/80/160um, perturbation g0=1e-6,eta0=0,initial velocity0.

For Fourier exp(ikx), represent even horizontal velocity with even Legendre
polynomials phi_j(z), z[-h0/2,h0/2]. Set Psi_j=integral_0^z phi_j. Then
u=sum a_j phi_j, w=-ik sum a_j Psi_j exactly obey incompressibility and symmetry.
Full linear Navier-Stokes weak form (no slender-gap pressure assumption):
M_ij=rho*integral(phi_i phi_j+k^2 Psi_i Psi_j),
D_ij=mu*integral(4k^2 phi_i phi_j+(phi_i'+k^2 Psi_i)(phi_j'+k^2 Psi_j)).
Let B=-ik*[integral(phi_j);phi_j(top)], y=(eta,g),
K=diag(gamma0*k^2/2,2E). Evolve M*a_dot=-D*a-B^*K*y and y_dot=B*a.
This includes the normal and tangential natural stress conditions through
surface energy and kinematics. Complex-amplitude energy .5a^*Ma+.5y^*Ky
is twice spatially averaged real-field energy, with the same convention for
dissipation. No diagonal/dropped bulk inertia approximation.

Use energy coordinates v=M^(1/2)*a,e=K^(1/2)*y. The generator has blocks
[-L,-C^*;C,0], with L=M^(-1/2)D M^(-1/2), C=K^(1/2)B M^(-1/2).
Propagate the finite linear system by matrix exponential; do not claim that
this replaces actual native temporal convergence. Retain all modal states
and observable time histories, reference source/material hashes and failures.

Independent controls: inviscid clean symmetric capillary wave
omega^2=(gamma0/rho)*k^3*tanh(k*h0/2), independently derived from cosh/sinh
potential flow and the normal stress condition. Require relative frequency
error<1e-8 for N8/16/32 at wavelengths40/80um. At k0, require first3 nonzero
even stress-free shear decay rates match(mu/rho)*(2pi*n/h0)^2 within1e-6
at N32; constant-velocity mode has zero decay within normalized1e-10.

Reference runs N8/16/32, times0,.5,1,2,5,10,20,50,100,200us. Require Hermitian
generator energy identity relative1e-10, D positive semidefinite to roundoff,
energy nonincrease to1e-8 normalized initial energy. Require N16->32 RMS
observable differences(eta,g,u_surface,u_mean) <1%, and smaller than N8->16.
After t0, reconstruct top stress independently from polynomial derivatives:
tau=mu*(u_z+ik*w), expected=-ik*E*g;
p=i/k*(rho*u_t-mu*(u_zz-k^2*u)), normal condition
p+2mu*ik*u_surface=gamma0*k^2*eta/2. Normalize RMS boundary residuals by
sum of component magnitudes, gate1% at N32. Initial g/velocity are not stress
compatible at t0, so exclude t0 from strong traction tests; keep its initial
energy and subsequent boundary-layer response. Failed gates remain failures.

Classical linear NS, free-surface modes, Legendre Galerkin and energy coordinates
are prior art. Literature context Lalli/Giusti2024 and vanNierop/Scheid/Stone2008,
local PDF hashes already retained. The capillary control is derived here;
do not claim an inaccessible primary source was read. No invented novelty.

Budget300s, isolated Python with one BLAS thread; no native/training task and
no changes to captured running uniform5 sources/binary. PNG/vectorPDF white
Nature red/blue. Original C1 physical flux/traction on the actual curved edge,
liquid/surfactant/work coupling,10/50/200us time/space,exchange coefficient,
thin-element ratio and global-fringe/heldout experiments all remain required.
