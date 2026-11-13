# Cut work: definitions, evidence and limits

For native piecewise-affine pyramid reconstruction H(u), constant viscosity
coefficient c=2 mu h, strain bilinear product B, nodal left-region ownership w_i
and pyramid fraction alpha, the native discrete port is

    P = dt sum c |T| [alpha B(Hu,Hu) - B(Hu,H(wu))]
      = dt sum c |T| (alpha-w_i) B(Hu,Hu)
        + dt sum c |T| B(Hu,w_i Hu-H(wu)).

This is a bilinear reconstruction commutator, not a residual invented to close
the ledger. Separately, exact midpoint boundary quadrature on each clipped
affine pyramid gives D_left = P_cut + P_wall + P_internal. Therefore

    P_discrete - P_cut = P_internal + P_wall - W_nodal.

The 64 saved mechanical steps pass both identities, with worst normalized error
1.74e-14 (all components checked).48random affine fields verify quadrature and
constitutive factors to5.28e-16; these are not native nonzero-wall affine tests.
Internal terms are sums of broken-field one-sided boundary powers. Calling them
physical dissipation, or adding them as a physical closure, would be unjustified.
Prior cut work was left trace; this analysis uses the same physical orientation.

The independent smooth-field problem uses SI L320um,W160um,h6um,mu.0015Pa s,
rho1045kg/m3,U1mm/s. Manufactured cut is x=.63L=201.6um, a numerical
choice to avoid grid-aligned cuts at every tested resolution. It is NOT the
original moving-film physical cut (still240um, unchanged). Streamfunction gives
u=U sin^2(pi x/L) sin(2pi y/W),
v=-U W/L sin(2pi x/L) sin^2(pi y/W).
It is divergence-free and zero on all walls. Thus f=-mu h Laplacian(u,v).
Exact force, strain dissipation and cut work are evaluated independently in
Python with analytic derivatives and high-accuracy quadrature. Native operator
patch cases apply unchanged C to known velocities; separate steady solve cases
solve C u = A f(x_i) by diagonal-preconditioned CG to relative residual1e-10.
Area-weighted sampled body loads are a discretization choice, not exact dual
basis integrals. No body source is inserted in the physical drainage solver.
Metric_Hess_Vec has its known mass term subtracted, leaving actual native C.
No time evolution or geometry derivative is exercised in these two campaigns.

Both8-case campaigns include n8/16/24/32across short axis, regular and8%jitter.
All original analytic-field node files are byte-exact across operator and solve
campaigns. Solver residuals are checked natively at all8cases; independent
SI/GEOS operator residual is checked at the two n8guards only. Actual velocity
L2 converges14.11->4.30->2.25->1.44% regular,13.55->4.13->2.22->1.39% jitter.
Pointwise force per area fails global convergence and jittered interior errors
are not monotone, so integral agreement must not be called pointwise traction
validation. Manufactured velocity convergence does not qualify curved-height,
moving-geometry, split drainage, physical normal flux or experimental black films.

Public primary reference context: HMM (Droniou et al.2010),
https://arxiv.org/abs/0812.2097 ; SUSHI (Eymard et al.2010),
https://arxiv.org/abs/0801.1430 . These are scalar diffusion frameworks; equivalence
to our vector moving-film viscosity is not claimed. Generic energy accounting
and gradient recovery are prior art, not the proposed film-method novelty.

Next: retain this verified steady-operator baseline; isolate the actual moving
mechanics/drainage split and BE numerical inertia. Do not expand simulation
horizon or tune physical viscosity until the current temporal motion criterion
is met. Physical normal flux and dynamic spatial convergence remain open.
