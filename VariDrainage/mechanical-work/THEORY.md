# Regional mechanical observer: definitions and limits

The observer is const and changes no state, force, parameter, tolerance or accepted step. Six old trajectory files match the frozen film-side campaign byte for byte for all four long runs (off plus three on).

For displacement q = dt u1, old/new regional ownership weights wa/wb, node kinetic energy K and curved-area energy e:

    dK = sum(wb K1 - wa K0); ownK = sum((wb-wa) K1)
    dE = sum(wb e1 - wa e0); ownE = sum((wb-wa) e1)
    WI = sum(wa gI dot q); WV = sum(wa gV dot q); WC = sum(wa gC dot q)
    C_local = grad(sum(wa e)) dot q  (weights held fixed in differentiation)
    BE = sum(wa rho V |u1-u0|^2 / 2)
    remainder = C_local - (dE-ownE)
    capillary_port = C_local-WC
    viscous_port = D_region-WV
    dK+dE+BE+D_region+remainder-ownK-ownE-capillary_port-viscous_port
       = sum(wa (gI+gV+gC) dot q)

D_region independently integrates the actual native corrected viscous strain on old face pyramids clipped to the region. Both complementary regions have opposite discrete ports. Neither port is defined as a leftover balance residual. The capillary remainder is not asserted positive or physical dissipation. BE is numerical backward-Euler inertial loss, distinct from viscosity.

The independent Python SI/GEOS audit reconstructs inertia, native viscosity including wall midpoint samples and face corrections, clipped geometry and regional work formulas at 64 SAVED mechanical steps. Curved scalar/energy formulas are independently audited at 85 SAVED state steps. Regional weighted capillary AD is tested against finite differences in six new native unit checks; Python does not independently reimplement capillary AD. All-step native statistics are distinct from these saved-step independent audits.

Recovered cut stress integrates the local corrected-gradient Cauchy-like stress times reconstructed velocity over the fixed cut, using left-trace half-open triangle ownership. This nonconforming recovery need not equal the numerical traction of the discrete energy. Its 42.6-52.7% discrepancy is an unresolved diagnostic, not proof of invalid physical viscosity. No independent capillary Cauchy traction or exact physical normal flux is validated here.

Primary prior art: Droniou, Eymard, Gallouet and Herbin (2010), A unified approach to mimetic finite difference, hybrid finite volume and mixed finite volume methods. https://arxiv.org/abs/0812.2097 . Eymard, Gallouet and Herbin (2010), SUSHI, https://arxiv.org/abs/0801.1430 . These scalar diffusion methods motivate distinguishing stabilized numerical fluxes from recovered gradients; equivalence to this moving-film vector-viscosity discretization has NOT been established. Generic energy balance or gradient recovery is not claimed novel.

Next bounded hypothesis: derive a viscosity port directly from the native bilinear form, including face-correction contributions, and compare on affine patch and fixed-cut manufactured tests before another physical run. Then test temporal splitting/BE inertia without fitted drainage rates. Full spatial/time convergence, continuum flux/traction and experimental drainage/black-film validation remain open.
