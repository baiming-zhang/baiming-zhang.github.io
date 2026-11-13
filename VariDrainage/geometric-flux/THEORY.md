# Geometric material ownership vs continuum transport

No native model or frozen trajectory changed. The geometry-only CLI records
the actual native Moving_Voronoi_Area_Jacobian times an input velocity.
Independent GEOS faces use Springel2010eq32-34: d=xj-xi,n=d/|d|,m=(xi+xj)/2,
    face_normal_speed = .5(vi+vj).n - (vj-vi).(facepoint-m)/|d|.
Fixed box walls have zero normal speed. Integrating full and right-clipped
face segments gives Adot_i and Ardot_i. Fractional ownership w_i=Ar_i/A_i has
    wdot_i = (Ardot_i - w_i Adot_i)/A_i.
Native liquid and material monolayer rates are sum_i V_i wdot_i and S_i wdot_i.
Three/five-node Gauss quadrature integrates along actual saved linear nodal
displacement; compares directly to native endpoint ownership, never a ledger
leftover. Constant nodal cut velocity is a different reconstruction; equality
to continuum fluid transport is NOT established by this geometric identity.

Primary8guards+60actual saved-step checks pass. Native Jv, independent face
derivative and stabilized central differences are separate comparisons. FD step
selection uses only two-FD-level differences<=2.5e-7, at most9levels, and original
native/FD1e-6 gates remain unchanged. Direct GEOS regular Voronoi polygons were
invalid in the first attempt; convexity repair is bounded<1e-10relative area
per raw Voronoi polygon before box clipping, no site jitter or snap.
Measured primary maximum5.15e-16; clipped native areas/weights are independently
checked against1e-8 gates. First two
oracle attempts and all raw observations remain archived.

Manufactured continuum transport uses zero-wall divergence-free smooth stream
velocity, flat-film affine positive h and Gamma (Gamma0=Gamma_inf/2, density
range0.4-0.6Gamma_inf), cut240um. Odd short-axis n7/15/23/31 avoids cut/regular
face coincidence. Continuity gives exact cut rate integral h*u_x dy and Gamma*u_x
dy, evaluated by independent analytic field quadrature. Geometry and simple
nodal-velocity traces are BOTH compared against that exact target, with all
intermediate values retained. Finest5% and overall-decrease flux gates pass for
both quantities/layouts, but the FULL manufactured protocol is NOT passed:
four regular-layout auxiliary area FD/stability guards fail. Their ownership FD
checks pass; all four jittered full geometry guards pass. The actual native Jv
and face formula agree in all eight. These facts do not establish why area FD
fails; no failed criterion is loosened or silently discarded.

This is instantaneous prescribed smooth transport, not time-integrated
manufactured advection, coupled curved-height drainage spatial convergence,
resolved dynamic meniscus, physical Stokes normal flux or experiment validation.
Current coupled-film spatial differences and physical flux remain open C1 work.
No empirical drainage rate or new transport closure was installed.

Primary references: Springel2010AREPOsec3.3eq32-34,
https://arxiv.org/abs/0901.4107 ; https://libgeos.org/usage/faq/#robustness .
Author PDF is saved locally, not rehosted. Moving-Voronoi/ALE kinematics is prior
art, not claimed as the innovative film-edge coupling contribution.
