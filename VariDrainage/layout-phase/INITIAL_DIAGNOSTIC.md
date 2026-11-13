# Read-only initial meniscus diagnostic

2026-11-13, registered before computing this diagnostic. Does not change the
running phase matrix, preparation, native trajectories or acceptance criteria.

For constant coverage, the native preparation uses constant one-surface
tension gamma=gamma_clean+RT Gamma_inf log(1-theta0). Its one-dimensional
continuum counterpart per transverse width is

  integral[2 gamma (sqrt(1+(h'/2)^2)-1) - mu h] dx.

Euler-Lagrange pressure is

  p=-(gamma/2) h''/(1+(h'/2)^2)^(3/2).

With mu=-gamma/R, the graph

  h(x)=h0+2 (x-cut)^2/(R+sqrt(R^2-(x-cut)^2))

has p=mu for cut<x<W and matches the native endpoint formula. It is h0 for
x<=cut. This tests the constant-activity INITIAL constrained preparation only:
not the released trace, not coupled time evolution, not an imposed native
interior profile and not a measured experimental meniscus. Boundary/trace
constraint forces are excluded from the constant-pressure ODE assertion.

Read guarded split20 and frozen background40/30/25/uniform20 initial states.
Verify original raw hashes, constant coverage/tension, h0 and farH endpoints,
and analytic ODE relation. Report A-weighted nodal height discrepancy versus
this continuum graph, separately from the already preserved cell-overlap
initial-field discrepancy. Sampling a smooth graph at different points and
comparing piecewise-constant cells are different operations: neither diagnostic
replaces original moving-overlap dh/displacement gates.

This is NOT an independent in-plane force or flux test. It cannot establish
that initial-shape error causes the final motion discrepancy. It supplies a
specific follow-up before any proposed change to geometry or preparation.
