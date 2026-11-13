# Actual moving-film temporal refinement

Only dt changes, with unchanged original model/binary, constants, geometry,
variational preparation and solver options. Fine patch starts200um; physical
cut240um. Horizon200us. All runs observe native mechanics without state mutation.
Nine frozen trajectory/observer files of the dt.25us control are byte-exact.
Material surfactant S remains insoluble and attached to cells;6mM SDS equilibrium
sets initial coverage only. No C2 exchange rate is introduced or fitted.

Original comparison: area-weighted overlap of actual moving Voronoi cells, exclude
bath intersections; relative L2 of dh=h(t)-h(0) and x(t)-x(0), normalized by the
finer member. Threshold2% for both at10/50/200us. Matched material-ID and absolute
errors supplement, never replace, this criterion. Adjacent pair agreement is not
an error estimator against a continuum exact solution or an asymptotic proof.

All accepted native steps record conservation, work-subtracted energy and force
target ratios, plus actual viscosity and BE inertial numerical loss. Independent
SI/GEOS formula and mechanics reconstruction audits saved steps only (first,
every10us and final). Finite-difference checks of weighted capillary AD are native
unit guards; no independent capillary Cauchy traction is newly validated here.
No decomposition causally isolates temporal splitting versus BE. A smaller ratio
of numerical inertia to viscosity is not an experimental viscosity correction.

No physical normal-flux or spatial-convergence threshold is relaxed by any time
pair pass. Previous local/uniform displacement space differences remain open.
Do not enlarge horizon or fit rates to fringe exponents on this evidence alone.

Supplemental pressure-gradient/Darcy and recorded-velocity physical-cut recovery
uses the inherited independent SI/GEOS procedure at20saved times10-200us. These
nonconforming estimates are compared to native fractional-ownership ledgers,
not asserted equal to an independently resolved continuum normal flux. Native
flux rate temporal differences are diagnostic; no retrospective2%flux gate is
added or substituted for the original dh/displacement criteria.
