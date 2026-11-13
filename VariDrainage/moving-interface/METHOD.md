# Moving film / fractional interface experiment (2026-11-13)

Purpose: carry the stationary-reference capillary work channel into the actual
C++ moving-particle / relative-drainage split. This is a diagnostic integration,
not yet a resolved-meniscus or localized reduced solver.

The existing SI SDS material, 1.2 x 0.8 mm planar film, 6 um initial thickness,
Langmuir surface energy and actual-wall viscous discretization are unchanged.
`tf_moving_voronoi_bench ... output [ledger_cut_um]` enables an observer only.
The default eight-argument output remains unchanged. S is amount on ONE face.
The observer reconstructs the mechanics endpoint from accepted positions and
pre-drain V (mechanics leaves V/S attached); the drainage endpoint is saved too.

## Frozen matrix before calculation

T=10 us. Closed local100: dt2.5/1.25/.625 us; closed uniform100: dt1.25 us;
closed local50/uniform50: dt.625 us. Open local100: dt1.25us/radius500um.
Equilibrium local100 dt2.5us/amplitude0; other amplitudes .005.
Control: local100 dt1.25us with observer disabled, compared byte-exact frames
and numeric histories to its enabled twin. Cut800um; extra cut775um observer
control on local100 dt1.25us. Cut is bookkeeping only, no new forcing.
10 total executions including control and alternative cut.

Each native attempt <=180s, whole matrix <=1800s. No >2h training job.
Preserve all failed cases; no solver tolerance repair in this experiment.
Native unit checks first. Observer V/S closure normalized by whole-domain
initial inventories <=1e-10. Capillary and surface energy product identities
normalized by sum of absolute reported work plus endpoint change <=1e-7;
roundoff floors are 1e-30 J. Pair-exchange cancellation uses the same scale.
At least one driven case must activate nonzero geometric S exchange (>1e-8
of initial S) and capillary-shape work (>1e-8 of its absolute work scale).
Observer-independent state output must be exact.
Time/grid relative L2 comparisons of final actual h-h0 and displacement must
be reported; gates 2% time and10% space. Closed and open cases never mixed.

## Work identity and interpretation

Region right ownership w_i=area(cell_i intersect {x>=cut}) / A_i. The same
positive fraction assigns V and S. Mechanics regional transfer=sum dw_i V_i
and sum dw_i S_i; drainage=sum (w_j-w_i) F_ij with positive F i->j.
This is a piecewise-constant reconstructed regional inventory, NOT a separately
resolved physical q dot n at the cut. Keeping V/S material-cell-attached also
does not by itself establish a continuum ALE geometric-conservation law.

For e_ij=k c (h_i-h_j)^2, k=sigma/8, assign each endpoint half face energy:
E_R=sum alpha e, alpha=(w_i+w_j)/2. Arithmetic bars denote endpoint averages.
Exactly dE_R=alpha_bar [t dh_i -t dh_j +g dc]+e_bar d_alpha, with
t=2k c_bar (h_i-h_j)_bar and g=k ((h_i-h_j)^2)_bar.
Split height increments exactly: dh_i=(1/A_i)_bar dV_i + V_i_bar d(1/A_i).
Weighted nodal pressure work leaves the paired shape term
t (w_j_bar-w_i_bar)(dh_i+dh_j)/2. Area, metric and ownership terms complete
the capillary identity; shape and ownership cancel between complementary
regions, whereas geometry work is generally globally nonzero.
All terms are evaluated from explicit states/geometry; no residual source.
Monolayer energy uses the solver's same conserved uniform-reference relative
entropy; its regional area and ownership terms are recorded separately.
This checks free-energy bookkeeping, NOT regional mechanical traction/kinetic
power or an independently verified local meniscus constitutive response.

## Literature and novelty boundary

[Springel, moving Voronoi hydrodynamics](https://arxiv.org/abs/1109.2218)
and [Barrett/Garcke/Nurnberg, stable insoluble-surfactant approximation](https://numdam.org/articles/10.1051/m2an/2014039/)
establish moving conservative geometry and surfactant stability as prior art.
[Thermodynamically consistent moving contact lines, 2021](https://www.sciencedirect.com/science/article/pii/S0045782521003649)
also treats moving boundaries with soluble surfactant. No component priority
claim. Candidate future contribution: a stateful localized edge response with
verified shared material/pressure/shape ports on a moving film, including
validity/convergence limits. This observer is a prerequisite, not that method.
