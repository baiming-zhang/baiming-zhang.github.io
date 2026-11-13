# Static prepared-edge feasibility, not moving validation

Command: python thin_film/scripts/edge_preparation/initial_edge.py
Execution script frozen byte-exact as run_at_execution.py. No native dynamics.
Regular sites plus exact initial h traces at x240/320um; fixed6um bulk h and
500um radius. Constant equilibrium activity/SDS6mM capillary grand potential,
same directional dual-face area quadrature as C1. Linear Newton starting guess
is not imposed on interior solution. Initialization constraints must be released
in any future dynamics. New initial S would equal theta0*Gamma_inf*curved_area;
it is NOT the same total S as the flat-jump preparation. No kinetics or fit.

All4static pressure residuals <=9.47e-14 in3Newton updates. Bath graph-flux
estimates are near zero, consistent with prepared-edge equilibrium, not proof
of sustained drainage. Releasing trace constraints yields nonzero pressure
reactions26-46Pa. Local-vs-uniform graph flux differs42.97%/42.92%, and20-to10um
changes are about48.1% using the fine result as denominator: all fail an
EXPLORATORY10percent sensitivity check added after execution. The graph faces
adjoining trace sites are not located at x=cut; they are not independent normal
flux validation. No moving assay or larger calculation is authorized by these
static results. A spatially consistent release/physical interface flux is needed.

Saved summary_at_execution.json preserves original labels; summary.json clarifies
free-node inventory and graph-flux labels and adds postprocessing comparisons.
All original fields, protocol, figure PNG/vector PDF and driver bytes retained.
No >2h job. Full C1 conservation/normal-flux/mechanical-power/time-space gates
and experimental validation remain open. No empirical drainage law inserted.

## Physical interpretation guard

Do not demand mesh-independent instantaneous t=0 flux from this piecewise
flat/constant-curvature preparation, or tune an arbitrary smoothing width until
it passes. A curvature/pressure transition can generate a short-time relaxation
layer. This is an inference, not a proven explanation of this discretization.
Aradian et al2001 (https://arxiv.org/abs/cond-mat/0105318) distinguish initial
flat film, transition and fixed-curvature border, and examine positive-time
pinch dynamics. Salez et al2012 (https://arxiv.org/abs/1210.5905) analyze capillary
relaxation of a step; their polymer-film mobility is not our soap-film closure.
Next bounded evidence must compare fixed POSITIVE physical times and integrated
exchange/profiles, with reference resolution chosen before seeing outcomes.
No analytical profile/time law goes in the solver. Instantaneous graph-flux
checks remain diagnostics, not a new blocker replacing C1's full requirements.
