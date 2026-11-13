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
