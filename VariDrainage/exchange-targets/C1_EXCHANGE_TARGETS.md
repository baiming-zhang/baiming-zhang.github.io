# C1 physical exchange: what must be predicted next

2026-11-13. Research specification, not a completed physical benchmark.

## Primary-source check

[Vigna-Brummer et al., arXiv2605.27108v1](https://arxiv.org/html/2605.27108v1)
Eq.(1) defines net volume flux per contact length as
`q = ke * (gamma / eta) * h^(5/2) / r^(3/2)`.
The reported combined estimate is `ke=(2.4 +/- 0.3)*10^-2` (discussion).
The experiments use fed vertical films with Dreft or SLES/CAPB solutions,
not the present insoluble SDS assay (Sec.2.1). Appendix A Eq.(24)-(25)
contains an empirical lambda, about1.8 in the horizontal-film comparison;
Appendix B reports a different early-contact fit, lambda=1.3+/-0.2.
These are external observations/model ingredients, not parameter-free predictions.
The thin-element thickness in Eq.(24) is a model relation, not an independent
measurement of our generated elements. Section2.2's Surface Evolver solves
meniscus equilibrium; it does not resolve film-to-meniscus exchange dynamics.

The public PDF was already downloaded, not downloaded again:
`03_DRIFT/literature/edge_exchange_20261113/vigna_2026_flux.pdf`,
8101281bytes, SHA256
`8b9f30dffd7b6cb3dab8a950bfcae95263a51e8ba5ab27708b41774c0054489e`.
Keep it local; link the primary source on the public site.

## Present eligibility decision from current native code

**NOT ELIGIBLE for an experimental exchange-coefficient claim.**
`tf_curved_edge_bench.cpp` explicitly sets insoluble SDS, zero gravity and
zero disjoining pressure. Its320x160um geometry has a prescribed far bath
height/potential from a500um radius, and the native record ends at10us.
The x=240um regional cut is an accounting/observation boundary; a flux across
it is not automatically the net physical film-meniscus exchange.
No coefficient or generated-TFE thickness ratio is extracted from this run.
Its numerical conservation, work and convergence tasks remain necessary.
No conclusion about experimental ke follows even if both time/space gates pass.

## Proposed resolved-edge measurements, before a new physical run

These are project design choices, not formulas claimed from the paper.

1. Declare the actual film-meniscus contact and its normal/sign convention,
   separate from a convenient interior accounting cut. Move the observation
   surface with the geometry when required. For a validated bulk velocity,
   relative crossing is `integral h*(u_bulk-w_cut).n ds`; use the corresponding
   model decomposition only after establishing which velocity each native
   variable represents. Do not double count node motion and relative drainage.
2. Export liquid and surfactant crossings on BOTH sides, the edge inventory,
   all exterior sources and boundary work. Reconcile edge inventory change
   with contact plus exterior fluxes, and reconcile film/edge interface work
   with opposite orientation. Pressure times flux alone is insufficient if
   tangential stress/shape work is present. Preserve the original global and
   regional gates; an improved postprocessing recovery is not a replacement.
3. Export contact h, local curvature radius, tension and viscosity with units.
   First check that these observables have physical meaning at the chosen
   contact. If uniform/quasisteady assumptions hold, evaluate ke as an OUTPUT.
   For nonuniform contacts, define the reported effective diagnostic beforehand:
   `ke_eff = integral Q(t) dt / integral integral (gamma/eta)*h^(5/2)/r^(3/2) ds dt`.
   Report signed net transfer, gross bidirectional crossings and the averaging
   window. This diagnostic need not be an intrinsic constant. Report failure
   to admit a scalar closure instead of averaging away intermittent dynamics.
4. Detect genuinely generated thin elements in the two-dimensional thickness
   field. Declare segmentation, reference thick-film region, minimum resolved
   feature width and treatment of merging BEFORE observing validation outputs.
   Report thickness-ratio distributions and sensitivity to observation scale,
   not just the thinnest node or a ratio to a prescribed initial thickness.
   A1D pinch cannot establish2D element nucleation/growth/merging.
5. Freeze physical inputs and extraction choices before validation. First
   establish actual local time/space accuracy and contact-placement robustness;
   then assess whether a reduced edge retains inventory and interface power.
   Compare a resolved edge, a conservative reduced edge and the same-accuracy
   whole-domain reference. Retain a full edge if reduction is not justified.
6. Keep coefficient measurements and the generated thickness ratios as external
   tests; do not prescribe ke or empirical lambda to manufacture them. Use
   independently justified material/geometry inputs and a heldout geometry or
   boundary condition. C2 rates require independent physical evidence. The
   final check remains global fringes plus a main and heldout physical case.

## Sequence and stopping decisions

Finish the ongoing uniform5 time pair, then the original same-dt spatial
comparison. Preserve failures and the original10/50/200us checkpoints. If
10us fails, diagnose its source before expanding the trajectory horizon.
If10us passes, it licenses only that checkpoint; decide the next bounded
run from the measured cost and unresolved physical requirements.

Do not launch a new experiment-matching run until geometry, material mapping,
surface/bulk velocity meaning and the extraction protocol are reviewable.
There is no new native run, coefficient estimate, TFE statistic, model
adoption or experiment validation in this note. Existing long-job monitor
continues every30minutes. Full C1 remains active and incomplete.
