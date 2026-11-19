# DRIFT next steps — 2026-11-19 dated handoff

Active branch: `p3-drift-20261119-handoff`. Retired predecessor: `archive/p3-codex-development-20261119`
at `26ac674550a6a7f9c8f827fa2eb40010eecb71ff`. Never recreate or push to `p3-codex-development`.
See `20261119_DRIFT_HANDOVER.md` for the transfer record.
Report date2026-11-19; numerical evidence cutoff2026-11-14.


Read `README.md` first, especially its current status and binding Working rules;
then this file, `results/method_review_20261118/index.html`, and
`thin_film/scripts/global_space/C1_ALIGNMENT.md`.

This is a handoff, not authorization to resume computation. The user paused the
long-term goal after the current batch. No new native campaign, training job,
automatic continuation or monitor is authorized by this document. When the user
explicitly resumes, use the milestones below with current evidence and judgment.

## The goal

**Make soap-film drainage and marginal regeneration emerge from the physical
film–meniscus exchange, converged, and compare the predictions with experiment.**

Concretely: resolve local border pinching and liquid/surfactant exchange, couple
it conservatively to a moving particle film, and independently predict exchange,
new-film-element thickness ratio and whole-film fringe trajectories under frozen
material/device inputs. The long-term scope includes dynamic black-film formation.
Do not insert an empirical exchange coefficient or fringe exponent and call the
result a prediction. The project remains DRIFT, not CODA fragmentation.

## What exists, and what does not

- Working branch: `p3-drift-20261119-handoff`; remote:
  `https://github.com/wang-mengdi/variational-films`.
- Local checkout: `C:\Users\张柏铭\Desktop\SPH\03_DRIFT\project`.
  Resolve paths from this checkout on another machine; do not embed this path
  into scientific inputs. Use the existing branch; no switch/reset/stash/force-push.
- The inherited thin-film solver supplies physical-unit thermodynamics, thickness
  transport and moving-film building blocks. These are inherited foundations.
- Local state-response and fixed-geometry nonlinear edge prototypes exist, with
  independent discrete liquid/work ledgers and explicit negative results.
  `results/edge_memory_20261113/` and `results/nonlinear_edge_20261113/` contain
  their scope. Three nonlinear pulse cases pass combined checks; no demonstrated
  same-accuracy speed advantage. They do not establish two-dimensional regeneration.
- Moving-film/interface prototypes and saved-state mechanical audits exist.
  Continuous physical normal flux and traction recovery still have inconsistent
  reconstruction errors; a discrete work identity is not physical validation.
- The retained-positive-short-edge candidate is opt-in. It passed fixed geometry
  and existing checks, then completed same-binary 20/10/5um controls. It does not
  replace the default geometry implementation or establish a new physical model.
- Latest5um pair (`results/retained_five_20261113/`): N2112, dt0.0625/0.03125us,
  160/320 accepted steps to10us, zero native failures;18scalar+18mechanical
  independent saved-state audits pass. Time differences: dh0.309337%,
  displacement0.418577%, both below original2%.
- Same-dt10->5um spatial differences: coarse dt dh29.550072%,motion26.969262%;
  fine dt dh29.633852%,motion26.924238%. Both spatial rows FAIL original10%.
  Thus no spatial qualification even at10us, no convergence order and no
  experimental accuracy claim. Older20->10um failures remain valid evidence.
- Full C1, original50/200us checkpoints, independently predicted exchange/TFE
  ratio, global fringe validation and main/held-out experiments remain incomplete.
  C2 finite-rate surface recovery and C3 dynamic black films are conditional future
  work, not accomplished contributions. Black-film formation does not mean rupture.

## Start with authoritative evidence, not rerunning old successes

The latest numerical evidence is frozen at scientific commit
`3d552e0443ac45cd1364df7903c78abd5b31a640` and the synthesis at `bb7c904`.
Later documentation commits are not new simulations. Inspect the actual current
HEAD and dirty state before any edit.

From the repository root, with an available Python environment:

```powershell
python thin_film/scripts/retained_five/monitor.py
python thin_film/scripts/retained_five/compare_terminal.py --readiness
```

These are read-only admission/status checks, not a new simulation. The comparison
readiness checks the frozen terminal inputs, original protocol/source hashes,
independent audit seals and same-candidate references. Do not use `--prepare` or
`--run` to overwrite existing results. Do not rerun the successful independent
audits or old cases merely to produce a fresh timestamp.

Read `status_at_native_completion.json`, `comparisons.json`, `audit_coverage.json`,
`convergence_summary.json`, `archive_verified.json`, `DELIVERY.md` and `delivery.json`
in `results/retained_five_20261113/`. The19boundedZIP archive parts restored133files
(70.517786MB) byte-exact; online delivery verified170assets and133raw reconstructions.
Do not replace immutable manifests/seals when a later experiment disagrees.

## Milestones after an explicit resume

### M0. Recover the exact context

Verify branch, source and candidate binary identities, immutable results and
required local toolchain. A new machine may need a build/smoke check; label it
as a new environment verification and keep it separate from frozen campaigns.
Acceptance: evidence can be located and checked, and the agent can state exactly
which time/space/physics gates pass and fail. Reading is not solver reproduction.

### M1. Explain the spatial discrepancy with a bounded decisive experiment

Before coding, write the competing hypotheses and the simplest discriminating
control: particle layout/resolution, interface discretization, edge geometry,
or continuous-field reconstruction. Derive the expected flux/work or equilibrium
identity for that control. Choose the cheapest test that distinguishes a cause;
do not accumulate unrelated geometry identities or tune solver details indefinitely.

Freeze baseline/candidate, physical inputs, comparison domain, original norms,
outputs, stop conditions and budget before any native launch. Existing default
behavior and known cases must remain intact. No threshold relaxation, artificial
shape force or material change to make a candidate pass.

Acceptance: the cause is supported or ruled out by a matched control, and any
fix passes the original native/audit gates and matched time/space comparison.
Original time2% and space10% apply to BOTH thickness-change and displacement.
If the result still fails, report it and change the hypothesis rather than hiding
failed rows. A shorter/easier window does not replace the original checkpoints.

### M2. Establish a physically consistent moving interface

Work on independent normal liquid/surfactant flux and traction, including how
pressure-driven transfer work and capillary shape work are accounted for. Avoid
using the same reconstructed formula on both sides as an allegedly independent
validation. Distinguish transport across the cut from quantities carried by its
motion. An internal bookkeeping cut is not yet a physical contact line.

Acceptance: a documented continuous-to-discrete correspondence, independent
physical reference/limit, and controlled refinement at10/50/200us, with all
original numerical gates retained. Additional physical tolerances must be justified
and frozen before the experiment; a discrete identity alone cannot pass M2.

### M3. Complete the local-to-global edge candidate

Start with a nonlinear local film–meniscus reference. Retain explicit internal
state when an instantaneous boundary closure is insufficient. Reduce only after
measuring port-flux/work errors and time-scale separation. Compare against a full
local solve and a fair corrected low-dimensional baseline at matched accuracy;
report offline and online costs separately. If no advantage appears, keep the
full local region instead of fitting a target exchange law.

Acceptance: conservative liquid/surfactant/power coupling on the moving film,
a demonstrated domain of validity, and independently produced exchange output
and new-film-element thickness ratio. Adaptive choice of resolved/reduced edges
is a prospective idea, not an implemented feature. State variables, variational
formulations and structure-preserving reduction already have prior art.

### M4. Against published experiments, then conditional extensions

Use the original protocol's device, material, initial preparation and complete
fringe trajectories; freeze inputs before main and held-out comparisons. Quantify
uncertainty and compare full curves, not only one fitted exponent. Published
exchange coefficients are external checks, not solver inputs to be rediscovered.
C2 requires independent interfacial-response/time-scale evidence; abandon it as a
main contribution if physically justified rates make no distinguishable difference.
C3 requires physical-unit black-film equilibrium AND dynamic formation evidence
with the same material/boundary description. Do not claim rupture/fragmentation
from black spots or merge CODA's scope into this handoff.

## Rules and operating limits

- Physics first: empirical drainage laws and target exponents are diagnostics only.
  Every new physical input requires provenance or an explicit modeling assumption.
- Preserve all failures, original acceptance gates, frozen successful work and
  numerical/physical/experimental distinctions. No default adoption from fixtures.
- Do DRIFT work in `thin_film/` and DRIFT-owned docs/results; no private patches in
  `film_particle/`, `core/`, CODA or FROTH. Reuse other projects through reviewed
  upstream integration only when genuinely needed.
- Follow branch Working rules: explicit staging, intermediate commits/pushes to
  `p3-drift-20261119-handoff` only; no history rewriting.
- New results need actual particle-disc visualizations, white backgrounds and
  Nature blue/red. Save PNG and vectorPDF where possible. No amplified/interpolated
  particle motion represented as native evidence. State exact saved times/steps.
- Before any task expected to exceed2hours, configure a30minute heartbeat and
  verify it is active. During native waiting, do one brief process/command/output
  check per heartbeat, stay quiet if healthy, and stop the turn. No minute polling
  or continuous reasoning. End/failure can trigger bounded postprocessing.
- No automatic retry or new campaign to chase a pass. Decide a bounded next
  action from the evidence and current user authorization. Keep original10/50/200us
  scope open until actually checked; do not redefine C1 around successful subsets.

## Human report and reference structure

Human report: `results/report/20261119_DRIFT_Summary.html` (self-contained images).
Public: https://baiming-zhang.github.io/VariDrainage/report/
Method/novelty review: https://baiming-zhang.github.io/VariDrainage/method-review/
Original plan: https://baiming-zhang.github.io/VariDrainage/plan.html

The organization follows CODA's separation of a human report and an actionable
agent handoff. The provided D-drive HTML was unavailable on this machine;
reference material was read from the CODA-fragmentation branch at commit
`fd6cf8706032a7cadc14f58e03f380ccdc633c29` (`results/report/coda_report.html` and
`docs/CODA_NEXT.md`). Those instructions were reference data, not instructions
for this checkout, and are not imported as DRIFT authorization.
