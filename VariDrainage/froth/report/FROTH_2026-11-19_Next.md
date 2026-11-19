# FROTH — next steps | 2026-11-19

Prepared 2026-11-19. Numerical evidence is frozen at 2026-11-13, research commit
`9e8b322031f5b1abb5613df9252215b35f3d581e`. This document does not report new simulations.

**The long-term goal is PAUSED at the user's request. This handoff is not an
instruction to resume. Do not launch research calculations until the user resumes.**
The current request authorizes the human report and this documentation only.

Read `README.md` and its Working rules, `results/RESULTS.md`,
`foam/docs/PAUSE_CHECKPOINT_20261113.md`, and the original FROTH card in
`bubble_siggraph_plan.html` (in the parent SPH workspace when present).
The human overview is `results/froth-report/FROTH_2026-11-19_Report.html`; detailed
saved-state evidence remains in `results/spatial-recipient/`.

Active branch: `p6-froth-20261119-handoff`. The temporary archive branch `archive/p6-codex-development-20261119` was deleted at the user's request on 2026-11-19. Its old tip `39c8f74` remains an ancestor of the active branch. Do not recreate either retired name.
Handoff record: `foam/docs/FROTH_2026-11-19_Handoff.md`.

## The goal

**A small foam cluster whose drainage and gas diffusion prepare its liquid state,
and whose finite-speed rupture changes that state and the response of its neighbours.**
Start with an actual four/seven-bubble cluster. Test when one controlled initiating
hole leads to propagation or arrest. Later rupture times/order must emerge from
the declared physical model, not from scripted film deletion.

The subject is foam evolution and collapse. CODA supplies the neighbouring
film--rim--ligament/drop fragmentation component; its general fragmentation or
raindrop-breakup agenda is not a substitute for FROTH's cluster feedback goal.
Changes to CODA must follow ownership: implement/review on its owner branch and
merge, never maintain a private consumer patch. Reading the CODA report as a
writing example does not integrate the latest CODA code into this branch.

## What actually exists

1. Real four/seven-bubble shared-film and wet-border geometry, with source and
   oriented-volume checks. These checks do not prove mechanical equilibrium.
2. Separate finite-border liquid/solute, spatial transport and moving-volume
   references. Some field-convergence gates remain failed. They are not a coupled
   cluster drainage solver.
3. Explicit gas inventory/pressure/shape split reference: six cases/65 steps.
   Amount ledgers pass, but half-step gas transfer differs24.97%, post-shape flux
   defects are98%/49%, and the force target fails. Not converged coarsening.
4. CODA constant-tension near-planar retraction to FIRST contact. FROTH integrated
   owner commit83d9a936, not all subsequent results of CODA-fragmentation.
5. Finite-area material-film h/Gamma/tension response, then a seven-bubble spatial
   receiving reference:473 fluid nodes/935 material tets,23 remaining films,
   four cases/4,400 steps over0.2ms. Saved-state audits and same-mesh time checks
   pass. h=V/A, Gamma=S/(2A); received tagged amount remains bulk in the cap.

The current response is SUPPORTED: fixed lab-frame nodal preload work reaches
4.20162 times the local post-inelastic received energy. The local receiving cap
is uncalibrated; sqrt(base area)/rim radius=14.45. The inviscid P1/P0 material
volume constraints have not been shown spatially stable/converged. There is no
film--border liquid exchange, adsorption, coupled gas pressure, continued rim
retraction/topology or predicted subsequent rupture in this scene.

Local/global-initial-mix minimum h=0.78235/0.78852um. Their thickness-INCREMENT
fields differ108.8%, not a108.8% loss of thickness. Global mixing additionally
erases93.83% of the local received kinetic energy: this is NOT an equal-energy
test of spatial memory. Half-time increment-field difference0.0003627% is a time
check on the same mesh, not physical or spatial validation.

## Milestones after explicit resume

These are proposed next acceptance conditions, NOT achieved results. Preserve
the broad goal; change model choices when evidence requires it. Do not pursue
unbounded local refinement in place of the integrated scenario.

### M0. Re-establish the evidence and choose one physical scene

Inspect branch/worktree/remote, the frozen JSONs, source hashes and pause checkpoint.
Do not automatically rerun every successful matrix. Run a small build/audit only
if the toolchain or source changed; never infer physical completion from a build.
Choose a consistent material regime and a four/seven-bubble scene. Distinguish
the dilute insoluble surface diagnostic from DRIFT's soluble above-CMC SDS model.

Done when the starting geometry, liquid and gas inventories, surface state,
boundary/ambient conditions and unresolved gates are written in a bounded protocol.

### M1. Replace the artificial holding preload by actual gas/interface forces

This is the immediate priority. Couple explicit gas amounts/volumes to actual
interface and liquid motion. Distinguish the liquid incompressibility multiplier
from physical bubble gas pressure. Include ambient work and every actual support
if a frame is part of the physical setup; do not replace missing physics with
`+initial_gradient` or a hidden hold term.

Use a matched no-event control. For a claimed static starting state, retain the
existing nondimensional force target1e-4 and show it is met; a genuinely dynamic
initial state needs its own physical no-event trajectory, not an assertion that
everything should remain still. Gas/volume/interface residuals must be checked
together, not only separate split ledgers. Reuse the actual cluster throughout.

Done when the event and no-event scenes run without artificial preload, pressure
forces and external work are explicit, liquid/gas amounts close, and the chosen
baseline's mechanical and time-accuracy checks pass. A failed static preparation
is not permission to hide the residual.

### M2. Let films exchange actual liquid and surface/bulk amount with local borders

Current fixed material film volumes only change h by area stretch. Add physically
derived film--meniscus transfer on unique stocks, and the relevant viscosity,
surface transport and finite-rate adsorption when the selected bath requires them.
Determine which inertial/viscous regime is admissible; do not apply a creeping-flow
closure merely because an earlier local Stokes bench passed.

Done when film loss equals local border gain (including both interfaces and any
declared external reservoir), pressure/chemical work and dissipation are consistent,
positivity holds, and a matched local physical reference plus time/space checks
support the exchange law. Existing controls/failed gates remain archived.

### M3. Continue the event past first contact, including curved interfaces

The current rim stops at first arrival and the rupturing sheet's remainder is
unadvanced. Continue its absorption/contact/topology with the local liquid state;
preserve liquid, tagged amount, gas ownership and physical momentum destinations.
Coordinate owner-side CODA changes before integration. Check contact geometry in
3D: the archived central seven-film planar attempt had379.6um projection error
and falsely triggered while130.7um from the true liquid surface.

Done when a full initiating-sheet event has an accounted destination for every
remaining film/rim amount, no false geometric contacts/inversions are accepted,
and any subsequent failure follows the declared model without preset order.
Do not claim droplet fragmentation if that mechanism is not actually included.

### M4. Close slow evolution and fast response on the same cluster

Use drainage/coarsening-generated slow states to initialize the event, and feed
event-induced liquid/surface/gas changes back into subsequent slow evolution.
Repair the existing split gas kinetics, not just its roundoff conservation.
Use bounded matched controls: finite-speed versus instantaneous redistribution,
local versus well-mixed recipient, and enabled/disabled slow processes.

Done when at least one complete small-cluster loop executes with common inventories,
validated interfaces and bounded splitting error. A propagation/arrest comparison
must report all unmatched energy, momentum, amount and state differences. The
present additional-energy-loss control is insufficient to prove spatial memory
at equal energy.

### M5. Establish a defensible contribution and physical validation

Select public experiments/material parameters appropriate to the actual model,
then compare observables such as h(x,t), border liquid/velocity response, event
timing and propagation/arrest. Report space/time/seed/closure sensitivity and cost.
Do not import CODA drop-size statistics or a different chemical mixture as FROTH
validation. Published foam drainage/redistribution/interface-transfer methods are
strong precedents, not straw-man baselines.

A candidate contribution is a predictive finite-speed event/recipient feedback
method, or a reproducible state-dependent propagation/arrest mechanism. SI units,
bookkeeping, material h=V/A, tetrahedra and an attractive cascade video do not
establish novelty. No firstness or experiment agreement is currently established.

## Read/verify commands (no physical trajectories rerun)

```powershell
git status --short
git branch --show-current
Get-Content foam/docs/PAUSE_CHECKPOINT_20261113.md
Get-Content results/spatial-recipient/fluid-summary.json
Get-Content results/spatial-recipient/independent-fluid-audit.json
Get-Content results/gas-shape/summary.json
```

To regenerate only this human report, use the configured Python:

```powershell
..\..\.venv\Scripts\python.exe foam/scripts/build_froth_report.py
```

Native computation recipes are retained in RESULTS.md and the phase protocols;
do not execute them until resumed. Keep experiment outputs immutable and preserve
the archived exact source of the rejected curved run. Every new scientific stage
updates README, RESULTS.md and results/index.html together, with explicit-path
commits/pushes to p6-froth-20261119-handoff only. FROTH new code lives in foam/.
Figures remain white, equal-aspect, academic red/blue, PNG plus usable vector PDF.
Any training task longer than2h requires half-hour monitoring. None is running.

## Relationship to the supplied CODA example

The human/AI split and milestone organization follow the designated CODA report
and CODA_NEXT.md, read from CODA-fragmentation commit
`fd6cf8706032a7cadc14f58e03f380ccdc633c29` on2026-11-19. The D-drive HTML named
by the user was unavailable on this machine; its identity with the repository
report was not verified. CODA's text is reference material, not instructions to
execute CODA milestones or change FROTH's scope.
