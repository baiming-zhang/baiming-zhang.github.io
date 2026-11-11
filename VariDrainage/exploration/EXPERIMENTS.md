# E1/E2 exploration protocol — 2026-11-11

This bounded batch advances the revised plan's local-border and uniform-response diagnostics.
It does not train a neural network, fit a drainage exponent, or validate a complete film/meniscus model.

## Frozen questions and cases

1. Can a fixed-cell, locally refined capillary reference reproduce a uniform fine reference while closing
   the two compartment ledgers? Twelve cases: the prior 1D regression, uniform dx=10/5/2.5/1.25 um,
   tangential dy=5/2.5 um, local patches of 30/60 um, bulk dx=10/5 um, and dt=5/2.5 us.
   Domain 200 x 80 um, initial film 2 um, a 1% tangential perturbation, reservoir -gamma/r.
   Fixed immobile interfaces; no surfactant, gravity, inertia, DLVO or resolved meniscus.
   A shared face has one liquid flux with opposite signs in the two compartment ledgers.
   This is one global finite-volume solve on an orthogonal mesh, not asynchronous subdomain coupling
   or the proposed particle-film interface. Depletion width is not a measured marginal-pinch width.
2. Does the existing SDS thermodynamics permit a distinguishable finite-rate uniform area response?
   Native C++ local cells: 2/6 um x area ratios 0.99/0.999/1.001/1.01 x four fixed-bulk response
   times 0.1/1/10/100 ms x two steps = 64 cases, each to 20 ms. Paired-salt totals and liquid
   volume stay fixed; both surfaces share the same coverage. The effective reversible pair-activity
   law obeys detailed balance for the inherited free energy. Its sole rate is a modelling choice,
   without an independently measured SDS kinetic value. Instantaneous bulk micellization remains.
   The rate sweep is sensitivity analysis and cannot establish H2 for a real material or a causal
   effect on drainage. Area jumps do external work; only subsequent fixed-area relaxation dissipates.

## Numerical gates

- Native tests: feasible coverage/bulk amount, total surfactant, fixed-area free-energy descent,
  zero-rate and large-rate limits, quiescent equilibrium, first-order convergence against independent
  explicit RK4, and equilibrium elasticity against the existing local-equilibrium energy derivative.
- Border: positive thickness, open/compartment liquid balances <1e-9 and grand-potential/implicit
  energy inequality increases <1e-9. Reproduce the independent prior 1D final minimum within 1e-8.
- Against the finest 2D case: field relative L2 <2%, integrated liquid loss, boundary minimum and
  half-depletion width changes <10%. Also report the relative error of the depletion itself so the
  nearly uniform film cannot conceal boundary errors. Report failures rather than tune to pass.
- Wall time and sampled process RSS are descriptive single-host measurements; source, data hashes,
  raw commands, native regression logs, and per-case output remain traceable.

## Run and resource bounds

Build `test_surface_exchange`, `tf_response_bench` with the repository's MSVC/xmake setup; rebuild
and check the inherited native tests after changing `Units.h`. From the repository root:

    python thin_film/scripts/phase2/run.py

The native subprocess limits plus border limit bound this batch to 55 minutes; expected runtime
is measured first with a small pilot. No task longer than two hours is launched. Native tasks use
two threads; BLAS thread limits are recorded. No caller VF_* variables are inherited or mutated.
Full response CSV and compressed border states stay in ignored `results/local-runs/exploration_20261111/`.
Public summaries and actual saved samples go in `results/exploration_20261111/`, logs in `results/logs/`.

## Decision after the batch

Conservation and resolution gates justify a subsequent moving-surface/particle local-edge prototype,
not a completed C1 claim. Finite-rate sensitivity without independent kinetics leaves C2 conditional.
Neither diagnostic justifies a Monier-device long run, an exponent claim, or SI black-film validation.

## Primary sources

- [Aradian, Raphael and de Gennes (2001)](https://arxiv.org/abs/cond-mat/0105318): physical motivation;
  the artificial fixed-potential reference is not a replication of their asymptotic calculation.
- [Thiele, Archer and Pismen (2016)](https://doi.org/10.1103/PhysRevFluids.1.083903): thermodynamic
  consistency of soluble-surfactant exchange; their substrate-film model is not directly our free film.
- [Blokhuis (2025)](https://arxiv.org/abs/2509.20110): the equilibrium model already used in `Solution.h`;
  it supplies equilibrium inputs, not the new effective kinetic coefficient.
