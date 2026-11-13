# Prepared moving edge: frozen diagnostic protocol

Original DRIFT purpose: predict edge pinching, exchange and experimental drainage
without empirical laws. This is C1 model validation, not experiment validation.

Build tf_curved_edge_bench and test_edge_preparation. Run `python
thin_film/scripts/prepared_motion/run.py`; inspect actual guards, then use
`--matrix` with unchanged source/binary hashes. No native overwrite or retry.

The new optional variational preparation fixes h=6um at physical x240um and
far h from independently specified radius500um at x320um. It minimizes the
constant-equilibrium-activity directional-area capillary grand potential; no
analytic interior shape is retained. Initial S=theta0*Gamma_inf*actual area,
using unchanged6mM SDS material. This is a NEW preparation with different liquid
and S inventories from old flat-jump cases. The physical domain remains320x160um,
inherited mesh placement jitter remains fixed, gravity/disjoining disabled.
Exact x240um trace sites are mobile and their h is freely evolved after init.
Only far bath is pinned/prescribed. No continuing constant-activity closure.

Four native checks: stationarity in actual variable-S energy; actual-area
coverage and mobile trace; signed-counter origin independence; released state
evolves conserving V/S. Separate existing14curved checks must still pass.
Legacy default-control fields/faces/history must remain byte-exact.

Four1us native guards cover local/uniform20/10um. Only if ALL guards pass may
the predeclared200us cases run: local/uniform20 dt1us, local20 dt.5us;
local/uniform10 dt1us, uniform10 dt.5us. Per guard180s, per matrix600s,
whole computation3300s. Guard wall times extrapolate the case cost before each
attempt; predicted>600s is skipped, not called a physics failure. No >2h job.

Observables at fixed POSITIVE10/50/200us: minimum thickness/location, evolving
thickness/displacement fields, liquid exchange into actual meniscus boundary,
fractional fixed-x240um regional relative/material transport and surfactant.
Keep cumulative quantities and signs; do not select windows after seeing flow.
Time dynamic relative L2<=2%, space<=10%, using exact final-cell overlap with
initial differences subtracted per mesh. Incomplete horizons never qualify.

All accepted steps retain scalar history. Full before/mechanics/drainage states
and faces save step1, every10us and final, with actual initial state separately.
Independent SI/GEOS audits apply ONLY at saved steps; do not claim all-step
independent formula verification. Native scalar ledgers cover all accepted steps.
Physical normal flux and regional mechanical/kinetic power remain separate C1
requirements; fractional-inventory identities do not satisfy them.

Instantaneous t0 graph flux is not a convergence gate. Aradian2001
https://arxiv.org/abs/cond-mat/0105318 distinguishes flat film/transition/border
and studies positive-time pinching. Salez2012 https://arxiv.org/abs/1210.5905
analyzes capillary relaxation of a step with a DIFFERENT polymer mobility.
Their analytical laws are diagnostic context only, not inserted in this model.
