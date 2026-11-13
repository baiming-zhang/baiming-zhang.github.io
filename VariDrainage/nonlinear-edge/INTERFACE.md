# Interface definition and remaining C1 work

This is the stationary1D nonlinear reference at source commits in the frozen
protocols. It is not the native moving solver, surfactant transport, a dynamic
deep meniscus or experimental validation. No second tolerance repair is applied
after the representation-only replay. The original48attempts stay archived.

## Defined physical/numerical terms

For cell liquid volume v_i=w dx_i h_i and fixed endpoint heights, the two
symmetric interfaces have capillary energy

`E=2 gamma w sum_faces d_f sqrt(1+((h_right-h_left)/(2d_f))^2)`.

The pressure is its derivative with respect to extensive liquid volume. Mobility
is the inherited immobile-interface h^3/(12mu), frozen at the old state within
each step. Exact area energy does not upgrade long-wave hydrodynamics to a full
free-boundary Navier-Stokes model. A radius specifies boundary heights and
baseline Laplace pressure; the internal discrete equilibrium is solved, not
reset during evolution. New forcing changes only the prescribed port pressure.

All reported energies and efforts are consistently relative to equilibrium:
per-face convex Bregman energy, pressure difference and capillary effort
difference. The sum of these face energies is total relative capillary free
energy. A region's assigned half-face Bregman energy is NOT in general the raw
regional surface energy minus mu0 times its own liquid change; the interface
reference traction is part of that gauge. This definition fixes the bookkeeping.

## Shared interface and two work channels

At a fixed cut between cells L/R, let q=k(mu_L-mu_R), F=dt q, and partition
the shared resistance by the two half-cell distances. The shared potential
mu_Gamma=mu_L-q R_L gives liquid transfer -F/+F and pressure work
`-mu_Gamma F / +mu_Gamma F`. Half the shared face energy belongs to each region.

For capillary face effort `tau=gamma w [z/sqrt(1+z^2)-z0/sqrt(1+z0^2)]`,
the left shape work is `tau_new*(Delta h_L+Delta h_R)/2`; the right is its
negative. It follows directly from the assigned face energy. It is not a
computed budget residual inserted back into the ledger.

Each region satisfies `Delta E+D+N=W_external+W_pressure+W_shape` up to the
recorded nonlinear residual. N is convex backward-Euler loss, separate from
physical viscous D. The two interface-work terms cancel pairwise. The
pressure-only ablation retains its missing-work residual explicitly.

## Localized reduction and honest conservation scope

The full left region and reduced right region have block-local trial/test
spaces V/W with W^T V=I. Each region's constant liquid covector lies in the
test space. Boundary volume directions are retained in the trial space, making
the reduced boundary pressure trace equal to the nonlinear variational pressure.
The true nonlinear energy is evaluated on the reconstructed state. Port and
standard DC-corrected modal bases use the same retained directions and cap17;
no transient snapshots are used for basis construction.

These properties preserve regional liquid amounts and interface power in the
projected equations. They do NOT imply that every reconstructed edge cell obeys
the original full finite-volume continuity equation. The pointwise discrepancy
is exported. Budget closure alone is not a field-accuracy or performance claim.

## What the evidence currently permits

- 12 checks pass, including independent derivatives, extended-precision weak
 response, regional power and an unreduced two-block coordinate control.
- Original48 attempts:1strict pass. Replay48 plus stress10:56/58horizons complete,
 21/58strict passes. Failures are retained rather than replaced by terminal-only
 or lower-accuracy acceptance rules.
- Three pulse input pairs have BOTH valid full/reduced budgets and matched
 time/grid references. Other small flux differences are descriptive, unvalidated.
- The local dense projection costs more than the banded full reference here.
 Use the resolved local region until native coupling demonstrates a real cost
 benefit. State-count reduction is not the C1 contribution by itself.

Next required: port these pressure/shape ledgers into the actual moving
particle-film solver; include changing face geometry, material relative motion
and surfactant amount; verify the coupled local field, flux, full cost and
time/space convergence. A resolved nonlinear2D film-meniscus region, independent
physical constraints and experimental holdouts remain necessary. No completed
C1 or novelty-priority claim is supported yet.

The reused principles have precedents: nonlinear passive projection in
[Chaturantabut, Beattie & Gugercin2016](https://doi.org/10.1137/15M1055085),
and variational dynamic-boundary energy/force consistency in
[Liu & Wu2019](https://doi.org/10.1007/s00205-019-01356-x).
