# Same-candidate 5 um pair at the original 10 us checkpoint

This is a NEW pair, not a restart of the original failed uniform5 campaign.
Reuse the exact retained_motion binary fae0ce667a1456afd81850fe2c833fb95287b291c0319589d355df37c4e34c2c,
without rebuilding or changing native sources, tolerances, materials or defaults.
The four retained_reference trajectories must be terminal, all native gates and
36 scalar + 36 mechanical saved audits verified, both original 2% time pairs
passing. The original failed 20->10 um spatial rows are preserved.

Run sequentially: u5_dt00625 (160 steps, dt=.0625 us), then u5_dt003125
(320 steps, dt=.03125 us), each N=2112, T=10 us, radius=500 um, amplitude=0,
insoluble SDS=6 mM, variational preparation, refine_start=200 um, observer=true,
internal accounting cut=240 um. Save first step and every 1.25 us (stride20/40).
No prescribed physical exchange coefficient or empirical drainage law is added.

Original native gates: complete all steps with failures=0, force/target<=1.000001,
all global/regional V/S residuals<=1e-10, energy-budget excess<=1e-7,
preparation pressure residual<=1e-8. On any native failure skip the later case.
No automatic retry, new output directory, source edit or tolerance relaxation.

Reuse the frozen reference runner's direct Popen ownership and kill/wait of the
exact held process on timeout. Per-case caps18000/28000s, whole native cap46000s.
Expected combined cost is several hours, estimated from prior actual 5 um costs;
caps are termination limits, not predictions. A verified ACTIVE app heartbeat
every30min for THIS new campaign is REQUIRED before --run. Capture its official
configuration and protocol hash; do not hand-edit scheduler configuration.
Parent additionally writes actual handle/progress checks every1800s.

After terminal cases, independently audit original SI/GEOS scalar/mechanical
saved formulas. Terminal immutable subsets may be audited while the disjoint
second native case runs, never read its live frames. Preserve successful seals.
Only after both cases and their audits pass, use original moving-overlap norms:
5 um time pair2%; matching-dt retained_reference10->new5 um space10%, both dh
and displacement, retain every failed row. Do not mix old original5 um with this
new candidate or use first-step equivalence as full-trajectory qualification.

This is only the original10us checkpoint. Original50/200us, continuous physical
normal flux/traction, independent exchange and thin-element predictions, global
fringes and held-out physical experiments remain required for fullC1. No claimed
innovation from the generic geometry correction itself. Figures: actual particle
discs, white background, Nature red/blue PNG plus vectorPDF. Archive closed logs,
exact raw/source/audit bytes in boundedZIPs with an actual restoration check.
