# Same nominal resolution, different deterministic site phase

Frozen diagnostic,2026-11-13. C1 is not complete. Background-space campaign
showed uniform20 time-pair passes, but local40/30/25 vs uniform20 motion fails.

Native code restarts i=0 in append(refineStart,W,dx), whereas uniform mode
continues the global column index. Set background_ratio=1: both nominal meshes
are20um everywhere; only fine-side deterministic1.5%-spacing jitter differs.
Exact trace240um and fixed bath320um sites stay identical. The physical
preparation procedure is identical but its discrete prepared state can differ.
This isolates the combined layout/preparation perturbation, NOT dynamics alone.

No native source/model/options changes. SHA256 of binary and five native
headers/bench must match the frozen background-space execution protocol.
Before long runs:24native checks, one-step local ratio1 guards atdt.0625 and
.03125us; N=144 and all native balance<=1e-10,energy<=1e-7,force ratio<=1.000001,
preparation pressure<=1e-8. Source and binary fixed during matrix.
Then200us local ratio1 runs at both dt, same original material/preparation.
Uniform20 comparator trajectories at matching dt are read byte-exact from
background_space_20261113, NEVER rerun or overwritten. Original10%SPACE and
2%TIME dh/motion overlap gates remain at10/50/200us. Treat same-dx comparison
as layout sensitivity, not spatial convergence or physical validation.

Before interpreting, verify site matching by known logical index: local's
first80bulk sites map to uniform's first10columns of each16-column row;
local's next48sites map to columns10..15; trace8and bath8are identical.
Domain-cell counts equal; undeformed nominal sites equal. Record actual
position offsets and exact unchanged left/trace/bath comparisons. Do not
silently align different particle ids or move native sites.

Independent saved-step SI/GEOS scalar and mechanical audits. New diagnostics
retain full overlap L2 and its numerator/denominator in nm; optionally partition
by fixed x bands and Cartesian components. Such diagnostics never replace the
frozen acceptance metrics. Report initial-state field differences separately.

Budget: frozen uniform native cost574.859/1143.781s; same N,margin1.3 =>
expected747.3/1486.9s. Each<=2400s,total native compute<=4000s (<2h). On budget
or guard failure retain partial output and skip; never restart a timeout.
No job over2h, so no half-hour automation needed for this bounded campaign.
Every raw native byte retained in bounded ZIP parts; plotsPNG+vectorPDF,
white Nature red/blue. No fitted physics, measured-elasticity substitution,
new adsorption rate, or claim that jitter sensitivity explains its own cause.
