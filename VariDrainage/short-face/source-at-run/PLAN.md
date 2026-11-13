# Fixed short-face diagnostic

This executes the previously committed GEOMETRY_FAILURE_NEXT.md proposal.
Original native headers, binaries, long-run source and tolerance stay unchanged.
Use a separate xmake project under this directory to avoid changing archived
thin_film/xmake.lua or root project configuration. No new library adoption.

Three9-site observations in a64x32spacing box, grid(-1,0,1)^2, perturb only
the(1,1)site's y coordinate downward by0,8e-9,8e-7spacing. Native geometry
cutoff uses the ORIGINAL tolerance1e-10; closure bound remains1e-9. Original
constructor outcomes are data, not forced to pass. Export sites, central
polygon edges, original exceptions and closure after omitting short edges.
The all-edge polygon vector sum is an explanatory telescoping identity ONLY;
it is not a repaired solver, independent physics test or new conservation gate.

Expected falsifiable pattern: unperturbed/resolved controls admitted, short
case rejected specifically by original normal closure; central filtered closure
exceeds1e-9 while full vector sum is below1e-12. Preserve all three outcomes
if this pattern fails. No adaptive perturbation or threshold search, no retry.
No dynamic C1 acceptance or exact replay of unavailable step78 geometry.

Build/config caps300s each; native cap30s. White Nature red/blue PNG/vectorPDF
and full raw bytes/metadata. No long computation or monitor needed. Original
10/50/200us, actual time/space, interface physical flux/work and experiment
scope retained; fullC1 remains incomplete.
