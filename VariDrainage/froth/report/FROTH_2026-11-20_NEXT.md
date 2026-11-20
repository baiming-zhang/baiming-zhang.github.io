# FROTH implementation handoff —2026-11-20 UTC

Active branch p6-froth-foam-drainage-v2; retired branch names must not return.
User explicitly requested actual topology implementation AND DRIFT continuation.
This supersedes the previous reporting-only pause for the work recorded here.
Goal automation was not resumed; no background native process remains.

READ: results/froth-report/FROTH_2026-11-20_Report.html (human),
results/topology-native-20261120/trajectory-audit.json and final-release/summary.json.

Implemented: pinned BSD MultiTracker/Los Topos headless adapter under foam.
Assertions always on, /fp:strict, collision safety on. One active vendor change:
const comparator for MSVC. No volume/normal-inversion bypass. Rejected cleanup
candidate was removed after a failed exact-input diagnostic. MPL2 SoapFilm3D
reference source remains clearly separate and licensed, not merged into solver.

Evidence:
- Final release8/8 fixture/control/rejection tests pass; independent crossing
  control is caught. Three triangulations14/15 pass; finest T2 remains timeout.
- Four prescribed-geometry trajectories pass saved-state independent audit:
  merge23states2->1; neck29states1->2; T1-local-thin24states removes2-3/adds0-1
  through native collapse2/VERTEX_POP2; T2 25states removesregion1/retains2.
- Actual seven-bubble saved dry mesh imports, all8 labels preserved.
- Merge surgery changes geometric enclosed volume +1.039%; neck -0.02389%.
  Not a conservative gas/film-liquid physical simulation. T1 open-boundary
  fixture may move outside the original coordinate box: no hard-wall claim.
- Original MCF T1 did not transition; insufficiently thin compression failed
  its unchanged volume-change bound. These failures are retained. Further
  geometric approach, not a relaxed volume gate, produces the successful T1.

Build xmake build froth_topology_native (portable Windows MSVC).
Fresh outputs only. Python foam/scripts/run_native_topology.py t1-local-thin
reproduces successful T1 inputs. Other runners and independent auditor are under
foam/scripts. Original output locations refuse reuse; set a NEW output root for
new experiments. Each process cap<=120s; no >2h job. PNG+vectorPDF and MP4 derive
from saved states; render_native_topology.py/build_topology_report.py reproduce.
Large failed-MCF state JSONs are losslessly zipped with hash manifest, not deleted.

NEXT ORDER:
1. Diagnose refined T2 at exact last collapse21->8. Broadly compare established
   geometry implementations; do not disable collision, volume or normal gates.
   Strict float, BLAS1, altered operation order and native flap cleanup did not
   fix it. Preserve all attempts and distinguish a budget exit from a crash.
2. T1 boundary/resolution tests, then apply operations to the actual free foam.
3. Couple REAL liquid/surfactant/gas inventories using event-local conservative
   transfer and new-film state initialization. No automatic material deletion
   on T2 or unaccounted new interface. Verify energy and physical flux too.
4. Run short coupled scenario/controls, then original lifecycle; no prescribed
   rupture order, no physical or novelty claim from geometric fixtures alone.

DRIFT separately remains below original space gate. Its boundary native probe
and conservative P1 remap are numerical components, not a corrected physical
solver. See its owner branch; do not patch thin_film privately in FROTH.
