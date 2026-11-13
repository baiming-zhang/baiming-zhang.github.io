## Native prepared edge at positive times (2026-11-13)

[White report, actual particles, PNG/vector PDF and raw archive](index.html).
9/9executions complete: byte-exact legacy control,4native1us guards and4conditional
200us matrix cases.4new+14existing curved unit checks pass. Independent SI/GEOS
formula checks pass at168SAVED steps; all accepted steps have native scalar
history. This is sampled independent auditing, not all-step formula validation.

Minimum h falls from6um to5.80-5.83um. Thin points x221-230um lie OUTSIDE the
current fine patch beginning x240um. Film-to-edge relative transfer is positive;
far-bath counter at200us is negative (still net replenishment), so this is not
system-wide draining or experiment validation. Time dh1.73%/displacement6.75%
fails2%; local/uniform20 space31.23%/144.08% fails10%. Local20/local10 also fails.
Two uniform10time cases are skipped by the predeclared600s cost estimate gate.

A separate recovered-pressure-gradient/recorded-velocity physical-cut estimate
uses exact GEOS cut lengths. Recovery-vs-ledger relative/material differences
8.46/29.47% at local20 become2.02/8.61% at local10. This is a nonconforming
reconstruction, not exact normal-flux validation or regional mechanical power.
Next falsification: keep physical x240um cut, extend fine coverage toward the
film-side thin points and compare matched conditions. No material/rate tuning.
60raw files17.94MB restore byte-exact from6ZIP parts each<5MB;4PNG/PDF pairs.
Largest native case519.42s. Full C1 remains open; no >2h job was launched.

## Reproduce and inspect

Build tf_curved_edge_bench / test_edge_preparation, then run
python thin_film/scripts/prepared_motion/run.py followed by --matrix in a FRESH
output location. Existing protocol refuses overwrite. Frozen driver/parameters,
source hashes, binary hashes and native environment in protocol.json; source
matches50f23a3 allowing explicitly verified Windows source line-endings.
Use analyze.py, package.py, report.py for saved-step analysis/archive/report.
The archive builder refuses overwriting an existing manifest. Original raw
files stay in results/local-runs/prepared_motion_20261113; >5MB never staged.
restore_raw.py was verified into a new directory against every original SHA256.
See CUT_FLUX.md for added-analysis scope and source hashes in analysis_sources.

No additional solver tolerance or physical material change. Frozen guards.json
and METHOD.md preserve startup protocol; completed status.json supersedes the
previous in-progress readme. The initial unit build used a nonexistent Summary
helper; corrected to repository Report API before executing the four unit checks.
