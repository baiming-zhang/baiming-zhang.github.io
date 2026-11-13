# Actual full-domain refinement pilot

Completed source/results: https://github.com/wang-mengdi/variational-films/tree/220e6b1
Scripts: thin_film/scripts/global_space/; native execution recipe1a72891.
Native source/binary/materials unchanged; hashes/arguments in protocol.json.
No source or solver modification during the native campaign.

4new runs,482accepted steps,6saved scalar and6saved mechanical audits pass.
The existing24unit checks were verified via log/binary hashes, NOT rerun.
Independent case audits used byte-identical immutable copies, excluding live
native files. Per-case verified.json seals inputs, core auditor sources and
outputs; final coverage includes all4terminal cases and no pending cases.

ONLY10us checkpoint: uniform10time differences0.508%dh/0.500%motion pass2%.
Uniform20vs10fine-dt space differences49.562%dh/76.118%motion FAIL10%.
Both original50/200us checkpoints remain required. No full time/spatial/C1
qualification, no continuum truth or experimental validation claimed.
One-step nonconforming physical-cut recovery is not an independent Stokes
normal-flux/traction validation. C1_ALIGNMENT.md retains original goals.

Two bounded ZIP parts restore44new raw files byte-exact. Namespace native/
maps to results/local-runs/global_space_20261113/.
python restore_raw.py REPORT_DIRECTORY NEW_OUTPUT_DIRECTORY
Uniform20reference data: ../background-space/, hashes in protocol.json and
verified_inputs.json. Never overwrite frozen data or invent omitted states.
Cost projections are not launched jobs and are not controlled performance
benchmarks. Longer jobs estimated>2h require half-hour monitoring before launch.
