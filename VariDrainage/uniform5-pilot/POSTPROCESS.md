# Independent terminal-case audit plan

Prepared while the actual uniform5 native batch is running. No native source,
binary, run driver, captured README, parameters or thresholds are changed.

Run `python thin_film/scripts/uniform5_pilot/audit_terminal.py --prepare` once
to freeze this adapter, inherited physical auditors and original native protocol.
`--preflight` verifies unchanged native source/binary hashes and refusal of
unrecorded pending cases, without reading their native frames or calculating
physical audits. This is infrastructure verification, not a new simulation result.

After a case is terminal in the parent's status record, use `--ready` to audit
eligible completed cases. The parent must have recorded returncode zero, all
native gates passing, N2112, full step count and a complete raw file inventory.
The matching case cannot be marked running in live.json. Check each raw SHA256
before copying to a separate immutable audit directory; the inherited auditor
checks the copy again. A fine native case may continue in a different directory.

Reuse the unchanged independent SI/GEOS scalar and mechanical formulas from
global_space/audit_ready.py and mechanical_work/analyze.py, including all saved
steps, regional liquid/surfactant ledgers, inertia, viscous dissipation and work.
Adapt only snapshot path and reused-check metadata. Earlier 56 candidate checks
were executed in the control campaign; none are rerun here. All physical
assertions and tolerances stay unchanged. Record adapter/source/input/output
hashes, a separate subprocess log, return code and runtime. Cap each case at
900 seconds, preserve failures and partial outputs, and do not automatically retry.

Repeated ready invocations reuse only complete source/input/output-sealed audits.
An unsealed attempt must be inspected. Terminal failed or skipped cases remain
failures/skips and are not converted to passing coverage. No live raw files are
archived or published by this entry point.

Passing these saved discrete identities does not establish physical normal flux
or traction, nor temporal/spatial convergence. The original next gates remain:
uniform5 time pair <=2% at10us, then matching-dt uniform10->5 spatial comparisons
<=10%, using the existing moving-cell-overlap thickness-change/displacement
definitions. Additional 1.25us snapshots are supplemental. Original50/200us,
independent physical flux/traction, exchange coefficient/thin-element ratio and
global-fringe/heldout experiments remain open. The C1 goal remains active.
