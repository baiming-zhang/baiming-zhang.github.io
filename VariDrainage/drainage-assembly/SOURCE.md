Results/source: https://github.com/wang-mengdi/variational-films/tree/e86001d66e2e89421bc6bfe06d9affb59166cfc7

# Drainage cost and local assembly provenance

Original plan/probe:a57b92b; original driver/measurement:d755b45.
Conditional controls:fdb5c1f; local assembly/header/tests:5d00f7e;
candidate actual/operator driver:e098cb3; independent audit wrappers:1da4c76.
Code in thin_film/scripts/drainage_cost/,scripts/local_assembly/ and
include/vf/thinfilm/CurvedEdgeLocalAssembly.h. Default native class unchanged.
Existing generic solver/core/backend untouched. Macro selection only inside
separate candidate translation units. No geometry candidate dispatch here.

Observed build tool output summary(not a full captured raw build log):
tf_drainage_cost_probe succeeds10.234s;6local assembly targets succeed28.25s.
Native executable hashes/full arguments and measurement source commits in
protocol.json and candidate/protocol.json. Original bench remains
76e49ed253e2341df37f93d02803272d020b8500ed9265da6f7d0ec124b6c5c3.
All source subsets/candidate binaries fixed throughout their native campaigns.
Native stdout logs preserved verbatim,all cases terminal before packaging.

Original/ namespace restores results/local-runs/drainage_cost_20261113/;
operator/ and candidate/ restore local_assembly_20261113/ under local-runs.
Capitalization in manifest is authoritative(original/ is lowercase).
Full exact nonzero matrix entries are raw native CSV; zeros not serialized,
no magnitude cutoff. All4operator inputs+comparisons and3actualcases included.
Byte flags for all11actual files retained; numerical-equivalence criteria
preregistered separately, including explicit exclusion of wall_s diagnostics.
This does not relabel the earlier geometry-candidate all-byte failure.

Existing independent scalar/mechanical algorithms reused on byte-exact
immutable snapshots; path/test-count metadata only adapted. Audit protocol
freezes wrappers/core source hashes before audits; case seals record raw and
derived-output hashes. Processed text may undergo Git newline normalization;
native ZIP bytes are exact. 56NEW unit checks executed once, not rerun by audits.
No new task>2h in this campaign; future long pilot needs a half-hour monitor.
Original C1 scientific acceptance and unresolved evidence remain unchanged.
