# Build observations

Initial profiling-target build at6d9ff71failed before any native measurement:
missing explicit `<filesystem>` include (MSVC C2079/C2039 at line19). Also
reported NOMINMAX redefinition. Existing tf_curved_edge_bench binary remains
76e49ed253e2341df37f93d02803272d020b8500ed9265da6f7d0ec124b6c5c3.

Repair: explicit standard header; guard Windows macro definitions and API
include. Non-Windows builds use std::clock CPU fallback, not tested here.
This is a tooling compile failure, not a simulation/physical-model failure.
No native result was produced or overwritten by the failed build.

Repaired build succeeds (16.0s); build.log preserved. Original native binary SHA256 verified unchanged again. Windows CPU timer tested in the next profile; non-Windows fallback remains unverified.
