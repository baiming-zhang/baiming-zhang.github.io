# Conditional candidate integration controls (registered before execution)

First four frozen states all byte-exact. Measured N544 clipping medians:
64.683->27.856ms(initial),66.130->28.557ms(t10). Both exceed preregistered1.25x
screen. Full original geometry71.487/75.312ms, so integration is worth testing.
These are clipping timings, not a demonstrated whole-solver speedup.

Use test-only translation units: include unmodified OverlapTransfer primitives,
then select candidate at downstream native Transfer_Voronoi call sites with a
locally scoped macro. Include the original bench/test source, undefine macro.
This prevents copies of native physics/solver code and preserves default targets.
No default solver adoption in this campaign. Candidate binary distinct.

Before actual trajectory: run candidate versions of existing curved-edge,
preparation, mechanical-observer and transfer controls; each must return zero,
report zero failed and at least original14/4/6/14 passed. These are NEW candidate
executions, not claims of re-executing baseline tests. Retain any failures.
Run candidate frozen-state profile on the same4inputs. Require exact equality
of every non-timing scalar with the existing baseline profile, including energy
import error, Hessian/gradient-dependent sink, AD tape size, Newton/CG counts,
force ratio and mechanical endpoint difference. Import must meet old1e-9 gate.

Conditional on ALL controls passing: run ONE actual uniform10,N544,T10us,
dt.0625us candidate. Use EXACT arguments/save stride from completed global
u10_dt00625. Reference original raw files verified by original SHA256, not rerun.
Require exactly same set/names of11native files and byte-exact contents, including
history, actual states, mechanical outputs and summary. Native force/energy/
inventory gates therefore must agree. Any byte mismatch fails the exactness
claim; retain and investigate without loosening acceptance thresholds.

Pre-execution inventory correction: initial plan said9files; inspecting the
authoritative frozen reference shows11(including mechanical observer outputs).
All11must match, none excluded; corrected before candidate execution.

Per test120s; frozen probes120s each,total600s; actual trajectory1800s;
total native budget2400s hard cap(<2h). No task>2h, no monitor required for this
bounded experiment. Preserve raw bytes and all failed/partial attempts, build
logs and frozen source/binary hashes. Original production binary76e49... unchanged.
Timing versus an older run is NOT a controlled whole-trajectory speedup; only
state-equivalence and operational cost are established by this control.
Original actual space failures, long checkpoints, physical flux/traction and
full C1 predictions remain outstanding. No C2 kinetics or empirical force.

## Preserved first driver stop and count correction

Initial driver accidentally assigned curved/observer counts6/14 rather than
14/6. All actual candidate binaries returned zero failures(14/4/6/14);
the script falsely marked the observer incomplete and skipped all profile and
trajectory work. Original status.json/protocol.json remain unchanged.
Authoritative baseline layout_phase status lists observer,preparation,curved;
its original logs and current source establish6,4,14 respectively. Transfer14.
Continuation verifies existing candidate log/binary hashes and reuses these
four executions; no native test rerun. New continuation_protocol.json and
continued_status.json record the repaired driver, all subsequent probes and
conditional trajectory. The guard-inventory bug is not a physics test failure.
Caps apply to new execution as above; prior tests consumed<1s, total<2h.
