# Frozen-state native cost profile, not a new solver or convergence test

Registered2026-11-13before measurements. Motivation: actual uniform20->10
pilot passes10us time gate but fails space49.56%dh/76.12%motion. Original
50/200us and C1 goals remain outstanding. Longer jobs must not be launched
without measured cost/accuracy review and half-hour monitoring when>2h.

New standalone thin_film CLI, no existing native header/solver modification.
Read frozen uniform20/10, initial and10us stage2 states atdt.03125us.
SDS6mM,insoluble,gravity/DLVO off, all original SI conversions, boundary heights
and potentials from snapshot. Input round-trip area/energy<=1e-9. Import guards
and standalone/coupled single-step mechanical endpoint comparison check that
the probe uses the intended state/model; they are not full C1 validation.

Measure geometry, energy-gradient, metric construction plus one value,
full problem gradient, Hessian-vector action. One warm-up and3samples each.
Separately time ONE native mechanics step, ONE coupled step and ONE saved-step
observer (no warm-up for these). No full trajectory replay or original-file
overwrite. Standalone/coupled endpoints agree<=1e-10local spacing. Require
native solve convergence,force ratio<=1.000001,finite timings. Samples are not
additive: functions overlap and full steps invoke kernels multiple times.
No percentage-of-whole claim without actual call attribution.

Wall time uses std::chrono::steady_clock; process kernel+user CPU time uses
Windows GetProcessTimes (100ns units; actual accounting resolution can be much
coarser). Primary API reference:
https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-getprocesstimes
Small CPU samples can quantize to zero. Host load and warm caches affect data;
this is a diagnostic, not a controlled whole-trajectory speedup claim.

Per case120s,total native600s,4cases sequential, no >2h task. Preserve failure
and raw inputs/output/logs; single bounded ZIP<5MB. White Nature red/bluePNG
and vectorPDF. Freeze input/source/binary hashes before the profile. Verify
the existing tf_curved_edge_bench binary remains byte-identical before/after
building/running the new target. Any optimization is a separate later step,
requires original physical/conservation/accuracy controls, and cannot change
materials,shape forces,acceptance thresholds or the failed findings.
