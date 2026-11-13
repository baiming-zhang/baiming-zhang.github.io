# Conditional DRIFT-local nonzero Hessian assembly, before candidate execution

Original profile completed4/4 atd755b45. N2112 initial/mechanical endpoint:
Fields(false).0730/.0706s;Fields(true)16.423/16.566s;matrix densities.4605/.5490%.
CachedH*v~1.0ms. This establishes a costly assembly operation, not its exact
share in the actual103s step. It warrants testing ALL local nonzero terms.

Implement separate CurvedEdgeLocalAssembly in new thin_film header. Delegate
Value/Pressure/first fields to ORIGINAL CurvedEdgeEnergy. For second fields
reconstruct identical area derivatives on same native face order, assemble
outer products ONLY on exactly nonzero support, no drop threshold. Preserve
all diagonal disjoining terms,curved cross-face dilution/capillary terms,
signed bath counter convention and gravity. Keep full dense output/cache and
generic H1 solver unchanged: assembly change only, not a sparse solve backend.
Default native source/binary unchanged; candidate selected in separate wrappers.
Classical local/sparse assembly prior art, not new thin-film physics.

Operator comparison on ALL4original frozen inputs, not a selected easier set:
energy,area,pressure,tension,node energy exactly equal; ALL Hessian entries,
diagonal and3deterministic actions relativeLinf<=2e-13; symmetry<=2e-13.
Use common max|reference| normalization,floor1e-30. No omitted small entries.
Compare timings in same process,3samples,alternate order; original full matrix
pre-evaluation serves its warmup, candidate warmup likewise. Record all samples.
Per case240s,total original/candidate operator budget600s.

New manufactured test5x4 from existing curved-edge fixture; variants:
curved bath,flat closed,signed bath,nonzero disjoining+gravity,steep/high coverage.
Manufactured disjoining parameters(code units only):eps.007,hb.7,c8.003,c2.002,
beta.01,lambda.6,rho_g.04; not new physical parameters or production settings.
Each variant tests first-field equality,full Hessian,diagonal,3actions,bath
columns,cache invalidation after V perturbation. Also invalid physical V/domain
guards. At least32checks, zero failures. Candidate wrappers of existing
curved14,preparation4,mechanical-observer6 must pass; baseline tests not rerun.

Only if ALL operator/unit controls pass, run candidate actual cases:
- original uniform10,T10us,dt.0625us (global_space u10_dt00625);
- original uniform5,one-step dt.0625 and.03125us (global_five guards).
Use EXACT native arguments,full11-file set and baseline raw SHA verification.
Numerical equivalence is preregistered rather than demanding equal wall_s:
  mesh/static parameter fields exact; prep physical fields relative<=1e-10,
  counts/indices/masks and saved row structure exact; prep residual must pass
  original1e-8(its diagnostic iteration count recorded, not forced identical).
  For initial/states:coords maxdiff/dx<=1e-10;velocity vector component errors
  /common max reference speed<=1e-8;other physical columns relativeLinf<=1e-10.
  Faces:indices exact,c<=1e-10,F<=1e-8. Summary physical fields<=1e-10 except
  max displacement<=1e-8;exclude wall_s and numerical residual diagnostics from
  relative-equality test, but ALL original force/energy/inventory/prep gates
  enforced unchanged. No physical time2% or space10% gate changed.
  Keep raw byte-match flags for ALL11files separately, not a new all-byte claim.
  Force/region-work rows remain covered by original independent SI/GEOS audit.
Actual cap900s for10us,300s per5um start,total actual1800s. Candidate unit runs
120s each, no>2h batch. Independent saved-step audits cap300s per actualcase,
existing formulas/tolerances on immutable copies. Preserve all failures/skips,
no retries or gate tuning. All phases together hard declared native<3600s and
audits<900s(<2h); no>2h task launched. Further long jobs need half-hour monitor.

Original space failure49.56%/76.12%,10/50/200us scope,physical flux/traction,
exchange coefficient/thin-element ratio and experiment validation remain open.
Known four regular area-FD failures and prior11-byte timing-only failure retain
their original labels. Report measured candidate wall costs separately from
historical baseline; no controlled full-trajectory speedup ratio without rerun.
