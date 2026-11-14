# Same-candidate20/10um reference trajectories at the first original checkpoint

Admit ONLY after retained_motion's three actual single-step native/equivalence
and independent scalar/mechanical controls pass. Reuse its FROZEN compiled
tf_retained_face_bench; no rebuild, original/default change or old-trajectory
replay. New trajectories use the SAME retained-positive-face geometry, local
Hessian assembly, CHOLMOD/H1, SI material and original nonlinear tolerances.

Four fixed cases, in order: dx20um then10um, each dt.0625us and.03125us,
all T10us, N144/544, radius500um, amplitude0,SDS6mM insoluble,variational
preparation, internal cut240um, refine_start200um, observer true. Save first
step and every1.25us(stride20/40). Native gates unchanged: complete160/320
steps with failures0,force/target<=1.000001,allglobal/regionalV/S<=1e-10,
energy budget excess<=1e-7,preparation pressure<=1e-8. Stop further native
cases after a native failure; never restart or modify parameters automatically.

Per-native caps300/600/1200/3600s, whole native batch6600s; every native process
is held by Popen, timeout kills/waits that exact process. Native whole cap is
below2h. Parent emits its own half-hour handle/progress record; this is not a
new >2h training job. No5um trajectory launched here. Future >2h5um batch
requires a verified scheduled half-hour monitor BEFORE launch.

After all terminal, audit original independent SI/GEOS saved-step scalar and
mechanical formulas. Apply original moving-overlap metrics: time2% for EACH
resolution, spatial10% at BOTH matching dt. Preserve all failed rows. Existing
old20->10space failures remain historical evidence, not replaced or erased.
Spatial failure at these grids does not by itself eliminate the planned5um
refinement study; time/native failures must be diagnosed before escalation.

This establishes new SAME-implementation references before any new5um pair;
never mix old coarse/new fine implementations as a qualifying convergence pair.
Original10/50/200us, physical flux/traction, independent exchange/TFE predictions
and experiments remain open. T10us is the first checkpoint, not fullC1.
