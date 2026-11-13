# Actual moving-film background resolution: frozen plan

Keep fine spacing20um, numerical fine start200um, physical cut240um, horizon200us,
same SDS6mM variational preparation and all material/force/flux/Newton settings.
New optional background_ratio changes only numerical site spacing (default2).
Nominal background40->30->25um is rounded to integer cell counts: report actual
x/y spacings from mesh.json. Uniform20is a comparator, not continuum truth.

Guards:24native checks; rebuilt default dt.0625us200us control reproduces nine
frozen time_refinement/dt00625files byte-exact. Omitted/explicit ratio2one-step
files also byte-exact. One-step bg1.5/bg1.25/uniform coarse/fine guards.
Conditional200us matrix: uniform20dt.0625 and.03125us, bg1.5and1.25dt.0625us.
Read existing local-bg2dt.0625/.03125results; do not rerun other frozen batches.

Original overlap-weighted relative dh and displacement<10%SPACE at10/50/200us.
Uniform comparator must independently satisfy original2%TIME pair at all three
times; local-bg2time pair already passed. New background layouts have no own
time pair in this first pilot; identify this limitation even if space passes.
Independent SI/GEOS saved-step formula/mechanics audits and physical-cut
recovery remain separate from native all-step ledgers and force statistics.

Cost estimate uses frozen u20_dt1full200us54.657s, linear step scaling times1.3.
Background layouts have fewer sites than uniform20. Engineering estimate only;
each native case<=2400s, whole<=6600s (under2h). Stop/skip by frozen cost budget,
preserve every partial failure. No >2h job and no new kinetics/fitted rate.
Do not claim global dx convergence from background-ratio agreement; uniform10
and actual coupled physical normal flux/experiment remain open.
