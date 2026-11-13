# Film-side coverage, fixed physical interface

Previous200us native data have thinnest nodes at221-230um, upstream of the old
fine-patch start240um. Test whether covering this film-side zone improves C1
dynamic fields and physical-cut recovery. This is a numerical extent choice,
NOT a new fitted physical parameter or experimental validation.

Run `python thin_film/scripts/film_side/run.py`, inspect controls/guards, then
`--matrix` with unchanged source/binary hashes. Build tf_curved_edge_bench first.
The optional final CLI argument refine_start_um changes local mesh placement
only. Defaults are unchanged; uniform ignores it. Physical trace and bookkeeping
cut stay240um; far bath320um, radius500um, h_bulk6um, SDS6mM, initial activity,
curved energy, viscosity/inertia, force, flux and all tolerances stay fixed.
Initial discrete V/S can differ with quadrature; report them rather than claim
identical initial inventories. No imposed evolving thickness/kinetic rate law.

Two choices fixed before new runs: start200um (40um film-side coverage) and160um
(80um coverage); thin-point location from prior data motivates these numerical
diagnostics. Compare the same20um fine spacing. This is not a search for a
physical rate or exponent. The boundary trace is released after initialization.

Controls repeat original local240/uniform20, dt1us,T200us; their states, faces,
scalar history and initial files must be BYTE-EXACT.4preparation+14curved checks
and two1us extent guards must pass before the new matrix. Matrix:200um dt1us,
160um dt1us,200um dt.5us and200um dt.25us, all200us. Save full three-stage
states at first step, every10us and final; retain all native scalar rows.
Guard wall extrapolation must be<=600s per case, whole computation<=3300s;
no task>2h is planned. Preserve failed prefixes; no tolerance/material retries.

Use identical independent SI/GEOS saved-step formula audits and fixed positive
10/50/200us comparisons. Dynamic thickness-change/displacement relative L2:
space<=10%, time<=2%, exact final-cell overlap, bath cells excluded. Incomplete
horizons never qualify. Recovered gradient/recorded-velocity cut flux remains a
NONCONFORMING diagnostic, not exact normal-flux verification; include all parts.
Regional mechanical traction/kinetic power and experiment prediction remain open.
Do not equate film-to-edge transfer with net draining out of the far bath.
