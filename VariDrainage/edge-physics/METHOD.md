# C1 material / conserved preparation / refinement protocol (2026-11-13)

Motivation: previous20um20us runs stop after7.25us on a curved Langmuir/tension
evaluation. Last saved states are admissible; no proven continuum-domain exit.
Keep original failed cases. Do not clip Gamma/tension or change Newton tolerances.

## Independent physics and frozen inputs

[Blokhuis2025](https://arxiv.org/abs/2509.20110), Table3 (p18 of preprint), gives
the existing pure SDS/water isotherm/mass-action parameters Gamma_inf3.70e-6,
1/K2=2.16(mM)^2,r=.75,m=65,x0=2.7e-7. These are fixed; no new fit. Concentration
1.94/6/19.4mM are deliberately selected INPUT controls, not measurements of this
device or selected by a drainage exponent. Equilibrium is used ONLY for initial
coverage/tension; evolution stays insoluble. Mixture viscosity1.5mPa.s and
density1045kg/m3 remain inherited assumptions. Applying pure-water surface
parameters to the glycerol mixture is an unvalidated material approximation.
No adsorption/desorption or hidden kinetic coefficient is introduced.

Preparation projected: original S_i=theta0 Gamma_inf A_i (amp0 in main matrix).
Preparation balanced: weights=S_i*A^s_i/A_i, S'_i=(sum S)/(sum weights)*weights.
This positive redistribution preserves total S, V and x; actual curved coverage
is uniform for amp0. It changes INITIAL spatial distribution, not a transport
repair during evolution and not proof of preparation-independent prediction.
Both use the SAME flat initial interior h and same far boundary height/pressure.
The h jump at far boundary remains an initialization limitation; this test does
not by itself establish grid-independent initial surface geometry.

## Frozen executions before any new SI run

Native default control local40 dt.25us T20us c19.4 amp.005: frames/faces/history
must be byte-identical to archived long_l40_t025. Two original failed local20 /
uniform20 dt.125us T20us c19.4 amp.005 diagnostic replays add stage/node/theta
to failure messages and retain a native prefix summary; no numerical change.
Then 3concentrations x2preparations x2meshes(local20/uniform20), all dt.125us,
T20us, amp0, radius500um =12runs.15total primary executions; each<=180s,
whole primary<=3300s. Preserve all failures, stop on default-control mismatch.

Conditional refinement was chosen BEFORE primary results: if BOTH c6 balanced
local20/uniform20 complete and independent gates pass, run uniform10 dt.03125us
one-step guard. If guard passes, run c6 balanced at T5us dt.03125us on
local20/uniform20/local10/uniform10; additionally uniform20 T20us dt.0625us.
Each refinement<=600s, whole refinement<=3300s. No run longer than2h.
This is a predeclared material regime for diagnosis, not a replacement of the
original c19.4 failed condition. Publish all inputs, attempts and missing gates.

## Frozen gates

12curved native checks first. Global/regional V/S balances<=1e-10; initial
projected/balanced total S difference<=1e-12, all Spositive, actual balanced
theta uniform (amp0) within1e-10. Independent material isotherm/species audit
relative<=1e-10. Same curved-area/pressure/energy/GEOS/drainage-power audits and
1e-7 complete-history power gate as previous curved_edge/README.md; independent
geometry<=1e-8. Accepted prefixes never qualify incomplete requested horizons.
Time h-change/displacement relative L2<=2%; spatial<=10% using exact final-cell
overlap. Compare ONLY matched input, preparation, timestep and completed horizon.
Local-vs-uniform and uniform20-vs-uniform10 both reported; boundary and bulk
background h must not hide dynamic error. No tuning of viscosity or noise.

New diagnostics label motion gradient/Hessian +/- and failing node/coverage;
the native benchmark catches failure only to save metadata and return failure.
All original defaults/force/flux formulas unchanged. Failure diagnosis does not
assert that the last accepted state or true continuum solution was invalid.
Full C1 still needs spatial validation, independent physical normal flux and
regional mechanical traction/kinetic power; current regional S is reconstructed.

## Bounded diagnostic closeout after primary evidence

Primary15attempts: only default control completes. All12material/preparation
cases stop early. All15saved accepted-prefix formula audits pass. Independent
isotherm3/3 and conserved-preparation6/6 checks pass; frozen refinement condition
fails, so no10um refinement runs. No successful material/preparation regime was
selected to replace the baseline. First metadata serialization failed on a NumPy
bool; cast to native bool and reused unchanged saved formula audits, no native
replay. Root ALIGNMENT.md restores original physical evaluation priorities.

Exactly TWO observational replays, c6balanced local20/uniform20, same20us,
.125us dt and all materials/tolerances. New failed_states/failed_faces outputs
separate the rejected mechanical/drainage candidate from accepted histories.
Accepted histories must remain byte-exact. Candidate diagnostics measure per-node
implicit continuity<=1e-10, face constitutive flux<=1e-7 and energy excess<=1e-7;
these are OBSERVATIONS, not a replacement Newton convergence or accepted steps.
No additional solver-tolerance loop. Then inspect inherited framed-film evidence
and return to edge-pinching, exchange and fringe/drainage observables.

## Signed infinite-bath counter protocol, frozen before correction runs

The two observed rejected candidates never moved V (F=0), despite constitutive
flux needing about7.1e-18m3 per leading face; film-only continuity error.36%.
Their boundary counters were1.1-5.6e-33m3, whereas physical film cells stayed
>2.1e-15m3 and tension>40mN/m. Source inspection confirms a positivity condition
on the declared INFINITE bath's nonphysical exchange counter. It was therefore
acting like a finite, mesh-dependent inventory. Do not accept rejected states.

Correct only CurvedEdgeEnergy/coupling: bath counter may be signed, film V/S
remain positive, and nonphysical bath-bath links have zero conductance. Same
energy derivatives and physical parameters.14native checks include signed-offset
energy/pressure invariance and identical physical trajectories for offset counters.
Frozen5runs: existing default_control, c6balanced local/uniform20, and the two
original c19.4projected amp.005 failures; same dt and20us. Each<=180s, whole900s.
Default-control h difference<=1e-10, dh/displacement<=1e-6. All old failures stay.
No original10um refinement: its old prerequisites remain failed. This correction
is not proof that the original initial-boundary problem models drainage: ALL15
original accepted prefixes have net BATH-TO-FILM exchange. Physical flow direction
and a compatible prepared edge must be reviewed before expanding computation.

## Return to physical observable definitions

Counter-corrected5attempts: default and both c6balanced cases complete20us;
original c19.4Hessian probes still stop after7.25us. All5saved series/prefix
formula audits pass. Space c6local20/uniform20 dh20.62%, displacement61.56% FAIL.
All net flows still replenish the film. Do not enlarge this setup into a drainage
validation simply because bath bookkeeping no longer stops the run.

Inherited6s CSV inspection is not a solver rerun. Reproduce fixed-height speed
proxies, then separately attempt Monier's peak-speed/tmax operator with fixed
5/7/9point quadratic derivative windows, exclude censored end maxima, require
>=4valid fringe peaks, affine tmax-vmax R2>=.95 and collapse RMS<=.05H. These
operational gates are diagnostic choices, not fitted model inputs. None of the
windows supplies >=4uncensored peaks (3each); paper n is therefore unqualified.
Raw13.6MB particle frames and the named legacy stdout are absent from current
Git; no pretend reconstruction of side pinching or exchange from scalar CSV.
The legacy4x5mm/D250um setup has no exact tuple in Monier TableI.

One NEW GEOMETRY-MATCHED inherited baseline startup is frozen before execution:
Monier TableI W16mm,H5mm,D250um,mu1.5mPa.s; SDS19.4mM,dx100um,h0=2um at rest,
r_ref=.5mm,dt5ms,T10ms,two steps,report/save every5ms.180s hard timeout. It is
NOT the new curved model, a fringe-exponent measurement, resolved marginal
regeneration or a full material match (model gamma38.6mN/m vs approx35reported).
Purpose: collect actual physical setup profiles/edge observables and demonstrate
which measurements need finer normal resolution. No >2h job is launched.


The geometry-matched inherited startup completed two steps in93.55s, with no
mechanics/drainage Newton failures. One OBSERVER replay is authorized by this
protocol to collect the missing original-geometry interface ledger: same inputs,
only VF_INTERFACE_CSV added,180s maximum. Compare rounded particle frame bytes
to the unobserved run; report any difference. The first frames have no particle
areas, so liquid/solute conservation cannot be independently inferred from them;
this additional ledger is needed for the original physical comparison.
Actual discrete inner height is5.022945mm (hex-lattice rounding), not exactly5mm.
No third identical startup or long run. Interpret as NOMINAL geometry matching,
not full experimental matching, finite-resolution pinching or fringe-law proof.

## Completed observation and physical report

The single observer replay completed in93.079s; frames byte-identical, native
ledger residual1.7197187491850097e-13 and180material crossings. See matched_audit
summary for separate material2S/B and remesh-region changes. These are bookkeeping
shares, not experimental mechanism fractions. No additional native replay.

matched_analyze.py reads only true saved frames and existing native ledger;
package.py deduplicates raw byte chunks, then report.py makes the white report
and4PNG/vector-PDF pairs. restore_raw.py recovered103files byte-exact from56parts.
Raw total218623793bytes, preserved including rejected and failed observations.
Local desktop/mobile browser checks pass without horizontal overflow. Full C1
remains incomplete; PHYSICAL_GATES.md freezes next physical decision priorities.
