# Curved-area / full-state moving-edge coupling (2026-11-13)

This opt-in native prototype replaces the planar-surface + quadratic thickness
energy by a COMMON curved-area / Langmuir energy for mechanics and drainage.
The local edge is fully evolved on a conforming fine Voronoi patch; no reduced
basis, prescribed interior profile, empirical drainage law or adsorption rate.
It is not yet a deep-meniscus Stokes solver or a complete regional mechanical
traction closure. All previous modes/results are unchanged.

## Discrete energy and prior art

For each dual face, d=site distance, l=face length, z=(h_i-h_j)/(2d), and
e_ij=l*d*(sqrt(1+z*z)-1). Each end receives e_ij/2 additional surface area.
A^s_i=A_i+sum e_ij/2, Gamma_i=S_i/A^s_i (one face); two faces in free energy.
This is a DIRECTIONAL DUAL-FACE AREA QUADRATURE, exact for a Cartesian 1D
section and matching the inherited quadratic energy in the small-slope limit.
It is not the exact reconstructed 2D surface graph area at finite slope.

Energy =2sum[A^s_i*KL(Gamma_i;meanGamma)+gamma_ref*(A^s_i-A_i)]
+disjoining+gravity+sum_bath mu_bath V_bath. meanGamma=sum S/sum projected A
is constant on the fixed domain. gamma_ref=gamma_clean-Pi(meanGamma), so
dE/dA^s=2[gamma_clean-Pi(Gamma_i)]. Clean tension is the existing measured
FilmMaterial.gamma0; aspect H/L comes from FilmScales. No new material input.
Hessian includes cross-face surface-dilution terms. Reverse geometry derivative
uses the same energy. Physical inertia / actual-wall face-pyramid viscosity and
lagged h^3/(12mu) relative drainage are inherited long-wave approximations.

[Thiele/Archer/Pismen2016](https://arxiv.org/abs/1609.00946) distinguishes
curved-area surfactant concentration from the projected conserved quantity and
derives coupled free-energy dynamics. This is established physics, NOT a new
thermodynamic principle. Candidate future method value is the stateful local
edge / moving-film exchange with verified accuracy, powers and costs.

## Boundary and initialization

320x160um domain, local cut x=240um, initial INTERIOR h=6um. All interior
h evolve from their V/A; fixed boundary sites lie exactly at x=320um for ALL
meshes, avoiding a mesh-dependent pressure boundary location. Specified radius
r=500/1000um sets ONLY far boundary height h0+2[r-sqrt(r^2-Ledge^2)] and
mu=-gamma_ref/r. They are two boundary data for the fourth-order film model,
not an imposed evolving meniscus shape. Boundary node V is bath bookkeeping,
NOT a resolved finite-meniscus volume and NOT its height. No bulk S exchange.
Initial S=theta0*Gamma_inf*A_projected with optional .005 projected coverage
perturbation; hence actual curved coverage is initially diluted near the edge.
This deliberately nonequilibrium initialization is documented, not calibrated.

## Frozen SI matrix (before calculation)

All T=.5us except guard T=.25us. Cases:
1. guard local40um dt.125us r500um amplitude0;
2-4. local40um dt.125/.0625/.03125us r500 amplitude.005;
5. uniform40um dt.0625us;
6-7. local20um / uniform20um dt.03125us;
8. local40um dt.0625us r1000um.
Each <=180s; whole batch<=1800s. Stop on failed startup guard; preserve all
failures; no second tolerance repair to make matrix pass. No >2h training job.

11 native checks: pressure/Hessian/geometry derivatives, bath independence,
actual surface dilution, flat equilibrium, small-slope limit, actual coupled
motion, liquid/S conservation, shared face flux, global dissipative inequality.
Initial development test used an unjustified 0.1% dilution assertion on a small
fixture; it failed (actual dilution .0738%). Replaced by a resolvable nonzero
dilution check, before SI calculation. Original failure log retained.

## Frozen acceptance and independent audit

Native accepted horizons, positive V and 0<theta<1; total V (including bath
bookkeeping), material S and regional V/S reconstruction residual<=1e-10.
Grand-potential + kinetic + viscous/drainage/upwind-loss excess<=1e-7.
Independent GEOS projected geometry/ownership <=1e-8; independent curved area,
pressure and saved energy <=1e-7 of corresponding complete-history scales.
Drainage REGIONAL free-energy balance: endpoint pressure work + shape work +
physical dissipation + explicit convex implicit-step losses, <=1e-7 (scale=sum
absolute step energies/works, floor1e-26J). Complementary exchanges cancel.
This does NOT claim regional mechanics/traction power is already verified.

For fixed geometry drainage, endpoint shape term on face i->j is
(w_j-w_i)*e'_ij*(gamma_j dh_i+gamma_i dh_j), with e'=l*z/(2sqrt(1+z*z));
dh_bath=0. Regional dissipation uses (w_i+w_j)/2 * F_ij*(mu_i-mu_j).
Pressure port = (w_j-w_i)*(mu_i+mu_j)/2 * F_ij.
Implicit losses are the sum of separately positive convex area-energy and
surface-dilution Bregman terms, derived from the model, never a residual source.
Report pressure-only ablation and common-state normal-flux/node residual.

Time relative L2 final dh=h-h_initial and displacement<=2%; spatial local/full
<=10%, via final-cell intersection weights; bulk background h cannot hide
dynamic failure. Local window errors may later be exploratory, not replacing
these frozen whole-domain gates. No experimental validation claim.

## Post-evidence extension, frozen before execution

Initial8native cases and8independent regional audits pass; all4dynamic
comparisons fail. The first .5us horizon has very little developed edge motion.
Extend to20us with unchanged material, far geometry and Newton tolerances.
Metric-only access was added to avoid evaluating/subtracting a cancelled FLAT
Langmuir energy inside the curved model; that unnecessarily ruled out valid
curved concentrations with S/A_projected>Gamma_inf. New12th native check tests
that domain explicitly. Common energy and inherited physical metric unchanged.
Two matched original .5us controls must have h relative L2<=1e-10 and dynamic
h-change/displacement relative L2<=1e-6 before interpreting this extension.

Longer matrix: original l40_t0625 and l20 controls; guard local40 dt.5us
T20us amp0; driven local40 dt.5/.25/.125us, uniform40 dt.25us, local20/uniform20
dt.125us, radius1000 local40 dt.25us. All driven amp.005.10executions,
per-case180s / whole1800s. Same conservation/power/field gates as original.
Failure is retained; do not reduce Newton tolerances or change viscosity.
Large original CSV files stay in results/local-runs/curved_edge_longer_20261113;
byte-preserving split ZIP archives will accompany the report.

## Frozen outcomes and reproduction

Original8/8 horizons complete. Longer10attempts:8complete including2matched
controls;2fine20um force/tangent evaluations stop after7.25us accepted prefixes.
All18accepted saved series/prefixes pass independent formulas; this does not
qualify the2incomplete horizons. Maximum original/qualified power residual
1.18e-13. Longer dt.25->.125us thickness-change.736%, displacement1.273% pass2%;
coarser time pair displacement2.393% fails. Space40um26.33/39.77% fails10%;
20um20us comparison missing. No continuum/exact2D-area/experiment claim.
A default22moving-check executable was rebuilt and rerun after metric access;
it passes. Curved12checks pass;34relevant checks, not45distinct checks.

Both failed runs have native logs and accepted CSV but no native summary after
an exception. Independent partial metadata uses observed IDs/steps, initial
saved node energy and an identified completed run's SAME material/boundary
constants; see longer/partial_attempts.json. Last saved coverages .999410/.999492
still have positive tension3.81/2.43mN/m. Exception is subsequent force/tangent
admissibility; it does NOT prove the continuum state crossed the EOS domain.
No concentration/tension clipping, material change or tolerance replay applied.

Original driver bytes match protocol hash in run_at_execution.py. Current
run.py permits >=11checks since12th check was added, but refuses overwriting
frozen outputs. To regenerate scripted matrices, select fresh O/RAW/LOG paths
in a copy of the driver before executing; full native commands/parameters and
source/binary hashes are in protocol/status JSON. A direct fresh native replay:

    tf_curved_edge_bench local 40 .25 20 500 .005 ABSOLUTE_NEW_OUTPUT_DIRECTORY

Run binaries with local xmake DLL paths and VF_THREADS=2 as recorded by run.py.
Analyze/report source scripts never synthesize trajectories. Large raw CSV
remains local; package.py splits into bounded ZIP chunks and verifies all70file
hashes. restore_raw.py REPORT_DIRECTORY NEW_OUTPUT_DIRECTORY reconstructs every
byte and refuses overwriting. Report includes PNG and vector PDF, native logs,
original/longer protocols and all failed attempts.

Remaining: actual regional mechanical traction/kinetic power, independently
resolved physical normal flux (current regional inventories are reconstructed),
space convergence and finite-slope2D area/mobility validity. Next use independent
SDS thermodynamics to separate initial-state/material-domain effects from
near-domain force/tangent evaluation. C2adsorption needs independent kinetics;
no fitted drainage exponent or arbitrary rate to cure this failure.
