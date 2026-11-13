# Actual frozen-state full curved-energy directional response

Registered2026-11-13before measurements. Original C1 full scope unchanged:
actual20->10 at10us still fails49.56%dh/76.12%motion;50/200us, independent
physical normal flux/traction, exchange coefficient, thin-element ratio and
experimental/global-fringe predictions remain open. No fitted forces or C2 rate.
Previous area-spectrum omitted curved-surface stiffness and actual boundaries.
Question: does weak short-wave restoring response persist in the FULL native
curved-energy geometry, and is it present in actual saved velocity projections?
This is a falsifiable directional diagnostic, not full coupled eigenanalysis,
trajectory convergence, a cause attribution or an innovation claim.

Use ORIGINAL native headers/functions, no candidate fast dispatch. Four frozen
cost_profile_20261113 inputs:N144/544,initial/10us,SDS6mM,insoluble,no gravity or
DLVO,320x160um,actual trace and prescribed-potential bath. Preserve SI materials,
V,S and bath heights/potentials during position perturbations. Require native
area/energy import agreement<=1e-9. Default source and production binary unchanged.

Mode q=(sin(m*pi*x/W),0), bath displacement exactly zero. Uniform20 modes
1,4,8,12,16; uniform10 modes1,4,8,12,16,24,32. Same physical m can be compared
across resolutions; m=16/32 are each nominal two-particle wavelengths. No claim
these trial directions are eigenmodes or span y-motion. Initial/t10 snapshots
are not mechanical equilibria. Keep signed stiffness, including negative values.

Let a=epsilon*dx in code length. epsilon fixed1e-3,1e-4,1e-5; all3 are saved.
K_full=q dot[g(x+a q)-g(x-a q)]/(2a), full native AD gradients.
Independent value check K_E=[E(x+aq)+E(x-aq)-2E(x)]/a^2.
Surface-rate Js q uses native curved control-surface area differences; intrinsic
surface-elastic term K_surface=sum 2 E_i/As_i*(Js_i q)^2,
E_i=RT Gamma_inf theta_i/(1-theta_i),theta_i=S_i/(Gamma_inf As_i).
Projected-area surrogate uses native exact J_A q in the SAME coefficient;
this is not a claimed additive block of the curved model. Full-minus-intrinsic
is the geometric/prestress residual for these gravity/DLVO-off fixed-inventory
states. It is not automatically positive, and is not viscosity.

Normalize stiffness by Mq*(c0*k)^2, Mq=sum_nonbath rho V_i |q_i|^2,
c0=sqrt(2 E0/(rho*6um)) using ORIGINAL material activity. Report K inN/m,
Mq inkg, c0 inm/s and k in1/m. No absolute-value square-root frequency shortcut.
Native viscous coefficient D=q^T Cq via MovingVoronoiProblem(full gradient+wall
samples),SI conversion S*T; record D/Mq(1/s) separately. Viscosity is not added
to free energy. No linearized drainage block or full eigenvalue calculation.
Mass projection of saved velocity onto q gives a one-direction kinetic fraction;
directions not orthogonal, so fractions do not add. Initial zero kinetic energy
is unavailable, not a successful zero-error dynamic check.

Predeclared qualification PER direction:
- At epsilon1e-4/1e-5, native pair topology must match unperturbed topology.
- K gradient stability and intrinsic K_surface stability each<=1% normalized by
  max(abs(finest term),1e-4*continuum K scale).
- At1e-4, value-vs-gradient K difference<=1% using same normalization.
If any fail, mark unresolved; no smaller eps/tolerance tuning this campaign.
No physical/time/space gate is modified by these diagnostic controls.

Two flat regular native fixtures n16x8/n32x16,no bath/trace,V=h0*A,S=theta0*A;
test m1,8,n. Independently known J_A q=dy*sin(k dx)*cos(k x_i).
For uniform flat h0 and uniform activity, full discrete stiffness is
K_flat=2E0 sum(Jq^2/A)+gamma0_equilibrium*aspect^2/2 *
 sum_native_pairs (face_length/distance)*(h'_i-h'_j)^2,
h'_i=-h0*Jq_i/A_i. Compare finest native gradient K against this analytic
discrete reference with abs difference/continuum K<=1e-5; exact J expression
normalized bydy<=1e-10. m=n has zero first-area response, hence vanishing flat
full stiffness as well. This guards the probe, not actual drainage accuracy.

4actual snapshot cases+2flat fixtures,30directions,total90epsilon rows.
Per CLI120s,total native600s hard cap, no >2h job. No successful frozen
trajectory is rerun. All input/output/logs/failures preserved with hashes,
bounded raw ZIP, white Nature red/bluePNG+vectorPDF. Report all inconclusive
directions. Results update all4research entry docs and own VariDrainage page.

Relevant prior art (local PDFs/text already retained in literature/layout_20261113):
Hess&Springel2010,section5.3,https://arxiv.org/abs/0912.0629;
Duffell&MacFadyen2015,https://arxiv.org/abs/1407.7300.
Their mesh regularization forces/velocities are not adopted here.
