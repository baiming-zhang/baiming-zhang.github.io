# Flat-film area-response spectrum: frozen observer diagnostic

2026-11-13. No native evolution, model change or C1 gate replacement.
Motivation: Hess & Springel2010, section5.3 and Fig5, documents loss of the
shortest-wave pressure restoring response in plain Voronoi particle hydro:
https://arxiv.org/abs/0912.0629 . PDF saved locally, page11visually checked.
This is prior art, not an innovation claim or proof of our moving-film cause.

Use the existing tf_voronoi_kinematics read-only CLI, native rectangular box
320x160um. Cell-centred n x n/2 grids, n16/32; regular and deterministic1.5%dx
jitter matching the native assay function. No trace or bath nodes: this is a
flat geometric/linearized bulk diagnostic, not the actual drainage geometry.
Prescribe u_x=U sin(m pi x/W),u_y=0,U=1mm/s, zero normal velocity at box walls.
For each n use unique modes1,2,4,n/4,n/2,3n/4,n-1,n.

On regular grids verify native Jv against the exact cell-area expression
Adot_i=U dy sin(k dx) cos(k x_i), k=m pi/W. Normalized maximum discrepancy
<=1e-10, normalization U dy (not near-zero expected rate at Nyquist).
Area partition and summed Adot checked separately, unchanged native sites.
No cancellation-sensitive finite-difference area oracle is added here.

For flat constant coverage, the area-thermodynamic Hessian block is
J^T diag(2E/A)J, E=RT Gamma_inf theta/(1-theta). Its Rayleigh frequency scale
is sqrt[sum(2E Adot^2/A)/sum(rho h A |u|^2)]. Compare to c*k where
c=sqrt(2E/(rho h)). This is ONLY the area-elastic block: omits curved-surface
gradient stiffness, viscosity and actual evolving boundary/shape fields.
Do not call it a measured full-mode frequency or a drainage stability proof.
Long-wave agreement and near-Nyquist loss are diagnostics, not pass criteria
for the full physical model. On regular grids m=n gives zero first area
response; uniform surface tension and fixed total area also remove the
baseline linear-tension Hessian term. Report finite-precision residuals.

Use material constants exclusively from frozen background-space uniform20
summary, keep source SHA256. Read raw prescribed inputs/outputs, retain full
logs and hashes. Native geometry CLI timeout15s/case,total120s. No >2h task.
Keep every unsuccessful observation, bounded ZIP<5MB for raw bytes, white
Nature red/blue PNG+vectorPDF. Do not add the literature's shape-correction
forces, change elasticity/viscosity, fit dispersion or loosen actual-film gates.
