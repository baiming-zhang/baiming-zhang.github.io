# Prior art and next physical diagnostic

Primary sources checked2026-11-13. These studies do not validate our native model.

- Alharbi & Naire2017, [adaptive thin-film moving mesh](https://doi.org/10.1016/j.cam.2017.01.019), and2019 [two-dimensional extension](https://keele-repository.worktribe.com/output/412693), DOI10.1016/j.cam.2019.02.010. Their mesh-density/MMPDE approach resolves multiple structures and moving contact-line/fingering features in substrate-film settings. This is prior art against generic local/adaptive thin-film refinement. The2019primary repository abstract was available via search; direct full-page/PDF request returned403. Full-PDF review is not claimed.
- Winchenbach & Kolb, [Optimized Refinement for Spatially Adaptive SPH](https://doi.org/10.1145/3363555), TOG40(1),8,2021. The publisher abstract reports refinement-pattern optimization and residual refinement error. Its kernel-SPH framework differs from our moving Voronoi model. It motivates measuring operator/refinement consistency beyond covering the visible thin region. Author PDF front matter says2020; publication record says2021. Author download redirected to HTML; no PDF was saved.
- Chen et al2021, [Momentum-Conserving Implicit MPM for Surface Energies with Spatial Gradients](https://arxiv.org/abs/2101.12408). Surface-energy-based treatment already handles spatial energy variation/Marangoni effects; the study also addresses linear/angular momentum in particle resampling. Generic energy-based capillary coupling is not our novelty. Original26.94MB PDF and extracted text are local under03_DRIFT/literature/refinement_20261113; no full paper is rehosted here. pdftotext reported an xref-reconstruction warning; abstract text was inspected.

Our inference: extending coverage helps thickness, but remaining motion and
material-cut discrepancies need mechanical consistency/transport diagnosis,
without an unexplained material-rate or exponent fit. Candidate research value
still requires a verified stateful moving-film/local-edge coupling that predicts
physical observables at useful cost. No novel contribution is established yet.

Next bounded diagnostic should observe actual mechanical initial/final force
residuals, inertia/viscous actions and regional power. Check the observer leaves
accepted states byte-identical. Distinguish nonlinear solve accuracy, split-step
error and nonconforming cut reconstruction before changing any solver parameter.
Keep the original fixed physical interface, materials, normal-flux and dynamic
convergence requirements. C2 needs independent physical kinetic data.
