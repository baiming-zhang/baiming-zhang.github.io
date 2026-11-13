# Prior art and inference limits

Provenance clarification: the frozen plan's phrase "five native headers/bench"
is a counting typo. Five source paths in the reference protocol comprise
three native files (bench plus two headers), runner and README. The phase
runner checks those three native hashes, all four binary hashes, and its own
runner/README hashes. Binary SHA256 covers the actual executed native program;
the source commit identifies the full source tree. No additional five-header
audit was performed. Original protocol and frozen plan are retained verbatim.

Primary papers downloaded locally with SHA256 manifest in
03_DRIFT/literature/layout_20261113/; no paper PDF republished.

- Duffell & MacFadyen2015, [High-Frequency Voronoi Noise Reduced by Smoothed
  Mesh Motion](https://arxiv.org/abs/1407.7300): moving-mesh hydro can produce
  short-scale noise associated with mesh motion/topological changes. They
  smooth mesh-generator velocities. This is an ALE hydro technique; it does
  not authorize smoothing our material-particle velocities or forces.
- Mocz et al.2015, [Reducing noise in moving-grid codes with strongly-centroidal
  Lloyd mesh regularization](https://arxiv.org/abs/1503.07868): mesh centering
  improves their hydro cases; generic Lloyd regularization is established
  prior art, not our new contribution.
- Steinberg, Yalinewich & Sari2016, [Grid Noise in Moving Mesh Codes: Fixing The
  Volume Inconsistency Problem](https://arxiv.org/abs/1603.06753): distinguishes
  geometric volume changes from approximate hydro flux volumes. Their remedy
  belongs to their Godunov/ALE system. Exact global inventory or geometry
  identities alone do not establish continuum accuracy in another solver.

These papers motivate separating perturbation sensitivity, transport
consistency and reconstruction error; none proves our discrepancy has the same
cause. Native DRIFT site positions carry material inventories and momentum.
Changing them for mesh regularity needs consistent V/S/momentum transfer and
its own work/energy argument. Do not add a smoother or residual flux correction
to obtain a desired pattern.

Kim & Mandre2017 [Marangoni elasticity of flowing soap films](https://doi.org/10.1103/PhysRevFluids.2.082001) gives a related long-wave speed
relation. The current model EOS implies c~6.98m/s,320um crossing~45.8us;
this is a model scale, not measured data or proof of wave/reflection causation.
Their experimental elasticity is not substituted for this model's chemistry.

Current frozen evidence: background40/30/25 vs uniform20 displacement RMS
differences43.15/43.89/43.99nm against reference45.79nm at200us. Bulk x<200um
accounts for65.05/66.27/67.16%of the squared differences. These are exact
partitions of the original overlap diagnostic, not replacement acceptance
metrics. New same-dx phase runs must finish and pass their own time pair before
interpreting their layout sensitivity; even then global dx convergence and
physical normal flux remain separate requirements.
