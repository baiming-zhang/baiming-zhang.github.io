# FROTH process, topology and geometry addendum — 2026-11-20

This is a reporting and evidence addendum, not a new solver run or a goal resume.
Active branch: `p6-froth-foam-drainage-v2`. Continue to read the root README and
`FROTH_2026-11-19_Next.md`; the numerical evidence remains the frozen November13 stage.

## Actual process

1. Select the actual near-planar seven-bubble shared film [2,4], max plane error1.70um.
   The initial hole is prescribed; nucleation and conservation of the initial hole
   liquid are not demonstrated.
2. Run the existing CODA constant-tension planar first-contact diagnostic, stopping
   at270.716us after527steps. The133saved rim states are now rendered, without
   interpolation. Gas and cluster geometry in this movie are frozen references.
3. Attach a volume-exact but uncalibrated tetrahedral receiving cap on a selected
   liquid-skin facet. This is a declared preparation operation, not a resolved
   collision/topology evolution. Its footprint length is14.45rim radii.
4. Fit incoming linear/angular impulse on local support, then velocity-project.
   Keep inelastic/projection loss separate from cap surface-energy estimates.
5. Independently advance the supported inviscid fixed-connectivity P1/P0 recipient
   and finite-area material films for0.2ms;41saved states per main comparison.
   The two movies have distinct time origins. Never present them as a continuous
   coupled collision, free-foam motion or a rupture cascade.

## Topology claim boundary

Moving mesh vertices are not mesh-topology updates. The recipient keeps935tets
and1472film triangles; no general merge, pinch-off, T1/T2, conservative remeshing
or subsequent rupture is implemented in this response. One added cap and the
prepared exclusion of the initiating sheet do not satisfy this gate.
Same-mesh time checks are not spatial or topology convergence. Surface motion
alone does not validate geometric forces or material transfer.

## Prior art to reuse or compare

- El Topo (2009): https://www.cs.ubc.ca/labs/imager/tr/2009/eltopo/eltopo.html
- Multimaterial Mesh-Based Surface Tracking / Los Topos (2014):
  https://www.cs.columbia.edu/cg/multitracker/
- Double Bubbles (2015): https://www.cs.columbia.edu/cg/doublebubbles/
- Algebraic Point Set Surfaces (2007): https://cgl.ethz.ch/research/past_projects/apss/
- Partial Optimal Transport (2021preprint/2022journal): https://arxiv.org/abs/2106.03936
- More Power to the Particles (2026preprint,v2): https://arxiv.org/abs/2601.05765v2

The screenshot's proposed smooth point-cloud/free-surface reconstruction is not
implemented by the current tetrahedral recipient. Do not conflate the two.
Geometry reconstruction, topology changes and analytic Laguerre/sphere clipping
have clear prior art. None automatically establishes our force consistency or
the proposed liquid/surfactant/gas/event feedback.

## Acceptance after explicit resume

G0: Independent strong deformation, merging, pinching and triple-junction T1/T2
tests of a selected geometry module; track intersections, element validity and
region labels, and report real before/after connectivity.
G1: Match geometry/energy derivatives and forces; flat-interface spurious flow,
Laplace pressure and capillary-oscillation controls on regular/irregular spatial
refinement. A visually smoothed normal is not a corrected physical force.
G2: Per-event liquid/surfactant, linear/angular momentum and energy transfer
ledgers, with explained physical losses and frozen tolerances. Preserve failures.
G3: Remove the artificial support via a physically justified gas/interface model,
connect drainage and continued rupture, then compare complete trajectories and
propagation/arrest against fair baselines and published experimental evidence.
Retain all original scientific gates. This list is proposed work, not passed tests.

## Rendering/reproduction

`python foam/scripts/render_report_process.py` reads frozen files under
`results/spatial-recipient/`, produces PNG/vectorPDF key states and two MP4s.
`python foam/scripts/build_froth_report.py` builds the standalone November20 report.
Neither command runs a physical simulation. See process-20261120/provenance.json
for source hashes, exact frame times and explicit interpolation/amplification flags.
The report embeds movies and posters; optical interference is not modeled.
The persistent user pause remains in effect for new physics calculations.
