# Research decision: stateful edge exchange

Date:2026-11-13. User steering: avoid unbounded small solver fixes; consult
existing work broadly and make the innovation testable. This assay follows C1's
original requirement to retain an evolving local region when scale separation
fails. It does not replace the moving-film failures with easier success claims.

## Closest work and consequences

| Primary source | Existing contribution | Consequence for this project |
|---|---|---|
| [Lalli & Giusti2024, JFM](https://doi.org/10.1017/jfm.2024.335) | Thin-film equations with interfacial mobility/surfactant and qualitative marginal-regeneration simulations | Simulating regeneration or adding surfactant alone is insufficient novelty. Need quantitative output and fair same-model references. |
| [Vigna-Brummer et al.,2025 preprint](https://arxiv.org/abs/2503.00630) | Measured meniscus shapes and coupled film supply; flowing regime and a gravito-exchange length | A changing meniscus radius or finite liquid supply is already motivated and studied. Distinguish prediction from imposed empirical flux. |
| [Vigna-Brummer et al.,2026 preprint](https://arxiv.org/abs/2605.27108) | Vertical/inclined plate experiments, steady/transient growth and independent estimates of the exchange coefficient | An external flux/meniscus benchmark is available. Its exchange law must not be inserted and then reported as predicted. |
| [Grassia2022](https://doi.org/10.1098/rspa.2022.0133) | Reduced film-to-film surfactant transport, driving and relaxation, explicit strain/internal states | Stateful meniscus-related transport itself has precedent. Need a different resolved-to-global coupling and validation. |
| [Cantat et al.,2026 review](https://doi.org/10.1146/annurev-fluid-100224-111138) | Capillary drainage, interfacial constraints, elasticity and symmetry breaking | Use the immobile-surface capillary limit only as a diagnostic. It cannot establish mobile-film or adsorption physics. |
| [Monier et al.,2024](https://arxiv.org/abs/2401.03931) | Experimental fringe trajectories and collapse/self-similarity | Keep the full trajectories and material/device provenance as final targets. No fitted fringe exponent in a closure. |
| [Gugercin et al.,2012](https://arxiv.org/abs/1101.3485) | Passive structure-preserving rational interpolation/model reduction | Projection, passivity and rational response bases are established methods to reuse, not new mathematical contributions. |
| [Steinberg et al.,2016](https://arxiv.org/abs/1603.06753) | Moving Voronoi grid noise and geometric-volume inconsistency | Keep the earlier motion-consistency failures visible; conservation does not establish accurate dynamics. |

SciSpace was used for broad discovery; technical interpretations above rely on
the primary papers/publisher or author records. The related Lian et al.2024 paper
(DOI10.1063/5.0235150) was discovered but full text was inaccessible in this pass;
its detailed method is not used as verified evidence. Review PDF download from
the author host timed out; publisher landing metadata is separate from a PDF.
Downloaded public PDFs and exact hashes are in the sibling local folder
`03_DRIFT/literature/edge_exchange_20261113/`; do not rehost the literature PDFs
on the public project site. Manifest records failed downloads explicitly.

## Candidate contribution and kill criteria

Working name: **Stateful Edge Exchange**. A prospective C1 method would attach a
small pressure/flow and surface-chemical state to the moving local edge, preserve
liquid/surfactant and interface power, and select a reduced or resolved edge from
measured port error rather than thickness error alone. The potentially distinctive
result is a validated resolved-to-global exchange operator with a reliable
validity domain; the components above are prior art. No priority claim follows
from this literature pass. More targeted graphics/ALE/reduced-boundary searches
and nonlinear physical validation remain necessary.

First falsifiable question: for an internally resolved thin strip, can identical
current boundary data and stored volume hide different outgoing flux? Then ask
whether a priori port-response modes preserve that flux in untrained pulse/chirp
inputs and several geometries. The controlled linear pilot deliberately omits
meniscus curvature, surfactant and moving-film coupling. It tests a numerical
architecture, not the mechanism of experimental marginal regeneration.

Execute one frozen matrix:18 geometry/input conditions,5/9/17-state port and
slow-mode baselines, instantaneous baseline,160/320-cell full reference, and
80-cell/time-halved checks on one geometry. Do not optimize the basis against
the validation results. If the static closure is adequate, drop the memory claim;
if17 states are inadequate, retain the full local region instead of fitting an
exchange law. Report offline and online cost separately.

Only a successful nonlinear local film-meniscus reference, coupled moving-film
ledger, same-accuracy cost comparison and independent physical holdout could
support the intended C1 contribution. The old17-41% motion spatial differences
are still unresolved and are not cured by this pilot.

## Frozen fair-baseline extension

The18-case run finishes. Add DC-response-corrected slow modes at the same
state caps before making an efficiency/accuracy comparison. Plain low modes
can otherwise have a known port-feedthrough cancellation error. This extension
does not change any original result, signal, physical parameter or gate.
