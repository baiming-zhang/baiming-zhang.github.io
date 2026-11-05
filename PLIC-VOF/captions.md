# Supplementary movie captions

All movies are silent, 30 fps, and 10–14 seconds long. Both 3840 × 2160 and
1920 × 1080 H.264 versions are supplied. Neural + AMGCL and Neural + compact
use the same independently rebuilt, frozen original-course model. This model
must not be identified as the unavailable historical paper checkpoint.

**Movie S1 — Rising bubble (14 s).** A 256 × 64 reflected half-domain is rotated
so that the rise direction appears vertical. All three panels use the same
physical camera. The plotted centroid traces retain their actual displacement
differences. The common visualization interval is t = 0–2.9. Aphros is an
independent numerical baseline, not an exact evolving interface.

**Movie S2 — Capillary wave (14 s).** Three 128² calculations are shown over
the complete common saved interval (approximately t = 0–2.24262). The enlarged
interface height is obtained by column integration of cell-average liquid
fraction. The lower plot shows the fundamental positive amplitude with the
documented cell-average sinc correction. The reference provides amplitude
only; no exact phase or evolving theoretical interface is asserted.

**Movie S3 — Falling drop (14 s).** Three 128² trajectories are compared over
t = 0–0.6 using an identical physical view and a common VOF color mapping.
The inset shows the absolute liquid-fraction difference between the AMGCL and
compact runs on a fixed 0–1 scale. This backend difference is not an error
against a known exact trajectory. No backend-specific recentring is applied.

**Movie S4 — Dam break (14 s).** The corrected, independently checked initial
geometry is used for all panels. Three 128² trajectories are compared over
t = 0–0.6, together with the absolute backend field-difference inset. Invalid
earlier all-liquid driver outputs are excluded. The full domain and fixed
color scales are retained throughout.

**Movie S5 — Static drop (10 s).** Three 128² pressure fields are shown over
t = 0–0.5. The mean pressure in nearly pure gas cells is removed independently
to fix the arbitrary gauge. A single 0–4.5 pressure scale is used in all panels
and frames. White contours indicate f = 0.5. The lower plot compares the
pressure jump with the analytic Laplace value, 4.

**Movie S6 — Reversible vortex (14 s).** Neural and Aphros cell-average liquid
fractions are shown over t = 0–8 at 128². Insets display min(f, 1−f) on the
same logarithmic scale, clipped between 1e−9 and 1e−2, with smaller values
white. Mixed-cell counts use 1e−6 < f < 1−1e−6. The exact circle is drawn only
at the initial and return times. Final MAE values use closed-form circle/cell
intersection averages. The initial circle is not an intermediate exact solution.

**Movie S7 — Reversed advection (14 s).** The same definitions and color limits
as Movie S6 are used for t = 0–5 at 128². Tail location, tail volume, mixed-cell
count, and actual neural face requests are distinct quantities. Prescribed-
velocity advection does not call the pressure solver.

## Temporal rendering and interpretation

Movies S1–S5 linearly interpolate the underlying scalar fields only inside
their recorded time intervals and then freshly render each contour/field at
the common physical time. The interpolation is for visualization, not an
additional simulated trajectory or an accuracy measurement. Movies S6–S7
show actual synchronized saved fields held until the next saved sample; their
tail maps and counts are not temporally interpolated. No crossfade, adjacent
RGB-frame blending, motion-trail overlay, or image mixing is used in any movie.
Quantitative report values come from the accepted original measurements.

Aphros retains its native CPU configuration, particle curvature, sharpening,
pressure tolerance, and documented initialization choices. Cross-code timings
are deployment observations rather than equal-algorithm, equal-resource
speedups. These movies do not establish a universal neural accuracy advantage.
