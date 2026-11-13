# Geometry/transport source and raw-data mapping

Completed source/results: https://github.com/wang-mengdi/variational-films/tree/6ffb3f5/thin_film/scripts/kinematic_flux

The primary diagnostic passes8motion guards and60actual saved-state checks.
The manufactured transport comparison is NOT fully qualified:4regular area-FD
checks fail. All four direct flux/layout convergence conditions pass, and all
four jitter geometry guards pass. Original criteria were not relaxed.

restore_raw.py restores172CSVfiles from native-observations.zip, preserving
SHA256 and rejecting overwrite. Namespaces map to source campaign native/:
primary -> results/kinematic_flux_adaptive_20261113/native/
transport -> results/kinematic_transport_20261113/native/
failed-direct -> results/kinematic_flux_20261113/native/
failed-fd -> results/kinematic_flux_convex_20261113/native/

Original frozen native film states are referenced by inputs.json. Download them
from https://baiming-zhang.github.io/VariDrainage/time-refinement/ and
https://baiming-zhang.github.io/VariDrainage/mechanical-work/ . No original
film trajectory was recomputed or overwritten here. The geometry CLI observes
midpoint positions and displacement-derived site velocities only.

Use a NEW campaign name/isolated checkout for reproducibility; retain every
frozen output. analyze.py and transport.py refuse existing protocol output;
report.py packages once, with --figures-only for later visual QA. Exact source
and binary SHA256 at execution are recorded by each protocol. The published
ZIP is the canonical original native byte archive, including earlier failures.

Primary reference: Springel2010AREPO, https://arxiv.org/abs/0901.4107 ,
section3.3equations32-34 (author PDF locally stored, not rehosted).
