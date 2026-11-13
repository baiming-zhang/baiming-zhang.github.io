# Native temporal refinement source and references

Completed research source/results: https://github.com/wang-mengdi/variational-films/tree/0efdaa6/thin_film/scripts/time_refinement

Protocol records exact native source and executable hashes. The native binary
is unchanged from mechanical_work_fixed_20261113; the new dt.25us control has
nine byte-exact trajectory/observer files. Only dt and output sampling stride
vary across the three new200us runs. No material/rate/tolerance fit.

Older dt1/.5us curves are frozen inputs, NOT recomputed here. Their raw archives:
https://baiming-zhang.github.io/VariDrainage/mechanical-work/
The source analysis checks every original input SHA256 before comparisons.

restore_raw.py creates NEW_OUTPUT_DIRECTORY/native/CASE. The native namespace
corresponds to original results/local-runs/time_refinement_20261113/CASE.
Keep all frozen directories. To reproduce computation, use an isolated checkout
and a NEW campaign name for O, RAW and LOG in run.py; record the resulting new
source hashes. run.py --matrix requires the SAME source/binary hashes as its
startup protocol. analyze.py, package.py and report.py then process that campaign.
No overwrite, reset or silent replay of existing completed native batches.

Native all-step statistics and saved-step independent formula audits are
separate. No full spatial, physical normal-flux or experimental qualification.
