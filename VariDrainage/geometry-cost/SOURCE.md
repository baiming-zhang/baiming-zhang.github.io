Results/source commit: https://github.com/wang-mengdi/variational-films/tree/ae50c2708be992d2b1ee14e34efb393a6b135913

# Geometry candidate provenance

Native geometry probe/driver:a76282d/a24f640; candidate controls:fd50ebc/ff87297.
Corrected guard-count continuation:c1be16d. Protocols include all source/binary
hashes. Original stop remains intact; no solver or default header changed.
Sources: thin_film/scripts/geometry_cost/, VoronoiClipProbe.h and named wrappers.

Observed tool build output (summary, not a full captured raw build log):
tf_geometry_cost_probe build succeeds7.531s; candidate profile/bench and four
test targets build succeeds30.375s, Windows/MSVC release+xmake.
Commands: xmake build tf_geometry_cost_probe; xmake build
tf_geometry_candidate_profile tf_geometry_candidate_bench
test_geometry_candidate_curved test_geometry_candidate_preparation
test_geometry_candidate_observer test_geometry_candidate_transfer.
No failed native compilation in this campaign. Native execution stdout logs
are preserved in logs/, separately from these tool-observed build summaries.

Original production tf_curved_edge_bench remains SHA256
76e49ed253e2341df37f93d02803272d020b8500ed9265da6f7d0ec124b6c5c3.
Baseline actual raw reference: ../global-space/, u10_dt00625; its native bytes
verified before comparison and not rerun. Native namespace probe/ restores to
geometry_cost_20261113/; candidate/ to geometry_candidate_20261113/ under
results/local-runs/. All snapshots and partial/failed attempts retained.
Frozen native campaigns terminal before report packaging. No new physical
convergence claimed. Full C1 remains open.
All11bytes gate failed because native summary includes wall_s. Ten other raw
files are byte-exact, and summary differs only in wall_s. The original failure
is retained; posthoc exclusion of this timing field is separately reported.
