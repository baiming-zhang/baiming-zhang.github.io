## Film-side coverage completed (2026-11-13)

[White report, actual layouts, PNG/PDF and raw data](index.html).
8/8native runs complete:2byte-exact200us controls,2one-step guards,4new200us
matrix cases;18native checks pass. All128SAVED-step independent formula audits
pass; native scalar history covers every accepted step. Numerical fine-patch
start changes240->200/160um; physical cut/trace remains240um, far bath320um.
Materials, forces, fluxes, initialization law and tolerances are unchanged.

Coverage now includes the thin points.200us thickness-change errors vs uniform20
fall31.23->19.93->9.59%; displacement errors144.08->87.34->91.21% still fail10%.
Start200 time pairs1/.5 and.5/.25us have displacement errors10.74/7.58%, fail2%.
Matched-ID supplementary time diagnostic also fails (~7.41% last pair); frozen
moving-cell-overlap criteria are not replaced. Uniform20is a numerical comparator,
not converged truth. Extent improvement does not validate full dynamics.

Recovered relative cut-flux discrepancy improves8.46->2.47%, but material flux
recovery still differs41.43% at start160. Nonconforming recovery is not exact
normal-flux validation. Far-bath cumulative transfer remains negative at200us.
Regional mechanical/kinetic power and full C1/experiment requirements stay open.
Next: observe mechanical force/action and transport consistency before expansion.
56raw files15.67MB restore from6bounded ZIP parts;4PNG/vector-PDF pairs. No >2h job.


Frozen source/binary hashes and commands: protocol.json. Source implementation ff6bde0 with byte/EOL verification in native_sources.json. Initial discrete inventories are reported in summary.json coverage.

Build tf_curved_edge_bench; run thin_film/scripts/film_side/run.py and --matrix in a fresh output namespace. Then analyze.py, package.py, report.py. Protocol and archive builder refuse overwrite. Raw files stay in results/local-runs/film_side_20261113. METHOD preserves the execution-time plan. Matched-ID supplementary diagnostics reuse saved states and leave acceptance unchanged.

LITERATURE.md records primary-source scope and the next bounded mechanical/transport diagnostic. No new C2 rate is introduced.
