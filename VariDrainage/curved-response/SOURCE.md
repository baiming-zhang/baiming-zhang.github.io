Results/source: https://github.com/wang-mengdi/variational-films/tree/b2710e0a5860c627b003840a5c670f0091aeb71c

# Original full curved-energy response provenance

Preregistered plan:f11ff30. Native probe/driver:eca745a.
Source: thin_film/scripts/curved_response/, bench/tf_curved_response_probe.cpp.
All ORIGINAL physical/native header hashes and binary hashes in protocol.json.
No experimental geometry dispatch; default solver source/binary unchanged.
Production SHA25676e49ed253e2341df37f93d02803272d020b8500ed9265da6f7d0ec124b6c5c3.

Observed tool build output summary (not a full captured raw build log):
xmake build tf_curved_response_probe, Windows/MSVC release, succeeds12.703s.
Execution: python thin_film/scripts/curved_response/run.py, one-shot/refuses
existing output. Native stdout logs in logs/ retained verbatim. No native
compilation/runtime failure in this batch. No prior successful trajectory rerun.

Original4actual inputs from cost_profile_20261113; preserved source SHA256 and
exact copied input.csv in raw archive. They came from background_space/global_space
uniform20/10 atdt.03125us, initial and10us stage2. Flat fixtures generated only
inside the read-only probe; they are analytic controls, not a prescribed dynamic
initial profile. input.json records descriptor/dx/source. response.nodes.csv
retains all native positions, inventories, surface areas, velocities and gradients.

Archive native/ maps to results/local-runs/curved_response_20261113/.
All3amplitudes and90rows retained. Exact raw byte hashes in download_manifest.
Scientific meaning: directional free-energy response+separate frozen physical
viscosity projections. Not a full coupled eigenproblem, actual space convergence,
physical normal traction/flux or experiment validation. Original full C1 remains.
