"""Accept only the documented 10-to-30 utilization-only solver revision."""
from pathlib import Path
import hashlib,json

B = Path(__file__).resolve().parents[1] / 'data/copv/c04_repair_v3'

def assert_same_numerical_solver(runs, summaries):
    proof = json.loads((B / 'gpu_utilization_30_amendment_v1.json').read_text())
    canonical = (B / 'gpu_solver_sources_v3/solver.py').read_bytes()
    old = b'or util>10:raise RuntimeError'
    new = b'or util>30:raise RuntimeError'
    assert canonical.count(old) == 1
    for run, summary in zip(runs, summaries):
        source = (run / 'solver_executed_snapshot.py').read_bytes()
        digest = hashlib.sha256(source).hexdigest()
        assert digest == summary['solver_sha256']
        assert digest in [proof['old_solver_sha256'], proof['new_solver_sha256']]
        assert source == canonical or source == canonical.replace(old, new)
    return dict(solver_sha256=[s['solver_sha256'] for s in summaries],
        resource_amendment='gpu_utilization_30_amendment_v1.json',
        scope='Full executed sources match except the authorized GPU utilization guard 10% to 30%; numerical code is byte-identical.')
