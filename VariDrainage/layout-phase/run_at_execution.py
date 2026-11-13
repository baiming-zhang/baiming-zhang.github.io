"""Bounded same-resolution actual film layout sensitivity experiment."""
from pathlib import Path
import importlib.util,json,subprocess,sys,csv,time
spec=importlib.util.spec_from_file_location('background_runner',Path(__file__).parents[1]/'background_space/run.py');bg=importlib.util.module_from_spec(spec);spec.loader.exec_module(bg)
R=bg.R;O=R/'results/layout_phase_20261113';RAW=R/'results/local-runs/layout_phase_20261113';LOG=R/'results/logs/layout_phase_20261113';save,sha=bg.save,bg.sha
bg.O=O;bg.RAW=RAW;bg.LOG=LOG;bg.base.LOG=LOG

def make(name,dt,T):
    c=bg.case(name,ratio=1.,dt=dt,T=T);c['observer']=True;return c
def frozen():
    old=json.loads((R/'results/background_space_20261113/protocol.json').read_text())
    for name,digest in old['source_hashes'].items():
        if name.endswith(('.cpp','.h')):assert sha(R/name)==digest,name
    for name,digest in old['binary_hashes'].items():assert sha(bg.base.BIN/(name+'.exe'))==digest,name
    return old
def main():
    for p in [O,RAW,LOG]:p.mkdir(exist_ok=True)
    prior=frozen();cases=[make('guard_dt00625',.0625,.0625),make('guard_dt003125',.03125,.03125),make('split20_dt00625',.0625,200),make('split20_dt003125',.03125,200)]
    if '--matrix' not in sys.argv:
        assert not (O/'protocol.json').exists()
        refs=json.loads((R/'results/background_space_20261113/status.json').read_text());refs=[c for c in refs['cases'] if c['name'] in ['u20_dt00625','u20_dt003125']];assert len(refs)==2 and all(c['native_gate'] for c in refs)
        for c in refs:
            for f in c['files']:assert sha(R/f['path'])==f['sha256']
        protocol=dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),native_reference_protocol='results/background_space_20261113/protocol.json',native_reference_sha256=sha(R/'results/background_space_20261113/protocol.json'),native_source_hashes={k:v for k,v in prior['source_hashes'].items() if k.endswith(('.cpp','.h'))},binary_hashes=prior['binary_hashes'],runner_sha256=sha(Path(__file__)),plan_sha256=sha(Path(__file__).with_name('README.md')),cases=cases,reference_cases=refs,budget_s=4000,per_case_timeout_s=2400,gates=dict(balance=1e-10,energy=1e-7,force_ratio=1.000001,preparation=1e-8,time=.02,space=.1),scope='Same nominal dx20, altered fine-site deterministic jitter phase. Combined geometry/preparation sensitivity, not isolated physical causation or global spatial convergence.')
        save(O/'protocol.json',protocol);(O/'run_at_execution.py').write_bytes(Path(__file__).read_bytes());state=dict(state='guards_running',tests=[],cases=[]);save(O/'status.json',state)
        for name,count in [('test_curved_mechanical_observer',6),('test_edge_preparation',4),('test_curved_edge',14)]:
            rec=bg.base.execute(name,[str(bg.base.BIN/(name+'.exe'))],120);state['tests'].append(rec);save(O/'status.json',state)
            assert rec['returncode']==0 and f'{count} passed, 0 failed' in (R/rec['log']).read_text()
        for c in cases[:2]:
            rec=bg.execute(c,120);rec['native_gate'] &= rec['summary']['N']==144;state['cases'].append(rec);save(O/'status.json',state)
            if not rec['native_gate']:break
        state['state']='guards_completed';save(O/'status.json',state);return
    p=json.loads((O/'protocol.json').read_text());state=json.loads((O/'status.json').read_text());assert state['state']=='guards_completed';assert sha(Path(__file__))==p['runner_sha256'];assert sha(Path(__file__).with_name('README.md'))==p['plan_sha256']
    good=len(state['cases'])==2 and all(c['native_gate'] for c in state['cases']);used=sum(c['wall_s'] for c in state['cases']+state['tests']);state.update(state='matrix_running',skipped=[]);save(O/'status.json',state)
    for c,estimate in zip(cases[2:],[747.4,1487.0]):
        frozen()
        if not good or used+estimate>4000:state['skipped'].append(dict(name=c['name'],reason='guard/cost gate',estimate_s=estimate));save(O/'status.json',state);continue
        rec=bg.execute(c,min(2400,4000-used));state['cases'].append(rec);used+=rec['wall_s'];save(O/'status.json',state)
    state.update(state='completed',compute_wall_s=used);save(O/'status.json',state)
if __name__=='__main__':main()
