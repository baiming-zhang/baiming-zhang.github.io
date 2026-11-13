"""Bounded actual uniform10 refinement, frozen references and no solver edit."""
from pathlib import Path
import importlib.util,json,subprocess,sys
spec=importlib.util.spec_from_file_location('background_runner',Path(__file__).parents[1]/'background_space/run.py');bg=importlib.util.module_from_spec(spec);spec.loader.exec_module(bg)
R=bg.R;O=R/'results/global_space_20261113';RAW=R/'results/local-runs/global_space_20261113';LOG=R/'results/logs/global_space_20261113';save,sha=bg.save,bg.sha
bg.O=O;bg.RAW=RAW;bg.LOG=LOG;bg.base.LOG=LOG
def make(name,dt,T):
    c=bg.case(name,mode='uniform',dt=dt,T=T);c.update(dx_um=10,observer=True);return c
def validate(p):
    for name,digest in p['native_source_hashes'].items():assert sha(R/name)==digest,name
    for name,digest in p['binary_hashes'].items():assert sha(bg.base.BIN/(name+'.exe'))==digest,name
    for name,digest in p['driver_hashes'].items():assert sha(R/name)==digest,name
def execute(c,timeout):
    rec=bg.execute(c,timeout);rec['native_gate'] &= rec.get('summary',{}).get('N')==544 and rec.get('summary',{}).get('failures')==0;return rec
def main():
    for p in [O,RAW,LOG]:p.mkdir(exist_ok=True)
    guards=[make('guard_dt00625',.0625,.0625),make('guard_dt003125',.03125,.03125)];matrix=[make('u10_dt00625',.0625,10),make('u10_dt003125',.03125,10)]
    if '--matrix' not in sys.argv:
        assert not (O/'protocol.json').exists();previous=json.loads((R/'results/layout_phase_20261113/protocol.json').read_text());checks=json.loads((R/'results/layout_phase_20261113/status.json').read_text())['tests'];assert len(checks)==3
        for rec,count in zip(checks,[6,4,14]):assert rec['returncode']==0 and sha(R/rec['log'])==rec['log_sha256'] and f'{count} passed, 0 failed' in (R/rec['log']).read_text(encoding='utf-8')
        refs=json.loads((R/'results/background_space_20261113/status.json').read_text());refs=[c for c in refs['cases'] if c['name'] in ['u20_dt00625','u20_dt003125']];assert len(refs)==2 and all(c['native_gate'] for c in refs)
        for c in refs:
            for f in c['files']:assert sha(R/f['path'])==f['sha256']
        historic=json.loads((R/'results/prepared_motion_20261113/guards.json').read_text());historic=next(c for c in historic['cases'] if c['name']=='guard_uniform10');assert historic['native_gate'] and historic['summary']['N']==544
        native=[p for p in (R/'thin_film').rglob('*') if p.is_file() and p.suffix in ['.h','.cpp']];drivers=[Path(__file__),Path(__file__).with_name('README.md'),Path(__file__).parents[1]/'background_space/run.py',Path(__file__).parents[1]/'moving_voronoi/run.py']
        p=dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),native_source_hashes={str(p.relative_to(R)):sha(p) for p in native},binary_hashes=previous['binary_hashes'],driver_hashes={str(p.relative_to(R)):sha(p) for p in drivers},reused_verified_checks=checks,reused_check_count=24,reused_checks_scope='Previous execution verified, not rerun in this campaign',reference_cases=refs,historical_cost=dict(source='results/prepared_motion_20261113/guards.json:guard_uniform10',wall_s=historic['wall_s'],scope='one-step engineering predictor, before observer-era actual runtime'),guards=guards,conditional_matrix=matrix,cost_margin=1.3,whole_compute_budget_s=6600,per_case_budget_s=5400,gates=dict(balance=1e-10,energy=1e-7,force_ratio=1.000001,preparation=1e-8,time=.02,space=.1),scope='Actual full-domain dx20->10 pilot at10us only. Original50/200us remain required. Neither full time/space qualification nor physical flux/experiment validation.',environment=dict(VF_THREADS='2',inherited_VF='excluded'))
        validate(p);save(O/'protocol.json',p);(O/'run_at_execution.py').write_bytes(Path(__file__).read_bytes());state=dict(state='guards_running',tests=[],reused_native_checks=24,cases=[]);save(O/'status.json',state)
        for c in guards:
            state['cases'].append(execute(c,120));save(O/'status.json',state)
        state['state']='guards_completed';save(O/'status.json',state);return
    p=json.loads((O/'protocol.json').read_text());state=json.loads((O/'status.json').read_text());assert state['state']=='guards_completed';validate(p);good=len(state['cases'])==2 and all(c['native_gate'] for c in state['cases']);used=sum(c['wall_s'] for c in state['cases']);perstep=max([p['historical_cost']['wall_s']]+[c['wall_s'] for c in state['cases']]);estimates=[perstep*round(c['T_us']/c['dt_us'])*p['cost_margin'] for c in matrix];cost_ok=used+sum(estimates)<=6600 and max(estimates)<=5400;state.update(state='matrix_running',skipped=[],cost_prediction=dict(per_step_s=perstep,predictions_s=estimates,pair_fits_budget=cost_ok));save(O/'status.json',state)
    for c,estimate in zip(matrix,estimates):
        validate(p)
        if not good or not cost_ok or used+estimate>6600:state['skipped'].append(dict(name=c['name'],reason='frozen guard or joint/remaining cost gate',estimated_s=estimate));save(O/'status.json',state);continue
        rec=execute(c,min(5400,6600-used));rec['estimated_s']=estimate;state['cases'].append(rec);used+=rec['wall_s'];save(O/'status.json',state)
    state.update(state='completed',compute_wall_s=used);save(O/'status.json',state)
if __name__=='__main__':main()
