"""Frozen native prepared-edge guards then conditional positive-time matrix."""
from pathlib import Path
import importlib.util,json,subprocess,time,sys
spec=importlib.util.spec_from_file_location('native_base',Path(__file__).resolve().parents[1]/'moving_voronoi/run.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
R=base.R;O=R/'results/prepared_motion_20261113';RAW=R/'results/local-runs/prepared_motion_20261113';LOG=R/'results/logs/prepared_motion_20261113'
sha,save=base.sha,base.save
def case(name,mode='local',dx=20,dt=1,T=200,prep='variational',stride=10,c=6,amp=0):return dict(name=name,mode=mode,dx_um=dx,dt_us=dt,T_us=T,preparation=prep,save_stride=stride,c_bulk_mM=c,amplitude=amp,radius_um=500)
def execute(c,timeout):
    dest=RAW/c['name'];assert not dest.exists();dest.mkdir()
    command=[str(base.BIN/'tf_curved_edge_bench.exe'),c['mode']]+[str(c[k]) for k in ['dx_um','dt_us','T_us','radius_um','amplitude']]+[str(dest),str(c['c_bulk_mM']),c['preparation'],str(c['save_stride'])]
    rec=base.execute(c['name'],command,timeout);rec['parameters']=c
    for f,key in [('summary.json','summary'),('preparation.json','preparation'),('failure.json','failure')]:
        if (dest/f).exists():rec[key]=json.loads((dest/f).read_text())
    s=rec.get('summary',{});rec['completed']=rec['returncode']==0 and s.get('steps')==round(c['T_us']/c['dt_us']);rec['native_gate']=rec['completed'] and max(s[k] for k in ['max_V_res','max_S_res','max_regionV_res','max_regionS_res'])<=1e-10 and s['max_budget_excess']<=1e-7
    if c['preparation']=='variational':rec['native_gate'] &= rec.get('preparation',{}).get('pressure_residual',1)>=0 and rec.get('preparation',{}).get('pressure_residual',1)<=1e-8
    if c['name']=='legacy_control':
        original=R/'results/local-runs/edge_counter_20261113/default_control';rec['legacy_bytes']={n:sha(dest/n)==sha(original/n) for n in ['states.csv','faces.csv','history.csv']};rec['native_gate'] &= all(rec['legacy_bytes'].values())
    rec['files']=[dict(path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in dest.iterdir() if p.is_file()]
    print(c['name'],rec['returncode'],'gate',rec['native_gate'],'wall',round(rec['wall_s'],3),rec.get('failure',{}),flush=True);return rec
def main():
    for p in [O,RAW,LOG]:p.mkdir(parents=True,exist_ok=True)
    base.LOG=LOG
    if '--matrix' not in sys.argv:
        assert not (O/'protocol.json').exists()
        guards=[case('legacy_control',dx=40,dt=.25,T=20,prep='projected',stride=1,c=19.4,amp=.005)]+[case(f'guard_{mode}{dx}',mode,dx,T=1,stride=1) for mode,dx in [('local',20),('uniform',20),('local',10),('uniform',10)]]
        matrix=[case('l20_dt1'),case('u20_dt1','uniform'),case('l20_dt05',dt=.5,stride=20),case('l10_dt1',dx=10),case('u10_dt1','uniform',10),case('u10_dt05','uniform',10,.5,stride=20)]
        paths=[Path(__file__),Path(__file__).with_name('README.md')]+[R/'thin_film'/p for p in ['include/vf/thinfilm/CurvedEdgePreparation.h','include/vf/thinfilm/CurvedEdgeMotion.h','include/vf/thinfilm/CurvedEdge.h','bench/tf_curved_edge_bench.cpp','tests/test_edge_preparation.cpp']]
        save(O/'protocol.json',dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),source_hashes={p.relative_to(R).as_posix():sha(p) for p in paths},binaries={n:sha(base.BIN/(n+'.exe')) for n in ['tf_curved_edge_bench','test_edge_preparation','test_curved_edge']},guards=guards,conditional_matrix=matrix,raw_directory=RAW.relative_to(R).as_posix(),per_guard_timeout_s=180,per_matrix_timeout_s=600,whole_compute_budget_s=3300,positive_observation_times_s=[10e-6,50e-6,200e-6],scope='Actual moving variable-S curved-area model with variational initial edge; trace heights constrained ONLY in initialization. No kinetic closure or empirical drainage law.',condition='all native tests and guards pass; guard wall * requested steps <=600s for each matrix case; no parameter/tolerance retry',gates={'V_S_region':1e-10,'energy':1e-7,'independent_formula':1e-7,'time_dynamic':.02,'space_dynamic':.1},sampling='All accepted history rows; three full states/faces for step1, every10us and final. Independent formulas audited at saved steps only.',environment={'VF_THREADS':'2','caller_VF':'excluded'}))
        (O/'run_at_execution.py').write_bytes(Path(__file__).read_bytes())
        state=dict(state='guards_running',tests=[],cases=[]);save(O/'status.json',state)
        for name,expected in [('test_edge_preparation','4 passed, 0 failed'),('test_curved_edge','14 passed, 0 failed')]:
            rec=base.execute(name,[str(base.BIN/(name+'.exe'))],120);state['tests'].append(rec);save(O/'status.json',state)
            if rec['returncode']!=0 or expected not in (R/rec['log']).read_text():state.update(state='stopped',reason='native unit check failed');save(O/'status.json',state);return
        for c in guards:
            state['cases'].append(execute(c,180));save(O/'status.json',state)
            if c['name']=='legacy_control' and not state['cases'][-1]['native_gate']:break
        state['state']='guards_completed';save(O/'status.json',state);return
    protocol=json.loads((O/'protocol.json').read_text());state=json.loads((O/'status.json').read_text());assert state['state']=='guards_completed'
    for name,h in protocol['source_hashes'].items():assert sha(R/name)==h,name
    for name,h in protocol['binaries'].items():assert sha(base.BIN/(name+'.exe'))==h,name
    guards={c['name']:c for c in state['cases']};used=sum(c['wall_s'] for c in state['cases'])+sum(c['wall_s'] for c in state['tests']);state['skipped']=[]
    allguards=len(guards)==5 and all(c['native_gate'] for c in guards.values())
    for c in protocol['conditional_matrix']:
        guard=guards.get(f"guard_{c['mode']}{c['dx_um']}");estimate=(guard or {}).get('wall_s',1e9)*round(c['T_us']/c['dt_us'])
        if not allguards or estimate>600 or used>=protocol['whole_compute_budget_s']:
            state['skipped'].append(dict(name=c['name'],reason='native guard failed' if not allguards else 'predeclared cost gate',estimated_s=estimate));continue
        rec=execute(c,min(600,protocol['whole_compute_budget_s']-used));state['cases'].append(rec);used+=rec['wall_s'];save(O/'status.json',state)
    state.update(state='completed',compute_wall_s=used);save(O/'status.json',state)
if __name__=='__main__':main()
