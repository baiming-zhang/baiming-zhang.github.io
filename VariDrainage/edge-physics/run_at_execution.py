"""Frozen material/preparation experiments; no overwrite or tolerance tuning."""
from pathlib import Path
import importlib.util,json,datetime,subprocess,time,numpy as np
spec=importlib.util.spec_from_file_location('prep_native_base',Path(__file__).resolve().parents[1]/'moving_voronoi/run.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
R=base.R;O=R/'results/edge_preparation_20261113';RAW=R/'results/local-runs/edge_preparation_20261113';LOG=R/'results/logs/edge_preparation_20261113';save,sha=base.save,base.sha
def case(name,mode='local',dx=20,dt=.125,T=20,c=6,prep='balanced',amp=0):return dict(name=name,mode=mode,dx_um=dx,dt_us=dt,T_us=T,c_bulk_mM=c,preparation=prep,amplitude=amp,radius_um=500)
def execute(c,raw=RAW,timeout=180):
    dest=raw/c['name'];assert not dest.exists();dest.mkdir(parents=True)
    command=[str(base.BIN/'tf_curved_edge_bench.exe'),c['mode']]+[str(c[k]) for k in ['dx_um','dt_us','T_us','radius_um','amplitude']]+[str(dest),str(c['c_bulk_mM']),c['preparation']]
    r=base.execute(c['name'],command,timeout);r['parameters']=c
    if (dest/'summary.json').exists():r['summary']=json.loads((dest/'summary.json').read_text())
    if (dest/'failure.json').exists():r['failure']=json.loads((dest/'failure.json').read_text())
    d=r.get('summary',{});r['completed']=r['returncode']==0 and d.get('steps')==round(c['T_us']/c['dt_us']);r['native_gate']=r['completed'] and max(d[k] for k in ['max_V_res','max_S_res','max_regionV_res','max_regionS_res'])<=1e-10 and d['max_budget_excess']<=1e-7
    r['files']=[dict(path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in dest.iterdir() if p.is_file()];return r
def main():
    O.mkdir(exist_ok=True);RAW.mkdir(exist_ok=True);LOG.mkdir(exist_ok=True);base.LOG=LOG;assert not (O/'protocol.json').exists();start=time.monotonic()
    cases=[case('default_control',dx=40,dt=.25,c=19.4,prep='projected',amp=.005),case('original_failure_local',c=19.4,prep='projected',amp=.005),case('original_failure_uniform','uniform',c=19.4,prep='projected',amp=.005)]
    for c in [1.94,6,19.4]:
        for prep in ['projected','balanced']:
            for mode in ['local','uniform']:cases.append(case(f'c{c}_{prep}_{mode}',mode,c=c,prep=prep))
    planned=[case('guard10','uniform',10,.03125,.03125)]
    planned += [case(f'{mode}{dx}_T5',mode,dx,.03125,5) for mode,dx in [('local',20),('uniform',20),('local',10),('uniform',10)]]
    planned += [case('uniform20_time','uniform',20,.0625,20)]
    paths=[Path(__file__),Path(__file__).with_name('README.md')]+[R/'thin_film'/p for p in ['bench/tf_curved_edge_bench.cpp','include/vf/thinfilm/CurvedEdge.h','include/vf/thinfilm/CurvedEdgeMotion.h','include/vf/thinfilm/Solution.h','include/vf/thinfilm/Units.h']]
    save(O/'protocol.json',dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),source_hashes={p.relative_to(R).as_posix():sha(p) for p in paths},binary_sha256=sha(base.BIN/'tf_curved_edge_bench.exe'),cases=cases,conditional_refinement=planned,condition='both c6 balanced local20/uniform20 complete and independent gates pass; then guard10 must pass',raw_directory=RAW.relative_to(R).as_posix(),whole_budget_s=3300,per_case_timeout_s=180,environment={'VF_THREADS':'2','caller_VF':'excluded'},started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    state=dict(state='running',tests=[],cases=[]);save(O/'status.json',state)
    r=base.execute('test_curved_edge',[str(base.BIN/'test_curved_edge.exe')],120);state['tests'].append(r);save(O/'status.json',state);assert r['returncode']==0 and '12 passed, 0 failed' in (R/r['log']).read_text()
    previous=R/'results/local-runs/curved_edge_longer_20261113/long_l40_t025'
    for c in cases:
        left=3300-(time.monotonic()-start)
        if left<1:state['skipped']='whole budget';break
        rec=execute(c,timeout=min(left,180))
        if c['name']=='default_control':
            rec['default_byte_exact']={name:sha(previous/name)==sha(RAW/c['name']/name) for name in ['states.csv','faces.csv','history.csv']};rec['native_gate'] &= all(rec['default_byte_exact'].values())
        state['cases'].append(rec);save(O/'status.json',state);print(c['name'],rec['returncode'],rec['native_gate'],f"{rec['wall_s']:.2f}s",rec.get('failure',{}).get('message',''),flush=True)
        if c['name']=='default_control' and not rec['native_gate']:state['skipped']='default control failed';break
    state.update(state='completed',wall_s=time.monotonic()-start);save(O/'status.json',state)
if __name__=='__main__':main()
