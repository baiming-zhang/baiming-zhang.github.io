"""Observer-only force/action and regional-work diagnostic; no solver tuning."""
from pathlib import Path
import importlib.util,json,subprocess,sys,time
spec=importlib.util.spec_from_file_location('native_base',Path(__file__).resolve().parents[1]/'moving_voronoi/run.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
R=base.R;O=R/'results/mechanical_work_20261113';RAW=R/'results/local-runs/mechanical_work_20261113';LOG=R/'results/logs/mechanical_work_20261113';save,sha=base.save,base.sha
def case(name,dt=.5,T=200,observer=True):return dict(name=name,mode='local',dx_um=20,dt_us=dt,T_us=T,radius_um=500,amplitude=0,c_bulk_mM=6,preparation='variational',save_stride=round(10/dt) if T>1 else 1,refine_start_um=200,observer=observer)
def execute(c,timeout):
    dest=RAW/c['name'];assert not dest.exists();dest.mkdir();cmd=[str(base.BIN/'tf_curved_edge_bench.exe'),c['mode']]+[str(c[k]) for k in ['dx_um','dt_us','T_us','radius_um','amplitude']]+[str(dest),str(c['c_bulk_mM']),c['preparation'],str(c['save_stride']),str(c['refine_start_um']),str(int(c['observer']))]
    rec=base.execute(c['name'],cmd,timeout);rec['parameters']=c
    for fn,key in [('summary.json','summary'),('preparation.json','preparation'),('failure.json','failure')]:
        if (dest/fn).exists():rec[key]=json.loads((dest/fn).read_text())
    d=rec.get('summary',{});rec['completed']=rec['returncode']==0 and d.get('steps')==round(c['T_us']/c['dt_us']);rec['native_gate']=rec['completed'] and max(d[k] for k in ['max_V_res','max_S_res','max_regionV_res','max_regionS_res'])<=1e-10 and d['max_budget_excess']<=1e-7
    if c['T_us']==200 and rec['completed']:
        tag={1:'f200_dt1',.5:'f200_dt05',.25:'f200_dt025'}[c['dt_us']];prior=R/'results/local-runs/film_side_20261113'/tag;rec['frozen_trajectory_bytes']={n:sha(dest/n)==sha(prior/n) for n in ['states.csv','faces.csv','history.csv','initial_states.csv','initial_faces.csv','preparation.json']};rec['native_gate'] &= all(rec['frozen_trajectory_bytes'].values())
    rec['files']=[dict(path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in dest.iterdir() if p.is_file()];print(c['name'],rec['returncode'],rec['native_gate'],round(rec['wall_s'],3),rec.get('failure',{}),flush=True);return rec
def main():
    for p in [O,RAW,LOG]:p.mkdir(parents=True,exist_ok=True)
    base.LOG=LOG
    if '--matrix' not in sys.argv:
        assert not (O/'protocol.json').exists();guards=[case('observer_off',observer=False),case('observer_guard',dt=1,T=1)];matrix=[case('obs_dt1',1),case('obs_dt05'),case('obs_dt025',.25)]
        paths=[Path(__file__),Path(__file__).with_name('README.md')]+[R/'thin_film'/p for p in ['bench/tf_curved_edge_bench.cpp','include/vf/thinfilm/CurvedEdgeMechanicalObserver.h','include/vf/thinfilm/CurvedEdgeRegional.h','include/vf/thinfilm/CurvedEdgeMotion.h','include/vf/thinfilm/MovingVoronoiFilm.h','tests/test_curved_mechanical_observer.cpp']]
        save(O/'protocol.json',dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),source_hashes={p.relative_to(R).as_posix():sha(p) for p in paths},binary_hashes={n:sha(base.BIN/(n+'.exe')) for n in ['tf_curved_edge_bench','test_curved_mechanical_observer','test_edge_preparation','test_curved_edge']},guards=guards,conditional_matrix=matrix,raw_directory=RAW.relative_to(R).as_posix(),scope='Observer only, actual inertial/viscous/curved-area model unchanged. All accepted Newton statistics; sampled regional energies, forces, positive clipped viscous dissipation and discrete ports; recovered physical-cut viscous traction is separate.',gates={'frozen_trajectory':'byte-exact6files','regional_derivative_unit':1e-7,'independent_work_formulas':1e-7,'physical_cut_geometry':1e-8,'native_balance':1e-10,'native_energy':1e-7},per_case_timeout_s=600,whole_compute_budget_s=3300,condition='all24checks and observer-off plus1us observer guard pass; estimated guard-wall*steps<=600s per case',sampling='All native Newton/kinetic/viscous statistics. Regional work and nodal forces only at first, every10us and final.',environment={'VF_THREADS':'2','caller_VF':'excluded'},limitations='Endpoint capillary remainder is not assumed positive physical dissipation. Discrete port identities do not independently validate continuum traction or normal flux. No new kinetic rate or changed convergence target.'))
        (O/'run_at_execution.py').write_bytes(Path(__file__).read_bytes());state=dict(state='guards_running',tests=[],cases=[]);save(O/'status.json',state)
        for name,count in [('test_curved_mechanical_observer',6),('test_edge_preparation',4),('test_curved_edge',14)]:
            rec=base.execute(name,[str(base.BIN/(name+'.exe'))],120);state['tests'].append(rec);save(O/'status.json',state)
            if rec['returncode']!=0 or f'{count} passed, 0 failed' not in (R/rec['log']).read_text():state.update(state='stopped',reason='native unit failure');save(O/'status.json',state);return
        for c in guards:
            rec=execute(c,600 if c['T_us']==200 else 180);state['cases'].append(rec);save(O/'status.json',state)
            if not rec['native_gate']:break
        state['state']='guards_completed';save(O/'status.json',state);return
    protocol=json.loads((O/'protocol.json').read_text());state=json.loads((O/'status.json').read_text());assert state['state']=='guards_completed'
    for n,h in protocol['source_hashes'].items():assert sha(R/n)==h,n
    for n,h in protocol['binary_hashes'].items():assert sha(base.BIN/(n+'.exe'))==h,n
    good=len(state['cases'])==2 and all(c['native_gate'] for c in state['cases']);used=sum(r['wall_s'] for r in state['tests']+state['cases']);state['skipped']=[];guard=next((r for r in state['cases'] if r['name']=='observer_guard'),None)
    for c in protocol['conditional_matrix']:
        estimate=(guard or {}).get('wall_s',1e9)*round(c['T_us']/c['dt_us'])
        if not good or estimate>600 or used>=3300:state['skipped'].append(dict(name=c['name'],reason='startup gate failed' if not good else 'frozen cost gate',estimated_s=estimate));continue
        rec=execute(c,min(600,3300-used));state['cases'].append(rec);used+=rec['wall_s'];save(O/'status.json',state)
    state.update(state='completed',compute_wall_s=used);save(O/'status.json',state)
if __name__=='__main__':main()
