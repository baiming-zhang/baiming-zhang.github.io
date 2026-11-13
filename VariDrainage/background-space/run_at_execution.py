"""Actual native spatial-background study, original metrics and bounded costs."""
from pathlib import Path
import importlib.util,json,subprocess,sys,csv
spec=importlib.util.spec_from_file_location('native_base',Path(__file__).parents[1]/'moving_voronoi/run.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
R=base.R;O=R/'results/background_space_20261113';RAW=R/'results/local-runs/background_space_20261113';LOG=R/'results/logs/background_space_20261113';save,sha=base.save,base.sha
CANONICAL=['states.csv','faces.csv','history.csv','initial_states.csv','initial_faces.csv','preparation.json','mechanics.csv','regional_mechanics.csv','mechanical_nodes.csv']
def case(name,mode='local',ratio=None,dt=.0625,T=200):return dict(name=name,mode=mode,ratio=ratio,dt_us=dt,T_us=T,dx_um=20,radius_um=500,amplitude=0,c_bulk_mM=6,preparation='variational',save_stride=round(10/dt) if T>1 else 1,refine_start_um=200)
def execute(c,timeout):
    dest=RAW/c['name'];assert not dest.exists();dest.mkdir();cmd=[str(base.BIN/'tf_curved_edge_bench.exe'),c['mode']]+[str(c[k]) for k in ['dx_um','dt_us','T_us','radius_um','amplitude']]+[str(dest),str(c['c_bulk_mM']),c['preparation'],str(c['save_stride']),str(c['refine_start_um']),'1']
    if c['ratio'] is not None:cmd.append(str(c['ratio']))
    rec=base.execute(c['name'],cmd,timeout);rec['parameters']=c
    for filename,key in [('summary.json','summary'),('preparation.json','preparation'),('mesh.json','mesh'),('failure.json','failure')]:
        if (dest/filename).exists():rec[key]=json.loads((dest/filename).read_text())
    d=rec.get('summary',{});rec['completed']=rec['returncode']==0 and d.get('steps')==round(c['T_us']/c['dt_us']);rec['native_gate']=False
    if rec['completed']:
        stats=list(csv.DictReader((dest/'mechanics.csv').open()));ratio=max(float(r['g1_N'])/float(r['target_N']) for r in stats);rec['maximum_force_target_ratio']=ratio;rec['native_gate']=len(stats)==d['steps'] and ratio<=1.000001 and max(d[k] for k in ['max_V_res','max_S_res','max_regionV_res','max_regionS_res'])<=1e-10 and d['max_budget_excess']<=1e-7 and rec['preparation']['pressure_residual']<=1e-8
    reference=None
    if c['name']=='control_default':reference=R/'results/local-runs/time_refinement_20261113/dt00625'
    if c['name']=='guard_explicit2':reference=RAW/'guard_omitted'
    if reference and rec['completed']:
        rec['byte_exact_canonical']={name:sha(dest/name)==sha(reference/name) for name in CANONICAL};rec['native_gate'] &= all(rec['byte_exact_canonical'].values())
    rec['files']=[dict(path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in dest.iterdir() if p.is_file()];print(c['name'],rec['returncode'],rec['native_gate'],round(rec['wall_s'],3),rec.get('failure',{}),flush=True);return rec
def main():
    for p in [O,RAW,LOG]:p.mkdir(exist_ok=True)
    base.LOG=LOG
    if '--matrix' not in sys.argv:
        assert not (O/'protocol.json').exists();prior=json.loads((R/'results/prepared_motion_20261113/status.json').read_text());cost=next(r for r in prior['cases'] if r['name']=='u20_dt1');assert cost['completed'] and cost['summary']['c_bulk_mM']==6
        guards=[case('control_default'),case('guard_omitted',T=.0625),case('guard_explicit2',ratio=2,T=.0625),case('guard_bg30',ratio=1.5,T=.0625),case('guard_bg25',ratio=1.25,T=.0625),case('guard_uniform','uniform',T=.0625),case('guard_uniform_fine','uniform',dt=.03125,T=.03125)]
        matrix=[case('u20_dt00625','uniform'),case('u20_dt003125','uniform',dt=.03125),case('bg30_dt00625',ratio=1.5),case('bg25_dt00625',ratio=1.25)]
        sources=[Path(__file__),Path(__file__).with_name('README.md'),R/'thin_film/bench/tf_curved_edge_bench.cpp',R/'thin_film/include/vf/thinfilm/CurvedEdgeMotion.h',R/'thin_film/include/vf/thinfilm/MovingVoronoiFilm.h']
        protocol=dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),source_hashes={p.relative_to(R).as_posix():sha(p) for p in sources},binary_hashes={name:sha(base.BIN/(name+'.exe')) for name in ['tf_curved_edge_bench','test_curved_mechanical_observer','test_edge_preparation','test_curved_edge']},guards=guards,conditional_matrix=matrix,cost=dict(reference='prepared_motion_20261113/u20_dt1',wall_s=cost['wall_s'],dt_us=1,margin=1.3),per_case_timeout_s=2400,whole_compute_budget_s=6600,gates=dict(balance=1e-10,energy=1e-7,force_ratio=1.000001,time=.02,space=.1),observation_times_s=[10e-6,50e-6,200e-6],environment=dict(VF_THREADS='2',inherited_VF='excluded'),scope='Only numerical background spacing changes; fine20um, physical cut240um and all model parameters fixed. Uniform20 comparator time pair required; no claim of full dx convergence or experimental validation.')
        save(O/'protocol.json',protocol);(O/'run_at_execution.py').write_bytes(Path(__file__).read_bytes());state=dict(state='guards_running',tests=[],cases=[]);save(O/'status.json',state)
        for name,count in [('test_curved_mechanical_observer',6),('test_edge_preparation',4),('test_curved_edge',14)]:
            rec=base.execute(name,[str(base.BIN/(name+'.exe'))],120);state['tests'].append(rec);save(O/'status.json',state)
            if rec['returncode']!=0 or f'{count} passed, 0 failed' not in (R/rec['log']).read_text():state.update(state='stopped',reason='unit failure');save(O/'status.json',state);return
        for c in guards:
            rec=execute(c,600);state['cases'].append(rec);save(O/'status.json',state)
            if not rec['native_gate']:break
        state['state']='guards_completed';save(O/'status.json',state);return
    p=json.loads((O/'protocol.json').read_text());state=json.loads((O/'status.json').read_text());assert state['state']=='guards_completed'
    for name,digest in p['source_hashes'].items():assert sha(R/name)==digest
    for name,digest in p['binary_hashes'].items():assert sha(base.BIN/(name+'.exe'))==digest
    good=len(state['cases'])==len(p['guards']) and all(r['native_gate'] for r in state['cases']);used=sum(r['wall_s'] for r in state['tests']+state['cases']);state['skipped']=[];state['state']='matrix_running';save(O/'status.json',state)
    for c in p['conditional_matrix']:
        estimate=p['cost']['wall_s']/c['dt_us']*p['cost']['margin']
        if not good or estimate>2400 or used+estimate>6600:state['skipped'].append(dict(name=c['name'],reason='frozen startup/cost gate',estimated_s=estimate));save(O/'status.json',state);continue
        rec=execute(c,min(2400,6600-used));rec['estimated_s']=estimate;state['cases'].append(rec);used+=rec['wall_s'];save(O/'status.json',state)
    state.update(state='completed',compute_wall_s=used);save(O/'status.json',state)
if __name__=='__main__':main()
