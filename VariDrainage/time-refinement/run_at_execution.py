"""Bounded native dt refinement, unchanged model, binary and frozen criteria."""
from pathlib import Path
import importlib.util,json,subprocess,sys,csv,time
spec=importlib.util.spec_from_file_location('native_base',Path(__file__).parents[1]/'moving_voronoi/run.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
R=base.R;O=R/'results/time_refinement_20261113';RAW=R/'results/local-runs/time_refinement_20261113';LOG=R/'results/logs/time_refinement_20261113';save,sha=base.save,base.sha
def case(name,dt,T=200):return dict(name=name,mode='local',dx_um=20,dt_us=dt,T_us=T,radius_um=500,amplitude=0,c_bulk_mM=6,preparation='variational',save_stride=round(10/dt) if T>1 else 1,refine_start_um=200,observer=True)
def execute(c,timeout):
    dest=RAW/c['name'];assert not dest.exists();dest.mkdir();cmd=[str(base.BIN/'tf_curved_edge_bench.exe'),c['mode']]+[str(c[k]) for k in ['dx_um','dt_us','T_us','radius_um','amplitude']]+[str(dest),str(c['c_bulk_mM']),c['preparation'],str(c['save_stride']),str(c['refine_start_um']),'1']
    rec=base.execute(c['name'],cmd,timeout);rec['parameters']=c
    for fn,key in [('summary.json','summary'),('preparation.json','preparation'),('failure.json','failure')]:
        if (dest/fn).exists():rec[key]=json.loads((dest/fn).read_text())
    d=rec.get('summary',{});rec['completed']=rec['returncode']==0 and d.get('steps')==round(c['T_us']/c['dt_us']);rec['native_gate']=False
    if rec['completed']:
        stats=list(csv.DictReader((dest/'mechanics.csv').open()));rec['maximum_force_target_ratio']=max(float(s['g1_N'])/float(s['target_N']) for s in stats);rec['native_gate']=len(stats)==d['steps'] and rec['maximum_force_target_ratio']<=1.000001 and max(d[k] for k in ['max_V_res','max_S_res','max_regionV_res','max_regionS_res'])<=1e-10 and d['max_budget_excess']<=1e-7 and rec['preparation']['pressure_residual']<=1e-8
    if c['name']=='control_dt025' and rec['completed']:
        prior=R/'results/local-runs/mechanical_work_fixed_20261113/obs_dt025';rec['byte_exact_to_frozen']={n:sha(dest/n)==sha(prior/n) for n in ['states.csv','faces.csv','history.csv','initial_states.csv','initial_faces.csv','preparation.json','mechanics.csv','regional_mechanics.csv','mechanical_nodes.csv']};rec['native_gate'] &= all(rec['byte_exact_to_frozen'].values())
    rec['files']=[dict(path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in dest.iterdir() if p.is_file()];print(c['name'],rec['returncode'],rec['native_gate'],round(rec['wall_s'],3),rec.get('failure',{}),flush=True);return rec
def main():
    for p in [O,RAW,LOG]:p.mkdir(exist_ok=True)
    base.LOG=LOG
    if '--matrix' not in sys.argv:
        assert not (O/'protocol.json').exists();prior=R/'results/mechanical_work_fixed_20261113';p=json.loads((prior/'protocol.json').read_text());s=json.loads((prior/'status.json').read_text());ref=next(r for r in s['cases'] if r['name']=='obs_dt025');assert sha(base.BIN/'tf_curved_edge_bench.exe')==p['binary_hashes']['tf_curved_edge_bench']
        for f in ref['files']:assert sha(R/f['path'])==f['sha256']
        guards=[case('control_dt025',.25),case('guard_dt0125',.125,.125)];matrix=[case('dt0125',.125),case('dt00625',.0625),case('dt003125',.03125)];sources=[Path(__file__),Path(__file__).with_name('README.md')]+[R/x for x in ['thin_film/bench/tf_curved_edge_bench.cpp','thin_film/include/vf/thinfilm/CurvedEdgeMotion.h','thin_film/include/vf/thinfilm/MovingVoronoiFilm.h','thin_film/include/vf/thinfilm/CurvedEdgeMechanicalObserver.h']]
        save(O/'protocol.json',dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),source_hashes={p.relative_to(R).as_posix():sha(p) for p in sources},binary_hashes={n:sha(base.BIN/(n+'.exe')) for n in ['tf_curved_edge_bench','test_curved_mechanical_observer','test_edge_preparation','test_curved_edge']},guards=guards,conditional_matrix=matrix,physical_cut_m=240e-6,fine_start_m=200e-6,per_case_timeout_s=600,whole_compute_budget_s=3300,cost_estimate=dict(source='mechanical_work_fixed_20261113/obs_dt025',wall_s=ref['wall_s'],dt_us=.25,margin=1.15),gates=dict(balance=1e-10,energy=1e-7,force_ratio=1.000001,time_dh_and_displacement=.02),observation_times_s=[10e-6,50e-6,200e-6],environment=dict(VF_THREADS='2',caller_VF='excluded'),scope='Only dt changes; unchanged model, physical cut, materials, preparation and solver. All native scalar/force statistics, sampled independent geometry/constitutive/mechanical audits. Time pair pass does not imply spatial or experimental validation.'))
        (O/'run_at_execution.py').write_bytes(Path(__file__).read_bytes());state=dict(state='guards_running',tests=[],cases=[]);save(O/'status.json',state)
        for name,count in [('test_curved_mechanical_observer',6),('test_edge_preparation',4),('test_curved_edge',14)]:
            rec=base.execute(name,[str(base.BIN/(name+'.exe'))],120);state['tests'].append(rec);save(O/'status.json',state)
            if rec['returncode']!=0 or f'{count} passed, 0 failed' not in (R/rec['log']).read_text():state.update(state='stopped',reason='native unit failure');save(O/'status.json',state);return
        for c in guards:
            rec=execute(c,600);state['cases'].append(rec);save(O/'status.json',state)
            if not rec['native_gate']:break
        state['state']='guards_completed';save(O/'status.json',state);return
    p=json.loads((O/'protocol.json').read_text());state=json.loads((O/'status.json').read_text());assert state['state']=='guards_completed'
    for n,h in p['source_hashes'].items():assert sha(R/n)==h,n
    for n,h in p['binary_hashes'].items():assert sha(base.BIN/(n+'.exe'))==h,n
    good=len(state['cases'])==2 and all(c['native_gate'] for c in state['cases']);used=sum(r['wall_s'] for r in state['tests']+state['cases']);state['skipped']=[]
    for c in p['conditional_matrix']:
        estimate=p['cost_estimate']['wall_s']*.25/c['dt_us']*p['cost_estimate']['margin']
        if not good or estimate>600 or used+estimate>3300:state['skipped'].append(dict(name=c['name'],reason='guard or frozen cost gate',estimate_s=estimate));continue
        rec=execute(c,min(600,3300-used));state['cases'].append(rec);used+=rec['wall_s'];good &= rec['native_gate'];save(O/'status.json',state)
    state.update(state='completed',compute_wall_s=used);save(O/'status.json',state)
if __name__=='__main__':main()
