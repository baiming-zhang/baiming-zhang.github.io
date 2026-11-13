"""Film-side refinement extent test. Physical interface remains x240um."""
from pathlib import Path
import importlib.util,json,subprocess,time,sys
spec=importlib.util.spec_from_file_location('native_base',Path(__file__).resolve().parents[1]/'moving_voronoi/run.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
R=base.R;O=R/'results/film_side_20261113';RAW=R/'results/local-runs/film_side_20261113';LOG=R/'results/logs/film_side_20261113';save,sha=base.save,base.sha
def case(name,start=200,mode='local',dt=1,T=200):return dict(name=name,mode=mode,dx_um=20,dt_us=dt,T_us=T,radius_um=500,amplitude=0,c_bulk_mM=6,preparation='variational',save_stride=round(10/dt) if T>1 else 1,refine_start_um=start)
def execute(c,timeout):
    dest=RAW/c['name'];assert not dest.exists();dest.mkdir()
    cmd=[str(base.BIN/'tf_curved_edge_bench.exe'),c['mode']]+[str(c[k]) for k in ['dx_um','dt_us','T_us','radius_um','amplitude']]+[str(dest),str(c['c_bulk_mM']),c['preparation'],str(c['save_stride'])]
    if c['name'] not in ['control_l240','control_uniform']:cmd.append(str(c['refine_start_um']))
    rec=base.execute(c['name'],cmd,timeout);rec['parameters']=c
    for f,key in [('summary.json','summary'),('preparation.json','preparation'),('failure.json','failure')]:
        if (dest/f).exists():rec[key]=json.loads((dest/f).read_text())
    m=rec.get('summary',{});rec['completed']=rec['returncode']==0 and m.get('steps')==round(c['T_us']/c['dt_us']);rec['native_gate']=rec['completed'] and max(m[k] for k in ['max_V_res','max_S_res','max_regionV_res','max_regionS_res'])<=1e-10 and m['max_budget_excess']<=1e-7 and rec.get('preparation',{}).get('pressure_residual',1)<=1e-8
    if c['name'] in ['control_l240','control_uniform']:
        prior='l20_dt1' if c['name']=='control_l240' else 'u20_dt1';source=R/'results/local-runs/prepared_motion_20261113'/prior
        rec['byte_exact_to_frozen_reference']={n:sha(dest/n)==sha(source/n) for n in ['states.csv','faces.csv','history.csv','initial_states.csv','initial_faces.csv','preparation.json']};rec['native_gate'] &= all(rec['byte_exact_to_frozen_reference'].values())
    rec['files']=[dict(path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in dest.iterdir() if p.is_file()];print(c['name'],rec['returncode'],rec['native_gate'],round(rec['wall_s'],3),rec.get('failure',{}),flush=True);return rec
def main():
    for p in [O,RAW,LOG]:p.mkdir(parents=True,exist_ok=True)
    base.LOG=LOG
    if '--matrix' not in sys.argv:
        assert not (O/'protocol.json').exists()
        guards=[case('control_l240',240),case('control_uniform',240,'uniform'),case('guard_f200',T=1),case('guard_f160',160,T=1)]
        matrix=[case('f200_dt1'),case('f160_dt1',160),case('f200_dt05',dt=.5),case('f200_dt025',dt=.25)]
        paths=[Path(__file__),Path(__file__).with_name('README.md'),R/'thin_film/bench/tf_curved_edge_bench.cpp',R/'thin_film/include/vf/thinfilm/CurvedEdgePreparation.h',R/'thin_film/include/vf/thinfilm/CurvedEdgeMotion.h',R/'thin_film/include/vf/thinfilm/CurvedEdge.h']
        protocol=dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),source_hashes={p.relative_to(R).as_posix():sha(p) for p in paths},binary_hashes={n:sha(base.BIN/(n+'.exe')) for n in ['tf_curved_edge_bench','test_edge_preparation','test_curved_edge']},guards=guards,conditional_matrix=matrix,physical_cut_m=240e-6,scope='Numerical fine-patch start200/160um. Original radius, physical initial trace/cut240um, far bath320um, materials, energy, force, mobility and tolerances unchanged.',prior_evidence='prepared_motion_20261113: minima221-230um outside old fine patch240um',condition='all18checks, both byte-exact controls and both1us guards pass; guard-wall*steps<=600s per case',gates={'native_balance':1e-10,'energy':1e-7,'saved_formula':1e-7,'space_dynamic':.1,'time_dynamic':.02},raw_directory=RAW.relative_to(R).as_posix(),observation_times_s=[10e-6,50e-6,200e-6],per_case_timeout_s=600,whole_compute_budget_s=3300,environment={'VF_THREADS':'2','caller_VF':'excluded'},sampling='All native scalar rows; three full states/faces at first step, each10us and final; initial state separately')
        save(O/'protocol.json',protocol);(O/'run_at_execution.py').write_bytes(Path(__file__).read_bytes());state=dict(state='guards_running',tests=[],cases=[]);save(O/'status.json',state)
        for name,expected in [('test_edge_preparation','4 passed, 0 failed'),('test_curved_edge','14 passed, 0 failed')]:
            rec=base.execute(name,[str(base.BIN/(name+'.exe'))],120);state['tests'].append(rec);save(O/'status.json',state)
            if rec['returncode']!=0 or expected not in (R/rec['log']).read_text():state.update(state='stopped',reason='unit failure');save(O/'status.json',state);return
        for c in guards:
            rec=execute(c,180 if c['T_us']==1 else 600);state['cases'].append(rec);save(O/'status.json',state)
            if not rec['native_gate']:break
        state['state']='guards_completed';save(O/'status.json',state);return
    protocol=json.loads((O/'protocol.json').read_text());state=json.loads((O/'status.json').read_text());assert state['state']=='guards_completed'
    for n,h in protocol['source_hashes'].items():assert sha(R/n)==h,n
    for n,h in protocol['binary_hashes'].items():assert sha(base.BIN/(n+'.exe'))==h,n
    guards={r['name']:r for r in state['cases']};used=sum(r['wall_s'] for r in state['cases']+state['tests']);state['skipped']=[];good=len(guards)==4 and all(r['native_gate'] for r in guards.values())
    for c in protocol['conditional_matrix']:
        guard=guards.get('guard_f'+str(c['refine_start_um']));estimate=(guard or {}).get('wall_s',1e9)*round(c['T_us']/c['dt_us'])
        if not good or estimate>600 or used>=3300:state['skipped'].append(dict(name=c['name'],reason='guard failed' if not good else 'frozen cost gate',estimated_s=estimate));continue
        rec=execute(c,min(600,3300-used));state['cases'].append(rec);used+=rec['wall_s'];save(O/'status.json',state)
    state.update(state='completed',compute_wall_s=used);save(O/'status.json',state)
if __name__=='__main__':main()
