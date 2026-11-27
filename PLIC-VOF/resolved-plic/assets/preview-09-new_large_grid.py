"""Independent true-grid fixtures, generated only after frozen flux deployment."""
from pathlib import Path
import datetime, hashlib, importlib.util, json, os, time
import numpy as np

R=Path(__file__).resolve().parent
W=R.parent
OUT=R/'flux_blind_test/large_grid'
SEED=202611270400
NS=[512,1024,2048,4096,8192]
SCENES=['single_drop','many_drops']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=-7))).isoformat()
def save(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(obj,indent=2));temp.replace(path)

def circles(scene):
    rng=np.random.default_rng(SEED+(0 if scene=='single_drop' else 1))
    if scene=='single_drop':
        return np.array([[*rng.uniform(.49,.51,2),rng.uniform(.20,.23)]])
    xx,yy=np.meshgrid((np.arange(16)+.5)/16,(np.arange(16)+.5)/16)
    centers=np.c_[xx.ravel(),yy.ravel()]+rng.uniform(-.002,.002,(256,2))
    return np.c_[centers,rng.uniform(.014,.018,256)]

def main():
    lock=R/'deployment_selection_locked.json'
    assert lock.exists(),'No new grid before deployment lock'
    choice=json.loads(lock.read_text())
    assert choice['operators']==['flux'] and choice['all_seeds_retained']
    source=W/'agent_experiments/large_grid_operator_fixtures.py'
    identity=dict(seed=SEED,resolutions=NS,scenes=SCENES,deployment_lock_sha256=sha(lock),generator_helper_sha256=sha(source),generator_sha256=sha(__file__))
    registration=OUT/'registration.json'
    if registration.exists():assert json.loads(registration.read_text())['identity']==identity
    else:save(registration,dict(identity=identity,registered=now(),new_after_lock=True,description='Actual full Cartesian float32 VOF fields with unique interface queries;fixed physical scenes refined over all N;no repeated validation rows.',scope='R/h>=4 for every circle;center VOF gate [.002,.998];Courant0.5;single local flux inference only',accuracy='Primary transported-volume MAE normalized by cell area,without multiplying by h^2;other flux metrics auxiliary.',selection='No architecture,runtime,weight or seed selection from these fields.'))
    spec=importlib.util.spec_from_file_location('_new_true_grid_helper',source)
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    helper.ROOT=W;helper.R=OUT;helper.scene_circles=circles
    oldcircles=np.load(W/'gpu-cpu-large-grid-20261125/fixtures/many_drops_N4096.npz')['circles']
    assert not np.array_equal(oldcircles,circles('many_drops'))
    for scene in SCENES:
        assert np.min(circles(scene)[:,2])*min(NS)>=4
    start=time.monotonic();records=[]
    for scene in SCENES:
        for n in NS:
            save(OUT/'status.json',dict(state='generating',pid=os.getpid(),scene=scene,N=n,completed=len(records),elapsed_seconds=time.monotonic()-start,updated=now()))
            rec=helper.build(scene,n)
            rec['minimum_radius_cells']=float(np.min(circles(scene)[:,2])*n)
            rec['path_sha256']=sha(W/rec['source'])
            records.append(rec)
            print(json.dumps(dict(scene=scene,N=n,queries=rec['active_interface_cells'],elapsed_seconds=rec['elapsed_seconds'])),flush=True)
    save(OUT/'manifest.json',dict(state='complete',identity=identity,cases=records,elapsed_seconds=time.monotonic()-start,finished=now()))
    save(OUT/'status.json',dict(state='complete',pid=os.getpid(),cases=len(records),elapsed_seconds=time.monotonic()-start,updated=now()))

if __name__=='__main__':main()
