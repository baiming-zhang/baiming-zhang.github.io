"""Registered stage2: same-source singleC/16C and raw-q objectives.

Only training/validation. No deployment selection, blind data, timing or model
routing. records() reads plans only and has no CUDA/model-construction side effect.
"""
from pathlib import Path
import os,sys
sys.dont_write_bytecode=True
for key in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
import argparse,copy,fcntl,hashlib,importlib.util,json,math,subprocess,time,traceback
import numpy as np
import torch
from torch import nn
from torch.nn import functional as F
from torch.func import functional_call,stack_module_state
R=Path(__file__).resolve().parent;ROOT=R.parent;W=ROOT.parent
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(W/'operator-ablation-20261122/shared'))
from gpu_lease import gpu_lease
SEEDS=[4101,4102,4103];STEPS=12000;CHECK=200

def now():return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False))
    for i in range(40):
        try:tmp.replace(p);return
        except PermissionError:
            if i==39:raise
            time.sleep(.05)
def save(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=Path(str(p)+'.tmp');torch.save(obj,tmp);tmp.replace(p)
def state(**kw):write(R/'status.json',dict(pid=os.getpid(),updated_utc=now(),**kw))

def records():
    p=read(R/'protocol.json');plans={v['group']:v for v in read(ROOT/'training_protocol.json')['groups']};out=[]
    assert p['planned_new_models']==24 and p['flux']['steps']==STEPS and p['curvature']['steps']==STEPS
    for source in p['flux']['source_groups']:
        for condition in p['flux']['conditions']:
            out.append(dict(group=source+'__'+condition,operator='flux',condition=condition,source_group=source,config=plans[source]['config'],sources=[str(ROOT/'groups'/source/f'seed{s}/best.pt') for s in SEEDS],init='stage1_checkpoint',guard_reference_path=str(ROOT/'groups'/source/'reference_validation.json')))
    source=p['curvature']['source_group']
    for condition in p['curvature']['conditions']:
        out.append(dict(group='curvature__'+condition,operator='curvature',condition=condition,source_group=source,config=plans[source]['config'],sources=[str(ROOT/'groups'/source/f'seed{s}/best.pt') for s in SEEDS],init='stage1_checkpoint',guard_reference_path=str(ROOT/'groups'/source/'reference_validation.json')))
    assert len(out)==8 and len({v['group'] for v in out})==8;return out

def build(config):
    if config['architecture']=='compact_monotone_flux_v1':
        from compact_models import build_model
        return build_model(config)
    name='_stage2_original_models'
    if name not in sys.modules:
        path=W/'operator-ablation-20261122/03-data-stencil/models.py';spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m)
    return sys.modules[name].OperatorMLP(config)

def load_source(rec,seed,device='cpu'):
    path=Path(rec['sources'][SEEDS.index(seed)]);ck=torch.load(path,map_location='cpu',weights_only=False);assert ck['config']==rec['config'] and ck['seed']==seed
    model=build(ck['config']);model.load_state_dict(ck['model'],strict=True);return model.to(device),ck

class FluxCurveTraining(nn.Module):
    """Training-only shared mass evaluation. State keys equal original model keys."""
    def __init__(self,model,config):
        super().__init__();self.net=model.net;self.k=int(config['stencil_size']);self.num_bins=int(config['num_bins']);self.register_buffer('bin_starts',model.bin_starts.detach().clone())
    def forward(self,x,c):
        n=x.shape[0];z=self.net(torch.cat((x,x.flip(1)),0).reshape(2*n,self.k*self.k)).reshape(2,n,self.num_bins)
        masses=(F.softplus(z)+1e-6).mean(0);cc=c.clamp(0,1)
        coverage=(self.num_bins*cc[:,:,None]-self.bin_starts[None,:,:]).clamp(0,1)
        ratio=(coverage*masses[:,None,:]).sum(-1)/masses.sum(-1)[:,None]
        f=x[:,self.k//2,self.k//2].clamp(0,1)[:,None];safe=f.clamp_min(1e-12);lo=((f+cc-1)/safe).clamp(0,1);hi=(cc/safe).clamp(0,1)
        return torch.where(f>0,torch.minimum(torch.maximum(ratio,lo),hi),torch.zeros_like(ratio))

def cpu_preflight(rec):
    target=R/'groups'/rec['group']/'preflight.json';identity=dict(sources=[sha(p) for p in rec['sources']],config=rec['config'],operator=rec['operator'],condition=rec['condition'])
    if target.exists() and read(target).get('identity')==identity and read(target).get('passed'):return
    torch.set_num_threads(1);models=[]
    for seed in SEEDS:
        m,_=load_source(rec,seed);m=m.double().eval();models.append(FluxCurveTraining(m,rec['config']) if rec['operator']=='flux' else m)
    k=rec['config']['stencil_size'];gen=torch.Generator().manual_seed(202611270521);x=torch.rand((3,7,k,k),generator=gen,dtype=torch.float64);c=torch.linspace(.03125,.96875,16,dtype=torch.float64)[None,:].expand(7,-1).clone();target_y=torch.rand((3,7,16),generator=gen,dtype=torch.float64)
    p,b=stack_module_state(models);base=copy.deepcopy(models[0]).to('meta');op=rec['operator'];diff=0.;grad_delta=0.
    if op=='flux':
        # Exact keys are exported without a training wrapper/core prefix.
        original,_=load_source(rec,SEEDS[0]);assert set(models[0].state_dict())==set(original.state_dict());ordinary=[]
        for seed in SEEDS:
            mm,_=load_source(rec,seed);ordinary.append(mm.double().eval())
        fn=torch.vmap(lambda pp,bb,xx:functional_call(base,(pp,bb),(xx,c)),in_dims=(0,0,0));fast=fn(p,b,x);f=x[:,:,k//2,k//2,None];ls=(f*(fast-target_y)).abs().mean((-1,-2));bg=torch.autograd.grad(ls.sum(),list(p.values()))
        for i,model in enumerate(ordinary):
            xx=x[i,:,None,:,:].expand(-1,16,-1,-1).reshape(-1,k,k);cc=c.reshape(-1);slow=model(xx,cc).reshape(7,16);diff=max(diff,float((slow-fast[i]).abs().max()));loss=(f[i]*(slow-target_y[i])).abs().mean();sg=torch.autograd.grad(loss,list(model.parameters()))
            grad_delta=max(grad_delta,max(float((a[i]-z).abs().max())/max(1.,float(z.abs().max())) for a,z in zip(bg,sg)))
    else:
        fn=torch.vmap(lambda pp,bb,xx:functional_call(base,(pp,bb),(xx,)),in_dims=(0,0,0));fast=fn(p,b,x);y=target_y[:,:,0];errors=fast-y;ls=errors.square().mean(-1) if rec['condition']=='unweighted_q_MSE' else errors.abs().mean(-1);bg=torch.autograd.grad(ls.sum(),list(p.values()))
        for i,model in enumerate(models):
            slow=model(x[i]);diff=max(diff,float(abs(slow-fast[i]).max()));e=slow-y[i];loss=e.square().mean() if rec['condition']=='unweighted_q_MSE' else e.abs().mean();sg=torch.autograd.grad(loss,list(model.parameters()));grad_delta=max(grad_delta,max(float((a[i]-z).abs().max())/max(1.,float(z.abs().max())) for a,z in zip(bg,sg)))
    assert torch.isfinite(fast).all() and all(torch.isfinite(z).all() for z in bg)
    assert diff<1e-10 and grad_delta<1e-9,(rec['group'],diff,grad_delta)
    write(target,dict(passed=True,identity=identity,device='CPU_float64',multic_vs_repeated_forward_max_abs=diff,per_seed_sum_batched_vs_serial_gradient_relative_max=grad_delta,training_state_keys='identical original net.* and bin_starts;no core prefix in saved weights',inference_architecture_changed=False,created_utc=now()))

def load_data():
    source=W/'operator-ablation-20261122/03-data-stencil/data';out={};identity={}
    for split in ['train','validation']:
        oldpath=source/f'enhanced_{split}.npz';newpath=ROOT/'data'/f'{split}.npz';multipath=ROOT/'data'/f'{split}_multic.npz'
        with np.load(oldpath) as z:
            keep=abs(z['curvature'])<=.25;old={k:torch.as_tensor(z[k][keep],device='cuda') for k in ['stencils','courants','flux','curvature']}
        with np.load(newpath) as z:
            baseid=z['base_id'].copy();new={k:torch.as_tensor(z[k],device='cuda') for k in ['stencils','courants','flux','curvature']}
        with np.load(multipath) as z:
            assert np.array_equal(baseid,z['base_id']),('base_id misalignment',split)
            assert z['flux_cdf'].shape==(len(baseid),16) and z['courant_grid'].shape==(16,)
            new['multi_c']=torch.as_tensor(z['courant_grid'],device='cuda');new['multi_y']=torch.as_tensor(z['flux_cdf'],device='cuda')
        out[split]=(old,new);identity[split]=dict(old_sha256=sha(oldpath),new_sha256=sha(newpath),multi_sha256=sha(multipath),old_bases=len(old['flux']),new_bases=len(new['flux']),exact_base_id_alignment=True)
    return out,identity

def crop(x,k):
    a=(7-k)//2;return x[:,a:a+k,a:a+k]

def augmented(data,count,g,op,multi=False):
    idx=torch.randint(len(data['flux']),(count,),device='cuda',generator=g);x=data['stencils'][idx];f=x[:,3,3];q=data['curvature'][idx]
    if op=='flux':
        c=data['multi_c'][None,:].expand(count,-1) if multi else data['courants'][idx,None];y=data['multi_y'][idx] if multi else data['flux'][idx,None]
    else:c=data['courants'][idx,None];y=q
    comp=torch.rand(count,device='cuda',generator=g)<.5
    if op=='flux':y=torch.where(comp[:,None],(c-f[:,None]*y)/(1-f[:,None]).clamp_min(1e-8),y)
    else:y=torch.where(comp,-y,y)
    x=torch.where(comp[:,None,None],1-x,x);flip=torch.rand(count,device='cuda',generator=g)<.5;x=torch.where(flip[:,None,None],x.flip(1),x)
    if op=='curvature':
        trans=torch.rand(count,device='cuda',generator=g)<.5;x=torch.where(trans[:,None,None],x.transpose(-1,-2),x)
    return x,c,y

@torch.no_grad()
def ordinary_metrics(params,buffers,config,data,op):
    normal=build(config).to('meta');k=config['stencil_size']
    if op=='flux':fn=torch.vmap(lambda p,b,x,c:functional_call(normal,(p,b),(x,c)),in_dims=(0,0,None,None))
    else:fn=torch.vmap(lambda p,b,x,c:functional_call(normal,(p,b),(x,)),in_dims=(0,0,None,None))
    pp=[]
    for first in range(0,len(data['flux']),2048):pp.append(fn(params,buffers,crop(data['stencils'][first:first+2048],k),data['courants'][first:first+2048]))
    predictions=torch.cat(pp,1);truth=data[op];f=data['stencils'][:,3,3];weak=abs(data['curvature'])<.001;plane=data['curvature']==0;out=[]
    for pred in predictions:
        valid=torch.isfinite(pred);count=int(valid.sum());n=len(pred);e=abs(pred[valid]-truth[valid]);weight=f[valid] if op=='flux' else torch.ones_like(e);w=weak[valid];p=plane[valid]
        out.append(dict(primary_mae=float((weight*e).mean()) if count else None,label_mae=float(e.mean()) if count else None,weak_mae=float(e[w].mean()) if w.any() else None,plane_mae=float(e[p].mean()) if p.any() else None,coverage=count/n,valid_count=count,samples=n,full_coverage=(count==n),max_error=float(e.max()) if count else None))
    return out

def validation(params,buffers,rec,oldval,newval,guard):
    old=ordinary_metrics(params,buffers,rec['config'],oldval,rec['operator']);new=ordinary_metrics(params,buffers,rec['config'],newval,rec['operator']);out=[]
    for i,(a,b) in enumerate(zip(old,new)):
        ref=guard[i];eligible=a['full_coverage'] and b['full_coverage'] and a['primary_mae']<=1.1*ref['primary_mae']+1e-7
        if rec['operator']=='curvature':
            for field in ['weak_mae','plane_mae']:
                if ref[field] is not None:eligible=eligible and a[field] is not None and a[field]<=1.1*ref[field]+1e-7
        out.append(dict(eligible=bool(eligible),circle_mae=b['primary_mae'] if b['full_coverage'] else None,old_mae=a['primary_mae'] if a['full_coverage'] else None,old=a,new=b,guard_reference='original_enhanced7_stage1_reference_not_cumulative'))
    return out

def selection_key(v):return (not v['eligible'],float('inf') if v['circle_mae'] is None else v['circle_mae'],float('inf') if v['old_mae'] is None else v['old_mae'])

def group(rec,data,data_identity):
    dest=R/'groups'/rec['group'];dest.mkdir(parents=True,exist_ok=True);status_path=dest/'status.json'
    if status_path.exists() and read(status_path).get('state')=='complete':return
    identity=dict(sources=[sha(p) for p in rec['sources']],data=data_identity,config=rec['config'],condition=rec['condition'],protocol_sha256=sha(R/'protocol.json'))
    guard=read(rec['guard_reference_path'])['old'];assert len(guard)==3
    for z in guard:assert z.get('valid_count',z['samples'])==z['samples'] and z['coverage']==1
    oldval,newval=data['validation'];models=[]
    for seed in SEEDS:
        model,_=load_source(rec,seed,'cuda');models.append(FluxCurveTraining(model,rec['config']) if rec['operator']=='flux' else model)
    params,buffers=stack_module_state(models);base=copy.deepcopy(models[0]).to('meta');del models
    optimizer=torch.optim.AdamW(params.values(),lr=3e-4,weight_decay=1e-5);generators=[torch.Generator(device='cuda').manual_seed(seed+270200) for seed in SEEDS]
    curves=[[] for _ in SEEDS];best=[(True,float('inf'),float('inf')) for _ in SEEDS];step=0;prior=0.;checkpoint_path=dest/'last.pt'
    if checkpoint_path.exists():
        ck=torch.load(checkpoint_path,map_location='cuda',weights_only=False);assert ck['identity']==identity
        for k,v in params.items():v.data.copy_(ck['params'][k])
        for k,v in buffers.items():v.data.copy_(ck['buffers'][k])
        optimizer.load_state_dict(ck['optimizer']);step=ck['step'];best=ck['bests'];curves=ck['curves'];prior=ck['elapsed_seconds']
        for g,s in zip(generators,ck['rng']):g.set_state(s.cpu())
    if not (dest/'source_start_validation.json').exists():write(dest/'source_start_validation.json',dict(identity=identity,source_selected_validation=validation(params,buffers,rec,oldval,newval,guard),original_enhanced7_guard=guard))
    start=time.monotonic();oldtrain,newtrain=data['train'];op=rec['operator'];k=rec['config']['stencil_size']
    if op=='flux':fn=torch.vmap(lambda p,b,x,c:functional_call(base,(p,b),(x,c)),in_dims=(0,0,0,0))
    else:fn=torch.vmap(lambda p,b,x,c:functional_call(base,(p,b),(x,)),in_dims=(0,0,0,0))
    def checkpoint():
        values=validation(params,buffers,rec,oldval,newval,guard);elapsed=prior+time.monotonic()-start
        for i,seed in enumerate(SEEDS):
            curves[i].append(dict(step=step,elapsed_seconds=elapsed,lr=optimizer.param_groups[0]['lr'],validation=values[i]))
            if selection_key(values[i])<tuple(best[i]):
                best[i]=selection_key(values[i]);weights={name:value[i].detach().cpu().clone() for name,value in params.items()}|{name:value[i].detach().cpu().clone() for name,value in buffers.items()}
                assert set(weights)==set(build(rec['config']).state_dict())
                save(dest/f'seed{seed}/best.pt',dict(config=rec['config'],seed=seed,model=weights,step=step,selection=values[i],source=rec,stage='stage2_courant_objective',identity=identity))
            write(dest/f'seed{seed}/curve.json',curves[i])
        save(checkpoint_path,dict(identity=identity,step=step,params={name:v.detach().clone() for name,v in params.items()},buffers={name:v.detach().clone() for name,v in buffers.items()},optimizer=optimizer.state_dict(),bests=best,curves=curves,rng=[g.get_state() for g in generators],elapsed_seconds=elapsed))
        st=dict(state='complete' if step==STEPS else 'training',group=rec['group'],operator=op,condition=rec['condition'],step=step,total=STEPS,elapsed_seconds=elapsed,estimated_remaining_seconds=elapsed/max(step,1)*(STEPS-step),validation=values,pid=os.getpid(),updated_utc=now());write(status_path,st);state(**{a:b for a,b in st.items() if a not in ['pid','updated_utc']});print(json.dumps(dict(group=rec['group'],step=step,circle_mae=[v['circle_mae'] for v in values],eligible=[v['eligible'] for v in values],elapsed_seconds=elapsed)),flush=True)
    if step==0:checkpoint()
    while step<STEPS:
        old_batches=[augmented(oldtrain,512,g,op) for g in generators];new_batches=[augmented(newtrain,256 if op=='flux' else 512,g,op,op=='flux' and rec['condition']=='sixteen_courants') for g in generators]
        def tensors(batches):return torch.stack([crop(b[0],k) for b in batches]),torch.stack([b[1] for b in batches]),torch.stack([b[2] for b in batches])
        xo,co,yo=tensors(old_batches);xn,cn,yn=tensors(new_batches);lr=3e-6+.5*(3e-4-3e-6)*(1+math.cos(math.pi*step/(STEPS-1)))
        for pg in optimizer.param_groups:pg['lr']=lr
        po=fn(params,buffers,xo,co);pn=fn(params,buffers,xn,cn)
        if op=='flux':
            fo=xo[:,:,k//2,k//2,None];fnew=xn[:,:,k//2,k//2,None];loss=.5*(fo*(po-yo)).abs().mean((-1,-2))+.5*(fnew*(pn-yn)).abs().mean((-1,-2))
        elif rec['condition']=='unweighted_q_MSE':loss=.5*(po-yo).square().mean(-1)+.5*(pn-yn).square().mean(-1)
        else:loss=.5*(po-yo).abs().mean(-1)+.5*(pn-yn).abs().mean(-1)
        assert torch.isfinite(loss).all(),(rec['group'],step,loss)
        optimizer.zero_grad(set_to_none=True);loss.sum().backward();optimizer.step();step+=1
        if step%CHECK==0 or step==STEPS:checkpoint()
    assert [sha(p) for p in rec['sources']]==identity['sources'],'Stage1 source changed'
    del params,buffers,base,optimizer;torch.cuda.empty_cache()

def main(args):
    if args.launch:
        folder=R/'groups';folder.mkdir(parents=True,exist_ok=True);pidpath=folder/'controller.pid'
        if pidpath.exists():
            pid=int(pidpath.read_text());proc=Path(f'/proc/{pid}/cmdline')
            if proc.exists() and 'stage2_courant_objective/train_stage2.py' in proc.read_text():print('already_running',pid);return
        stream=(folder/'controller.log').open('a');p=subprocess.Popen([sys.executable,'-B','-u',str(Path(__file__).resolve())],stdout=stream,stderr=stream,stdin=subprocess.DEVNULL,start_new_session=True,cwd=W);pidpath.write_text(str(p.pid));print(p.pid);return
    guard=(R/'training.lock').open('a+')
    try:fcntl.flock(guard,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:print('Existingstage2trainingretained');return
    if (R/'selection_locked.json').exists():state(state='training_complete',cached=True,models=24);return
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    try:os.sched_setaffinity(0,{max(os.sched_getaffinity(0))})
    except (OSError,AttributeError):pass
    plans=records()
    for rec in plans:cpu_preflight(rec)
    state(state='cpu_preflight_complete',groups=8,models=24,gpu_started=False)
    if args.preflight_only:return
    while not (ROOT/'data/multic_manifest.json').exists() or read(ROOT/'data/multic_manifest.json').get('state')!='complete':state(state='waiting_for_multic_labels',gpu_started=False);time.sleep(10)
    with gpu_lease('stage2-data-resident',formal_timing=False):
        torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False;torch.cuda.set_per_process_memory_fraction(.15);data,identity=load_data();torch.cuda.synchronize()
    write(R/'training_protocol.json',dict(registered_protocol_sha256=sha(R/'protocol.json'),groups=plans,models=24,steps=STEPS,checkpoint_every=CHECK,data=identity,optimizer='AdamW reset;lr3e-4->3e-6 cosine,wd1e-5',per_seed_loss='sum over3;16C mean overcurve targets perbase,not16x weight',augmentation_rng='seed+270200;identical draw counts/sequence forsingleCvs16C andqMSEvsqL1',guard='Originalenhanced7 stage1 reference,10%+1e-7;notcumulative',validation='Ordinary original model singleCourant forward only;no accelerateddecoder usedforselection',inference='Originalmodelstate keys,no trainingwrapper or metadata',blind_test_generated=False,formal_timing=False))
    for rec in plans:
        state(state='waiting_for_gpu_lease',group=rec['group'])
        with gpu_lease('stage2-train-'+rec['group'],formal_timing=False):group(rec,data,identity);torch.cuda.synchronize()
    ranking=[]
    for rec in plans:
        cks=[torch.load(R/'groups'/rec['group']/f'seed{s}/best.pt',map_location='cpu',weights_only=False) for s in SEEDS]
        ranking.append(dict(group=rec['group'],operator=rec['operator'],condition=rec['condition'],source_group=rec['source_group'],all_seed_retention_pass=all(c['selection']['eligible'] for c in cks),circle_mae=float(np.mean([c['selection']['circle_mae'] for c in cks])),selected_steps=[c['step'] for c in cks]))
    write(R/'selection_locked.json',dict(created=now(),models=24,ranking=ranking,selection='Validation only;all24 checkpoint choices retained;deployment Pareto lock/anyblinddata notcreatedhere'))
    state(state='training_complete',models=24,groups=8,selection_locked=True,blind_test_generated=False)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--launch',action='store_true');parser.add_argument('--preflight-only',action='store_true');args=parser.parse_args()
    try:main(args)
    except BaseException as error:
        failure=dict(state='failed',pid=os.getpid(),updated_utc=now(),error=repr(error),traceback=traceback.format_exc());write(R/'groups'/f'failure_{int(time.time())}.json',failure);state(state='failed',error=repr(error));raise
