"""Resolved full-circle training/validation geometry. Never creates blind tests."""
from pathlib import Path
import sys,os,time,json,hashlib,fcntl,traceback
import numpy as np
import torch

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'operator-ablation-20261122'/'shared'))
from gpu_lease import gpu_lease

DATA=HERE/'data';SPLITS=[('train',131072,202611270301,1000000000),('validation',16384,202611270302,2000000000)]
RADIUS_BINS=[4,8,16,32,64,128,256,512]
X32,W32=np.polynomial.legendre.leggauss(32);X64,W64=np.polynomial.legendre.leggauss(64)
CONFIG=dict(generator='resolved_complete_circle_v1',family='circle',stencil=7,radius_cells=[4,512],
    radius_distribution='96%loguniform4..512,2%at4,2%at512',angle='uniform0..2pi',offset='independent uniform[-.5,.5]^2',
    central_liquid_fraction=[.002,.998],phase='independent50%complement;one variant per base',
    courant_distribution='86%uniform[0,1],10%loguniform[1e-4,.05],2%exact0,2%exact1;storedfloat32 fixed before truth integration',
    flux_definition='right-face strip[.5-C,.5]x[-.5,.5],swept area/central liquid area;cellwidth1',
    curvature_definition='signedkappa*h=+1/R for liquid disk,-1/R for phase complement',
    inference_inputs='stencils only forcurvature;stencils+courants forflux;all geometry metadata offline only',
    training=[dict(split=s,count=n,seed=seed,base_id_start=base) for s,n,seed,base in SPLITS],
    precision='analyticGPU float64 integrals;float32 neural inputs/labels;independentCPU adaptive quadrature',
    device='CUDA;per-process allocation fraction .15;max2048 proposals',split_rule='independentseeds;no base or phase/Courant variant crosses splits',
    no_test_generated=True,ellipse='not included;circle only to retain fully verified truth')

def now():return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);temporary=path.with_suffix(path.suffix+'.tmp');temporary.write_text(json.dumps(value,indent=2));temporary.replace(path)
def status(**value):write(HERE/'data_status.json',dict(updated_utc=now(),pid=os.getpid(),**value))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def disk_rectangle(cx,cy,r,xlo,xhi,ylo,yhi):
    """Complete disk intersection with arbitrary per-cell/per-sample rectangles.

    Split at disk support endpoints and intersections with horizontal rectangle
    edges, then integrate the clipped top/bottom circular arcs analytically.
    No local graph approximation and no padding of the outer stencil.
    """
    cx=cx[:,None];cy=cy[:,None];r=r[:,None]
    xl,xh,yl,yh=torch.broadcast_tensors(xlo,xhi,ylo,yhi,cx)[0:4]
    rlo=((r-(yl-cy).abs())*(r+(yl-cy).abs())).clamp_min(0).sqrt()
    rhi=((r-(yh-cy).abs())*(r+(yh-cy).abs())).clamp_min(0).sqrt()
    cuts=torch.stack((xl,xh,(cx-r).expand_as(xl),(cx+r).expand_as(xl),cx-rlo,cx+rlo,cx-rhi,cx+rhi),-1)
    cuts=cuts.maximum(xl[...,None]).minimum(xh[...,None]).sort(-1).values
    left=cuts[...,:-1];right=cuts[...,1:];mid=(left+right)/2;width=right-left
    C=cx[...,None];Y=cy[...,None];RR=r[...,None]
    def primitive(x):
        v=(x-C).maximum(-RR).minimum(RR)
        radial=((RR-v.abs())*(RR+v.abs())).clamp_min(0).sqrt()
        return .5*(v*radial+RR*RR*torch.asin((v/RR).clamp(-1,1)))
    ds=primitive(right)-primitive(left)
    v=(mid-C).abs();radial=((RR-v)*(RR+v)).clamp_min(0).sqrt();top=Y+radial;bottom=Y-radial
    upper=torch.where(top<yh[...,None],Y*width+ds,yh[...,None]*width)
    lower=torch.where(bottom>yl[...,None],Y*width-ds,yl[...,None]*width)
    valid=(v<RR)&(torch.minimum(top,yh[...,None])>torch.maximum(bottom,yl[...,None]))
    return torch.where(valid,upper-lower,0.).sum(-1).clamp_min(0).minimum((xh-xl)*(yh-yl))

def cpu_oracle(cx,cy,r,xlo,xhi,ylo,yhi):
    """Independent numerical integration of the clipped vertical chord."""
    if xhi<=xlo:return 0.
    knots=[xlo,xhi,cx-r,cx+r]
    for y in (ylo,yhi):
        d=abs(y-cy)
        if d<r:
            q=np.sqrt((r-d)*(r+d));knots.extend((cx-q,cx+q))
    knots=np.unique(np.clip(knots,xlo,xhi));out=0.;reported=0.
    def integrand(x):
        v=abs(x-cx);a=np.sqrt(np.maximum(0,(r-v)*(r+v)))
        return np.where(v<r,np.maximum(0,np.minimum(yhi,cy+a)-np.maximum(ylo,cy-a)),0.)
    def adaptive(left,right,tolerance=2e-12,depth=0):
        midpoint=(left+right)/2;halfwidth=(right-left)/2
        low=halfwidth*np.dot(W32,integrand(midpoint+halfwidth*X32))
        high=halfwidth*np.dot(W64,integrand(midpoint+halfwidth*X64))
        if abs(high-low)<=tolerance or depth>=20:return high,abs(high-low)
        a,ea=adaptive(left,midpoint,tolerance/2,depth+1);b,eb=adaptive(midpoint,right,tolerance/2,depth+1)
        return a+b,ea+eb
    for left,right in zip(knots[:-1],knots[1:]):
        if right==left:continue
        part,error=adaptive(left,right)
        out+=part;reported+=error
    return out

def proposal(rng,n):
    radius=np.exp(rng.uniform(np.log(4),np.log(512),n));mode=rng.random(n)
    radius=np.where(mode<.02,4.,np.where(mode>.98,512.,radius))
    angle=rng.uniform(0,2*np.pi,n);ox=rng.uniform(-.5,.5,n);oy=rng.uniform(-.5,.5,n)
    return dict(radius_cells=radius,angle=angle,offset_x=ox,offset_y=oy,cx=radius*np.cos(angle)+ox,cy=radius*np.sin(angle)+oy)

def generate(split,count,seed,idstart):
    final=DATA/(split+'.npz');report=DATA/(split+'.json')
    if final.exists() and report.exists():return json.loads(report.read_text())
    partdir=DATA/(split+'_parts');partdir.mkdir(parents=True,exist_ok=True);resume_path=partdir/'state.json'
    rng=np.random.default_rng(seed);accepted=0;part_number=0
    if resume_path.exists():
        checkpoint=json.loads(resume_path.read_text());assert checkpoint['seed']==seed and checkpoint['count']==count
        accepted=checkpoint['accepted'];part_number=checkpoint['parts'];rng.bit_generator.state=checkpoint['rng']
    axis=torch.arange(-3,4,device='cuda',dtype=torch.float64);yy,xx=torch.meshgrid(axis.flip(0),axis,indexing='ij')
    xl=xx.flatten()[None,:]-.5;yl=yy.flatten()[None,:]-.5;started=time.monotonic()
    while accepted<count:
        n=min(2048,max(256,int((count-accepted)*1.1)));p=proposal(rng,n)
        cx,cy,r=[torch.from_numpy(p[key]).cuda() for key in ('cx','cy','radius_cells')]
        centre=disk_rectangle(cx,cy,r,torch.tensor([[-.5]],device='cuda',dtype=torch.float64),torch.tensor([[.5]],device='cuda',dtype=torch.float64),torch.tensor([[-.5]],device='cuda',dtype=torch.float64),torch.tensor([[.5]],device='cuda',dtype=torch.float64))[:,0]
        keep=((centre>=.002)&(centre<=.998)).cpu().numpy();indices=np.flatnonzero(keep)[:count-accepted]
        if len(indices)==0:continue
        p={key:value[indices] for key,value in p.items()};cx=cx[indices];cy=cy[indices];r=r[indices];inside_f=centre[indices];m=len(indices)
        comp=rng.random(m)<.5;mode=rng.random(m);c=rng.random(m).astype('f4')
        low=mode<.10;c[low]=np.exp(rng.uniform(np.log(1e-4),np.log(.05),int(low.sum()))).astype('f4')
        c[(mode>=.10)&(mode<.12)]=0.;c[(mode>=.12)&(mode<.14)]=1.
        cd=torch.from_numpy(c.astype('f8')).cuda();complement=torch.from_numpy(comp).cuda()
        x=disk_rectangle(cx,cy,r,xl,xl+1,yl,yl+1).reshape(-1,7,7)
        strip=disk_rectangle(cx,cy,r,(.5-cd)[:,None],torch.full_like(cd,.5)[:,None],torch.full_like(cd,-.5)[:,None],torch.full_like(cd,.5)[:,None])[:,0]
        f=torch.where(complement,1-inside_f,inside_f);swept=torch.where(complement,cd-strip,strip)
        x=torch.where(complement[:,None,None],1-x,x)
        # Endpoints use their exact geometric limits, not a learned/model output.
        swept=torch.where(cd==0,torch.zeros_like(swept),torch.where(cd==1,f,swept))
        ratio=(swept/f).clamp(0,1);q=torch.where(complement,-1/r,1/r)
        values=dict(stencils=x.float().cpu().numpy(),curvature=q.float().cpu().numpy(),courants=c,flux=ratio.float().cpu().numpy(),
            base_id=np.arange(idstart+accepted,idstart+accepted+m,dtype='i8'),family=np.full(m,'circle',dtype='U6'),complement=comp,
            center_fraction=f.cpu().numpy(),swept_volume=swept.cpu().numpy(),**p)
        assert np.isfinite(values['stencils']).all() and np.isfinite(values['flux']).all()
        temporary=partdir/f'part_{part_number:05d}.tmp.npz';np.savez(temporary,**values);temporary.replace(partdir/f'part_{part_number:05d}.npz')
        accepted+=m;part_number+=1;write(resume_path,dict(seed=seed,count=count,accepted=accepted,parts=part_number,rng=rng.bit_generator.state))
        status(state='generating',split=split,accepted=accepted,total=count,parts=part_number,elapsed_seconds=time.monotonic()-started)
    parts=[]
    for number in range(part_number):
        with np.load(partdir/f'part_{number:05d}.npz') as z:parts.append({key:z[key] for key in z.files})
    data={key:np.concatenate([part[key] for part in parts]) for key in parts[0]}
    assert len(data['base_id'])==count and len(np.unique(data['base_id']))==count
    checks=verify(data,split)
    temporary=DATA/(split+'.tmp.npz');np.savez(temporary,**data);temporary.replace(final)
    record=dict(passed=True,split=split,count=count,seed=seed,path=str(final.relative_to(ROOT)),sha256=sha(final),dtype='float32 neural tensors',
        radius_histogram=dict(edges=RADIUS_BINS,counts=np.histogram(data['radius_cells'],RADIUS_BINS)[0].tolist()),
        radius_min=float(data['radius_cells'].min()),radius_max=float(data['radius_cells'].max()),complement_count=int(data['complement'].sum()),
        courant_zero_count=int((data['courants']==0).sum()),courant_one_count=int((data['courants']==1).sum()),
        centre_fraction_min=float(data['center_fraction'].min()),centre_fraction_max=float(data['center_fraction'].max()),checks=checks,finished_utc=now())
    write(report,record);return record

def verify(data,split):
    radius=data['radius_cells'];f=data['center_fraction'];c=data['courants'];q=data['curvature']
    assert radius.min()>=4 and radius.max()<=512
    assert f.min()>=.002-1e-12 and f.max()<=.998+1e-12
    assert data['stencils'].shape==(len(f),7,7) and data['stencils'].dtype==np.float32
    assert np.max(abs(abs(q.astype('f8'))-1/radius))<2e-8
    ids=list(np.unique(np.linspace(0,len(f)-1,16).astype(int)))
    for lower,upper in zip(RADIUS_BINS[:-1],RADIUS_BINS[1:]):
        choices=np.flatnonzero((radius>=lower)&(radius<=upper))
        if len(choices):ids.extend(choices[np.linspace(0,len(choices)-1,3).astype(int)].tolist())
    for mask in (c==0,c==1,(c>0)&(c<.001)):
        choices=np.flatnonzero(mask)
        if len(choices):ids.extend(choices[np.linspace(0,len(choices)-1,3).astype(int)].tolist())
    ids.extend(np.argsort(f)[:4]);ids.extend(np.argsort(f)[-4:]);ids=np.unique(ids)
    max_fraction=max_swept=max_flux=0.;cases=[]
    for j,i in enumerate(ids):
        cx,cy,r=[float(data[name][i]) for name in ('cx','cy','radius_cells')];comp=bool(data['complement'][i]);courant=float(c[i])
        ref=np.array([[cpu_oracle(cx,cy,r,x-.5,x+.5,y-.5,y+.5) for x in range(-3,4)] for y in range(3,-4,-1)])
        inside_centre=ref[3,3];v=cpu_oracle(cx,cy,r,.5-courant,.5,-.5,.5)
        if comp:ref=1-ref;v=courant-v
        true_f=ref[3,3];expected=0. if courant==0 else (1. if courant==1 else v/true_f)
        e=float(np.max(abs(ref-data['stencils'][i])));s=abs(v-float(data['swept_volume'][i]));fe=abs(expected-float(data['flux'][i]))
        max_fraction=max(max_fraction,e);max_swept=max(max_swept,s);max_flux=max(max_flux,fe)
        cases.append(dict(index=int(i),base_id=int(data['base_id'][i]),radius_cells=r,courant=courant,complement=comp,fraction_max_abs=e,swept_abs=s,flux_ratio_abs=fe))
        status(state='independent_cpu_truth_check',split=split,checked=j+1,total=len(ids),fraction_max_abs=max_fraction,flux_ratio_max_abs=max_flux)
    assert max_fraction<2e-7 and max_flux<5e-7,(split,max_fraction,max_flux,max_swept)
    checks=dict(passed=True,method='independent adaptive32/64pointGauss numerical chord integration with analytic clipping breakpoints;notGPU primitive',geometries=len(ids),cell_integrals=49*len(ids),
        fraction_max_abs=max_fraction,swept_volume_max_abs=max_swept,flux_ratio_max_abs=max_flux,curvature_formula_max_abs=float(np.max(abs(abs(q.astype('f8'))-1/radius))),cases=cases)
    return checks

def geometry_keys(data):
    return np.ascontiguousarray(np.stack([data[key] for key in ('cx','cy','radius_cells')],1)).view(np.dtype((np.void,24))).reshape(-1)

def main():
    DATA.mkdir(parents=True,exist_ok=True);lock=(DATA/'generation.lock').open('a+')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:
        print(json.dumps(dict(state='already_running_no_duplicate')),flush=True);return
    original=json.loads((HERE/'protocol.json').read_text());assert original['scope']['radius_cells_min']==4
    config_path=DATA/'generation_protocol.json'
    if config_path.exists():assert json.loads(config_path.read_text())==CONFIG
    else:write(config_path,CONFIG)
    manifest=DATA/'manifest.json'
    if manifest.exists() and json.loads(manifest.read_text()).get('complete'):
        print(json.dumps(dict(state='already_complete_no_duplicate')),flush=True);return
    torch.set_num_threads(2);torch.backends.cuda.matmul.allow_tf32=False;torch.cuda.set_per_process_memory_fraction(.15)
    began=time.monotonic();status(state='preflight',gpu=torch.cuda.get_device_name(),memory_fraction=.15,command=sys.argv,config_path=str(config_path.relative_to(ROOT)))
    reports=[]
    for split,count,seed,idstart in SPLITS:
        with gpu_lease('03-plic-resolved-data-'+split,formal_timing=False):
            reports.append(generate(split,count,seed,idstart));torch.cuda.synchronize()
    with np.load(DATA/'train.npz') as source:train={name:source[name] for name in ('cx','cy','radius_cells')}
    with np.load(DATA/'validation.npz') as source:val={name:source[name] for name in ('cx','cy','radius_cells')}
    tk=geometry_keys(train);vk=geometry_keys(val);overlap=int(len(np.intersect1d(tk,vk)))
    assert overlap==0 and len(np.unique(tk))==len(tk) and len(np.unique(vk))==len(vk)
    note='\n'.join(['Resolved full-circle training data; no blind test created.',
        'train.npz:131072 bases,seed202611270301;validation.npz:16384 bases,seed202611270302. Variants stay within split;this version has one phase/Courant instance per base. Exact full-circle(cx,cy,R) keys have no duplicates or cross-split intersection.',
        'Radius R/h4..512,logarithmic coverage plus endpoint radii;uniform angle and independent subcell offset;central liquid fraction .002..998. No subgrid geometry,planes or ellipses included.',
        'Neural inputs:stencils float32[N,7,7],top row y=+3,x=-3..3;courants float32[N] forflux. Targets:curvature signedkappa*h=float32±1/R;flux float32=swept_liquid_area/central_liquid_area at right-face strip[.5-C,.5]x[-.5,.5].',
        'Metadata only:base_id,family,cx,cy,radius_cells,angle,offset_x,offset_y,complement,center_fraction,swept_volume. These must never be fed to the network.',
        'True full disk integrals in GPUfloat64;Courant rounded to storedfloat32 BEFORE strip integration;independent adaptive CPU chord integration checks full7x7,strip volume and ratio across radius/low-C/endpoint/boundary strata.',
        *[f"{r['split']}:fraction max_abs={r['checks']['fraction_max_abs']:.9g},swept volume={r['checks']['swept_volume_max_abs']:.9g},flux ratio={r['checks']['flux_ratio_max_abs']:.9g};{r['checks']['geometries']}CPU-checked geometries." for r in reports]])
    (DATA/'README.txt').write_text(note+'\n',encoding='utf-8')
    write(manifest,dict(complete=True,created_utc=now(),generator=CONFIG,files=reports,no_duplicate_full_geometry=True,cross_split_geometry_overlap=overlap,
        neural_input_keys=dict(flux=['stencils','courants'],curvature=['stencils']),target_keys=dict(flux='flux',curvature='curvature'),blind_test_generated=False))
    status(state='complete',files=[r['path'] for r in reports],train_count=131072,validation_count=16384,cross_split_geometry_overlap=overlap,elapsed_seconds=time.monotonic()-began)
    print(json.dumps(dict(state='complete',seconds=time.monotonic()-began)),flush=True)

if __name__=='__main__':
    try:main()
    except BaseException as error:
        failure=dict(state='failed',error=repr(error),traceback=traceback.format_exc(),updated_utc=now(),pid=os.getpid())
        DATA.mkdir(parents=True,exist_ok=True)
        with (DATA/'failures.jsonl').open('a') as stream:stream.write(json.dumps(failure)+'\n')
        write(HERE/'data_status.json',failure);raise
