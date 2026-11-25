"""Actual Cartesian VOF fields and exact local circular-interface targets."""
from pathlib import Path
import json,time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
R=ROOT/'gpu-cpu-large-grid-20261125'
NS=[512,1024,2048,4096]

def save(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix('.tmp');temp.write_text(json.dumps(obj,indent=2))
    for attempt in range(30):
        try:temp.replace(path);return
        except PermissionError:
            if attempt==29:raise
            time.sleep(.05)

def circle_area(x0,x1,y0,y1,r):
    """Exact disk/rectangle intersection, split at every branch change."""
    x0,x1,y0,y1,r=np.broadcast_arrays(*[np.asarray(v,dtype=np.float64) for v in (x0,x1,y0,y1,r)])
    lo=np.clip(x0,-r,r);hi=np.maximum(lo,np.clip(x1,-r,r))
    roots0=np.sqrt(np.maximum(r*r-y0*y0,0));roots1=np.sqrt(np.maximum(r*r-y1*y1,0))
    nodes=np.stack([lo,hi,-roots0,roots0,-roots1,roots1],-1)
    nodes=np.sort(np.maximum(lo[...,None],np.minimum(hi[...,None],nodes)),axis=-1)
    a=nodes[...,:-1];b=nodes[...,1:];mid=(a+b)*.5;rr=r[...,None]
    h=np.sqrt(np.maximum(rr*rr-mid*mid,0));dy=b-a
    def primitive(x):
        return .5*(x*np.sqrt(np.maximum(rr*rr-x*x,0))+rr*rr*np.arcsin(np.clip(x/rr,-1,1)))
    integral=primitive(b)-primitive(a)
    upper=np.where(y1[...,None]<h,y1[...,None]*dy,integral)
    lower=np.where(y0[...,None]>-h,y0[...,None]*dy,-integral)
    area=np.where(np.minimum(y1[...,None],h)>np.maximum(y0[...,None],-h),upper-lower,0).sum(-1)
    return np.maximum(area,0)

def scene_circles(scene):
    if scene=='single_drop':return np.array([[.5013,.4987,.22]])
    rng=np.random.default_rng(202611250401)
    xx,yy=np.meshgrid((np.arange(16)+.5)/16,(np.arange(16)+.5)/16)
    centers=np.c_[xx.ravel(),yy.ravel()]+rng.uniform(-.002,.002,(256,2))
    return np.c_[centers,rng.uniform(.014,.018,256)]

def build(scene,n):
    path=R/'fixtures'/f'{scene}_N{n}.npz';meta=path.with_suffix('.json')
    if path.exists() and meta.exists():return json.loads(meta.read_text())
    begin=time.perf_counter();circles=scene_circles(scene);field=np.zeros((n,n),dtype=np.float32);owner=np.full((n,n),-1,dtype=np.int16)
    for ci,(cx,cy,radius) in enumerate(circles):
        cx*=n;cy*=n;radius*=n
        xa=max(0,int(np.floor(cx-radius)));xb=min(n,int(np.ceil(cx+radius)))
        ya=max(0,int(np.floor(cy-radius)));yb=min(n,int(np.ceil(cy+radius)))
        xx,yy=np.meshgrid(np.arange(xa,xb)+.5-cx,np.arange(ya,yb)+.5-cy)
        nearest=np.maximum(abs(xx)-.5,0)**2+np.maximum(abs(yy)-.5,0)**2
        farthest=(abs(xx)+.5)**2+(abs(yy)+.5)**2
        inside=farthest<=radius*radius;partial=(nearest<radius*radius)&~inside
        patch=np.zeros(xx.shape,dtype=np.float64);patch[inside]=1
        patch[partial]=circle_area(xx[partial]-.5,xx[partial]+.5,yy[partial]-.5,yy[partial]+.5,radius)
        patch=np.clip(patch,0,1);view=field[ya:yb,xa:xb]
        assert not np.any((view>0)&(patch>0)),'Disks overlap'
        view[:]=np.maximum(view,patch.astype('f4'));own=owner[ya:yb,xa:xb];own[patch>0]=ci
    active=(field>=.002)&(field<=.998);iy,ix=np.nonzero(active)
    assert min(iy.min(),ix.min())>=5 and max(iy.max(),ix.max())<n-5
    indices=owner[iy,ix];c=circles[indices];x=ix+.5-c[:,0]*n;y=iy+.5-c[:,1]*n;radius=c[:,2]*n
    area=circle_area(x-.5,x+.5,y-.5,y+.5,radius)
    swept=circle_area(x,x+.5,y-.5,y+.5,radius)
    fractions=field[iy,ix].astype(np.float64)
    assert np.max(abs(fractions-area))<5e-7
    # An independent split Gauss integral checks the analytical primitive.
    nodes,weights=np.polynomial.legendre.leggauss(256);rng=np.random.default_rng(202611250405+n)
    oracle_error=0.
    for i in rng.choice(len(x),min(32,len(x)),replace=False):
        rr=radius[i];lo=max(x[i]-.5,-rr);hi=min(x[i]+.5,rr)
        cuts=[lo,hi]
        for yy in [y[i]-.5,y[i]+.5]:
            if abs(yy)<rr:
                xx=np.sqrt(rr*rr-yy*yy);cuts.extend([max(lo,min(hi,-xx)),max(lo,min(hi,xx))])
        cuts=np.sort(cuts);approx=0.
        for a,b in zip(cuts[:-1],cuts[1:]):
            xx=(a+b)*.5+(b-a)*.5*nodes;hh=np.sqrt(np.maximum(rr*rr-xx*xx,0))
            approx+=(b-a)*.5*np.dot(weights,np.maximum(np.minimum(y[i]+.5,hh)-np.maximum(y[i]-.5,-hh),0))
        oracle_error=max(oracle_error,abs(approx-area[i]))
    assert oracle_error<2e-6,oracle_error
    rows=iy[:,None,None]+np.arange(5,-6,-1)[None,:,None]
    cols=ix[:,None,None]+np.arange(-5,6)[None,None,:]
    stencil=np.ascontiguousarray(field[rows,cols],dtype=np.float32)
    assert np.array_equal(stencil[:,5,5],field[iy,ix])
    path.parent.mkdir(exist_ok=True)
    np.savez_compressed(path,field=field,stencils=stencil,ix=ix,iy=iy,circle_id=indices,circles=circles,
        center_fraction=fractions,oracle_center_fraction=area,flux_ratio_truth=swept/area,
        flux_density_truth=swept/.5,q_truth=1/radius,kappa_truth=1/c[:,2])
    result=dict(scene=scene,N=n,total_grid_cells=n*n,active_interface_cells=len(ix),circle_count=len(circles),
        source=str(path.relative_to(ROOT)),genuine_field=True,stencil_width=11,subcrops=[3,5,7,9,11],
        center_gate=[.002,.998],CFL=.5,truth='Analytical circle/rectangle integrals;known signed circle curvature',
        physical_domain=[0,1,0,1],water='inside disjoint circles',row_order='high-y to low-y',
        area_mass_error=abs(float(field.sum(dtype=np.float64))/n**2-float(np.pi*(circles[:,2]**2).sum())),
        independent_quadrature_checks=32,independent_quadrature_max_abs=oracle_error,
        elapsed_seconds=time.perf_counter()-begin)
    save(meta,result);return result

def main():
    R.mkdir(exist_ok=True);cases=[];start=time.perf_counter()
    for scene in ['single_drop','many_drops']:
        for n in NS:
            save(R/'fixture_status.json',dict(state='running',scene=scene,N=n,completed=len(cases)))
            m=build(scene,n);cases.append(m);print(json.dumps(m),flush=True)
    save(R/'fixture_manifest.json',dict(state='complete',cases=cases,elapsed_seconds=time.perf_counter()-start))
    save(R/'fixture_status.json',dict(state='complete',cases=len(cases),elapsed_seconds=time.perf_counter()-start))

if __name__=='__main__':main()
