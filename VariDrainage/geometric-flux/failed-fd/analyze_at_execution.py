"""Native geometry observer vs independent face kinematics, GEOS FD and quadrature."""
from pathlib import Path
import sys,json,csv,time,subprocess,hashlib
import numpy as np,shapely
from shapely import box,MultiPoint,STRtree
sys.path.insert(0,str(Path(__file__).parents[1]/'moving_voronoi'))
import run as base
R=base.R;O=R/'results/kinematic_flux_convex_20261113';RAW=O/'native';DOMAIN=box(0,0,320e-6,160e-6);RIGHT=box(240e-6,0,320e-6,160e-6)
POLYGON_STATS=dict(calls=0,invalid_raw_polygons=0,max_hull_area_change_over_domain=0.,max_hull_relative_cell_area_change=0.)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,data):p.write_text(json.dumps(data,indent=2,allow_nan=False),encoding='utf-8')
def polys(xy):
    raw=np.array(shapely.get_parts(shapely.voronoi_polygons(MultiPoint(xy),extend_to=DOMAIN,ordered=True)),dtype=object);hull=shapely.convex_hull(raw);difference=abs(shapely.area(hull)-shapely.area(raw));relative_change=float(max(difference/np.maximum(shapely.area(hull),1e-60)))
    assert relative_change<1e-10,'Convex representation repair exceeds preregistered cell-area bound'
    POLYGON_STATS['calls']+=1;POLYGON_STATS['invalid_raw_polygons']+=int(sum(~shapely.is_valid(raw)));POLYGON_STATS['max_hull_area_change_over_domain']=max(POLYGON_STATS['max_hull_area_change_over_domain'],float(max(difference)/DOMAIN.area));POLYGON_STATS['max_hull_relative_cell_area_change']=max(POLYGON_STATS['max_hull_relative_cell_area_change'],relative_change)
    return shapely.intersection(hull,DOMAIN)
def geometry(xy,velocity):
    cells=polys(xy);area=shapely.area(cells);right=shapely.area(shapely.intersection(cells,RIGHT));weight=right/area;i,j=STRtree(cells).query(cells,predicate='intersects');ok=i<j;i,j=i[ok],j[ok];face=shapely.intersection(cells[i],cells[j]);length=shapely.length(face);ok=length>1e-14;i,j,face,length=i[ok],j[ok],face[ok],length[ok]
    mid=shapely.get_coordinates(shapely.centroid(face));d=xy[j]-xy[i];dist=np.linalg.norm(d,axis=1);normal=d/dist[:,None];skew=mid-.5*(xy[i]+xy[j]);speed=np.sum(.5*(velocity[i]+velocity[j])*normal,axis=1)-np.sum((velocity[j]-velocity[i])*skew,axis=1)/dist
    clipped=shapely.intersection(face,RIGHT);rlength=shapely.length(clipped);rmid=mid.copy();active=rlength>1e-14;rmid[active]=shapely.get_coordinates(shapely.centroid(clipped[active]));rspeed=np.sum(.5*(velocity[i]+velocity[j])*normal,axis=1)-np.sum((velocity[j]-velocity[i])*(rmid-.5*(xy[i]+xy[j])),axis=1)/dist
    Adot=np.zeros(len(xy));Rdot=Adot.copy();np.add.at(Adot,i,length*speed);np.add.at(Adot,j,-length*speed);np.add.at(Rdot,i,rlength*rspeed);np.add.at(Rdot,j,-rlength*rspeed)
    wdot=(Rdot-weight*Adot)/area
    return dict(area=area,weight=weight,area_rate=Adot,right_rate=Rdot,weight_rate=wdot,faces=len(i))
def fd(xy,velocity,epsilon):
    plus=polys(xy+epsilon*velocity);minus=polys(xy-epsilon*velocity);ap=shapely.area(plus);am=shapely.area(minus);wp=shapely.area(shapely.intersection(plus,RIGHT))/ap;wm=shapely.area(shapely.intersection(minus,RIGHT))/am
    return (ap-am)/(2*epsilon),(wp-wm)/(2*epsilon)
def relative(a,b):return float(np.max(np.abs(a-b))/max(1e-60,float(np.max(np.abs(b)))))
def native(name,xy,velocity):
    source=RAW/(name+'-input.csv');target=RAW/(name+'-native.csv');assert not source.exists() and not target.exists()
    with source.open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['id','x_m','y_m','ux_m_s','uy_m_s']);writer.writerows([[i,*xy[i],*velocity[i]] for i in range(len(xy))])
    rec=base.execute(name,[str(base.BIN/'tf_voronoi_kinematics.exe'),str(source),'0.00024',str(target)],5);assert rec['returncode']==0,rec
    return np.atleast_1d(np.genfromtxt(target,delimiter=',',names=True)),rec
def verify(name,xy,velocity):
    geo=geometry(xy,velocity);obs,record=native(name,xy,velocity);eps=1e-5*np.sqrt(float(min(geo['area'])))/max(1e-30,float(max(np.linalg.norm(velocity,axis=1))));a1,w1=fd(xy,velocity,eps);a2,w2=fd(xy,velocity,eps/2)
    errors=dict(native_area=relative(geo['area'],obs['A_m2']),native_weight=float(max(abs(geo['weight']-obs['w']))),native_Jv=relative(geo['area_rate'],obs['Jv_m2_s']),area_FD=relative(a2,geo['area_rate']),weight_FD=relative(w2,geo['weight_rate']),FD_level_change=relative(a1,a2),global_area_sum=float(abs(sum(geo['area_rate']))/max(1e-60,np.sum(abs(geo['area_rate'])))))
    passed=max(errors['native_area'],errors['native_weight'])<1e-8 and max(errors[k] for k in ['native_Jv','area_FD','weight_FD'])<1e-6
    return dict(name=name,errors=errors,passed=passed,FD_epsilon_s=eps/2,native=record),geo
def integrated(x0,x1,dt,count):
    nodes,weights=np.polynomial.legendre.leggauss(count);answer=np.zeros(len(x0));velocity=(x1-x0)/dt
    for node,weight in zip(nodes,weights):answer+=.5*weight*geometry(x0+.5*(node+1)*(x1-x0),velocity)['weight_rate']
    return answer
def main():
    start=time.monotonic();O.mkdir(exist_ok=True);RAW.mkdir(exist_ok=True);(O/'logs').mkdir(exist_ok=True);base.LOG=O/'logs';assert not (O/'protocol.json').exists()
    sources=[Path(__file__),Path(__file__).with_name('README.md'),R/'thin_film/bench/tf_voronoi_kinematics.cpp',R/'thin_film/include/vf/thinfilm/MovingVoronoiFilm.h'];save(O/'protocol.json',dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),source_hashes={p.relative_to(R).as_posix():sha(p) for p in sources},binary_sha256=sha(base.BIN/'tf_voronoi_kinematics.exe'),cases=['obs_dt1','control_dt025','dt003125'],samples_us=list(range(10,201,10)),fixed_cut_m=240e-6,whole_budget_s=600,native_timeout_s=5,gates=dict(geometry=1e-8,area_and_weight_derivatives=1e-6,liquid_surface_and_quadrature=1e-5),scope='Native geometry-only observations of frozen states. Independent face-speed shape derivatives and GEOS FD. Not continuum normal-flux validation.',reference='Springel2010AREPOsec3.3eq32-34, author PDF locally preserved and page17 visually checked.'))
    (O/'analyze_at_execution.py').write_bytes(Path(__file__).read_bytes());state=dict(state='running',guards=[],samples=[]);save(O/'status.json',state)
    grid=np.array([[(i+.5)*320e-6/9,(j+.5)*160e-6/5] for j in range(5) for i in range(9)])
    for jitter in [0.,.08]:
        xy=grid.copy();k=np.arange(len(xy));xy[:,0]+=jitter*320e-6/9*np.sin(k*1.7);xy[:,1]+=jitter*160e-6/5*np.cos(k*2.3);x,y=(xy-np.array([160e-6,80e-6])).T;X,Y=np.pi*xy[:,0]/320e-6,np.pi*xy[:,1]/160e-6
        fields=dict(translation=np.tile([1e-3,.4e-3],(len(xy),1)),rotation=20*np.c_[-y,x],strain=np.c_[10*x+7*y,-3*x-5*y],zero_wall=1e-3*np.c_[np.sin(X)**2*np.sin(2*Y),-.5*np.sin(2*X)*np.sin(Y)**2])
        for name,velocity in fields.items():
            result,_=verify(f'guard_{name}_j{int(jitter*100)}',xy,velocity);state['guards'].append(result);save(O/'status.json',state);print(result['name'],result['passed'],flush=True)
    if not all(r['passed'] for r in state['guards']):state.update(state='guard_failed');save(O/'status.json',state);return
    sources=[('obs_dt1','mechanical_work_fixed_20261113'),('control_dt025','time_refinement_20261113'),('dt003125','time_refinement_20261113')];inputs=[]
    for name,campaign in sources:
        status=json.loads((R/'results'/campaign/'status.json').read_text());rec=next(r for r in status['cases'] if r['name']==name);assert rec['completed']
        for f in rec['files']:assert sha(R/f['path'])==f['sha256']
        inputs.append(dict(campaign=campaign,name=name,files=rec['files']));folder=R/'results/local-runs'/campaign/name;states=np.atleast_1d(np.genfromtxt(folder/'states.csv',delimiter=',',names=True));dt=rec['summary']['dt_s']
        for step in sorted(set(map(int,states['step']))):
            if step*dt<9.999e-6:continue
            assert time.monotonic()-start<600
            a,b=[states[(states['step']==step)&(states['stage']==stage)] for stage in [0,1]];x0=np.c_[a['x_m'],a['y_m']];x1=np.c_[b['x_m'],b['y_m']];velocity=(x1-x0)/dt;test,geo=verify(f'{name}_s{step}',.5*(x0+x1),velocity);q3=integrated(x0,x1,dt,3);q5=integrated(x0,x1,dt,5);endpoint=(b['w']-a['w'])/dt
            item=dict(case=name,step=step,t_us=round(step*dt*1e6,8),dt_s=dt,geometry_check=test,coordinate_velocity_roundtrip_relative=relative(velocity,np.c_[b['ux_m_s'],b['uy_m_s']]))
            for label,amount in [('liquid',a['V_m3']),('surface',a['S_mol'])]:
                native_rate=float(amount@endpoint);r3=float(amount@q3);r5=float(amount@q5);scale=max(abs(native_rate),float(np.sum(abs(amount*q5))),1e-60)
                item[label]=dict(native_endpoint_rate=native_rate,face_quadrature_rate=r5,gauss3_rate=r3,normalized_endpoint_error=abs(r5-native_rate)/scale,normalized_quadrature_change=abs(r5-r3)/scale,geometric_boundary_term=float(np.sum(amount/geo['area']*geo['right_rate'])),density_area_term=float(-np.sum(amount/geo['area']*geo['weight']*geo['area_rate'])))
            item['passed']=test['passed'] and max(item[k][v] for k in ['liquid','surface'] for v in ['normalized_endpoint_error','normalized_quadrature_change'])<1e-5;state['samples'].append(item);save(O/'status.json',state)
        print(name,'saved mechanics samples complete',flush=True)
    save(O/'inputs.json',inputs);series={}
    for name,_ in sources:
        rows=[r for r in state['samples'] if r['case']==name];metrics={}
        for key in ['liquid','surface']:
            a=np.array([r[key]['native_endpoint_rate'] for r in rows]);b=np.array([r[key]['face_quadrature_rate'] for r in rows]);metrics[key+'_relative_L2']=float(np.linalg.norm(a-b)/max(1e-60,np.linalg.norm(a)))
        series[name]=dict(samples=len(rows),passed=sum(r['passed'] for r in rows),**metrics)
    passed=all(r['passed'] for r in state['samples']) and all(r[key]<1e-5 for r in series.values() for key in ['liquid_relative_L2','surface_relative_L2']);save(O/'summary.json',dict(series=series,guards=len(state['guards']),guard_passed=sum(r['passed'] for r in state['guards']),passed=passed,polygon_representation=POLYGON_STATS,full_C1_complete=False,scope='Analytical geometric ownership rate verified; not continuum physical normal flux or an added closure.'))
    state.update(state='completed',passed=passed,wall_s=time.monotonic()-start);save(O/'status.json',state);print(json.dumps(series,indent=2))
if __name__=='__main__':
    try:main()
    except Exception as exc:
        if (O/'status.json').exists():
            state=json.loads((O/'status.json').read_text());state.update(state='failed',exception=repr(exc));save(O/'status.json',state)
        raise
