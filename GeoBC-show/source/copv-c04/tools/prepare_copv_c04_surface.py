"""Original TET10 face extraction and independent closed-pressure load check."""
from pathlib import Path
import argparse,json,time
import numpy as np
from geobc.copv_structured_mesh import analytic_segments
from geobc.copv_geometry import DEFAULT_PARAMETERS

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--mesh',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--parameters',type=Path,required=True);args=ap.parse_args();params=json.loads(args.parameters.read_text(encoding='utf-8'));end=params['straight_length_mm']/2+params['dome_depth_mm']+params['neck_length_mm'];liner_end=end-3.;upper_z=end+22.;lower_z=-(end+15.);upper_bore=.68*params['boss_bore_radius_mm'];lower_bore=.55*params['boss_bore_radius_mm'];start=time.perf_counter();out=args.output;out.mkdir(parents=True,exist_ok=True);m=np.load(args.mesh);x=m['xyz'];c=m['cells'];pids=m['cell_part_ids'];local=np.array([[1,2,3,5,8,9],[0,3,2,7,8,6],[0,1,3,4,9,7],[0,2,1,6,5,4]])
    faces=c[:,local].reshape(-1,6);owners=np.repeat(np.arange(len(c)),4);key=np.ascontiguousarray(np.sort(faces[:,:3],axis=1)).view(np.dtype((np.void,24))).ravel();order=np.argsort(key);sortedkey=key[order];startidx=np.r_[0,np.flatnonzero(sortedkey[1:]!=sortedkey[:-1])+1];counts=np.diff(np.r_[startidx,len(order)])
    if counts.max()>2:raise RuntimeError('Nonmanifold original tetra faces')
    ext=order[startidx[counts==1]];shared=startidx[counts==2];pair0=order[shared];pair1=order[shared+1];different=pids[owners[pair0]]!=pids[owners[pair1]];interfaces=np.r_[pair0[different],pair1[different]];selected=np.r_[ext,interfaces]
    centroid=x[faces[ext,:3]].mean(1);area=.5*np.cross(x[faces[ext,1]]-x[faces[ext,0]],x[faces[ext,2]]-x[faces[ext,0]]);opposite=x[c[owners[ext],np.tile(np.arange(4),len(c))[ext]]];flip=np.sum(area*(opposite-centroid),axis=1)>0;area[flip]*=-1
    rz,_=analytic_segments(params,(4,2500,16,2500,4));z=centroid[:,2];r=np.linalg.norm(centroid[:,:2],axis=1);inner=np.interp(abs(z),rz[rz[:,1]>=0,1],rz[rz[:,1]>=0,0]);inner=np.where(z>liner_end,upper_bore,np.where(z<-liner_end,lower_bore,inner));tags=np.full(len(ext),4);tags[r<=inner+.30]=1;tags[np.isclose(z,upper_z,atol=1e-5)]=2;tags[np.isclose(z,lower_z,atol=1e-5)]=3
    pressure=-area[tags==1];pressure_center=centroid[tags==1];corner_ids=np.unique(c[:,:4]);xx=x[corner_ids];caps=[]
    for zend,bore in [(upper_z,upper_bore),(lower_z,lower_bore)]:
        pp=xx[np.isclose(xx[:,2],zend,atol=1e-7)&np.isclose(np.linalg.norm(xx[:,:2],axis=1),bore,atol=1e-7),:2];assert len(pp)>=3;pp=pp[np.argsort(np.arctan2(pp[:,1],pp[:,0]))];caps.append(float(.5*abs(np.sum(pp[:,0]*np.roll(pp[:,1],-1)-pp[:,1]*np.roll(pp[:,0],-1)))))
    cap_top,cap_bottom=caps;force=pressure.sum(0)+[0,0,cap_top-cap_bottom];moment=np.cross(pressure_center,pressure).sum(0)
    for tag,value in [(2,cap_top),(3,-cap_bottom)]:
        mask=tags==tag;weight=np.linalg.norm(area[mask],axis=1);loads=np.c_[np.zeros((len(weight),2)),weight/weight.sum()*value];moment+=np.cross(centroid[mask],loads).sum(0)
    scale=np.pi*params['inner_radius_mm']**2;relative=float(np.linalg.norm(force)/scale);report={'pressure_boundary_facets':int((tags==1).sum()),'top_cap_facets':int((tags==2).sum()),'bottom_cap_facets':int((tags==3).sum()),'pressure_plus_caps_force_resultant_N':force.tolist(),'pressure_plus_caps_moment_resultant_Nmm':moment.tolist(),'force_relative_to_p_pi_R2':relative,'moment_relative_to_p_pi_R2_L':float(np.linalg.norm(moment)/(scale*(upper_z-lower_z))),'method':'Independent outward tetra face topology and exact planar triangle pressure area vectors; same idealized fluid boundary rule; caps each applied once','passed':relative<1e-4 and np.linalg.norm(moment)/(scale*(upper_z-lower_z))<1e-4,'exterior_faces':len(ext),'material_interface_faces':len(interfaces),'seconds':time.perf_counter()-start}
    np.savez_compressed(out/'surface.npz',faces=faces[selected],owner_cells=owners[selected],part_ids=pids[owners[selected]],is_exterior=np.r_[np.ones(len(ext),bool),np.zeros(len(interfaces),bool)])
    report['passed']=bool(report['passed']);(out/'load_closure_review.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    if not report['passed']:raise RuntimeError('Pressure/cap closure failed independently of gauge reactions')
if __name__=='__main__':main()
