"""Read-only export of actual completed-task BC/load inputs; no FE solves."""
from pathlib import Path
import base64,gzip,hashlib,json,math,re,time
import numpy as np
W=Path(__file__).resolve().parents[1];SITE=W/'showcase-public-site/GeoBC-show';OUT=SITE/'assets/fe-bc';OUT.mkdir(parents=True,exist_ok=True)
SUB=np.array([[0,3,5],[3,1,4],[5,4,2],[3,4,5]])
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(a,t='<f4'):return base64.b64encode(np.asarray(a,dtype=t).tobytes()).decode()
def save(name,d):
    p=OUT/(name+'.json.gz');p.write_bytes(gzip.compress(json.dumps(d,separators=(',',':'),ensure_ascii=False).encode(),mtime=0));return dict(asset=p.name,bytes=p.stat().st_size)
def fields(p,ids):
    a=np.load(p);assert np.array_equal(ids,a['node_ids']);return a['U_peak_mm']
def glyph(p,u,d,kind,value=None):return dict(p=np.asarray(p).tolist(),u=np.asarray(u).tolist(),d=np.asarray(d).tolist(),kind=kind,value=value)
def sample_nodes(x,ids,count=16):
    # Deterministic farthest-point sampling, always selecting actual nodes.
    ids=np.asarray(ids);q=x[ids];sel=[int(np.argmin(q[:,0]))];distance=np.full(len(ids),np.inf)
    for _ in range(min(count,len(ids))-1):
        distance=np.minimum(distance,np.sum((q-q[sel[-1]])**2,axis=1));sel.append(int(np.argmax(distance)))
    return ids[sel]
def triangulate(faces,remap):return remap[faces[:,SUB]].reshape(-1,3)

def rod_geometry(m,id):
    x=m['xyz'];faces=m['faces'];groups=m['face_group_ids'];support=np.unique(faces[groups==1]);pin=np.unique(faces[groups==2]);normal=m['radial_normal']
    members=np.flatnonzero(np.linalg.norm(normal,axis=1)>.5);assert np.array_equal(members,support)
    nodes=np.union1d(support,pin);lookup=np.full(len(x),-1);lookup[nodes]=np.arange(len(nodes))
    q=x[support,:2];circle=np.linalg.lstsq(np.c_[2*q,np.ones(len(q))],np.sum(q*q,axis=1),rcond=None)[0];L=float(circle[0]);assert abs(circle[1])<1e-7
    center=x[m['cells'][:,:4]].mean(1);valid=np.flatnonzero((center[:,0]>.16*L)&(center[:,0]<.84*L))
    cellids=valid[np.linspace(0,len(valid)-1,8,dtype=int)];N=np.r_[np.full(4,-.125),np.full(6,.25)];body_pos=np.einsum('j,ijk->ik',N,x[m['cells'][cellids]])
    d=dict(id=id,kind='rod',xyz=enc(x[nodes]),groups=[dict(kind='constraint',tri=enc(triangulate(faces[groups==1],lookup),'<u4')),dict(kind='pressure',tri=enc(triangulate(faces[groups==2],lookup),'<u4'))],support_nodes=len(support),support_facets=int((groups==1).sum()),pressure_facets=int((groups==2).sum()),boundary_text=['大头销孔各节点：nᵣ·U = 0（双向径向销支承）','大头孔面：〈uθ〉A = 0，〈uz〉A = 0（两个面积均值定位）','切向/轴向不是逐点固定；无整面夹死、墙面或额外铰链'],samples={})
    return d,dict(nodes=nodes,support=sample_nodes(x,support,16),pin=pin,L=L,body_cells=cellids,body_N=N,body_pos=body_pos)

def rod_sample(m,ctx,U,load,state,load_path,force_path):
    x=m['xyz'];n=m['radial_normal'];axis=math.atan2(*np.array(load['small_pin_target_force_N'])[::-1])+math.radians(load['pressure_axis_quadrature_correction_deg']);direction=np.array([math.cos(axis),math.sin(axis)])
    width=float(re.search(r'z/([0-9.eE+\-]+)',load['pressure_law'])[1]);p0=float(load['pressure_peak_coefficient_MPa']);skew=float(load['width_skew']);pin=ctx['pin'];er=x[pin,:2]/np.linalg.norm(x[pin,:2],axis=1)[:,None];p=p0*np.maximum(0,er@direction)*(1+skew*x[pin,2]/width)
    chosen=sample_nodes(x,pin[p>p0*.08],14);arrows=[]
    for i in ctx['support']:arrows.append(glyph(x[i],U[i],n[i]/np.linalg.norm(n[i]),'constraint'))
    for i in chosen:
        radial=np.r_[x[i,:2]/np.linalg.norm(x[i,:2]),0.];v=p0*max(0,np.dot(radial[:2],direction))*(1+skew*x[i,2]/width);arrows.append(glyph(x[i],U[i],radial,'pressure',float(v)))
    snap=load.get('snapshot');body=np.zeros(3)
    if snap:
        omega=snap['omega_rod_rad_s'];alpha=snap['alpha_rod_rad_s2'];a0=snap['a_small_local_mm_s2'];bu=np.einsum('j,ijk->ik',ctx['body_N'],U[m['cells'][ctx['body_cells']]])
        for point,up in zip(ctx['body_pos'],bu):
            density_force=-4.43e-9*np.array([a0[0]-omega**2*point[0]-alpha*point[1],a0[1]+alpha*point[0]-omega**2*point[1],0.]);length=np.linalg.norm(density_force)
            if length>0:arrows.append(glyph(point,up,density_force/length,'body',float(length)))
    F=np.load(force_path);fn=F['force_N'] if 'force_N' in F else F['F_N'];assert np.array_equal(m['node_ids'],F['node_ids'])
    pinforce=np.array(load['small_pin_assembled_force_N']);body=np.array(load['body_force_N']);assert np.linalg.norm(fn.sum(0)-pinforce-body)<3e-6
    radial_error=float(abs(np.sum(U*n,axis=1)).max());assert state=='input_only' or radial_error<5e-5
    means=[float(np.sum(m[key]*U)) for key in ['mean_tangent','mean_axial'] if key in m]
    assert state=='input_only' or not means or max(abs(v) for v in means)<1e-5
    return dict(state=state,u=enc(U[ctx['nodes']]),glyphs=arrows,p0_MPa=p0,pressure_direction=direction.tolist(),half_width_mm=width,skew=skew,pin_force_N=pinforce.tolist(),body_force_N=body.tolist(),snapshot=snap,radial_constraint_error_mm=radial_error,mean_constraint_error_mm=means,load_sha256=sha(load_path),nodal_force_sha256=sha(force_path))

def export_rods():
    manifest=read(SITE/'data/r04-viewer-manifest.json');protocol=read(W/'data/rod_three_solver_0006_v3/protocol.json');index={};samples=0
    for case in manifest['cases']:
        id=case['id'];root=W/'data/rod_fenicsx_r04_pilot_v1/cases'/id;gpath=root/'geometry.npz';m=dict(np.load(gpath));d,ctx=rod_geometry(m,id);d['geometry_sha256']=sha(gpath)
        for entry in [s for s in manifest['samples'] if s['case_id']==id]:
            scenario=entry['scenario'];folder=root/scenario
            if entry['outcome']=='reused':
                assert id=='rodv2_0006' and scenario=='combined_bending';load_path=W/'data/rod_remesh_comparison_0006_v3/solution/load.json';force_path=W/'data/rod_remesh_comparison_0006_v3/solution/input_force.npz';fpath=W/'data/rod_three_solver_0006_v3/fenicsx/fields.npz'
            else:load_path=folder/'load.json';force_path=folder/'load.npz';fpath=folder/'fenicsx/fields.npz'
            state='input_only' if entry['outcome']=='failed' else 'accepted';U=np.zeros_like(m['xyz']) if state=='input_only' else fields(fpath,m['node_ids'])
            item=rod_sample(m,ctx,U,read(load_path),state,load_path,force_path);item['outcome']=entry['outcome'];d['samples'][scenario]=item;samples+=1
        if id=='rodv2_0006':
            fine=np.load(W/'data/rod_three_solver_0006_v3/inputs.npz');assert np.array_equal(m['xyz'],fine['xyz']) and np.array_equal(m['cells'],fine['cells'])
            lp=W/'data/rod_remesh_comparison_0006_v3/solution/load.json';fp=W/'data/rod_remesh_comparison_0006_v3/solution/input_force.npz'
            for key,path in [('cal_peak',W/'publication/GeoBC-Neural-Operator/results/rod_cal_peak_snapshot_20261024/fine_peak_snapshot.npz'),('fen_peak',W/'data/rod_three_solver_0006_v3/fenicsx/fields.npz')]:d['samples'][key]=rod_sample(m,ctx,fields(path,m['node_ids']),read(lp),'accepted',lp,fp)
        index[id]=save(id,d)
    id='rodv2_0006-coarse';meshpath=W/'data/rod_v2_sixload_v1/cases/rodv2_0006/h3/mesh.npz';m=dict(np.load(meshpath));bc=np.load(meshpath.parent/'constraints.npz');m['radial_normal']=bc['radial_normal'];d,ctx=rod_geometry(m,id);folder=W/'data/rod_mixed_1000_v4_bending15/pilot/rodv2_0006/combined_bending/target';lp=folder/'load.json';fp=folder/'input_force.npz'
    # The frozen coarse deck contains the same two distributed equations.
    text=(folder/'rod.inp').read_text();large=[int(v) for v in re.findall(r'\*EQUATION\s*\n(\d+)',text) if int(v)>3];assert len(large)==2
    d['distributed_equation_terms']=large;d['geometry_sha256']=sha(meshpath);d['samples']['cal_peak']=rod_sample(m,ctx,fields(folder/'fields.npz',m['node_ids']),read(lp),'accepted',lp,fp);index[id]=save(id,d)
    return index,samples

def export_patch():
    root=W/'data/rod_three_solver_0006_v1/benchmark';p=root/'inputs.npz';m=np.load(p);x=m['xyz'];U=fields(root/'fenicsx/fields.npz',m['node_ids']);F=m['force_N'];normal=m['radial_normal'];arrows=[]
    for i in range(len(x)):
        arrows.append(glyph(x[i],U[i],normal[i],'constraint'))
        if np.linalg.norm(F[i])>1e-12:arrows.append(glyph(x[i],U[i],F[i]/np.linalg.norm(F[i]),'nodal_force',float(np.linalg.norm(F[i]))))
    fs=np.array([[1,2,3,5,8,9],[0,3,2,7,8,6],[0,1,3,4,9,7],[0,2,1,6,5,4]]);tri=fs[:,SUB].reshape(-1,3)
    d=dict(id='r02-reference-patch',kind='patch',title='R02 独立曲边单元核对试件',xyz=enc(x),base_tri=enc(tri,'<u4'),groups=[],support_nodes=10,boundary_text=['10 个节点：实际 nᵢ·Uᵢ=0 线性方向约束','另有两条均值线性方程；这里不是实体销孔/墙面/铰链','人为构造的节点力用于数值核对，不是发动机载荷'],constraint_rows=[m['mean_tangent'].tolist(),m['mean_axial'].tolist()],force_N=F.tolist(),geometry_sha256=sha(p),samples={'peak':dict(state='accepted',u=enc(U),glyphs=arrows,pin_force_N=F.sum(0).tolist(),body_force_N=[0.,0.,0.])})
    return save(d['id'],d)

def export_c01():
    root=W/'data/copv/c01_geometry_v1';p=root/'fenicsx_geometry.npz';m=np.load(p);xy=m['xy'];ij=m['node_ij'];src=np.load(root/'fenicsx_smoke/fields.npz');ur=src['U_rz_mm'];assert np.array_equal(src['xy'],xy)
    nt=32;theta=np.arange(nt)*2*np.pi/nt;c=np.cos(theta);s=np.sin(theta);x=np.stack([xy[:,0,None]*c,xy[:,0,None]*s,np.broadcast_to(xy[:,1,None],(len(xy),nt))],axis=-1).reshape(-1,3);U=np.stack([ur[:,0,None]*c,ur[:,0,None]*s,np.broadcast_to(ur[:,1,None],(len(xy),nt))],axis=-1).reshape(-1,3)
    def ring_surface(ids):
        tris=[]
        for a,b in zip(ids[:-1],ids[1:]):
            for j in range(nt):
                k=(j+1)%nt;tris.extend([[a*nt+j,b*nt+j,b*nt+k],[a*nt+j,b*nt+k,a*nt+k]])
        return np.asarray(tris,dtype=np.uint32)
    inner=np.flatnonzero(ij[:,0]==ij[:,0].min());inner=inner[np.argsort(ij[inner,1])];outer=np.flatnonzero(ij[:,0]==ij[:,0].max());outer=outer[np.argsort(ij[outer,1])];bottom=np.flatnonzero(ij[:,1]==0);bottom=bottom[np.argsort(ij[bottom,0])];top=np.flatnonzero(ij[:,1]==ij[:,1].max());top=top[np.argsort(ij[top,0])]
    pressure=ring_surface(inner);support=ring_surface(bottom);alltri=np.concatenate([pressure,support,ring_surface(outer),ring_surface(top)]);arrows=[]
    for j in range(0,nt,4):
        i=int(bottom[-1])*nt+j;arrows.append(glyph(x[i],U[i],[0,0,1],'constraint'))
    assert sha(p)==read(SITE/'data/copv-c01-gpu-fenicsx.json')['input_geometry_sha256']
    for k in np.linspace(2,len(inner)-3,11,dtype=int):
        edge=min(int(k)//2,(len(inner)-3)//2);xi=(int(k)-2*edge)/2;tangent=np.array([4*xi-3,4-8*xi,4*xi-1])@xy[inner[2*edge:2*edge+3]];normal=np.array([tangent[1],-tangent[0]]);normal/=np.linalg.norm(normal)
        for j in [0,8,16,24]:
            i=int(inner[k])*nt+j;d=[normal[0]*c[j],normal[0]*s[j],normal[1]];arrows.append(glyph(x[i],U[i],d,'pressure',1.))
    d=dict(id='c01-smoke',kind='axisym',title='C01 v1 历史 GPU 接入试件',xyz=enc(x),base_tri=enc(alltri,'<u4'),groups=[dict(kind='constraint',tri=enc(support,'<u4')),dict(kind='pressure',tri=enc(pressure,'<u4'))],geometry_sha256=sha(p),support_nodes=len(bottom),boundary_text=['底端整条经线边：uz=0；径向 ur 自由','上端自由；未施加闭口端合力','仅为 2 层各向同性轴对称接入试件；绕轴显示不是独立三维求解'],original_nodes=len(xy),original_cells=len(m['cells']),samples={'peak':dict(state='accepted',u=enc(U),glyphs=arrows,p0_MPa=1.,pin_force_N=None,body_force_N=[0.,0.,0.],snapshot=None)})
    return save(d['id'],d)

def main():
    start=time.perf_counter();index,count=export_rods();index['r02-reference-patch']=export_patch();index['c01-smoke']=export_c01();report=dict(schema='completed_fe_bc_loads_v1',cases=index,r04_inputs=count,no_fe_recomputed=True,seconds=time.perf_counter()-start)
    (OUT/'index.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(report))
if __name__=='__main__':main()
