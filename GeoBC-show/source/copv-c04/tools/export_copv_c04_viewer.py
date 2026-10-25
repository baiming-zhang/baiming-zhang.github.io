"""Create scientific display payloads from original mesh and solved fields.

Tri6 faces use four display triangles. Stress color is original-cell
volume-weighted GP mean, never averaged across a material/ply interface.
Cut faces evaluate the solved FE displacement. P1 midside output values are
affine interpolation and introduce no independent quadratic unknowns.
"""
from pathlib import Path
import argparse,base64,gzip,json,time
import numpy as np

EDGES=[(0,1),(1,2),(2,0),(0,3),(2,3),(1,3)]
def enc(a,dtype):return base64.b64encode(np.asarray(a,dtype=dtype).tobytes()).decode()
def save(p,value):p.write_bytes(gzip.compress(json.dumps(value,separators=(',',':'),ensure_ascii=False).encode(),mtime=0))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();root=args.root;out=args.output;out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();m=np.load(root/'mesh_h6_netgen/mesh.npz');a=np.load(root/'fenicsx/nodal_fields.npz');s=np.load(root/'fenicsx/cell_display_fields.npz');surface=np.load(root/'surface/surface.npz');report=json.loads((root/'fenicsx/summary.json').read_text());review=json.loads((root/'mesh_h6_netgen/mesh_review.json').read_text(encoding='utf-8'));parts={p['part_id']:p for p in review['parts'] if p['structural']};x=m['xyz'];u=a['U_mm'];cells=m['cells'];pids=m['cell_part_ids'];fiber=s['mean_stress_local_MPa'][:,0];vm=s['mean_GP_von_Mises_MPa'];materials={name:i for i,name in enumerate(['carbon','aluminum','titanium','steel','adhesive'])};sub=np.array([[0,3,5],[3,1,4],[5,4,2],[3,4,5]]);body=[]
    for pid,p in parts.items():
        take=(surface['part_ids']==pid)&surface['is_exterior'];faces=surface['faces'][take];owners=surface['owner_cells'][take];n=len(faces)
        if not n:continue
        tri=(np.arange(n)[:,None,None]*6+sub[None,:,:]).reshape(-1,3)
        body.append(dict(part_id=pid,name=p['name'],material=p['material'],kind=materials[p['material']],xyz=enc(x[faces].reshape(-1,3),'<f4'),u=enc(u[faces].reshape(-1,3),'<f4'),fiber=enc(np.repeat(fiber[owners],6),'<f4'),vm=enc(np.repeat(vm[owners],6),'<f4'),tri=enc(tri,'<u4')))
    save(out/'surface.json.gz',{'parts':body})
    # Display section of the FE domain, not the separate cosmetic C01 CAD.
    points={pid:[] for pid in parts};motion={pid:[] for pid in parts};fs={pid:[] for pid in parts};vs={pid:[] for pid in parts};tc={pid:[] for pid in parts};corner=x[cells[:,:4]]
    for axis,other,sign,normal in [(0,1,-1,np.array([1.,0,0])),(1,0,1,np.array([0.,-1,0]))]:
        values=corner[:,:,axis];candidates=np.flatnonzero((values.min(1)<=0)&(values.max(1)>=0))
        for eid in candidates:
            verts=corner[eid];poly=[]
            for i,j in EDGES:
                vi,vj=verts[i,axis],verts[j,axis]
                if vi*vj<0:
                    t=-vi/(vj-vi);lam=np.zeros(4);lam[i]=1-t;lam[j]=t;poly.append(lam)
                elif abs(vi)<1e-10:
                    lam=np.zeros(4);lam[i]=1;poly.append(lam)
                elif abs(vj)<1e-10:
                    lam=np.zeros(4);lam[j]=1;poly.append(lam)
            if len(poly)<3:continue
            poly=np.unique(np.round(poly,12),axis=0)
            if len(poly)<3:continue
            q=poly@verts;axes=[i for i in range(3) if i!=axis];center=q.mean(0);angle=np.arctan2(q[:,axes[1]]-center[axes[1]],q[:,axes[0]]-center[axes[0]]);poly=poly[np.argsort(angle)];clipped=[]
            for i,A in enumerate(poly):
                B=poly[(i+1)%len(poly)];da=sign*(A@verts)[other];db=sign*(B@verts)[other]
                if da>=-1e-10:clipped.append(A)
                if (da>0 and db<0) or (da<0 and db>0):clipped.append(A+da/(da-db)*(B-A))
            if len(clipped)<3:continue
            poly=np.array(clipped);q=poly@verts;N=np.c_[poly*(2*poly-1),np.stack([4*poly[:,i]*poly[:,j] for i,j in EDGES],axis=1)];U=N@u[cells[eid]];pid=int(pids[eid]);offset=len(points[pid]);points[pid].extend(q);motion[pid].extend(U);fs[pid].extend([float(fiber[eid])]*len(q));vs[pid].extend([float(vm[eid])]*len(q))
            for j in range(1,len(q)-1):
                idx=[0,j,j+1]
                if np.dot(np.cross(q[j]-q[0],q[j+1]-q[0]),normal)<0:idx=[0,j+1,j]
                tc[pid].append(np.array(idx)+offset)
    section=[]
    for pid,p in parts.items():
        if not points[pid]:continue
        section.append(dict(part_id=pid,name=p['name'],material=p['material'],kind=materials[p['material']],xyz=enc(points[pid],'<f4'),u=enc(motion[pid],'<f4'),fiber=enc(fs[pid],'<f4'),vm=enc(vs[pid],'<f4'),tri=enc(tc[pid],'<u4')))
    save(out/'section.json.gz',{'parts':section})
    carbon=np.isin(pids,[pid for pid,p in parts.items() if p['material']=='carbon']);metal=np.isin(pids,[pid for pid,p in parts.items() if p['material'] in ('aluminum','titanium','steel')]);manifest={'schema':'copv_c04_field_viewer_v1','parts':[dict(part_id=pid,name=p['name'],material=p['material']) for pid,p in parts.items()],'ranges':{'U_mm':float(np.linalg.norm(u,axis=1).max()),'fiber_MPa':float(abs(fiber[carbon]).max()),'metal_VM_MPa':float(vm[metal].max())},'states':json.loads((root/'fenicsx/load_unload_states.json').read_text()),'summary':report,'display_method':'Original TET10/Tri6 nodes, 4 sub-triangles per face; P2 nodal displacement; stress is volume-weighted GP mean per original cell, with no averaging across materials/plies. Quarter section is cut from the FE domain and evaluates P2 displacement at section vertices.','state_method':'Exact linear-elastic field superposition; no RGB frame blending or video smoothing','seconds':time.perf_counter()-start};manifest['display_method']=f"Straight TET10 geometry; solved P{report['displacement_degree']} displacement with true independent P2 midside unknowns only for degree 2. Tri6 faces use four display triangles; section vertices evaluate the actual solved P1/P2 field. Original-cell GP volume means, no material/ply mixing.";(out/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps({'surface_bytes':(out/'surface.json.gz').stat().st_size,'section_bytes':(out/'section.json.gz').stat().st_size,'seconds':manifest['seconds']}),flush=True)
if __name__=='__main__':main()
