"""Exact affine-TET10 field evaluation and stratified common physical probes.

Candidate cell searches never replace containment, native IDs or material checks.
Original GP statistics are reconstructed from saved nodal fields, not smoothed.
"""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from scipy.spatial import cKDTree
from geobc.copv_material import engineering_stiffness
from geobc.copv_structured_mesh import analytic_segments
from geobc.copv_geometry import DEFAULT_PARAMETERS
EDGES=[(0,1),(1,2),(2,0),(0,3),(2,3),(1,3)]
DL=np.array([[-1,-1,-1],[1,0,0],[0,1,0],[0,0,1.]])
FACTOR=np.array([1,1,1,2,2,2])
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def write(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def shape(L):return np.concatenate([L*(2*L-1),np.stack([4*L[...,i]*L[...,j] for i,j in EDGES],axis=-1)],axis=-1)
def derivatives(L):return np.concatenate([(4*L[...,None]-1)*DL,np.stack([4*(L[...,j,None]*DL[i]+L[...,i,None]*DL[j]) for i,j in EDGES],axis=-2)],axis=-2)
class Model:
    def __init__(self,id,meshroot,solution):
        self.id=id;self.meshroot=meshroot;self.solution=solution;self.report=read(solution/'summary.json');self.review=read(meshroot/'mesh_review.json')
        with np.load(meshroot/'mesh.npz') as m:self.x=m['xyz'];self.cells=m['cells'];self.pids=m['cell_part_ids'];self.nodeids=m['node_ids']
        with np.load(solution/'nodal_fields.npz') as a:
            assert np.array_equal(self.nodeids,a['node_ids']);self.U=a['U_mm'];self.gauge=a['gauge_coordinates_mm']
        assert digest(meshroot/'mesh.npz')==self.report['mesh_sha256']
        self.corner=self.x[self.cells[:,:4]];self.J=np.swapaxes(self.corner[:,1:]-self.corner[:,:1],1,2);self.det=np.linalg.det(self.J);assert self.det.min()>0
        for k,(a,b) in enumerate(EDGES):assert abs(self.x[self.cells[:,4+k]]-.5*(self.corner[:,a]+self.corner[:,b])).max()<1e-9
        self.inv=np.linalg.inv(self.J);self.volume=self.det/6;self.parts={p['part_id']:p for p in self.review['parts'] if p['structural']};self.trees={}
    def evaluate(self,points,part_ids):
        found=np.full(len(points),-1,np.int64);L=np.full((len(points),4),np.nan)
        for pid in np.unique(part_ids):
            rows=np.flatnonzero(part_ids==pid);ids=np.flatnonzero(self.pids==pid)
            if int(pid) not in self.trees:self.trees[int(pid)]=(ids,cKDTree(self.corner[ids].mean(1)))
            ids,tree=self.trees[int(pid)]
            for k in [24,128,768]:
                rows=rows[found[rows]<0]
                if not len(rows):break
                _,near=tree.query(points[rows],k=min(k,len(ids)));candidate=ids[np.atleast_2d(near)]
                if len(rows)==1:candidate=candidate.reshape(1,-1)
                q=np.einsum('nkij,nkj->nki',self.inv[candidate],points[rows,None]-self.corner[candidate,0]);bary=np.concatenate([1-q.sum(-1,keepdims=True),q],axis=-1)
                inside=(bary.min(-1)>=-1e-9)&(bary.max(-1)<=1+1e-9);ok=inside.any(1);which=inside.argmax(1);found[rows[ok]]=candidate[np.flatnonzero(ok),which[ok]];L[rows[ok]]=bary[np.flatnonzero(ok),which[ok]]
        valid=found>=0;out=np.full((len(points),3),np.nan);grad=np.full((len(points),3,3),np.nan)
        if valid.any():
            eid=found[valid];assert np.array_equal(self.pids[eid],part_ids[valid]);out[valid]=np.einsum('na,nai->ni',shape(L[valid]),self.U[self.cells[eid]])
            reference=np.einsum('nai,naj->nij',self.U[self.cells[eid]],derivatives(L[valid]));grad[valid]=reference@self.inv[eid]
            mapped=np.einsum('na,nai->ni',L[valid],self.corner[eid]);assert np.linalg.norm(mapped-points[valid],axis=1).max()<1e-8
        return out,grad,found
class Material:
    def __init__(self,contract,parts):
        self.parts=parts;self.matrices={}
        for name,m in contract['materials'].items():
            if 'E_MPa' in m:
                E=m['E_MPa'];nu=m['nu'];G=E/(2*(1+nu));self.matrices[name]=engineering_stiffness(E1=E,E2=E,E3=E,G12=G,G13=G,G23=G,nu12=nu,nu13=nu,nu23=nu)
            else:self.matrices[name]=engineering_stiffness(**{k:m[k] for k in ['E1','E2','E3','G12','G13','G23','nu12','nu13','nu23']})
        rz,n=analytic_segments(DEFAULT_PARAMETERS,(4,2500,16,2500,4));self.tree=cKDTree(rz);self.tangent=np.c_[-n[:,1],n[:,0]]
    def stress(self,positions,grad,pids):
        pids=np.broadcast_to(pids,positions.shape[:-1]);Q=np.broadcast_to(np.eye(3),positions.shape[:-1]+(3,3)).copy();carbon=np.isin(pids,[pid for pid,p in self.parts.items() if p['material']=='carbon'])
        if carbon.any():
            x=positions[carbon];r=np.linalg.norm(x[:,:2],axis=1);_,near=self.tree.query(np.c_[r,x[:,2]]);t=self.tangent[near];er=np.c_[x[:,:2]/r[:,None],np.zeros(len(x))];et=np.c_[-er[:,1],er[:,0],np.zeros(len(x))];m=-t[:,0,None]*er;m[:,2]-=t[:,1];normal=np.cross(m,et);angle=np.deg2rad([self.parts[int(pid)].get('angle',0) for pid in pids[carbon]])[:,None];ca=np.cos(angle);sa=np.sin(angle);Q[carbon]=np.stack([ca*m+sa*et,-sa*m+ca*et,normal],axis=-1)
        eps=.5*(grad+np.swapaxes(grad,-1,-2));e=np.swapaxes(Q,-1,-2)@eps@Q;e=np.stack([e[...,0,0],e[...,1,1],e[...,2,2],2*e[...,1,2],2*e[...,0,2],2*e[...,0,1]],axis=-1);s=np.empty_like(e)
        for pid,p in self.parts.items():
            chosen=pids==pid;s[chosen]=e[chosen]@self.matrices[p['material']].T
        return s,e
def bulk(model,material,L,qw,out):
    file=out/(model.id+'_GP_statistics.json')
    if file.exists():return read(file)
    dN=derivatives(L);N=shape(L);stats={pid:dict(name=p['name'],material=p['material'],volume_mm3=0.,sum=np.zeros(6),sum_square=np.zeros(6),energy=0.) for pid,p in model.parts.items()};mean_error=0.;mean_scale=0.;display=np.load(model.solution/'cell_display_fields.npz')['mean_stress_local_MPa']
    for first in range(0,len(model.cells),4000):
        ids=np.arange(first,min(first+4000,len(model.cells)));u=model.U[model.cells[ids]];G=np.einsum('eai,qaj->eqij',u,dN)@model.inv[ids,None];positions=np.einsum('qa,eai->eqi',L,model.corner[ids]);s,e=material.stress(positions,G,model.pids[ids,None]);w=model.det[ids,None]*qw;mean=np.sum(s*w[...,None],axis=1)/w.sum(1)[:,None];mean_error+=float(np.sum((mean-display[ids])**2));mean_scale+=float(np.sum(display[ids]**2))
        for pid in np.unique(model.pids[ids]):
            chosen=model.pids[ids]==pid;st=stats[int(pid)];ww=w[chosen];ss=s[chosen];ee=e[chosen];st['volume_mm3']+=float(ww.sum());st['sum']+=np.sum(ss*ww[...,None],axis=(0,1));st['sum_square']+=np.sum(ss*ss*ww[...,None],axis=(0,1));st['energy']+=float(.5*np.sum(ss*ee*ww[...,None]))
    energy=sum(v['energy'] for v in stats.values());error=abs(energy/model.report['strain_energy_Nmm']-1);check=float(np.sqrt(mean_error/max(mean_scale,1e-30)))
    assert error<1e-6,(model.id,'GP energy reconstruction',error)
    assert check<1e-6,(model.id,'raw GP cell mean reconstruction',check)
    for st in stats.values():
        st['mean_local_stress_MPa']=(st.pop('sum')/st['volume_mm3']).tolist();square=st.pop('sum_square')/st['volume_mm3'];st['RMS_components_MPa']=np.sqrt(square).tolist();st['RMS_tensor_MPa']=float(np.sqrt(np.dot(square,FACTOR)))
    result=dict(case=model.id,association='Actual degree-4 GP and per-material local frames; no cross-cell or material averaging',GP_energy_Nmm=energy,GP_energy_relative_reconstruction_error=error,GP_cell_mean_relative_reconstruction_error=check,parts=stats);write(file,result);return result
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workspace',type=Path,required=True);args=ap.parse_args();W=args.workspace;O=W/'data/copv/c03_mesh_v1';B=W/'data/copv/c02_structural_v1';out=O/'analysis';out.mkdir(exist_ok=True);start=time.perf_counter()
    definitions=[('h6_P1',B/'mesh_h6_netgen',B/'fenicsx'),('h4_P1',O/'mesh_h4',O/'h4_P1'),('h6_P2',B/'mesh_h6_netgen',O/'h6_P2'),('h4_P2',O/'mesh_h4',O/'h4_P2'),('h6_P1_discrete',B/'mesh_h6_netgen',O/'h6_P1_discrete')]
    base=Model(*definitions[0]);contract=read(W/'publication/GeoBC-Neural-Operator/configs/copv_c03_mesh_v1.json');mat=Material(contract,base.parts)
    pivot=int(np.argmax(base.det));qxyz=np.load(B/'fenicsx/GP_xyz_mm.npy',mmap_mode='r')[pivot];q=(qxyz-base.corner[pivot,0])@base.inv[pivot].T;L=np.c_[1-q.sum(1),q];qw=np.load(B/'fenicsx/GP_volume_weights_mm3.npy',mmap_mode='r')[pivot]/base.det[pivot];assert abs(qw.sum()-1/6)<1e-12
    rng=np.random.default_rng(20261026);pids=np.array(sorted(base.parts));counts=np.full(len(pids),2048//len(pids));counts[:2048-counts.sum()]+=1;points=[];labels=[];weights=[]
    for pid,count in zip(pids,counts):
        candidates=np.flatnonzero(base.pids==pid);v=base.volume[candidates];chosen=rng.choice(candidates,size=int(count),replace=False,p=v/v.sum());points.extend(base.corner[chosen].mean(1));labels.extend([pid]*int(count));weights.extend([v.sum()/count]*int(count))
    points=np.array(points);labels=np.array(labels);weights=np.array(weights);assert len(np.unique(points,axis=0))==2048;np.savez_compressed(out/'common_probes.npz',xyz_mm=points,part_ids=labels,part_volume_weights_mm3=weights)
    fields={};summaries={};gauges={};reports={};bulks={}
    for id,meshroot,solution in definitions:
        if not (solution/'summary.json').exists():continue
        model=base if id=='h6_P1' else Model(id,meshroot,solution);U,G,eid=model.evaluate(points,labels);valid=eid>=0;stress=np.full((len(points),6),np.nan);stress[valid],_=mat.stress(points[valid],G[valid],labels[valid]);fields[id]=dict(U=U,stress=stress,contained=valid);gauges[id]=model.gauge.tolist();reports[id]=model.report
        summaries[id]={k:model.report[k] for k in ['nodes','cells','dofs','displacement_degree','max_U_mm','strain_energy_Nmm','true_relative_residual','peak_solve_seconds','total_seconds']};summaries[id].update(mapped_probes=int(valid.sum()),probes_total=len(points),cap_convention=model.report.get('cap_convention','circular'))
        if id!='h6_P1_discrete':bulks[id]=bulk(model,mat,L,qw,out)
        print(json.dumps(dict(case=id,summary=summaries[id])),flush=True)
        if model is not base:del model
    comparisons={}
    for a,b,key in [('h6_P1','h4_P1','P1_mesh'),('h6_P1','h6_P2','degree_on_h6'),('h6_P2','h4_P2','P2_mesh'),('h4_P1','h4_P2','degree_on_h4'),('h6_P1','h6_P1_discrete','cap_closure')]:
        if a not in fields or b not in fields:continue
        aa,bb=fields[a],fields[b];valid=aa['contained']&bb['contained'];du=bb['U'][valid]-aa['U'][valid];ds=bb['stress'][valid]-aa['stress'][valid];den=aa['U'][valid];sigma=aa['stress'][valid];w=weights[valid]
        comparisons[key]=dict(reference=a,candidate=b,common_probes=int(valid.sum()),unmapped_explicitly_excluded=int((~valid).sum()),relative_probe_vector_L2=float(np.linalg.norm(du)/np.linalg.norm(den)),probe_vector_RMSE_mm=float(np.sqrt(np.mean(np.sum(du*du,axis=1)))),p95_vector_difference_mm=float(np.quantile(np.linalg.norm(du,axis=1),.95)),max_vector_difference_mm=float(np.linalg.norm(du,axis=1).max()),volume_weighted_probe_stress_Frobenius_difference=float(np.sqrt(np.sum(ds*ds*FACTOR*w[:,None])/np.sum(sigma*sigma*FACTOR*w[:,None]))),energy_relative_change=reports[b]['strain_energy_Nmm']/reports[a]['strain_energy_Nmm']-1,node_max_U_relative_change=reports[b]['max_U_mm']/reports[a]['max_U_mm']-1,gauge_coordinate_max_difference_mm=float(abs(np.array(gauges[a])-np.array(gauges[b])).max()),no_rigid_registration=True,by_part={})
        for pid in pids:
            mask=valid&(labels==pid);delta=bb['U'][mask]-aa['U'][mask];comparisons[key]['by_part'][str(pid)]=dict(mapped=int(mask.sum()),relative_U_L2=float(np.linalg.norm(delta)/np.linalg.norm(aa['U'][mask])) if mask.any() else None)
        if key.startswith('degree_'):
            with np.load((B/'fenicsx' if a=='h6_P1' else O/a)/'nodal_fields.npz') as ua, np.load(O/b/'nodal_fields.npz') as ub:
                assert np.array_equal(ua['node_ids'],ub['node_ids']);delta=ub['U_mm']-ua['U_mm'];comparisons[key]['same_mesh_all_stored_nodes_relative_U_L2']=float(np.linalg.norm(delta)/np.linalg.norm(ua['U_mm']))
    columns={'x_mm':points[:,0],'y_mm':points[:,1],'z_mm':points[:,2],'part_id':labels}
    for id,f in fields.items():
        columns[id+'_contained']=f['contained'].astype(int)
        for i,c in enumerate('xyz'):columns[id+'_U'+c+'_mm']=f['U'][:,i]
    np.savetxt(out/'common_2048_probe_fields.csv',np.column_stack(list(columns.values())),delimiter=',',header=','.join(columns),comments='',fmt='%.12g')
    result=dict(schema='copv_c03_comparison_v1',probe_count=len(points),method='2048 fixed material-stratified interior cell-centroid probes, volume-proportional cell sampling; affine tetra containment and P1/P2 shape functions; no nearest-node matching, extrapolation or rigid registration',cases=summaries,comparisons=comparisons,GP_statistics={k:k+'_GP_statistics.json' for k in bulks},limitations=['Probe norms are sampled differences, not exact global discretization errors','GP RMS is material-wise, from the actual solved degree and original quadrature','Sharp-corner raw peaks are diagnostic and not convergence certificates'],seconds=time.perf_counter()-start);write(out/'comparison.json',result);print(json.dumps(result))
if __name__=='__main__':main()
