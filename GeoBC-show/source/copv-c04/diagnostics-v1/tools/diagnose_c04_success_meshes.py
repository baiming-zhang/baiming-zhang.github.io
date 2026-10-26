"""Localize existing native C04 mesh defects without modifying FE results."""
from pathlib import Path
import argparse,csv,json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

W=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    review=read(W/'data/copv/c04_research_v1/combined_review.json');results=[];fig,axes=plt.subplots(2,2,figsize=(13,10),constrained_layout=True)
    for col,suffix in enumerate(['0002','0018']):
        tick=time.perf_counter();case=next(r for r in review['cases'] if r['parent_case_id']=='copv_c04_'+suffix);root=W/case['artifact_root'];m=np.load(root/'mesh/mesh.npz');qfile=np.load(root/'quality/cell_quality.npz');assert np.array_equal(m['cell_ids'],qfile['cell_ids'])
        xyz=m['xyz'];c=m['cells'];pids=m['cell_part_ids'];q=qfile['sicn'];cellids=m['cell_ids'];records=[];worst=[];plot=[];bad_all=[]
        for part in case['quality']['parts']:
            pid=part['part_id'];ix=np.flatnonzero(pids==pid);v=xyz[c[ix,:4]];center=v.mean(1);det=np.linalg.det(np.swapaxes(v[:,1:]-v[:,:1],1,2));volume=det/6;bad=q[ix]<.05
            area=np.stack([.5*np.linalg.norm(np.cross(v[:,b]-v[:,a],v[:,d]-v[:,a]),axis=1) for a,b,d in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]],1)
            altitude=det/2/area.max(1);edge=np.stack([np.linalg.norm(v[:,b]-v[:,a],axis=1) for a,b in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]],1);ratio=edge.max(1)/altitude
            zabs=abs(center[:,2]);r=np.linalg.norm(center[:,:2],axis=1);base=case['parameters_mm']['straight_length_mm']/2;depth=case['parameters_mm']['dome_depth_mm'];zeta=(zabs-base)/depth
            def stats(values):return dict(min=float(values.min()),p01=float(np.quantile(values,.01)),median=float(np.median(values)),p99=float(np.quantile(values,.99)),max=float(values.max())) if len(values) else None
            rows=dict(part_id=pid,name=part['name'],cells=len(ix),low_SICN_cells=int(bad.sum()),low_SICN_cell_fraction=float(bad.mean()),low_SICN_volume_fraction=float(volume[bad].sum()/volume.sum()),low_SICN_volume_mm3=float(volume[bad].sum()),signed_volume_error=(part['mesh_volume_mm3']/part['CAD_volume_mm3']-1),SICN=stats(q[ix]),all_min_altitude_mm=stats(altitude),bad_min_altitude_mm=stats(altitude[bad]),bad_longest_edge_over_min_altitude=stats(ratio[bad]),bad_abs_z_mm=stats(zabs[bad]),bad_radius_mm=stats(r[bad]),bad_zeta_dome=stats(zeta[bad]))
            records.append(rows)
            for k in np.argsort(q[ix])[:15]:worst.append(dict(part_id=pid,cell_id=int(cellids[ix[k]]),sicn=float(q[ix[k]]),x_mm=float(center[k,0]),y_mm=float(center[k,1]),z_mm=float(center[k,2]),r_mm=float(r[k]),min_altitude_mm=float(altitude[k]),max_edge_mm=float(edge[k].max()),max_edge_over_min_altitude=float(ratio[k]),volume_mm3=float(volume[k])))
            keep=np.linspace(0,len(ix)-1,min(1800,len(ix)),dtype=int);axes[0,col].scatter(center[keep,2],r[keep],s=.6,c='#8fa7ad',alpha=.24,rasterized=True)
            axes[0,col].scatter(center[bad,2],r[bad],s=1.2,c='#c33652',alpha=.5,rasterized=True)
            bad_all.extend(ix[bad]);plot.append((pid,part['name'],bad.sum(),volume[bad].sum()/volume.sum()))
        selected=np.asarray(bad_all,dtype=int);vbad=xyz[c[selected,:4]];np.savez_compressed(out/('copv_c04_'+suffix+'_low_SICN_native.npz'),cell_ids=cellids[selected],cell_part_ids=pids[selected],sicn=q[selected],corner_node_ids=m['node_ids'][c[selected,:4]],corner_xyz_mm=vbad)
        fields=list(worst[0]);f=out/('copv_c04_'+suffix+'_worst_native_cells.csv')
        with f.open('w',encoding='utf-8-sig',newline='') as stream:writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(worst)
        axes[0,col].set(title='Case '+suffix+': native cell centroids',xlabel='z (mm)',ylabel='radius (mm)');axes[0,col].text(.02,.02,'Red: SICN < 0.05; grey: sampled context',transform=axes[0,col].transAxes,fontsize=9)
        ids=[p['part_id'] for p in records];axes[1,col].bar(np.arange(len(ids)),[100*p['low_SICN_cell_fraction'] for p in records],color='#c33652',label='cell count %');axes[1,col].plot(np.arange(len(ids)),[100*p['low_SICN_volume_fraction'] for p in records],'o-',color='#176f7e',label='part volume %');axes[1,col].set(xticks=np.arange(len(ids)),xticklabels=ids,xlabel='Structural part ID',ylabel='Fraction below SICN 0.05 (%)');axes[1,col].legend()
        result=dict(case_id=case['parent_case_id'],artifact_root=case['artifact_root'],mesh_sha256=case['mesh_sha256'],scope='Existing straight-sided geometry; P2 solution unchanged. Diagnostic only, not a replacement acceptance protocol.',parts=records,worst_native_cells=worst,seconds=time.perf_counter()-tick)
        results.append(result);(out/(case['parent_case_id']+'.json')).write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(dict(case_id=result['case_id'],seconds=result['seconds'],low_SICN_cells=len(selected))),flush=True)
    fig.savefig(out/'defect_localization.png',dpi=180);(out/'summary.json').write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
