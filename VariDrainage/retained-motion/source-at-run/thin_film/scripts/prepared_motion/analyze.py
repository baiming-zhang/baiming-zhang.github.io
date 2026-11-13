"""Saved-step SI/GEOS audit and independent physical-cut constitutive estimate.

Pressure recovered from curved energy, not a balance residual. A least-squares
gradient and lagged cell mobility give a separate Darcy vector estimate; exact
GEOS cut lengths integrate it and the recorded particle velocity. This is a
nonconforming recovered field, not exact Stokes flux or a new solver closure.
"""
from pathlib import Path
import importlib.util,json,numpy as np,shapely
from shapely import LineString,STRtree
import run as cfg
R,O,RAW,save,sha=cfg.R,cfg.O,cfg.RAW,cfg.save,cfg.sha
spec=importlib.util.spec_from_file_location('independent_curved',Path(__file__).resolve().parents[1]/'curved_edge/analyze.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)

def physical_cut(a,b,c,faces,m,mean):
    q=old.quantities(c,faces,m,mean);xy=np.c_[c['x_m'],c['y_m']];i,j=q['i'],q['j'];grad=np.zeros_like(xy);neighbours=[[] for _ in xy]
    for ii,jj in zip(i,j):neighbours[ii].append(jj);neighbours[jj].append(ii)
    worst_condition=0
    for ii,ids in enumerate(neighbours):
        d=xy[ids]-xy[ii];length=np.linalg.norm(d,axis=1);matrix=d/length[:,None];dp=(q['pressure'][ids]-q['pressure'][ii])/length
        grad[ii],_,rank,_=np.linalg.lstsq(matrix,dp,rcond=None)
        assert rank==2;worst_condition=max(worst_condition,float(np.linalg.cond(matrix)))
    poly=old.polys(c);cut=LineString([(240e-6,0),(240e-6,160e-6)]);ell=shapely.length(shapely.intersection(poly,cut));assert abs(ell.sum()/160e-6-1)<1e-8
    h=b['V_m3']/b['A_m2'];bath=c['bath']>0;assert not any(bath&(ell>1e-15));mobility=h**3/(12*m['viscosity_Pa_s'])
    darcy=float(np.sum(ell*(-mobility*grad[:,0])))
    materialV=float(np.sum(ell*h*b['ux_m_s']));materialS=float(np.sum(ell*b['S_mol']/b['A_m2']*b['ux_m_s']))
    native_relative=float(np.sum((c['w'][j]-c['w'][i])*faces['F_m3'])/m['dt_s'])
    native_materialV=float(np.sum((c['w']-a['w'])*a['V_m3'])/m['dt_s']);native_materialS=float(np.sum((c['w']-a['w'])*a['S_mol'])/m['dt_s'])
    return dict(recovered_relative_m3_s=darcy,ledger_relative_m3_s=native_relative,recovered_material_m3_s=materialV,ledger_material_m3_s=native_materialV,recovered_surface_mol_s=materialS,ledger_surface_mol_s=native_materialS,cut_length_m=float(ell.sum()),cut_cells=int(sum(ell>1e-15)),max_gradient_condition=worst_condition)

def compare(a,b,kind,cache):
    ia,sa,ma=cache[a];ib,sb,mb=cache[b];answer=[]
    for time in [10e-6,50e-6,200e-6]:
        na,nb=round(time/ma['dt_s']),round(time/mb['dt_s']);ta=sa[(sa['step']==na)&(sa['stage']==2)];tb=sb[(sb['step']==nb)&(sb['stage']==2)]
        if not len(ta) or not len(tb):continue
        pa,pb=old.polys(ta),old.polys(tb);i,j=STRtree(pb).query(pa,predicate='intersects');w=shapely.area(shapely.intersection(pa[i],pb[j]));ok=(w>0)&(ta['bath'][i]==0)&(tb['bath'][j]==0);i,j,w=i[ok],j[ok],w[ok]
        x={'dh':ta['h_m']-ia['h_m'],'displacement':np.c_[ta['x_m']-ia['x_m'],ta['y_m']-ia['y_m']]};y={'dh':tb['h_m']-ib['h_m'],'displacement':np.c_[tb['x_m']-ib['x_m'],tb['y_m']-ib['y_m']]}
        errors={key:old.norm(x[key][i],y[key][j],w) for key in x};answer.append(dict(candidate=a,reference=b,kind=kind,t_s=time,**errors,passed=max(errors.values())<(.02 if kind=='time' else .1)))
    return answer

def main():
    status=json.loads((O/'status.json').read_text());assert status['state']=='completed';status.update(sampled_steps=True,native_checks=18);save(O/'status.json',status)
    old.O=O;old.RAW=RAW;old.main()
    # The inherited auditor's formulae are unchanged; only saved-step iteration
    # is enabled. No numerical rerun and no fabricated omitted states.
    audit=json.loads((O/'summary.json').read_text());audit['scope']='Independent formula audit at saved steps ONLY; native histories cover all accepted steps';save(O/'formula_audit.json',audit)
    cache={};observables={};cut_results={}
    for rec in status['cases']:
        root=RAW/rec['name'];m=rec.get('summary',{});states=old.read(root/'states.csv');faces=old.read(root/'faces.csv');hist=old.read(root/'history.csv')
        if not len(states) or not m.get('steps'):continue
        initial=states[(states['step']==1)&(states['stage']==0)];mean=initial['S_mol'].sum()/(m['Gamma_inf_mol_m2']*old.domain.area)
        assert len(hist)==m['steps'];rows=[];cuts=[]
        for step in sorted(set(map(int,states['step']))):
            a,b,c=[states[(states['step']==step)&(states['stage']==stage)] for stage in range(3)];f=faces[(faces['step']==step)&(faces['stage']==2)];film=c['bath']==0;mini=c[film][np.argmin(c['h_m'][film])];bath=c['bath']>0
            history=hist[:step];rows.append(dict(step=step,t_s=step*m['dt_s'],min_h_m=float(mini['h_m']),min_x_m=float(mini['x_m']),net_into_bath_m3=float(np.sum(c['V_m3'][bath]-initial['V_m3'][bath])),relative_into_edge_m3=float(history['relativeV_m3'].sum()),material_into_edge_m3=float(history['materialV_m3'].sum()),material_surface_into_edge_mol=float(history['materialS_mol'].sum())))
            if rec['parameters']['preparation']=='variational' and step*m['dt_s']>=9.999e-6:
                cuts.append(dict(step=step,t_s=step*m['dt_s'],**physical_cut(a,b,c,f,m,mean)))
        cut_results[rec['name']]=cuts;observables[rec['name']]=dict(records=rows,native_all_step_budget_max=float(max(hist['global_budget'])),native_all_step_balance_max=float(max(max(hist[k]) for k in ['V_res','S_res','regionV_res','regionS_res'])),native_history_rows=len(hist))
        if rec['completed']:cache[rec['name']]=(initial,states,m)
    comparisons=[]
    for a,b,kind in [('l20_dt1','l20_dt05','time'),('l20_dt1','u20_dt1','space'),('l20_dt1','l10_dt1','space'),('l10_dt1','u10_dt1','space'),('u10_dt1','u10_dt05','time')]:
        if a in cache and b in cache:comparisons+=compare(a,b,kind,cache)
    summary=dict(attempts=len(status['cases']),completed=sum(c['completed'] for c in status['cases']),checks=18,saved_step_formula_audits=audit['audits'],comparisons=comparisons,skipped=status.get('skipped',[]),observables=observables,normal_flux_scope='Independent recovered pressure-gradient/cell-mobility and recorded-velocity cut integrals; nonconforming reconstruction, NOT exact normal-flux or Stokes validation. Compare components and refine.',full_C1_complete=False)
    save(O/'summary.json',summary);save(O/'physical_cut.json',cut_results);print(json.dumps(comparisons,indent=2))
if __name__=='__main__':main()
