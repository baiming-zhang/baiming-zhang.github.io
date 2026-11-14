"""Independent SI energy/pressure, GEOS geometry and regional drainage power."""
from pathlib import Path
import json,numpy as np,shapely
from shapely import MultiPoint,box,STRtree
from run import R,O,RAW,save,sha
domain=box(0,0,320e-6,160e-6);right=box(240e-6,0,320e-6,160e-6)
def read(p):return np.atleast_1d(np.genfromtxt(p,delimiter=',',names=True))
def polys(t):return shapely.intersection(np.array(shapely.get_parts(shapely.voronoi_polygons(MultiPoint(np.c_[t['x_m'],t['y_m']]),extend_to=domain,ordered=True)),dtype=object),domain)
def quantities(t,faces,m,mean):
    i=faces['i'].astype(int);j=faces['j'].astype(int);xy=np.c_[t['x_m'],t['y_m']];d=np.linalg.norm(xy[i]-xy[j],axis=1);length=faces['c']*d
    z=(t['h_m'][i]-t['h_m'][j])/(2*d);f=np.hypot(1,z);extra=length*d*z*z/(f+1);asurf=t['A_m2'].copy();np.add.at(asurf,i,.5*extra);np.add.at(asurf,j,.5*extra)
    theta=t['S_mol']/(m['Gamma_inf_mol_m2']*asurf);delta=theta-mean;RTG=m['RT_J_mol']*m['Gamma_inf_mol_m2'];ref=m['clean_tension_N_m']+RTG*np.log1p(-mean)
    surf=2*(asurf*RTG*(theta*np.log1p(delta/mean)+(1-theta)*np.log1p(-delta/(1-mean)))+ref*(asurf-t['A_m2']))
    gamma=m['clean_tension_N_m']+RTG*np.log1p(-theta);first=length*z/(2*f);pressure=np.zeros(len(t));effort=(gamma[i]+gamma[j])*first
    np.add.at(pressure,i,effort/t['A_m2'][i]);np.add.at(pressure,j,-effort/t['A_m2'][j]);bath=t['bath']>0;pressure[bath]=m['bath_pressure_Pa'];energy=surf+bath*m['bath_pressure_Pa']*t['V_m3']
    return dict(i=i,j=j,d=d,length=length,z=z,first=first,area=asurf,theta=theta,gamma=gamma,surf=surf,energy=energy,pressure=pressure)
def norm(x,y,w):
    delta=(x-y)**2;base=y*y
    if x.ndim==2:delta=delta.sum(axis=1);base=base.sum(axis=1)
    return float(np.sqrt(np.sum(w*delta)/max(1e-60,np.sum(w*base))))
def main():
    status=json.loads((O/'status.json').read_text());audits={};cache={};regional={};partial={}
    for record in status['cases']:
        name=record['name'];root=RAW/name
        for f in record['files']:assert sha(R/f['path'])==f['sha256']
        if not (root/'states.csv').exists() or (root/'states.csv').stat().st_size<200:continue
        states=read(root/'states.csv');faces=read(root/'faces.csv');history=read(root/'history.csv');initial=states[(states['step']==1)&(states['stage']==0)]
        if not len(initial):continue
        if 'summary' in record:m=record['summary']
        else:
            reference=next(c for c in status['cases'] if c['completed'] and c['parameters']['radius_um']==record['parameters']['radius_um']);m=dict(reference['summary']);m.update(N=len(initial),steps=int(max(states['step'])),dt_s=record['parameters']['dt_us']*1e-6,initial_energy_J=float(initial['node_energy_J'].sum()))
            partial[name]=dict(material_and_boundary_source=reference['name'],accepted_steps=m['steps'],accepted_time_s=m['steps']*m['dt_s'],N=len(initial),initial_energy_from_saved_native_nodes_J=m['initial_energy_J'],completed_horizon=False,reason='native exception, no native summary; inferred prefix metadata only')
        n=m['N']
        mean=initial['S_mol'].sum()/(m['Gamma_inf_mol_m2']*domain.area);V0=initial['V_m3'].sum();S0=initial['S_mol'].sum()
        measurements={k:[] for k in ['geometry','curved_area','pressure','node_energy','V','S','node_flux','regionV','regionS','flux_law','global_energy','half_drain','global_budget','positive_theta']};budgets=[];energy_scale=1e-26;pressure_scale=1e-26
        audited_steps=sorted(set(map(int,states['step']))) if status.get('sampled_steps') else list(range(1,m['steps']+1))
        for step in audited_steps:
            ts=[states[(states['step']==step)&(states['stage']==z)] for z in range(3)];fs=[faces[(faces['step']==step)&(faces['stage']==z)] for z in range(3)];qs=[quantities(t,f,m,mean) for t,f in zip(ts,fs)];a,b,c=ts;qa,qb,qc=qs
            for stage in range(3):
                t,q=ts[stage],qs[stage];pressure_scale=max(pressure_scale,float(max(abs(t['pressure_Pa']))));energy_scale=max(energy_scale,float(sum(abs(t['node_energy_J']))))
                physical=t['bath']==0 if m.get('signed_bath_counter') else np.ones(len(t),dtype=bool)
                measurements['curved_area'].append(float(max(abs(q['area']/t['As_m2']-1))));measurements['pressure'].append(float(max(abs(q['pressure']-t['pressure_Pa']))));measurements['node_energy'].append(float(max(abs(q['energy']-t['node_energy_J']))));measurements['V'].append(abs(t['V_m3'].sum()-V0)/V0);measurements['S'].append(abs(t['S_mol'].sum()-S0)/S0);measurements['positive_theta'].append(bool(min(q['theta'])>0 and max(q['theta'])<1 and min(t['V_m3'][physical])>0 and np.all(np.isfinite(t['V_m3']))))
                if stage<2:
                    poly=polys(t);areas=shapely.area(poly);w=shapely.area(shapely.intersection(poly,right))/areas
                    measurements['geometry'] += [float(max(abs(areas/t['A_m2']-1))),float(max(abs(w-t['w'])))]
            i,j=qc['i'],qc['j'];F=fs[2]['F_m3'];dv=c['V_m3']-b['V_m3'];pred=np.zeros(n);np.add.at(pred,i,-F);np.add.at(pred,j,F);measurements['node_flux'].append(float(max(abs(dv-pred))/V0))
            dh=c['h_m']-b['h_m'];linarea=np.zeros(n);np.add.at(linarea,i,.5*qc['first']*(dh[i]-dh[j]));np.add.at(linarea,j,.5*qc['first']*(dh[i]-dh[j]))
            # Two explicit convex losses: monolayer free-energy in area and
            # convex geometric sqrt area. Neither is a budget residual repair.
            darea=qc['area']-qb['area'];dsurf=qc['surf']-qb['surf'];dilution_loss=2*qc['gamma']*darea-dsurf;geometry_loss=2*qc['gamma']*(linarea-darea);loss=dilution_loss+geometry_loss
            h=b['V_m3']/b['A_m2'];hf=.5*(h[i]+h[j]);bi=b['bath'][i]>0;bj=b['bath'][j]>0;hf=np.where(bi&~bj,h[j],np.where(bj&~bi,h[i],hf));k=hf**3/(12*m['viscosity_Pa_s'])*fs[2]['c'];dp=qc['pressure'][i]-qc['pressure'][j];dt=m['dt_s']
            if m.get('signed_bath_counter'):k[bi&bj]=0
            measurements['flux_law'].append(float(max(abs(F-dt*k*dp))/max(1e-40,max(abs(F)))));active=k>0;half=float(np.sum(F[active]*F[active]/(2*dt*k[active])))
            measurements['half_drain'].append(abs(half-(history['halfDrain_J'][step-1]-(history['halfDrain_J'][step-2] if step>1 else 0))))
            for region in [0,1]:
                w=c['w'] if region else 1-c['w'];oldw=a['w'] if region else 1-a['w'];wi,wj=w[i],w[j]
                shape=float(np.sum((wj-wi)*qc['first']*(qc['gamma'][j]*dh[i]+qc['gamma'][i]*dh[j])))
                diss=float(np.sum(.5*(wi+wj)*F*dp));port=float(np.sum((wj-wi)*.5*(qc['pressure'][i]+qc['pressure'][j])*F));de=float(np.sum(w*(qc['energy']-qb['energy'])));numloss=float(np.sum(w*loss));balance=de+numloss+diss-port-shape
                scale=abs(de)+abs(numloss)+abs(diss)+abs(port)+abs(shape)
                budgets.append(dict(step=step,t_s=step*dt,region=region,dE_J=de,convex_loss_J=numloss,dilution_loss_J=float(np.sum(w*dilution_loss)),geometry_loss_J=float(np.sum(w*geometry_loss)),dissipation_J=diss,pressure_work_J=port,shape_work_J=shape,residual_J=balance,pressure_only_J=balance+shape,scale_J=scale))
                materialV=np.sum((w-oldw)*a['V_m3']);materialS=np.sum((w-oldw)*a['S_mol']);relative=np.sum((wj-wi)*F)
                measurements['regionV'].append(float(abs(np.sum(w*c['V_m3']-oldw*a['V_m3'])-materialV-relative)/V0));measurements['regionS'].append(float(abs(np.sum(w*c['S_mol']-oldw*a['S_mol'])-materialS)/S0))
            totalE=float(qc['energy'].sum());kinetic=.5*m['rho_kg_m3']*np.sum(c['V_m3']*(c['ux_m_s']**2+c['uy_m_s']**2));row=history[step-1];measurements['global_energy'].append(abs(totalE-row['E_J']))
            e_initial=m['initial_energy_J'];budget=(totalE-e_initial+kinetic+row['visc_J']+row['halfDrain_J']+row['mix_J'])/max(1e-26,abs(e_initial)+abs(row['boundaryWork_J']));measurements['global_budget'].append(float(budget))
        scale=max(1e-26,sum(b['scale_J'] for b in budgets));power=max(abs(b['residual_J']) for b in budgets)/scale;omission=max(abs(b['pressure_only_J']) for b in budgets)/scale
        pair=max(abs(budgets[k]['shape_work_J']+budgets[k+1]['shape_work_J'])+abs(budgets[k]['pressure_work_J']+budgets[k+1]['pressure_work_J']) for k in range(0,len(budgets),2))/scale
        result={k:max(v) for k,v in measurements.items() if k!='positive_theta'};result['positive_theta']=all(measurements['positive_theta']);result['pressure']/=pressure_scale
        for key in ['node_energy','global_energy','half_drain']:result[key]/=energy_scale
        result.update(regional_power=power,pressure_only=omission,pair_exchange=pair,minimum_convex_loss=min(b['convex_loss_J'] for b in budgets)/scale,maximum_shape_work=max(abs(b['shape_work_J']) for b in budgets),complete_history_work_scale_J=scale)
        result['prefix_formulas_passed']=result['positive_theta'] and result['geometry']<=1e-8 and max(result[k] for k in ['V','S','node_flux','regionV','regionS'])<=1e-10 and max(result[k] for k in ['curved_area','pressure','node_energy','global_energy','half_drain','regional_power','pair_exchange','global_budget'])<=1e-7 and result['minimum_convex_loss']>=-1e-7
        result['passed']=record['native_gate'] and result['prefix_formulas_passed'];result['completed_horizon']=record['completed'];result['audited_saved_steps']=audited_steps;result['all_accepted_steps_audited']=len(audited_steps)==m['steps']
        audits[name]=result;regional[name]=budgets;last=states[(states['step']==m['steps'])&(states['stage']==2)]
        if record['completed']:cache[name]=(initial,last)
        print(name,'PASS' if result['passed'] else 'FAIL','regional',power,'omitted',omission,flush=True)
    comparisons=[]
    pairs=[('long_l40_t05','long_l40_t025','time'),('long_l40_t025','long_l40_t0125','time'),('long_l40_t025','long_u40','space'),('long_l20','long_u20','space')] if status.get('longer') else [('l40_t125','l40_t0625','time'),('l40_t0625','l40_t03125','time'),('l40_t0625','u40','space'),('l20','u20','space')]
    for a,b,kind in pairs:
        if a not in cache or b not in cache:continue
        sa,ta=cache[a];sb,tb=cache[b];va={'dh':ta['h_m']-sa['h_m'],'displacement':np.c_[ta['x_m']-sa['x_m'],ta['y_m']-sa['y_m']]};vb={'dh':tb['h_m']-sb['h_m'],'displacement':np.c_[tb['x_m']-sb['x_m'],tb['y_m']-sb['y_m']]}
        i=np.arange(len(ta));j=i;w=tb['A_m2'];overlap=None
        if kind=='space':
            pa,pb=polys(ta),polys(tb);i,j=STRtree(pb).query(pa,predicate='intersects');overlap=shapely.intersection(pa[i],pb[j]);w=shapely.area(overlap);mask=w>0;i,j,w,overlap=i[mask],j[mask],w[mask],overlap[mask]
        d=dict(candidate=a,reference=b,kind=kind,**{key:norm(va[key][i],vb[key][j],w) for key in va});d['passed']=max(d['dh'],d['displacement'])<(.02 if kind=='time' else .1)
        if overlap is not None:
            wr=shapely.area(shapely.intersection(overlap,right));d['exploratory_edge_window']={key:norm(va[key][i],vb[key][j],wr) for key in va}
        comparisons.append(d)
    save(O/'partial_attempts.json',partial);save(O/'regional_power.json',regional);save(O/'summary.json',dict(attempts=len(status['cases']),completed=sum(c['completed'] for c in status['cases']),audits=audits,comparisons=comparisons,native_checks=status.get('native_checks',34 if status.get('longer') else 11),shapely=shapely.__version__))
    save(O/'raw_manifest.json',[dict(path='data/'+p.relative_to(RAW).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in RAW.rglob('*') if p.is_file()]);print(json.dumps(comparisons,indent=2))
if __name__=='__main__':
    import sys
    if '--longer' in sys.argv:
        O=O/'longer';RAW=R/json.loads((O/'protocol.json').read_text())['raw_directory']
    main()
