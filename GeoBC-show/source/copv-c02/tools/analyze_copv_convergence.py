"""Compare solved vessel meshes in their structured parameters, never invF/probes.

No solve is launched. Displacements are tensor-quadratically interpolated from
all original coarse P2 nodes to fine nodes in periodic phi, segment-local
meridian parameter, and actual normal offset. Same-gauge raw differences are
the primary metrics; no rigid-motion fit or stress peak criterion is applied.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def meridian_positions(source_manifest,target_manifest,target_j):
    coarse=np.asarray(source_manifest['meridian_segment_elements'],int)
    fine=np.asarray(target_manifest['meridian_segment_elements'],int)
    if len(coarse)!=5 or len(fine)!=5:raise ValueError('Expected five meridian segments')
    bc=np.r_[0,2*np.cumsum(coarse)];bf=np.r_[0,2*np.cumsum(fine)]
    segment=np.minimum(np.searchsorted(bf,target_j,side='right')-1,4)
    local=(target_j-bf[segment])/(2*fine[segment])
    return bc[segment]+local*2*coarse[segment]


def quadratic_weights(u):
    return np.stack([(1-u)*(1-2*u),4*u*(1-u),u*(2*u-1)],axis=-1)


def interpolate(source,sm,target,tm,field):
    """Original source array [node,...] -> original target node order."""
    if sm['parameters_mm']!=tm['parameters_mm']:raise ValueError('Geometry parameters differ')
    ij=np.asarray(source['structured_node_ijk'],int);target_ij=np.asarray(target['structured_node_ijk'],int)
    shape=tuple(ij.max(0)+1);grid=np.full(shape,-1,np.int64)
    grid[tuple(ij.T)]=np.arange(len(ij))
    if np.any(grid<0) or np.prod(shape)!=len(ij):raise ValueError('Incomplete structured source grid')
    positions=[target_ij[:,0]*sm['nphi']/tm['nphi'],
               meridian_positions(sm,tm,target_ij[:,1])]
    base=[];weights=[]
    for axis,pos in enumerate(positions):
        elements=(shape[axis]-1)//2 if axis==1 else shape[0]//2
        cell=np.minimum(np.floor(pos/2).astype(int),elements-1)
        base.append(2*cell);weights.append(quadratic_weights((pos-2*cell)/2))
    source_offsets=np.asarray(source['thickness_offsets_mm']);target_offsets=np.asarray(target['thickness_offsets_mm'])[target_ij[:,2]]
    bounds=source_offsets[::2]
    if target_offsets.min()<bounds[0]-1e-12 or target_offsets.max()>bounds[-1]+1e-12:raise ValueError('Target offset outside source wall')
    rcell=np.clip(np.searchsorted(bounds,target_offsets,side='right')-1,0,len(bounds)-2)
    base.append(2*rcell);weights.append(quadratic_weights((target_offsets-bounds[rcell])/(bounds[rcell+1]-bounds[rcell])))
    field=np.asarray(field);out=np.zeros((len(target_ij),)+field.shape[1:])
    expand=(slice(None),)+(None,)*(field.ndim-1)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                ids=grid[(base[0]+i)%shape[0],base[1]+j,base[2]+k]
                weight=weights[0][:,i]*weights[1][:,j]*weights[2][:,k]
                out+=field[ids]*weight[expand]
    return out


def weighted_quantile(x,w,quantile):
    order=np.argsort(x);x=np.asarray(x)[order];w=np.asarray(w)[order]
    return float(np.interp(quantile*np.sum(w),np.cumsum(w),x))


def layer_statistics(data):
    sigma=data['stress_local'];epsilon=data['strain_local'];w=data['quadrature_weight'];ply=data['gp_ply_ids']
    result={}
    for pid in np.unique(ply):
        chosen=ply==pid;weights=w[chosen];s=sigma[chosen];e=epsilon[chosen];volume=float(weights.sum())
        norm=np.sqrt(np.sum(s[:,:3]**2,axis=1)+2*np.sum(s[:,3:]**2,axis=1))
        result[str(int(pid))]=dict(volume_mm3=volume,
            stress_local_weighted_mean_MPa=(np.sum(s*weights[:,None],axis=0)/volume).tolist(),
            stress_local_weighted_RMS_MPa=np.sqrt(np.sum(s*s*weights[:,None],axis=0)/volume).tolist(),
            stress_tensor_weighted_RMS_MPa=float(np.sqrt(np.sum(norm*norm*weights)/volume)),
            stress_tensor_weighted_p50_MPa=weighted_quantile(norm,weights,.5),
            stress_tensor_weighted_p95_MPa=weighted_quantile(norm,weights,.95),
            strain_energy_Nmm=float(.5*np.sum(np.sum(s*e,axis=1)*weights)),
            association='Each physical GP ply retained separately; no material/interface averaging')
    return result


def scalar_metrics(row):
    p=row['pressure_MPa'];work=row['pressure_work_Nmm']+row['cap_work_Nmm']
    return dict(strain_energy_Nmm=row['strain_energy_Nmm'],external_work_Nmm=work,
        pressure_compliance_Nmm_per_MPa2=work/(p*p),end_extension_mm=row['solver']['end_extension_mm'],
        central_cylinder_mean_radial_displacement_mm=row['central_cylinder_mean_radial_displacement_mm'],
        midsection_axial_force_N=row['mid_cylinder_stress_integral']['axial_force_from_one_sided_stress_N'],
        midsection_axial_force_relative_error=row['mid_cylinder_stress_integral']['relative_error'],
        reduced_relative_residual=row['solver']['reduced_relative_residual'])


def relative_change(coarse,fine):
    return float(abs(coarse-fine)/max(abs(fine),1e-30))


def geometry_diagnostics(coarse,fine,output,cm,fm):
    """Analytic meridian continuity plus sampled P2 geometry error, no FE."""
    from geobc.copv_structured_mesh import analytic_segments
    p=cm['parameters_mm']
    if p!=fm['parameters_mm']:raise ValueError('Different analytic geometries')
    R=p['inner_radius_mm'];b=p['boss_bore_radius_mm'];D=p['dome_depth_mm'];L=p['straight_length_mm']
    control=np.array([[R,L/2],[R,L/2+.4*D],[b,L/2+.6*D],[b,L/2+D]])
    t=np.linspace(0,1,10001)[:,None]
    d=3*((1-t)**2*(control[1]-control[0])+2*t*(1-t)*(control[2]-control[1])+t*t*(control[3]-control[2]))
    dd=6*((1-t)*(control[2]-2*control[1]+control[0])+t*(control[3]-2*control[2]+control[1]))
    speed=np.linalg.norm(d,axis=1);curvature=(d[:,0]*dd[:,1]-d[:,1]*dd[:,0])/speed**3
    third=6*(control[3]-3*control[2]+3*control[1]-control[0])
    cross=lambda a,b:a[:,0]*b[...,1]-a[:,1]*b[...,0]
    curvature_t=cross(d,third)/speed**3-3*cross(d,dd)*np.sum(d*dd,axis=1)/speed**5
    curvature_s=curvature_t/speed
    total=p['liner_thickness_mm']+p['overwrap_thickness_mm']
    parallel_factor=1+total*curvature
    # Dense analytic meridian sampling only at phi=0: error due to meridian
    # approximation is isolated from circumference polygon approximation.
    dense_counts=(16,192,64,192,16)
    rz,normal=analytic_segments(p,dense_counts);jj=np.arange(len(rz));labels=['left_neck','left_dome','cylinder','right_dome','right_neck']
    dense_manifest=dict(fm,meridian_segment_elements=list(dense_counts))
    boundaries=np.r_[0,2*np.cumsum(dense_counts)]
    segments=np.minimum(np.searchsorted(boundaries,jj,side='right')-1,4)
    samples={};traces={}
    for name,root,m in [('coarse',coarse,cm),('fine',fine,fm)]:
        with np.load(root/'mesh.npz') as a:
            maxk=len(a['thickness_offsets_mm'])-1
            positions=np.vstack([np.c_[np.zeros(len(jj),int),jj,np.zeros(len(jj),int)],
                                 np.c_[np.zeros(len(jj),int),jj,np.full(len(jj),maxk)]])
            target={'structured_node_ijk':positions,'thickness_offsets_mm':a['thickness_offsets_mm']}
            predicted=interpolate(a,m,target,dense_manifest,a['xyz']).reshape(2,len(jj),3)
            truth=np.stack([np.c_[rz[:,0],np.zeros(len(jj)),rz[:,1]],
                            np.c_[rz[:,0]+total*normal[:,0],np.zeros(len(jj)),rz[:,1]+total*normal[:,1]]])
            errors=np.linalg.norm(predicted-truth,axis=2)
            samples[name]={wall:{'RMS_mm':float(np.sqrt(np.mean(errors[i]**2))),'max_mm':float(errors[i].max()),
                'max_error_segment':labels[segments[np.argmax(errors[i])]],
                'max_error_segment_local_parameter':float((jj[np.argmax(errors[i])]-boundaries[segments[np.argmax(errors[i])]])/(2*dense_counts[segments[np.argmax(errors[i])]])),
                'by_segment':{label:{'RMS_mm':float(np.sqrt(np.mean(errors[i,segments==s]**2))),'max_mm':float(errors[i,segments==s].max())} for s,label in enumerate(labels)}} for i,wall in enumerate(['inner','outer'])}
            traces[name]=errors
    records={}
    for name,root,m in [('coarse',coarse,cm),('fine',fine,fm)]:
        with np.load(root/'mesh.npz') as a:
            if 'min_det_jac' in a: minimum=float(a['min_det_jac'].min())
            else:minimum=m.get('min_det_jac_mm3')
        records[name]=dict(min_gmsh_det_jac_mm3=minimum,min_sicn=m.get('min_sicn'),
            all_cells_positive_det_jac=bool(minimum is not None and minimum>0))
    diagnostic=dict(schema_version=1,geometry_only=True,parameters_mm=p,
        continuity=dict(tangent_continuous=True,curvature_continuous=False,
            cylinder_and_neck_curvature_per_mm=0.,right_dome_start_curvature_per_mm=float(curvature[0]),
            right_dome_end_curvature_per_mm=float(curvature[-1]),
            analytic_dome_endpoint_speed_mm=float(speed[0]),
            normal_derivative_jump_magnitude_per_mm=float(abs(curvature[0])),
            normal_second_arclength_derivative_dome_endpoint_per_mm2=np.sqrt(curvature_s[[0,-1]]**2+curvature[[0,-1]]**4).tolist(),
            curvature_arclength_derivative_dome_endpoints_per_mm2=curvature_s[[0,-1]].tolist(),
            straight_segment_normal_first_and_second_derivatives=[0.,0.],
            note='Piecewise analytic meridian is tangent continuous (C1 under arc-length parametrization), not C2. Each join is a mesh segment boundary; no P2 cell crosses a join.'),
        offset_risk=dict(total_normal_thickness_mm=total,min_sampled_parallel_curve_tangent_scale=float(parallel_factor.min()),
            sampled_max_abs_curvature_per_mm=float(abs(curvature).max()),samples=10001,
            actual_mesh_checks=records,
            note='Parallel-curve tangential factor 1+offset*kappa with n=(dz,-dr)/speed. Positive dense-sample factor and actual Gmsh minDetJac; not a proof for arbitrary parameter ranges.'),
        sampled_P2_meridian_geometry_error_mm=samples,
        error_interpretation='Normals stay continuous at cylinder/dome and dome/neck, but their first/second arc-length derivatives jump. Curvature and normal-offset high derivatives concentrate interpolation error in dome transition cells (outer-wall maximum near the cylinder join for this baseline). Meridian refinement can improve pressure normals, GP material frames and bending stresses even when global equilibrium was already satisfied.',
        convergence_claimed=False)
    output.mkdir(parents=True,exist_ok=True)
    (output/'geometry_diagnostics.json').write_text(json.dumps(diagnostic,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    for name,color in [('coarse','#527da4'),('fine','#cd8057')]:
        axs[0].plot(rz[:,1],1000*traces[name][0],color=color,label=name+' inner')
        axs[0].plot(rz[:,1],1000*traces[name][1],color=color,ls='--',label=name+' outer')
    axs[0].set(xlabel='Analytic inner meridian z (mm)',ylabel='P2 geometry error (micrometers)',title='Dense segment-local samples at phi=0')
    axs[0].legend(fontsize=8,frameon=False);axs[0].grid(alpha=.15)
    axs[1].plot(t[:,0],curvature,label='Analytic signed curvature (1/mm)')
    axs[1].plot(t[:,0],parallel_factor,label='Outer offset tangent scale')
    axs[1].axhline(0,color='#888',ls=':');axs[1].set(xlabel='Right dome parameter',title='Positive offset scale; finite curvature jumps');axs[1].legend(fontsize=8,frameon=False);axs[1].grid(alpha=.15)
    fig.savefig(output/'geometry_interpolation_error.png',dpi=180);plt.close(fig)
    return diagnostic


def analyze_case(coarse_root,fine_root,out,case,cm,fm):
    cp=coarse_root/f'solution_{case}.json';fp=fine_root/f'solution_{case}.json'
    cr=json.loads(cp.read_text());fr=json.loads(fp.read_text())
    for key in ['parameters_mm','pressure_MPa','wrap_ply_angles_degrees','boundary_model','material_case_note']:
        if cr[key]!=fr[key]:raise ValueError(f'Changed physical protocol: {key}')
    if cr['material_provenance']['material_json_sha256']!=fr['material_provenance']['material_json_sha256']:raise ValueError('Changed material source contract')
    with np.load(coarse_root/f'solution_{case}.npz') as c,np.load(fine_root/f'solution_{case}.npz') as f:
        uc=interpolate(c,cm,f,fm,c['displacement']);uf=f['displacement'];xyz=f['xyz'];delta=uc-uf
        er=xyz.copy();er[:,2]=0;er/=np.linalg.norm(er,axis=1)[:,None]
        urc=np.sum(uc*er,axis=1);urf=np.sum(uf*er,axis=1);radial_delta=urc-urf
        full=float(np.linalg.norm(delta)/np.linalg.norm(uf));radial=float(np.linalg.norm(radial_delta)/np.linalg.norm(urf))
        gc=interpolate(c,cm,f,fm,c['xyz']);geometry_error=np.linalg.norm(gc-xyz,axis=1)
        # Exact same-node reproduction protects indexing, periodic seam and all P2 node types.
        self_u=interpolate(c,cm,c,cm,c['displacement']);identity=float(np.max(np.abs(self_u-c['displacement'])))
        if identity>1e-11:raise ValueError('Source self-interpolation did not reproduce original quadratic nodes')
        cs,fs=layer_statistics(c),layer_statistics(f)
        if cs.keys()!=fs.keys():raise ValueError('Physical ply IDs differ')
        changes={pid:{key:relative_change(cs[pid][key],fs[pid][key]) for key in ['volume_mm3','stress_tensor_weighted_RMS_MPa','stress_tensor_weighted_p50_MPa','stress_tensor_weighted_p95_MPa','strain_energy_Nmm']} for pid in cs}
        np.savez_compressed(out/f'{case}_comparison_fields.npz',xyz=xyz,coarse_interpolated_displacement=uc,
            fine_displacement=uf,displacement_difference=delta,radial_difference=radial_delta,
            coarse_interpolated_geometry=gc,geometry_interpolation_error_mm=geometry_error,
            fine_structured_node_ijk=f['structured_node_ijk'])
    cscalar,fscalar=scalar_metrics(cr),scalar_metrics(fr)
    row=dict(case=case,n_coarse_nodes=cr['n_nodes'],n_fine_nodes=fr['n_nodes'],
        displacement=dict(full_vector_relative_L2_on_all_fine_nodes=full,radial_relative_L2_on_all_fine_nodes=radial,
            full_vector_RMS_difference_mm=float(np.sqrt(np.mean(np.sum(delta*delta,axis=1)))),radial_RMS_difference_mm=float(np.sqrt(np.mean(radial_delta**2))),
            radial_mean_signed_difference_mm=float(radial_delta.mean()),source_self_interpolation_max_abs_mm=identity,
            rigid_motion_removal=False,note='Raw same-gauge fields. Radial projection excludes axial/tangential displacement; no fitted rigid components removed.'),
        geometry_comparison=dict(coarse_to_fine_RMS_mm=float(np.sqrt(np.mean(geometry_error**2))),max_mm=float(geometry_error.max()),
            note='Different triquadratic approximations of the same analytic geometry; parameter interpolation includes this discretization difference.'),
        scalar_metrics=dict(coarse=cscalar,fine=fscalar,relative_changes={k:relative_change(cscalar[k],fscalar[k]) for k in cscalar
            if k not in ['reduced_relative_residual','midsection_axial_force_relative_error']}),
        physical_ply_statistics=dict(coarse=cs,fine=fs,relative_changes=changes),
        source_json_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [cp,fp]})
    return row


def plot(root,rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,3,figsize=(15,4.2),layout='constrained')
    for row,color in zip(rows,['#527da4','#cd8057']):
        case=row['case'];n=[row['n_coarse_nodes'],row['n_fine_nodes']];sc=row['scalar_metrics']
        axs[0].plot(n,[sc['coarse']['strain_energy_Nmm'],sc['fine']['strain_energy_Nmm']],'o-',label=case,color=color)
        axs[1].plot(n,[100*sc['coarse']['midsection_axial_force_relative_error'],100*sc['fine']['midsection_axial_force_relative_error']],'o-',label=case,color=color)
        ps=row['physical_ply_statistics']['relative_changes'];ids=list(ps)
        axs[2].plot(ids,[100*ps[i]['stress_tensor_weighted_RMS_MPa'] for i in ids],'o-',label=case,color=color)
    axs[0].set(xlabel='Original P2 nodes',ylabel='Strain energy (N mm)',title='Energy across two meshes')
    axs[1].set(xlabel='Original P2 nodes',ylabel='Midsection axial-force error (%)',title='Independent pressure resultant check')
    axs[2].set(xlabel='Physical ply ID',ylabel='Relative change (%)',title='Layerwise volume-RMS stress change')
    for ax in axs:ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False);ax.legend(frameon=False)
    fig.suptitle('COPV mesh comparison | same pressure, layup and gauge | two meshes are not a convergence proof',fontsize=11)
    fig.savefig(root/'convergence_summary.png',dpi=180);plt.close(fig)
    fig,axs=plt.subplots(len(rows),2,figsize=(11,4*len(rows)),squeeze=False,layout='constrained')
    for i,row in enumerate(rows):
        a=np.load(root/f"{row['case']}_comparison_fields.npz");x=a['xyz'];er=x.copy();er[:,2]=0;er/=np.linalg.norm(er,axis=1)[:,None]
        ijk=a['fine_structured_node_ijk'];maxk=ijk[:,2].max();sel=(ijk[:,0]==0)&(ijk[:,2]==maxk)
        sort=np.argsort(ijk[sel,1]);z=x[sel,2][sort]
        for key,label in [('coarse_interpolated_displacement','Coarse P2 interpolated'),('fine_displacement','Fine original nodes')]:
            ur=np.sum(a[key]*er,axis=1);axs[i,0].plot(z,ur[sel][sort],label=label)
        axs[i,1].plot(z,a['radial_difference'][sel][sort]);axs[i,0].set_title(row['case']+' outer meridian radial U')
        axs[i,0].set_ylabel('Radial U (mm)');axs[i,1].set_title('Raw coarse minus fine');axs[i,1].set_ylabel('Radial difference (mm)')
        for ax in axs[i]:ax.set_xlabel('z (mm)');ax.grid(alpha=.15)
        axs[i,0].legend(frameon=False)
    fig.savefig(root/'meridian_displacement_comparison.png',dpi=180);plt.close(fig)


def run(coarse,fine,output,cases=('isotropic','reported')):
    coarse=Path(coarse);fine=Path(fine);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    cm=json.loads((coarse/'mesh_manifest.json').read_text());fm=json.loads((fine/'mesh_manifest.json').read_text())
    geometry=geometry_diagnostics(coarse,fine,output,cm,fm)
    rows=[analyze_case(coarse,fine,output,case,cm,fm) for case in cases]
    report=dict(schema_version=1,status='two_mesh_comparison_not_production_approval',coarse_root=str(coarse.resolve()),fine_root=str(fine.resolve()),
        method='Tensor quadratic interpolation using original structured_node_ijk, five segment-local meridian parameters, periodic phi and actual normal-offset coordinates',
        physical_inverse_mapping_used=False,nearest_neighbor_used=False,rigid_motion_fit_used=False,
        stress_peak_convergence_claimed=False,production_ready=False,rows=rows,geometry_diagnostics=geometry,
        limitations=['Two-mesh differences alone do not prove asymptotic convergence','No single-point peak stress convergence criterion',
                     'Physical ply GP statistics remain separate; no cross-interface stress smoothing','Same-gauge raw full-vector and radial differences retained'])
    (output/'comparison.json').write_text(json.dumps(report,indent=2)+'\n');plot(output,rows)
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--coarse',type=Path,required=True);p.add_argument('--fine',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--cases',nargs='+',default=['isotropic','reported'])
    p.add_argument('--case',action='append',choices=['isotropic','reported'],help='Select a single case; repeat to select both')
    p.add_argument('--geometry-only',action='store_true',help='Inspect available meshes without requiring completed solutions')
    a=p.parse_args()
    if a.geometry_only:
        cm=json.loads((a.coarse/'mesh_manifest.json').read_text());fm=json.loads((a.fine/'mesh_manifest.json').read_text())
        r=geometry_diagnostics(a.coarse,a.fine,a.output,cm,fm)
        print(json.dumps({'geometry_only':True,'offset_risk':r['offset_risk'],'sampled_errors':r['sampled_P2_meridian_geometry_error_mm']},indent=2))
    else:
        r=run(a.coarse,a.fine,a.output,a.case or a.cases);print(json.dumps({'status':r['status'],'cases':[{'case':v['case'],'displacement':v['displacement'],'scalars':v['scalar_metrics']} for v in r['rows']]},indent=2))
