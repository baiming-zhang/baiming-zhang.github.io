"""Deterministic pre-label COPV geometry preparation; never launches FE.

The proposed labels originate from a separately solved three-component
axisymmetric reference field. Evaluation on a 3D Hex27 node grid is explicitly
sampling that field, not another independent 3D finite-element solve.
"""
from __future__ import annotations
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist
from scipy.stats import qmc
from geobc.copv_structured_mesh import analytic_segments

COORDINATES=('inner_radius_mm','straight_length_over_radius','dome_depth_over_radius',
    'boss_bore_radius_over_radius','neck_length_over_radius','liner_thickness_over_radius','overwrap_thickness_over_radius')
DIMENSIONS=('inner_radius_mm','straight_length_mm','dome_depth_mm','boss_bore_radius_mm',
    'neck_length_mm','liner_thickness_mm','overwrap_thickness_mm')


def physical_parameters(values):
    v=dict(zip(COORDINATES,np.asarray(values,float)))
    return {DIMENSIONS[i]:float(v[name] if i==0 else v[name]*v['inner_radius_mm']) for i,name in enumerate(COORDINATES)}


def independent_parameters(parameters):
    radius=parameters['inner_radius_mm']
    return np.array([parameters[key] if i==0 else parameters[key]/radius for i,key in enumerate(DIMENSIONS)])


def _cross(a,b):return a[...,0]*b[...,1]-a[...,1]*b[...,0]


def _point_segment_distance(point,start,end):
    delta=end-start
    t=np.clip(np.sum((point-start)*delta,axis=-1)/np.sum(delta*delta,axis=-1),0,1)
    return np.linalg.norm(point-start-t[...,None]*delta,axis=-1)


def polyline_clearance(a,b,known_upper_bound):
    """Continuous segment distances, including crossings, not point distances.

    Bounding balls exclude no segment pair that could improve the supplied
    upper bound: center distance <= gap + both half-segment lengths.  For our
    normal-offset profiles a paired endpoint gives the exact nominal gap as
    a known upper bound on the minimum distance between the polylines.
    """
    a0,a1=a[:-1],a[1:];b0,b1=b[:-1],b[1:]
    half_a=np.linalg.norm(a1-a0,axis=1)/2;half_b=np.linalg.norm(b1-b0,axis=1)/2
    neighbors=cKDTree((b0+b1)/2).query_ball_point((a0+a1)/2,
        known_upper_bound+half_a+half_b.max()+1e-10)
    counts=np.fromiter((len(v) for v in neighbors),int,count=len(neighbors))
    ia=np.repeat(np.arange(len(a0)),counts)
    ib=np.concatenate(neighbors).astype(int)
    if not len(ib):raise ValueError('No clearance segment candidates')
    p,q=a0[ia],a1[ia];r,s=b0[ib],b1[ib]
    distances=np.minimum.reduce([_point_segment_distance(p,r,s),_point_segment_distance(q,r,s),
                                _point_segment_distance(r,p,q),_point_segment_distance(s,p,q)])
    v,w=q-p,s-r;den=_cross(v,w);nonparallel=np.abs(den)>1e-14
    ta=np.divide(_cross(r-p,w),den,out=np.full_like(den,np.inf),where=nonparallel)
    tb=np.divide(_cross(r-p,v),den,out=np.full_like(den,np.inf),where=nonparallel)
    distances[nonparallel&(ta>=0)&(ta<=1)&(tb>=0)&(tb<=1)]=0.
    return float(distances.min())


def geometry_eligibility(parameters,gates):
    """Dense analytic local guards + continuous sampled-polyline clearance.

    Not a global analytic/CAD distance proof. Actual FE mesh Jacobian and
    field-resolution checks are separate, deferred until the reviewed pilot.
    """
    p={key:float(parameters[key]) for key in DIMENSIONS};reason=[];metrics={}
    if any(not np.isfinite(v) or v<=0 for v in p.values()):
        return dict(passed=False,reasons=['nonpositive_or_nonfinite_dimension'],metrics=metrics)
    R,L,D,b,neck,tl,tw=[p[k] for k in DIMENSIONS];wall=tl+tw
    if b>=.5*R:reason.append('boss_radius_not_less_than_half_inner_radius')
    if not .35*R<=D<=1.5*R:reason.append('dome_depth_outside_conservative_geometry_domain')
    if wall>=float(gates['total_wall_over_radius_strict_max'])*R:reason.append('wall_over_radius_excess')
    if reason:return dict(passed=False,reasons=reason,metrics=metrics)
    t=np.linspace(0,1,int(gates['dense_dome_curvature_samples']))[:,None]
    control=np.array([[R,L/2],[R,L/2+.4*D],[b,L/2+.6*D],[b,L/2+D]])
    d=3*((1-t)**2*(control[1]-control[0])+2*t*(1-t)*(control[2]-control[1])+t*t*(control[3]-control[2]))
    dd=6*((1-t)*(control[2]-2*control[1]+control[0])+t*(control[3]-2*control[2]+control[1]))
    speed=np.linalg.norm(d,axis=1)
    if np.any(speed<=0):return dict(passed=False,reasons=['zero_meridian_tangent'],metrics=metrics)
    curvature=_cross(d,dd)/speed**3
    ratio=float(np.max(np.abs(curvature))*wall);parallel=float(min(1.,(1+wall*curvature).min()))
    metrics.update(max_abs_curvature_times_wall=ratio,min_parallel_tangent_scale=parallel,
                   curvature_samples=len(t))
    if ratio>float(gates['sampled_max_abs_curvature_times_total_wall_max']):reason.append('offset_curvature_ratio_excess')
    if parallel<float(gates['sampled_min_parallel_tangent_scale_min']):reason.append('parallel_offset_tangent_scale_too_small')
    if reason:return dict(passed=False,reasons=reason,metrics=metrics)
    counts=(16,256,64,256,16)
    rz,normal=analytic_segments(p,counts)
    profiles=[rz,rz+tl*normal,rz+wall*normal]
    metrics.update(clearance_curve_samples=len(rz),min_profile_radius_mm=float(min(v[:,0].min() for v in profiles)),
        min_profile_axial_sample_increment_mm=float(min(np.diff(v[:,1]).min() for v in profiles)))
    if metrics['min_profile_radius_mm']<=0:reason.append('nonpositive_offset_radius')
    if metrics['min_profile_axial_sample_increment_mm']<=0:reason.append('folded_offset_profile')
    if reason:return dict(passed=False,reasons=reason,metrics=metrics)
    # Positive z progression precludes self-intersection of each individual
    # profile; closest continuous segments additionally detect layer crossings.
    gaps={name:polyline_clearance(profiles[i],profiles[j],gap)
        for name,i,j,gap in [('liner',0,1,tl),('wrap',1,2,tw),('total',0,2,wall)]}
    fraction=float(gates['sampled_nearest_layer_clearance_min_fraction_of_nominal_gap'])
    metrics['polyline_min_layer_clearance_mm']=gaps
    metrics['minimum_clearance_over_nominal']=float(min(gaps['liner']/tl,gaps['wrap']/tw,gaps['total']/wall))
    if metrics['minimum_clearance_over_nominal']<fraction:reason.append('nonlocal_layer_clearance_too_small')
    return dict(passed=not reason,reasons=reason,metrics=metrics,
        scope='Dense analytic local geometry and exact distances between sampled polyline segments; no FE/production claim')


def proposed_grids():
    export=dict(formulation='sampling_grid_only_not_independently_solved_3d',cell_type='Gmsh Hex27',
        nphi=16,meridian_segment_elements=[2,12,16,12,2],liner_elements=1,wrap_elements_per_physical_ply=1)
    reference=dict(formulation='independent_axisymmetric_3component_Quad9_FE',
        meridian_segment_elements=[2,96,64,96,2],liner_elements=4,wrap_elements_per_physical_ply=2)
    ne=sum(export['meridian_segment_elements']);nr=sum(reference['meridian_segment_elements'])
    te=export['liner_elements']+4*export['wrap_elements_per_physical_ply'];tr=reference['liner_elements']+4*reference['wrap_elements_per_physical_ply']
    export.update(nodes=2*export['nphi']*(2*ne+1)*(2*te+1),cells=export['nphi']*ne*te)
    reference.update(nodes=(2*nr+1)*(2*tr+1),cells=nr*tr)
    ratios=[]
    for target,source in zip(export['meridian_segment_elements'],reference['meridian_segment_elements']):
        if source%target:raise ValueError('Proposed meridian grids are not nested')
        ratios.append(source//target)
    if reference['liner_elements']%export['liner_elements'] or reference['wrap_elements_per_physical_ply']%export['wrap_elements_per_physical_ply']:
        raise ValueError('Proposed material/ply offset grids are not nested')
    # One meridian point maps to every phi node by cylindrical vector rotation.
    return dict(frozen=False,reference_field=reference,export_node_grid=export,
        nested_mapping=dict(verified_by_integer_reference_indices=True,meridian_subdivision_ratios=ratios,
            liner_subdivision_ratio=reference['liner_elements']//export['liner_elements'],
            wrap_ply_subdivision_ratio=reference['wrap_elements_per_physical_ply']//export['wrap_elements_per_physical_ply'],
            method='Exact segment-local meridian and physical-ply offset node indices, then cylindrical-to-global rotation at each original phi; no physical inverse or nearest-neighbor'),
        resolution_validation_pending=True,production_resolution_hash=None)


def _best_lhs(n,d,seed,trials=16):
    best=None;score=-1.
    for child in seed.spawn(trials):
        sample=qmc.LatinHypercube(d=d,seed=np.random.default_rng(child)).random(n)
        value=float(pdist(sample).min()) if n>1 else 0.
        if value>score:best,score=sample,value
    return best


def prepare_manifest(candidate_path,output_dir,seed=2026):
    """Write a reviewable draft; never freeze, mesh, solve or replace labels."""
    path=Path(candidate_path);out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    target=out/'draft_manifest.json'
    if target.exists():
        prior=json.loads(target.read_text())
        if prior.get('labels_generated') or prior.get('protocol_frozen'):
            raise RuntimeError('Refuse to replace a labelled or frozen manifest')
    if any(out.rglob('solution*.npz')):
        raise RuntimeError('Refuse to prepare over existing simulation labels')
    config=json.loads(path.read_text());ranges=config['geometry']['candidate_independent_ranges']
    lower=np.array([ranges[n][0] for n in COORDINATES]);upper=np.array([ranges[n][1] for n in COORDINATES])
    baseline=independent_parameters(config['geometry']['baseline_parameters_mm'])
    if np.any(baseline<lower) or np.any(baseline>upper):raise ValueError('Baseline outside candidate range')
    gates=config['geometry']['eligibility_gates_before_any_solve']
    stream=np.random.SeedSequence(seed);pilot_seed,main_seed,retry_seed=stream.spawn(3)
    candidates=[('pilot_baseline',baseline)]
    for j,name in enumerate(COORDINATES):
        for side,value in [('minimum',lower[j]),('maximum',upper[j])]:
            v=baseline.copy();v[j]=value;candidates.append((f'pilot_extreme_{name}_{side}',v))
    pilot_lhs=_best_lhs(9,len(COORDINATES),pilot_seed)
    candidates += [('pilot_lhs',lower+v*(upper-lower)) for v in pilot_lhs]
    main_lhs=qmc.LatinHypercube(d=len(COORDINATES),seed=np.random.default_rng(main_seed)).random(976)
    candidates += [('main_lhs',lower+v*(upper-lower)) for v in main_lhs]
    accepted=[];rejected=[];seen=set();attempt=0;next_candidate=0
    # A pre-label rejected LHS geometry may be deterministically resampled;
    # mandatory baseline/extrema are never disguised as random replacements.
    retry_rng=np.random.default_rng(retry_seed)
    while len(accepted)<1000:
        if attempt>=10000:raise RuntimeError('Too many rejected pre-label geometry candidates')
        if next_candidate<len(candidates):kind,values=candidates[next_candidate];next_candidate+=1
        else:kind='geometry_rejection_resample';values=lower+retry_rng.random(len(COORDINATES))*(upper-lower)
        parameters=physical_parameters(values)
        identity=json.dumps({k:round(parameters[k],12) for k in DIMENSIONS},sort_keys=True,separators=(',',':'))
        check=geometry_eligibility(parameters,gates)
        if identity in seen:check=dict(passed=False,reasons=['duplicate_geometry'],metrics={})
        attempt_id=f'candidate_{attempt:05d}';attempt+=1
        if not check['passed']:
            rejected.append(dict(attempt_id=attempt_id,design_kind=kind,parameters_mm=parameters,eligibility=check))
            if kind=='pilot_baseline' or kind.startswith('pilot_extreme'):
                (out/'draft_rejected_attempts.json').write_text(json.dumps(rejected,indent=2)+'\n')
                raise ValueError(f'Mandatory pilot design failed geometry gate: {kind}: {check["reasons"]}')
            # Preserve the first 24 accepted geometries as the pilot: retry a
            # rejected pilot LHS immediately, without taking a main-set point.
            if kind.startswith('pilot'):
                candidates.insert(next_candidate,('pilot_lhs_geometry_resample',lower+retry_rng.random(len(COORDINATES))*(upper-lower)))
            continue
        seen.add(identity);index=len(accepted)
        accepted.append(dict(case_id=f'copv_{index:04d}',accepted_index=index,attempt_id=attempt_id,
            design_kind=kind,pilot=index<24,parameters_mm=parameters,
            independent_coordinates={k:float(v) for k,v in zip(COORDINATES,values)},
            geometry_sha256=hashlib.sha256(identity.encode()).hexdigest(),eligibility=check,
            label_state='not_generated',solver_output=None))
    order=np.random.default_rng(seed).permutation(np.arange(24,1000))
    assignments={int(i):'train' for i in range(24)}
    for name,ids in [('train',order[:776]),('validation',order[776:876]),('test',order[876:])]:
        assignments.update({int(i):name for i in ids})
    for row in accepted:row['split']=assignments[row['accepted_index']]
    counts={name:sum(row['split']==name for row in accepted) for name in ('train','validation','test')}
    assert counts==dict(train=800,validation=100,test=100)
    assert len(seen)==1000 and sum(row['pilot'] for row in accepted)==24
    assert all(row['split']=='train' for row in accepted[:24])
    assert sum(row['design_kind'].startswith('pilot_extreme') for row in accepted)==14
    report=dict(schema='copv_prelabel_draft_manifest_v1',status='draft_for_root_review_not_production_authorization',
        protocol_frozen=False,production_protocol_hash=None,labels_generated=False,production_ready=False,
        candidate_config_path=str(path.resolve()),candidate_config_content_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        seed=seed,split_seed=seed,random_streams='SeedSequence(seed).spawn(3): pilot LHS, main LHS, pre-label retry; split independent default_rng(seed)',
        geometry_count=1000,pilot_count=24,split_counts=counts,attempt_count=attempt,rejected_attempts=rejected,
        sampling='Baseline + 14 one-coordinate extrema + 9 best-of-16 maximin LHS pilot; 976 LHS main candidates; deterministic pre-label retries logged',
        candidate_physics_and_material=config,grids=proposed_grids(),
        label_provenance='Independent fine three-component axisymmetric FE field sampled to ALL original 3D Hex27 nodes; NOT a direct FE solve on each export grid',
        post_label_policy='Once labels exist, case IDs/geometry/splits are immutable; failed FE cases retain IDs and failure records, never replacement geometries',
        preparation_checks=dict(unique_geometry_count=len(seen),all_geometry_gates_passed=all(r['eligibility']['passed'] for r in accepted),
            all_pilot_in_train=True,mandatory_extreme_count=14,FE_cases_started=0,
            max_offset_curvature_ratio=max(r['eligibility']['metrics']['max_abs_curvature_times_wall'] for r in accepted),
            min_parallel_tangent_scale=min(r['eligibility']['metrics']['min_parallel_tangent_scale'] for r in accepted),
            minimum_polyline_clearance_over_nominal=min(r['eligibility']['metrics']['minimum_clearance_over_nominal'] for r in accepted)),
        limitations=['Draft sampling protocol and grids require root review; no material/resolution/production hash freeze',
            'Eligibility is sampled analytic/polyline geometry; actual solver mesh Jacobians and range-wide field resolution remain necessary',
            'Axisymmetric labels cannot represent nonaxisymmetric supports, defects or loads'],cases=accepted)
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return report
