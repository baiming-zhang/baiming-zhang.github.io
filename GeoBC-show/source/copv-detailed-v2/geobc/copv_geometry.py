"""Geometry-only COPV test pieces; no material, winding or FE validation claims.

The dome is a C1 cubic-Bezier neck transition, not an isotensoid. Layers are
sampled along the analytical meridian normal, then interpolated by OCC splines.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import time

DEFAULT_PARAMETERS = dict(inner_radius_mm=100., straight_length_mm=400.,
                          dome_depth_mm=60., boss_bore_radius_mm=15.,
                          neck_length_mm=20., liner_thickness_mm=1., overwrap_thickness_mm=3.)


def generate_copv_geometry(output_dir, parameters=None, dome_samples=101):
    import gmsh
    import numpy as np
    started = time.perf_counter()
    p = dict(DEFAULT_PARAMETERS)
    if parameters:
        if set(parameters)-set(p):
            raise ValueError("Unknown COPV geometry parameters")
        p.update(parameters)
    if any(not np.isfinite(v) or v <= 0 for v in p.values()):
        raise ValueError("Positive finite dimensions required")
    R, length, depth, bore, neck, tl, tw = p.values()
    if not (bore < .5*R and .35*R <= depth <= 1.5*R and tl+tw < .1*R and dome_samples >= 25):
        raise ValueError("Geometry outside conservative pilot range")
    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    # (radius,z) control polygon, with z-parallel tangents at both ends.
    control = np.array([[R,length/2], [R,length/2+.4*depth],
                        [bore,length/2+.6*depth], [bore,length/2+depth]])
    t = np.linspace(0,1,dome_samples)[:,None]
    right = (1-t)**3*control[0]+3*(1-t)**2*t*control[1]+3*(1-t)*t*t*control[2]+t**3*control[3]
    deriv = 3*((1-t)**2*(control[1]-control[0])+2*t*(1-t)*(control[2]-control[1])+t*t*(control[3]-control[2]))
    second = 6*((1-t)*(control[2]-2*control[1]+control[0])+t*(control[3]-2*control[2]+control[1]))
    norm = np.stack([deriv[:,1],-deriv[:,0]],axis=1)
    norm /= np.linalg.norm(norm,axis=1)[:,None]
    curvature = (deriv[:,0]*second[:,1]-deriv[:,1]*second[:,0])/np.linalg.norm(deriv,axis=1)**3
    offset_curvature_ratio = float(np.max(np.abs(curvature))*(tl+tw))
    if offset_curvature_ratio >= .8:
        raise ValueError("Offset approaches local curvature radius; refine design")
    base = [np.array([[bore,-length/2-depth-neck],[bore,-length/2-depth]]),
            right[::-1]*[1,-1], np.array([[R,-length/2],[R,length/2]]),
            right, np.array([[bore,length/2+depth],[bore,length/2+depth+neck]])]
    normals = [np.array([[1.,0.],[1.,0.]]),norm[::-1]*[1,-1],
               np.array([[1.,0.],[1.,0.]]),norm,np.array([[1.,0.],[1.,0.]])]
    layers = [[v+n*offset for v,n in zip(base,normals)] for offset in [0.,tl,tl+tw]]
    all_meridians = [np.concatenate([s if i==0 else s[1:] for i,s in enumerate(layer)]) for layer in layers]
    all_normals = np.concatenate([s if i==0 else s[1:] for i,s in enumerate(normals)])
    max_normal_error = float(max(np.max(np.abs((all_meridians[j]-all_meridians[0])-o*all_normals))
                                 for j,o in [(1,tl),(2,tl+tw)]))
    for profile in all_meridians:
        if np.min(profile[:,0]) <= 0 or np.any(np.diff(profile[:,1]) <= 0):
            raise ValueError("Nonpositive radius or folded normal offset")
    np.savez_compressed(out/'meridian.npz', inner_rz=all_meridians[0],
                        interface_rz=all_meridians[1], outer_rz=all_meridians[2],
                        outward_normal_rz=all_normals, inner_dome_control_rz=control)
    gmsh.initialize()
    try:
        gmsh.model.add("copv_geometry_pilot")
        occ = gmsh.model.occ
        endpoint_cache = {}
        def point(rz):
            key = tuple(np.round(rz,10))
            if key not in endpoint_cache:
                endpoint_cache[key]=occ.addPoint(float(rz[0]),0,float(rz[1]))
            return endpoint_cache[key]
        curves = []
        for layer in layers:
            tags=[]
            for segment in layer:
                points = [point(v) for v in segment]
                tags.append(occ.addLine(*points) if len(points)==2 else occ.addSpline(points))
            curves.append(tags)
        material_faces=[]
        for j in [0,1]:
            end = occ.addLine(point(layers[j][-1][-1]), point(layers[j+1][-1][-1]))
            start = occ.addLine(point(layers[j+1][0][0]),point(layers[j][0][0]))
            loop = occ.addCurveLoop(curves[j]+[end]+[-v for v in curves[j+1][::-1]]+[start])
            material_faces.append((2,occ.addPlaneSurface([loop])))
        result = occ.revolve(material_faces,0,0,0,0,0,1,2*np.pi)
        volumes = [dt for dt in result if dt[0]==3]
        occ.synchronize()
        tm=(t[:-1]+t[1:])/2
        exact=(1-tm)**3*control[0]+3*(1-tm)**2*tm*control[1]+3*(1-tm)*tm*tm*control[2]+tm**3*control[3]
        dd=3*((1-tm)**2*(control[1]-control[0])+2*tm*(1-tm)*(control[2]-control[1])+tm*tm*(control[3]-control[2]))
        nn=np.stack([dd[:,1],-dd[:,0]],axis=1);nn/=np.linalg.norm(nn,axis=1)[:,None]
        cad_offset_error=[]
        for j,offset in enumerate([0.,tl,tl+tw]):
            expected=exact+offset*nn
            xyz=np.stack([expected[:,0],np.zeros(len(expected)),expected[:,1]],axis=1)
            closest,_=gmsh.model.getClosestPoint(1,curves[j][3],xyz.ravel())
            cad_offset_error.append(float(np.max(np.linalg.norm(np.array(closest).reshape(-1,3)-xyz,axis=1))))
        if len(volumes)!=2:
            raise RuntimeError(f"Expected two layered volumes, got {volumes}")
        frag,mapping = occ.fragment([volumes[0]],[volumes[1]])
        occ.synchronize()
        volumes=[dt for dt in frag if dt[0]==3]
        if len(volumes)!=2:
            raise RuntimeError("Fragment changed layer count")
        liner = [v for d,v in mapping[0] if d==3]
        wrap = [v for d,v in mapping[1] if d==3]
        def boundaries(tags):
            return {tag for d,tag in gmsh.model.getBoundary([(3,t) for t in tags],oriented=False) if d==2}
        bl,bw=boundaries(liner),boundaries(wrap)
        interface=bl&bw
        if not interface:
            raise RuntimeError("No conformal shared liner/overwrap interface")
        groups={"liner_volume":liner,"overwrap_volume":wrap,"bonded_interface":sorted(interface),
                "wetted_inner":[],"external_outer":[],"neck_end_liner":[],"neck_end_overwrap":[]}
        for tags,body,endname in [(bl-interface,"wetted_inner","neck_end_liner"),(bw-interface,"external_outer","neck_end_overwrap")]:
            for tag in tags:
                if gmsh.model.getType(2,tag)=="Plane": groups[endname].append(tag)
                else: groups[body].append(tag)
        if any(not tags for tags in groups.values()):
            raise RuntimeError(f"Missing semantic group: {groups}")
        for index,(name,tags) in enumerate(groups.items(),1):
            dim=3 if name.endswith("volume") else 2
            gmsh.model.addPhysicalGroup(dim,tags,index);gmsh.model.setPhysicalName(dim,index,name)
        volume_values={name:sum(occ.getMass(3,t) for t in groups[name]) for name in ['liner_volume','overwrap_volume']}
        if min(volume_values.values())<=0:
            raise RuntimeError("Nonpositive material volume")
        gmsh.write(str(out/'geometry.step'))
        gmsh.write(str(out/'geometry.brep'))
        gmsh.write(str(out/'geometry_groups.geo_unrolled'))
        # Gmsh serializes an XAO companion. Keep its reference portable.
        (out/'geometry_groups.geo_unrolled').write_text('Merge "geometry_groups.geo_unrolled.xao";\n',encoding='utf-8')
        checks=dict(positive_material_volumes_mm3=volume_values,shared_interface_face_count=len(interface),
                    max_sample_normal_offset_error_mm=max_normal_error,
                    max_offset_curvature_ratio=offset_curvature_ratio,
                    midpoint_CAD_profile_error_mm=cad_offset_error,
                    positive_radius=True,monotone_axial_profiles=True,cylinder_dome_neck_C1_analytic=True)
        report=dict(schema_version=1,status="geometry_only_unvalidated",production_ready=False,
                    parameters_mm=p,gmsh_version=gmsh.__version__,dome_samples=dome_samples,
                    profile="C1 cubic-Bezier transition, OCC interpolated meridian normal offsets",
                    checks=checks,semantic_entities=groups,
                    limitations=["Double-open neck test piece: end pressure closure/resultant not yet implemented",
                                 "No layup, orthotropic material constants, winding-path or manufacturability validation",
                                 "No FE solve, pressure selection or production mesh",
                                 "Constant normal thickness prescribed analytically; reported CAD interpolation error is sampled, not a global bound",
                                 "Bezier dome includes a neck inflection; not an isotensoid or certified vessel design"],
                    generation_seconds=time.perf_counter()-started)
        report['artifacts']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file() and p.suffix in {'.step','.brep','.geo_unrolled','.npz','.xao'}}
        (out/'geometry_manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        return report
    finally:
        gmsh.finalize()


def plot_copv_geometry(output_dir):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=Path(output_dir);a=np.load(out/'meridian.npz')
    fig,(left,right)=plt.subplots(1,2,figsize=(12,5),gridspec_kw={'width_ratios':[1.4,1]})
    for key,color,label in [('inner_rz','#7389a5','Inner cavity wall'),('interface_rz','#d79d43','Bonded interface'),('outer_rz','#4e9c89','Overwrap outer wall')]:
        rz=a[key];left.plot(rz[:,1],rz[:,0],color=color,label=label);left.plot(rz[:,1],-rz[:,0],color=color)
    inn,inter,outer=(a[k] for k in ['inner_rz','interface_rz','outer_rz'])
    left.fill_between(inn[:,1],inn[:,0],np.interp(inn[:,1],inter[:,1],inter[:,0]),color='#d79d43',alpha=.25)
    left.set_aspect('equal');left.set_xlabel('Axial z (mm)');left.set_ylabel('Radius r (mm)');left.legend(fontsize=8)
    # Actual 3D revolution of sampled outer meridian, no FE fields involved.
    fig.delaxes(right);right=fig.add_subplot(122,projection='3d',proj_type='ortho')
    theta=np.linspace(0,2*np.pi,65)
    rr=outer[:,0,None];zz=np.broadcast_to(outer[:,1,None],(len(outer),len(theta)))
    right.plot_surface(rr*np.cos(theta),rr*np.sin(theta),zz,color='#7caaa2',alpha=.9,linewidth=0)
    right.set_box_aspect((2*outer[:,0].max(),2*outer[:,0].max(),np.ptp(outer[:,1])))
    right.view_init(18,-60);right.set_axis_off();right.set_title('Open-neck geometry preview',fontsize=10)
    fig.suptitle('COPV geometry pilot | normal-offset layers | geometry only',fontsize=13)
    fig.tight_layout();fig.savefig(out/'geometry_preview.png',dpi=180,bbox_inches='tight');plt.close(fig)
