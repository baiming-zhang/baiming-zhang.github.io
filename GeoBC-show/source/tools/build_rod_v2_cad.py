"""Build the approved rod preview as a fused OpenCASCADE solid (mm).

The shoulder rho values remain loft shape controls, NOT constant-radius fillets.
Section internal/external corners are exact circular arcs. No FE solution runs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import gmsh
import numpy as np


APPROVED = dict(L=160., rb=20., rs=13.75, db=11.5, ds=6.25, Tb=27.5,
                Ts=18., H=18., B=11., tw=3., tf=3., kb=1.55, ks=1.25,
                lb=.22, ls=.18, rB=6., rS=4., rp=1., re=.4, e=0.)


def build(out: Path, mesh_size: float = 3., parameters=None, verification_mesh=True, threads=4,
          fusion_strategy='simultaneous', step_volume_rtol=1e-8, boolean_tolerance=0.,loft_parameterization='ChordLength'):
    out.mkdir(parents=True, exist_ok=True)
    p = APPROVED.copy()
    if parameters:
        if set(parameters) - set(p):
            raise ValueError('Unknown parameter')
        p.update(parameters)
    if any(not math.isfinite(v) for v in p.values()):
        raise ValueError('Nonfinite parameter')
    if p['e'] != 0:
        raise ValueError('This CAD implementation validates the approved straight rod only (e=0).')
    if p['B'] >= min(p['Tb'], p['Ts']) or p['H'] <= 2*p['tf']+2*p['rp']:
        raise ValueError('Incompatible cross section')
    start = time.perf_counter()
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.NumThreads', threads)
        gmsh.option.setNumber('Geometry.ToleranceBoolean',boolean_tolerance)
        gmsh.model.add('approved_rod_v2')
        occ = gmsh.model.occ
        Rs, Rb = p['rs']+p['ds'], p['rb']+p['db']

        def ring(cx, ri, ro, thick):
            outer = occ.addCylinder(cx, 0, -thick/2, 0, 0, thick, ro)
            cutter = occ.addCylinder(cx, 0, -thick, 0, 0, 2*thick, ri)
            solid, _ = occ.cut([(3, outer)], [(3, cutter)])
            occ.synchronize()
            # True rolling-ball fillets on circular ring rims only.
            edges = list({t for d, t in gmsh.model.getBoundary(
                gmsh.model.getBoundary(solid, oriented=False), oriented=False, combined=False)
                          if d == 1 and gmsh.model.getType(d, t) == 'Circle'})
            result = occ.fillet([t for d, t in solid], edges, [p['re']])
            return [(d, t) for d, t in result if d == 3]

        small = ring(0, p['rs'], Rs, p['Ts'])
        big = ring(p['L'], p['rb'], Rb, p['Tb'])
        # Ear rectangles overlap the ring wall; final boolean union removes seams.
        ear_solids = []
        half = .60*Rb
        base = math.sqrt(Rb*Rb-half*half)-.6
        for sign in [-1, 1]:
            low = base if sign == 1 else -1.09*Rb
            ear_solids.append((3, occ.addBox(p['L']-half, low,
                              -p['Tb']/2, 2*half, 1.09*Rb-base, p['Tb'])))

        def smooth(t, scale, length):
            t = max(0., min(1., t))
            a = 1.5+scale/max(2., length*.22)
            return 0. if t == 0 else 1. if t == 1 else t**a/(t**a+(1-t)**a)

        xs = math.sqrt(Rs**2-(p['ks']*p['H']/2)**2)
        xb = p['L']-math.sqrt(Rb**2-(p['kb']*p['H']/2)**2)
        if xs <= p['rs']+p['re'] or xb >= p['L']-p['rb']-p['re']:
            raise ValueError('Shoulder intrudes into bore')
        span = xb-xs

        def section(s):
            H, B = p['H'], p['B']
            if s < p['ls']:
                w = smooth(s/p['ls'], p['rS'], span*p['ls'])
                H *= p['ks']+(1-p['ks'])*w
                B = (p['Ts']-2*p['re'])+(p['B']-p['Ts']+2*p['re'])*w
            elif s > 1-p['lb']:
                w = smooth((s-1+p['lb'])/p['lb'], p['rB'], span*p['lb'])
                H *= 1+(p['kb']-1)*w
                B += (p['Tb']-2*p['re']-p['B'])*w
            x = xs+span*s
            offset = p['e']*.01*(p['L']-Rs-Rb)*math.sin(math.pi*s)**2
            tf, tw = p['tf'], p['tw']
            pts = np.array([[-H/2,-B/2],[-H/2,B/2],[-H/2+tf,B/2],
                            [-H/2+tf,tw/2],[H/2-tf,tw/2],[H/2-tf,B/2],
                            [H/2,B/2],[H/2,-B/2],[H/2-tf,-B/2],
                            [H/2-tf,-tw/2],[-H/2+tf,-tw/2],[-H/2+tf,-B/2]])
            arcs = []
            def point(v):
                return occ.addPoint(x, float(v[0]+offset), float(v[1]))
            for i, q in enumerate(pts):
                a, b = pts[(i-1) % 12]-q, pts[(i+1) % 12]-q
                a, b = a/np.linalg.norm(a), b/np.linalg.norm(b)
                r = p['rp'] if i in (3,4,9,10) else p['re']
                arcs.append((point(q+a*r), point(q+(a+b)*r), point(q+b*r)))
            edges = []
            for i, (a,c,b) in enumerate(arcs):
                edges.extend([occ.addCircleArc(a,c,b), occ.addLine(b,arcs[(i+1) % 12][0])])
            return occ.addWire(edges, checkClosed=True)

        # Split lofts at shoulder boundaries to retain a straight central shank.
        stations = sorted(set([0., p['ls'], 1-p['lb'], 1.] +
                              [p['ls']*i/6 for i in range(1,6)] +
                              [1-p['lb']+p['lb']*i/6 for i in range(1,6)]))
        wires = [section(s) for s in stations]
        # A single global spline would overshoot across the long straight span.
        # Independent shoulders + ruled center keep the center section constant.
        beam = []
        left = [w for s,w in zip(stations,wires) if s <= p['ls']+1e-10]
        right = [w for s,w in zip(stations,wires) if s >= 1-p['lb']-1e-10]
        center = [wires[stations.index(p['ls'])],wires[stations.index(1-p['lb'])]]
        for section_wires, ruled in [(left,False),(center,True),(right,False)]:
            beam.extend(occ.addThruSections(section_wires, makeSolid=True, makeRuled=ruled,
                                           maxDegree=3, continuity='C1', parametrization=loft_parameterization))
        if fusion_strategy=='simultaneous':
            solids, _ = occ.fuse(small, big+ear_solids+[(d,t) for d,t in beam if d == 3])
        elif fusion_strategy=='sequential':
            solids=small
            for tool in [(d,t) for d,t in beam if d==3]+big+ear_solids:
                solids,_=occ.fuse(solids,[tool])
        else:raise ValueError('Unknown fusion strategy')
        occ.removeAllDuplicates()
        occ.synchronize()
        solids = gmsh.model.getEntities(3)
        if len(solids) != 1:
            raise RuntimeError(f'Expected one fused solid, got {solids}')
        # Keep only the solid and its transitive boundary, not construction wires.
        parents = solids
        for dimension in [2,1,0]:
            used = set(gmsh.model.getBoundary(parents, oriented=False, combined=False))
            extra = [dt for dt in gmsh.model.getEntities(dimension) if dt not in used]
            if extra:
                occ.remove(extra,recursive=False)
                occ.synchronize()
            parents = list(used)
        solid = solids[0][1]
        faces = gmsh.model.getBoundary(solids, oriented=False)
        edge_incidence = {}
        for face in faces:
            for d,t in gmsh.model.getBoundary([face], oriented=False):
                edge_incidence[t] = edge_incidence.get(t, 0)+1
        # Periodic cylinder seam edges can occur twice within a single face.
        bad_edges = {str(t):n for t,n in edge_incidence.items() if n != 2}
        if bad_edges:
            raise RuntimeError(f'Nonclosed shell edge incidence: {bad_edges}')
        volume = occ.getMass(3, solid)
        if volume <= 0:
            raise RuntimeError('Nonpositive CAD volume')
        support, load = [], []
        face_report = []
        for d,t in faces:
            area, box = occ.getMass(d,t), occ.getBoundingBox(d,t)
            kind = gmsh.model.getType(d,t)
            if kind == 'Cylinder':
                if abs(area-2*math.pi*p['rb']*(p['Tb']-2*p['re'])) < 1e-3:
                    support.append(t)
                if abs(area-2*math.pi*p['rs']*(p['Ts']-2*p['re'])) < 1e-3:
                    load.append(t)
            face_report.append(dict(tag=t,type=kind,area_mm2=area,bbox_mm=box))
        if len(support) != 1 or len(load) != 1:
            raise RuntimeError(f'Bore selection failed: {support}, {load}')
        for name, tags in [('support_big_bore',support),('load_small_bore_full',load)]:
            tag = gmsh.model.addPhysicalGroup(2,tags)
            gmsh.model.setPhysicalName(2,tag,name)
        tag = gmsh.model.addPhysicalGroup(3,[solid])
        gmsh.model.setPhysicalName(3,tag,'rod_v2_solid')
        gmsh.write(str(out/'connecting_rod.step'))
        gmsh.write(str(out/'connecting_rod.brep'))
        gmsh.write(str(out/'geometry.geo_unrolled'))
        report = dict(parameters=p,units='mm',gmsh_version=gmsh.__version__,fusion_strategy=fusion_strategy,boolean_tolerance_mm=boolean_tolerance,step_volume_rtol=step_volume_rtol,
                      loft_parameterization=loft_parameterization,
                      solid_count=1,volume_mm3=volume,closed_shell_edge_incidence=True,
                      face_count=len(faces),faces=face_report,
                      boundary_semantics={'support_big_bore':support,'load_small_bore_full':load},
                      notes=['Shoulder rho values control loft curvature, not constant-radius fillets.',
                             'Section corners are exact circular arcs; ring rims are OCC fillets.',
                             'Ear-to-ring external edges remain sharp.',
                             'Whole bore semantic labels only; no physical load or FE solve assigned.'])
        # Validate serialized STEP independently; mesh the reopened file, not memory construction.
        original_area=sum(occ.getMass(d,t) for d,t in faces)
        original_bbox=np.asarray(occ.getBoundingBox(3,solid))
        gmsh.clear()
        gmsh.model.add('step_roundtrip')
        occ.importShapes(str(out/'connecting_rod.step'))
        occ.synchronize()
        entities = gmsh.model.getEntities(3)
        revolume = sum(occ.getMass(d,t) for d,t in entities)
        rel = abs(revolume-volume)/volume
        rearea=sum(occ.getMass(d,t) for d,t in gmsh.model.getBoundary(entities,oriented=False))
        area_error=abs(rearea-original_area)/original_area
        bbox_error=float(np.max(np.abs(np.asarray(occ.getBoundingBox(*entities[0]))-original_bbox))) if len(entities)==1 else float('inf')
        if len(entities) != 1 or rel > step_volume_rtol or area_error>1e-7 or bbox_error>5e-5:
            raise RuntimeError(f'STEP roundtrip mismatch: solids={len(entities)}, volume_relative={rel}, area_relative={area_error}, bbox_mm={bbox_error}')
        report['step_roundtrip'] = dict(solid_count=len(entities),volume_mm3=revolume,relative_volume_error=rel,
            relative_surface_area_error=area_error,max_bbox_difference_mm=bbox_error,area_rtol=1e-7,bbox_tolerance_mm=5e-5)
        reopened_groups = {'support_big_bore':[], 'load_small_bore_full':[], 'other_boundary':[]}
        for d,t in gmsh.model.getBoundary(entities, oriented=False):
            area = occ.getMass(d,t)
            name = 'other_boundary'
            if gmsh.model.getType(d,t) == 'Cylinder':
                if abs(area-2*math.pi*p['rb']*(p['Tb']-2*p['re'])) < 1e-3:
                    name = 'support_big_bore'
                elif abs(area-2*math.pi*p['rs']*(p['Ts']-2*p['re'])) < 1e-3:
                    name = 'load_small_bore_full'
            reopened_groups[name].append(t)
        for name, tags in reopened_groups.items():
            if not tags:
                raise RuntimeError('Missing reopened boundary: '+name)
            tag = gmsh.model.addPhysicalGroup(2,tags)
            gmsh.model.setPhysicalName(2,tag,name)
        tag = gmsh.model.addPhysicalGroup(3,[t for d,t in entities])
        gmsh.model.setPhysicalName(3,tag,'rod_v2_solid')
        report['step_reopened_boundary_semantics'] = reopened_groups
        if not verification_mesh:
            report['mesh_verification'] = dict(performed=False, purpose='Independent production TET10 mesher performs mesh validation.')
            report['elapsed_seconds'] = time.perf_counter()-start
            report['sha256'] = {name:hashlib.sha256((out/name).read_bytes()).hexdigest()
                                for name in ['connecting_rod.step','connecting_rod.brep']}
            (out/'cad_validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
            return report
        gmsh.option.setNumber('Mesh.MeshSizeMin', .7)
        gmsh.option.setNumber('Mesh.MeshSizeMax', mesh_size)
        gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 16)
        gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 1)
        gmsh.option.setNumber('Mesh.Algorithm3D', 1)
        gmsh.option.setNumber('Mesh.RandomSeed',2026)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize('Netgen')
        ids, cells = gmsh.model.mesh.getElementsByType(4)
        q = np.asarray(gmsh.model.mesh.getElementQualities(ids,'minSICN'))
        if not len(q) or not np.isfinite(q).all() or q.min() <= 0:
            raise RuntimeError('Invalid tetrahedral mesh')
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        xyz, tags = np.asarray(xyz).reshape(-1,3), np.asarray(tags)
        gmsh.write(str(out/'verification_mesh.msh'))
        gmsh.write(str(out/'surface_preview.stl'))
        np.savez_compressed(out/'verification_mesh.npz',node_ids=tags,xyz=xyz,
                            cell_ids=np.asarray(ids),cells_node_ids=np.asarray(cells).reshape(-1,4),min_sicn=q)
        _, triangles = gmsh.model.mesh.getElementsByType(2)
        sort = np.argsort(tags)
        np.savez_compressed(out/'surface_for_preview.npz',xyz=xyz[sort],
                            triangles=np.searchsorted(tags[sort],np.asarray(triangles).reshape(-1,3)))
        report['mesh_verification'] = dict(element='tetra4',nodes=len(tags),cells=len(ids),
                                           min_sicn=float(q.min()),median_sicn=float(np.median(q)),
                                           purpose='CAD meshability check; not FE convergence validation')
        report['mesh_verification']['solver_ready'] = False
        report['mesh_verification']['quality_note'] = 'Requires local quality checks and FE convergence before simulation use.'
        report['elapsed_seconds'] = time.perf_counter()-start
        report['sha256'] = {name:hashlib.sha256((out/name).read_bytes()).hexdigest()
                            for name in ['connecting_rod.step','connecting_rod.brep','verification_mesh.msh']}
        (out/'cad_validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
        print(json.dumps({k:report[k] for k in ['solid_count','volume_mm3','step_roundtrip','mesh_verification','elapsed_seconds']}))
    finally:
        gmsh.finalize()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path('results/rod_v2_cad/approved_160mm'))
    parser.add_argument('--mesh-size',type=float,default=3.)
    args = parser.parse_args()
    build(args.output.resolve(),args.mesh_size)
