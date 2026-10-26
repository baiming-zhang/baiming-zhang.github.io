"""Measure CAD corner tangent angles on already identified low-quality faces."""
from pathlib import Path
import argparse,itertools,json
import gmsh,numpy as np
W=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('candidate',type=Path);a=ap.parse_args();c=a.candidate.resolve();o=c/'port_surface_diagnosis_v1/corner_angles_v1.json';assert not o.exists();r=json.loads((c/'port_surface_diagnosis_v1/review.json').read_text());review=json.loads((c/'mesh/mesh_review.json').read_text());state=json.loads((W/review['parent_attempt']/'state.json').read_text());rows=[]
gmsh.initialize();gmsh.option.setNumber('General.NumThreads',1)
try:
 gmsh.model.occ.importShapes(state['partition_file']);gmsh.model.occ.synchronize()
 for f in r['bad_faces']:
  vertices={}
  for curve in {e['curve'] for e in f['curves']}:
   lo,hi=gmsh.model.getParametrizationBounds(1,curve);bounds=[float(lo[0]),float(hi[0])]
   ends=list({abs(v) for dim,v in gmsh.model.getBoundary([(1,curve)],combined=False)})
   for u,sign in [(bounds[0],1),(bounds[1],-1)]:
    pos=np.asarray(gmsh.model.getValue(1,curve,[u]));tangent=sign*np.asarray(gmsh.model.getDerivative(1,curve,[u]));tangent/=np.linalg.norm(tangent)
    for v in ends:
     xyz=np.asarray(gmsh.model.getValue(0,v,[]))
     if np.linalg.norm(xyz-pos)<1e-7:vertices.setdefault(v,[]).append(dict(curve=curve,tangent=tangent,xyz=xyz))
  corners=[]
  for vertex,edges in vertices.items():
   for p,q in itertools.combinations(edges,2):
    if p['curve']==q['curve']:continue
    corners.append(dict(vertex=vertex,curves=[p['curve'],q['curve']],xyz_mm=p['xyz'].tolist(),smaller_outgoing_tangent_angle_deg=float(np.degrees(np.arccos(np.clip(p['tangent']@q['tangent'],-1,1))))))
  rows.append(dict(surface=f['surface'],type=f['type'],part_ids=f['part_ids'],min_surface_SICN=f['min_surface_SICN'],corners=corners))
finally:gmsh.finalize()
o.write_text(json.dumps(dict(schema='c04_CAD_corner_tangent_angles_v1',faces=rows,scope='Exact boundary-curve endpoint derivatives. Smaller outgoing angle does not by itself distinguish convex interior from reflex exterior angle; use trimmed surface containment before concluding impossibility.'),indent=2)+'\n');print(json.dumps([dict(surface=f['surface'],angles_deg=[v['smaller_outgoing_tangent_angle_deg'] for v in f['corners']]) for f in rows]))
