"""Detailed C01 CAD concept assembly, not a qualified pressure-vessel design.

Ports/seal-stack dimensions are explicit design choices. Nominal materials are
identities for construction/rendering; no new FE material/pressure validation.
"""
import argparse,base64,gzip,hashlib,json,time
from pathlib import Path
import numpy as np
import gmsh
from geobc.copv_structured_mesh import analytic_segments
from geobc.copv_geometry import DEFAULT_PARAMETERS

MATERIALS={
 'aluminum':{'label':'Al 6061 内衬','color':[.66,.72,.76],'kind':1},
 'carbon':{'label':'碳纤维 / 环氧','color':[.075,.085,.092],'kind':0},
 'titanium':{'label':'Ti-6Al-4V boss','color':[.49,.52,.55],'kind':2},
 'steel':{'label':'316L 接口套','color':[.70,.73,.75],'kind':3},
 'elastomer':{'label':'弹性体 O 形圈','color':[.045,.052,.047],'kind':4},
 'peek':{'label':'PEEK 挡圈','color':[.64,.54,.35],'kind':5},
 'adhesive':{'label':'环氧粘接界面','color':[.53,.31,.13],'kind':6},
}
def encoded(a,dtype):return base64.b64encode(np.asarray(a,dtype=dtype).tobytes()).decode()
def write(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

def build(out,parameters):
    start=time.perf_counter();out.mkdir(parents=True,exist_ok=True);p=dict(DEFAULT_PARAMETERS);p.update(parameters)
    R,L,D,b,neck,tl,tw=p.values();shoulder=L/2+D;end=shoulder+neck;cutz=shoulder-8.;liner_end=end-3.;adh=.15
    rz,n=analytic_segments(p,(2,50,16,50,2));parts=[]
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal',1);gmsh.model.add('COPV_C01_detailed_construction_concept');occ=gmsh.model.occ
        def poly_solid(profile,sign=1):
            profile=np.array(profile,float)*[1,sign]
            pts=[occ.addPoint(float(r),0,float(z)) for r,z in profile]
            lines=[occ.addLine(pts[i],pts[(i+1)%len(pts)]) for i in range(len(pts))]
            face=occ.addPlaneSurface([occ.addCurveLoop(lines)])
            return [dt for dt in occ.revolve([(2,face)],0,0,0,0,0,1,2*np.pi) if dt[0]==3]
        def profile_curves(profile):
            # Separate straight segments and domes retain exact joins.
            cuts=np.cumsum([0,4,100,32,100,4]);curves=[];cache={}
            def point(v):
                key=tuple(np.round(v,11))
                if key not in cache:cache[key]=occ.addPoint(float(v[0]),0,float(v[1]))
                return cache[key]
            for i in range(5):
                segment=profile[cuts[i]:cuts[i+1]+1];pts=[point(q) for q in segment]
                curves.append(occ.addSpline(pts) if i in (1,3) else occ.addLine(pts[0],pts[-1]))
            return curves,point(profile[0]),point(profile[-1])
        def shell(inner,outer,limit):
            a,a0,a1=profile_curves(inner);c,c0,c1=profile_curves(outer)
            loop=occ.addCurveLoop(a+[occ.addLine(a1,c1)]+[-i for i in c[::-1]]+[occ.addLine(c0,a0)])
            face=occ.addPlaneSurface([loop]);vol=[dt for dt in occ.revolve([(2,face)],0,0,0,0,0,1,2*np.pi) if dt[0]==3]
            clip=occ.addCylinder(0,0,-limit,0,0,2*limit,2*R)
            result,_=occ.intersect(vol,[(3,clip)]);return [dt for dt in result if dt[0]==3]
        def enclosure(profile):
            c,a,z=profile_curves(profile);axis_top=occ.addPoint(0,0,float(profile[-1,1]));axis_bottom=occ.addPoint(0,0,float(profile[0,1]))
            loop=occ.addCurveLoop(c+[occ.addLine(z,axis_top),occ.addLine(axis_top,axis_bottom),occ.addLine(axis_bottom,a)])
            face=occ.addPlaneSurface([loop]);vol=[dt for dt in occ.revolve([(2,face)],0,0,0,0,0,1,2*np.pi) if dt[0]==3]
            clip=occ.addCylinder(0,0,-liner_end,0,0,2*liner_end,2*R);return [dt for dt in occ.intersect(vol,[(3,clip)])[0] if dt[0]==3]
        def add(name,material,vol,**kwargs):
            if not vol:raise RuntimeError('Empty part '+name)
            parts.append({'name':name,'material':material,'volumes':vol,**kwargs});return parts[-1]
        def cut(vol,tools):return [dt for dt in occ.cut(vol,tools,removeObject=True,removeTool=False)[0] if dt[0]==3]
        liner=add('金属气密内衬','aluminum',shell(rz,rz+tl*n,liner_end),group='liner')
        add('环氧粘接层 0.15 mm','adhesive',shell(rz+tl*n,rz+(tl+adh)*n,cutz),group='bond')
        for i,angle in enumerate([55,-55,-55,55]):
            a=tl+adh+i*tw/4;c=tl+adh+(i+1)*tw/4
            add(f'包覆铺层 {i+1} / {angle:+d}°','carbon',shell(rz+a*n,rz+c*n,cutz),group='ply',ply=i+1,angle=angle)
        # Separate polar reinforcement sleeves, with a tapered meridian edge.
        for sign,label in [(1,'上'),(-1,'下')]:
            select=(sign*rz[:,1]>=L/2+.45*D)&(sign*rz[:,1]<=cutz+3)
            ids=np.flatnonzero(select);a=rz[ids]+(tl+adh+tw)*n[ids]
            zsign=sign*rz[ids,1];q=np.clip((zsign-(L/2+.45*D))/(.38*D),0,1);extra=1.4*np.sin(q*np.pi/2)**2+.03
            outer=a+extra[:,None]*n[ids]
            # Piecewise polygon representation is recorded for these collars.
            vol=poly_solid(np.vstack([a,outer[::-1]]))
            add(label+'极区局部补强','carbon',vol,group='reinforcement',angle=25*sign)
        excluded=enclosure(rz+tl*n)
        boss_parts=[]
        for sign,label in [(1,'上端充装口'),(-1,'下端服务口')]:
            ext=10. if sign==1 else 6.
            outer=[(0,cutz-5),(39,cutz-5),(44,cutz-2),(44,cutz+2),(38,cutz+4),(26,shoulder+2),(26,end-1),(24,end+1),(24,end+ext-1),(22,end+ext),(0,end+ext)]
            raw=poly_solid(outer,sign)
            bore=occ.addCylinder(0,0,sign*(liner_end),0,0,sign*40,b+3.2)
            boss=cut(raw,excluded+[(3,bore)])
            occ.remove([(3,bore)],recursive=True)
            for z in [shoulder+5,shoulder+8]:
                groove=occ.addTorus(0,0,sign*z,26,.42)
                boss=cut(boss,[(3,groove)]);occ.remove([(3,groove)],recursive=True)
            bp=add(label+' · 嵌入法兰与车削 boss','titanium',boss,group='boss',end=sign);boss_parts+=boss
        # Relieve composite/adhesive against actual boss CAD; no cosmetic overlap.
        for part in parts:
            if part['group'] in ('bond','ply','reinforcement'):part['volumes']=cut(part['volumes'],boss_parts)
        occ.remove(excluded,recursive=True)
        for sign,label in [(1,'上端充装口'),(-1,'下端服务口')]:
            ext=22. if sign==1 else 15.;bore=.68*b if sign==1 else .55*b
            profile=[(bore,liner_end),(b+3.1,liner_end),(b+3.1,end+10),(b+1,end+12),(b+1,end+ext-1),(bore+.8,end+ext),(bore,end+ext)]
            sleeve=poly_solid(profile,sign)
            # Six-flat wrench land, connected to the insert; no decorative bolts.
            hexbottom=sign*(end+8 if sign==1 else end+6);height=sign*5
            pts=[occ.addPoint(22*np.cos(t),22*np.sin(t),hexbottom) for t in np.arange(6)*np.pi/3]
            lines=[occ.addLine(pts[i],pts[(i+1)%6]) for i in range(6)]
            face=occ.addPlaneSurface([occ.addCurveLoop(lines)]);hexvol=[dt for dt in occ.extrude([(2,face)],0,0,height) if dt[0]==3]
            hole=occ.addCylinder(0,0,hexbottom-sign,0,0,height+sign*2,bore)
            hexvol=cut(hexvol,[(3,hole)]);occ.remove([(3,hole)],recursive=True)
            sleeve=[dt for dt in occ.fuse(sleeve,hexvol)[0] if dt[0]==3]
            bp=next(a for a in parts if a['group']=='boss' and a['end']==sign)
            # Sleeve/hex receives an explicit mating pocket in boss.
            bp['volumes']=cut(bp['volumes'],sleeve)
            seals=[]
            for index,z in enumerate([end+2,end+5]):
                groove=occ.addTorus(0,0,sign*z,b+3.12,.82)
                bp['volumes']=cut(bp['volumes'],[(3,groove)]);sleeve=cut(sleeve,[(3,groove)])
                occ.remove([(3,groove)],recursive=True)
                tor=occ.addTorus(0,0,sign*z,b+3.12,.76)
                seals.append(add(label+f' · O 形圈 {index+1}','elastomer',[(3,tor)],group='seal',end=sign))
            z=end+3.35;backup=poly_solid([(b+2.3,z),(b+3.95,z),(b+3.95,z+.65),(b+2.3,z+.65)],sign)
            bp['volumes']=cut(bp['volumes'],backup);sleeve=cut(sleeve,backup)
            add(label+' · PEEK 挡圈','peek',backup,group='backup',end=sign)
            add(label+' · 中空接口套 / 六角扳手位','steel',sleeve,group='insert',end=sign)
        occ.synchronize()
        # Remove only residual construction faces/curves not bounding any solid.
        used=set(abs(t) for _,t in gmsh.model.getBoundary(gmsh.model.getEntities(3),combined=False,oriented=True) if _==2)
        orphan=[(2,t) for _,t in gmsh.model.getEntities(2) if t not in used]
        if orphan:occ.remove(orphan,recursive=True);occ.synchronize()
        for i,part in enumerate(parts,1):
            volumes=part['volumes'];mass=sum(occ.getMass(3,t) for _,t in volumes)
            if not np.isfinite(mass) or mass<=0:raise RuntimeError('Nonpositive solid '+part['name'])
            part['volume_mm3']=float(mass);part['solid_count']=len(volumes)
            gmsh.model.addPhysicalGroup(3,[t for _,t in volumes],i);gmsh.model.setPhysicalName(3,i,f'part_{i}_{part["material"]}')
        gmsh.write(str(out/'detailed_assembly.step'));gmsh.write(str(out/'detailed_assembly.brep'))
        for key,value in [('Mesh.MeshSizeMin',.65),('Mesh.MeshSizeMax',7.),('Mesh.MeshSizeFromCurvature',48),('Mesh.ElementOrder',1),('Mesh.Algorithm',6),('Mesh.MeshSizeExtendFromBoundary',0)]:gmsh.option.setNumber(key,value)
        def emit_model(filename,section=False):
            gmsh.model.mesh.generate(2);modelparts=[]
            for part in parts:
                positions=[];normals=[];triangles=[];count=0
                boundary=gmsh.model.getBoundary(part['volumes'],combined=True,oriented=True)
                for dim,face in boundary:
                    if dim!=2:continue
                    faceid=abs(face);tags,x,uv=gmsh.model.mesh.getNodes(2,faceid,includeBoundary=True,returnParametricCoord=True)
                    x=np.asarray(x).reshape(-1,3);normal=np.asarray(gmsh.model.getNormal(faceid,uv)).reshape(-1,3)
                    if face<0:normal=-normal
                    index={int(t):i+count for i,t in enumerate(tags)};types,_,nodes=gmsh.model.mesh.getElements(2,faceid)
                    for typ,ns in zip(types,nodes):
                        if typ!=2:raise RuntimeError('Expected Tri3 display tessellation')
                        triangles.append(np.array([index[int(t)] for t in ns],np.uint32).reshape(-1,3))
                    positions.append(x);normals.append(normal);count+=len(x)
                entry={k:v for k,v in part.items() if k!='volumes'}
                entry.update(xyz=encoded(np.concatenate(positions),'<f4'),normal=encoded(np.concatenate(normals),'<f4'),tri=encoded(np.concatenate(triangles),'<u4'))
                modelparts.append(entry)
            result={'schema':'copv_c04_detailed_assembly_v2_interface_alignment','case_id':out.parent.name,'parameters_mm':p,'materials':MATERIALS,'parts':modelparts,'section':section,'scope':'Parametric CAD construction concept. Port sizes, seal stack, metal identities and local reinforcement are illustrative design choices, not a qualified product.','display_tessellation':'Tri3 of actual component CAD surfaces, not a volume FE mesh. Section is a true CAD quarter cut, retaining solid cut faces.','cavity':'Empty storage cavity. Metallic liner, adhesive, four structural ply bands and polar reinforcement; no invented internal baffles.','thread_display':'Procedural helical thread texture on insert bore; not a CAD thread solid or a qualified thread form.','original_solver_scope':'C04 geometry input only until this case passes mesh gates and has its own GPU-FEniCSx result; no other geometry field is reused.','seconds':time.perf_counter()-start}
            (out/filename).write_bytes(gzip.compress(json.dumps(result,separators=(',',':'),ensure_ascii=False).encode(),mtime=0))
            return result,modelparts
        result,modelparts=emit_model('model.json.gz')
        manifest={k:v for k,v in result.items() if k!='parts'};manifest['parts']=[{k:v for k,v in a.items() if k not in ('xyz','normal','tri')} for a in modelparts];manifest['gmsh_version']=gmsh.__version__
        manifest['artifact_sha256']={f:hashlib.sha256((out/f).read_bytes()).hexdigest() for f in ['detailed_assembly.step','detailed_assembly.brep','model.json.gz']}
        write(out/'assembly_manifest.json',manifest)
        gmsh.model.mesh.clear();gmsh.model.removePhysicalGroups()
        box=occ.addBox(0,-2*R,-end-60,2*R,2*R,2*end+120)
        for part in parts:part['volumes']=cut(part['volumes'],[(3,box)])
        occ.remove([(3,box)],recursive=True);occ.synchronize()
        used=set(abs(t) for dim,t in gmsh.model.getBoundary(gmsh.model.getEntities(3),combined=False,oriented=True) if dim==2)
        orphan=[(2,t) for dim,t in gmsh.model.getEntities(2) if t not in used]
        if orphan:occ.remove(orphan,recursive=True);occ.synchronize()
        section,_=emit_model('section.json.gz',True)
        manifest['section_sha256']=hashlib.sha256((out/'section.json.gz').read_bytes()).hexdigest()
        manifest['total_build_seconds']=time.perf_counter()-start;write(out/'assembly_manifest.json',manifest)
        print(json.dumps({'parts':len(parts),'solids':sum(a['solid_count'] for a in parts),'seconds':manifest['total_build_seconds'],'model_bytes':(out/'model.json.gz').stat().st_size,'section_bytes':(out/'section.json.gz').stat().st_size}),flush=True)
    finally:gmsh.finalize()

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--parameters',type=Path,required=True);args=ap.parse_args();build(args.output,json.loads(args.parameters.read_text(encoding='utf-8')))
