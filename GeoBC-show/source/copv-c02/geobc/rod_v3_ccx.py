"""Independent CCX pilot adapter; original TET10 nodes and radial pin retained."""
import json,math,re,hashlib,time,subprocess,os
from pathlib import Path
import numpy as np
from geobc.rod_io import write_json

def material_contract(source):
    E=115100.;points=[]
    # NASA/TM2016-218864 TableII,20C,1e-3/s column. Offset-yield
    # engineering total strain = sigma/E+offset. Before necking:
    # true sigma=sigma_eng*(1+eps_eng); eps_pl=log(1+eps_eng)-sigma_true/E.
    # PL is a summary statistic, so the initial onset is an explicit approximation.
    points.append([781.6,0.])
    engineering=[]
    for stress,offset,total in [(881.1,.0002,None),(925.5,.002,None),(941.3,None,.018)]:
        strain=stress/E+offset if total is None else total
        true=stress*(1+strain);plastic=math.log1p(strain)-true/E
        assert plastic>points[-1][1] and true>points[-1][0]
        points.append([true,plastic]);engineering.append(dict(stress_MPa=stress,total_engineering_strain=strain,offset=offset))
    return dict(E_MPa=E,nu=.28,density_tonne_mm3=4.43e-9,plastic_true_stress_vs_equivalent_plastic_strain=points,
        source_url='https://ntrs.nasa.gov/api/citations/20160008952/downloads/20160008952.pdf',source_sha256=hashlib.sha256(Path(source).read_bytes()).hexdigest(),
        source_pages=dict(PDF_page=35,printed_page=31,table='II',material_pedigree_PDF_page=7),temperature_C=20,source_strain_rate_s_1=.001,
        source_engineering_points=engineering,visually_reviewed=True,plastic_strain_data_limit=points[-1][1],
        approximation='Piecewise isotropic rate-independent J2 summary-point approximation, not a digitized full curve or certified rod material. PL781.6MPa used as approximate initial yield; TableII combines different tests. Source plate is textured/orthotropic; isotropic reduction is declared.',
        inherited_assumptions='nu=.28 and density4430kg/m3 retained assumptions; neither inferred from TableII.',
        limitations='20C quasistatic reference only; no validated high-rate, temperature, cyclic Bauschinger, damage/fracture/fatigue behavior. Beyond last tabulated plastic strain is unsupported and flagged, never accepted silently.')

def ccx_order(cells):
    import gmsh
    gmsh.initialize()
    try:q=np.asarray(gmsh.model.mesh.getElementProperties(11)[4]).reshape(10,3)
    finally:gmsh.finalize()
    pairs=[(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)];order=list(range(4))
    for i,j in pairs:
        distance=np.linalg.norm(q-(q[i]+q[j])/2,axis=1);k=int(np.argmin(distance));assert distance[k]<1e-12;order.append(k)
    assert len(set(order))==10
    return cells[:,order],order

def read_dat(file):
    blocks=[];current=None
    for line in Path(file).read_text(errors='replace').splitlines():
        t=line.strip()
        if 'for set' in t and 'time' in t:
            kind='U' if t.startswith('displacements') else 'RF' if t.startswith('forces') else 'S' if t.startswith('stresses') else 'PEEQ' if t.startswith('equivalent plastic strain') else None
            current=None
            if kind:
                match=re.search(r'time\s+([0-9.Ee+\-]+)',t);current=dict(kind=kind,time=float(match.group(1)),rows=[]);blocks.append(current)
        elif current and t:
            fields=t.split()
            try:values=[float(x) for x in fields]
            except ValueError:current=None;continue
            expected={'U':4,'RF':4,'S':8,'PEEQ':3}[current['kind']]
            if len(values)==expected:current['rows'].append(values)
    for b in blocks:b['rows']=np.asarray(b['rows'],float)
    return blocks

def write_deck(file,m,bc,force,material,plastic=False,nlgeom=False,unload=False,increment=.1):
    cells,order=ccx_order(m['cells']);ids=m['node_ids'];anchor=int(bc['axial_gauge_original_index']);normals=bc['radial_normal'].astype(float)
    supported=np.flatnonzero(np.linalg.norm(normals,axis=1)>.5)
    assert anchor in supported and np.isfinite(force).all() and force.shape==m['xyz'].shape
    with Path(file).open('w',encoding='ascii') as f:
        f.write('*HEADING\nOriginal-node rod pilot, radial bilateral pin, reference loads\n*NODE,NSET=NALL\n')
        for nid,x in zip(ids,m['xyz']):f.write(f'{int(nid)},'+','.join(f'{v:.14g}' for v in x)+'\n')
        f.write('*ELEMENT,TYPE=C3D10,ELSET=EALL\n')
        for eid,c in zip(m['cell_ids'],cells):f.write(f'{int(eid)},'+','.join(str(int(v)) for v in ids[c])+'\n')
        f.write('*MATERIAL,NAME=REFERENCE\n*ELASTIC\n'+f"{material['E_MPa']:.12g},{material['nu']:.12g}\n")
        if plastic:
            f.write('*PLASTIC\n')
            for stress,eps in material['plastic_true_stress_vs_equivalent_plastic_strain']:f.write(f'{stress:.12g},{eps:.12g}\n')
        f.write('*SOLID SECTION,ELSET=EALL,MATERIAL=REFERENCE\n*BOUNDARY\n')
        f.write(f'{int(ids[anchor])},1,3,0\n')
        equations=[]
        for i in supported:
            if i==anchor:continue
            n=normals[i,:2];order_dof=np.argsort(-np.abs(n));active=[int(j) for j in order_dof if abs(n[j])>1e-9]
            if len(active)==1:f.write(f'{int(ids[i])},{active[0]+1},{active[0]+1},0\n')
            else:equations.append((int(ids[i]),n,active))
        for nid,n,active in equations:
            f.write('*EQUATION\n2\n'+','.join(f'{nid},{j+1},{n[j]:.14g}' for j in active)+'\n')
        f.write('*STEP,'+('NLGEOM,' if nlgeom else '')+'INC=500\n*STATIC\n'+f'{increment},1,1e-7,{increment}\n*CLOAD\n')
        for nid,force_row in zip(ids,force):
            for j,v in enumerate(force_row):
                if abs(v)>1e-14:f.write(f'{int(nid)},{j+1},{v:.14g}\n')
        output='*NODE PRINT,NSET=NALL,FREQUENCY=999999\nU,RF\n*EL PRINT,ELSET=EALL,FREQUENCY=999999\nS'+(',PEEQ' if plastic else '')+'\n*NODE FILE,FREQUENCY=999999\nU\n*EL FILE,FREQUENCY=999999\nS'+(',PEEQ' if plastic else '')+'\n*END STEP\n'
        f.write(output)
        if unload:
            f.write('*STEP,'+('NLGEOM,' if nlgeom else '')+'INC=500\n*STATIC\n'+f'{increment},1,1e-7,{increment}\n*CLOAD,OP=NEW\n');f.write(output)
    return dict(gmsh_to_CCX_order=order,anchor_original_index=anchor,anchor_node_id=int(ids[anchor]),support_nodes=len(supported),
        constraint='Same reference-normal u_r=0 linear equations and tangential/z anchor; no extra clamp/contact. With NLGEOM radial equations stay reference-linear, not exact finite sliding cylinder contact.',
        force='Consistent assembled nodal surface+body reference loads, dead-load directions in NLGEOM; no follower traction or explicit transient dynamics claim.')

def run_ccx(executable,deck,timeout=3600,threads=4):
    start=time.perf_counter();env=dict(os.environ,OMP_NUM_THREADS=str(threads),CCX_NPROC_RESULTS=str(threads),CCX_NPROC_STIFFNESS=str(threads),NUMBER_OF_CPUS=str(threads))
    with (deck.parent/'driver.log').open('w',encoding='utf-8') as log:
        child=subprocess.Popen([str(executable),'-i',deck.stem],cwd=deck.parent,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        write_json(deck.parent/'process.json',dict(pid=child.pid,started_unix=time.time(),deck=str(deck)))
        try:code=child.wait(timeout=timeout)
        except subprocess.TimeoutExpired:child.kill();child.wait();raise RuntimeError('CCX isolated pilot timeout; preserve incomplete result, not a physical failure label')
    seconds=time.perf_counter()-start;text=(deck.parent/'driver.log').read_text(errors='replace')
    if code or 'Job finished' not in text or '*ERROR' in text:raise RuntimeError(f'CCX failed exit{code}; inspect driver.log')
    return seconds

def export_result(deck,m,force,material,metadata,seconds):
    blocks=read_dat(deck.with_suffix('.dat'));times=sorted({x['time'] for x in blocks if x['kind']=='U'})
    assert times and abs(times[-1]-(2 if metadata.get('unload') else 1))<1e-6
    def at(kind,t):
        match=[x['rows'] for x in blocks if x['kind']==kind and abs(x['time']-t)<1e-7];assert match,(kind,t);return match[-1]
    ids=m['node_ids'];native={int(n):i for i,n in enumerate(ids)};field=[]
    for t in times:
        a=at('U',t);u=np.full((len(ids),3),np.nan)
        for row in a:u[native[int(row[0])]]=row[1:]
        assert np.isfinite(u).all() and len(a)==len(ids)
        field.append((t,u))
    peak=next((u for t,u in field if abs(t-1)<1e-7),field[-1][1]);s=at('S',1.)
    gpids=s[:,:2].astype('int64');stress=s[:,2:];von=np.sqrt(.5*((stress[:,0]-stress[:,1])**2+(stress[:,1]-stress[:,2])**2+(stress[:,2]-stress[:,0])**2)+3*np.square(stress[:,3:]).sum(1))
    if metadata.get('plastic'):
        pe=at('PEEQ',1.);assert np.array_equal(pe[:,:2].astype('int64'),gpids);plastic=pe[:,2]
    else:plastic=np.zeros(len(s))
    assert np.isfinite(stress).all() and np.isfinite(plastic).all() and plastic.min()>-1e-10
    result=dict(node_ids=ids,U_peak_mm=peak,U_final_mm=field[-1][1],force_peak_N=force,GP_element_ids=gpids[:,0],GP_indices=gpids[:,1],GP_stress_MPa=stress,GP_von_mises_MPa=von,GP_equivalent_plastic_strain=plastic)
    np.savez_compressed(deck.parent/'fields.npz',**result)
    report=dict(**metadata,completed=True,solver='CalculiX2.23',seconds=seconds,nodes=len(ids),elements=len(m['cells']),max_U_mm=float(np.linalg.norm(peak,axis=1).max()),
        max_residual_U_mm=float(np.linalg.norm(field[-1][1],axis=1).max()) if metadata.get('unload') else None,
        max_GP_von_Mises_MPa=float(von.max()),max_PEEQ=float(plastic.max()),plastic_GP_fraction=float(np.mean(plastic>1e-8)),
        material_table_range_exceeded=bool(plastic.max()>material.get('plastic_strain_data_limit',float('inf'))),
        stress_association='Actual integration point tensors retained; no nodal stress averaging; raw peak stress not converged',
        material=material,mesh_increment_convergence='PENDING; pilot solve completion is not validated dataset acceptance',native_node_identity_verified=True)
    write_json(deck.parent/'summary.json',report);return report
