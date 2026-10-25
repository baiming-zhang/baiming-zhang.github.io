"""Rotated orthotropic C3D10 traction patch, analytical and CalculiX evidence."""
from pathlib import Path
import argparse,json,subprocess,time
import numpy as np
import gmsh
from geobc.copv_material import engineering_stiffness,rotate_stiffness,stress_to_tensor
from geobc.rod_v3_ccx import read_dat

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--ccx',type=Path,required=True);args=ap.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True);cfg=json.loads((Path(__file__).resolve().parents[1]/'configs/copv_c02_structural_v1.json').read_text());m=cfg['materials']['carbon'];C=engineering_stiffness(**{k:m[k] for k in ['E1','E2','E3','G12','G13','G23','nu12','nu13','nu23']});az=np.deg2rad(37);el=np.deg2rad(22);Q=np.array([[np.cos(az),-np.sin(az),0],[np.sin(az),np.cos(az),0],[0,0,1]])@np.array([[np.cos(el),0,np.sin(el)],[0,1,0],[-np.sin(el),0,np.cos(el)]])
    B=np.array([[.001,.0002,-.0001],[.0002,-.0002,.0003],[-.0001,.0003,.0004]]);strain=np.array([B[0,0],B[1,1],B[2,2],2*B[1,2],2*B[0,2],2*B[0,1]]);Cg=rotate_stiffness(C,Q);stress=Cg@strain;sigma=stress_to_tensor(stress)
    corners=np.array([[0,0,0],[10,0,0],[0,10,0],[0,0,10]],float);edges=[(0,1),(1,2),(2,0),(0,3),(2,3),(1,3)];xyz=np.vstack([corners,[(corners[i]+corners[j])/2 for i,j in edges]]);U=xyz@B.T;F=np.zeros_like(xyz);faces=np.array([[1,2,3,5,8,9],[0,3,2,7,8,6],[0,1,3,4,9,7],[0,2,1,6,5,4]])
    gmsh.initialize()
    try:
        q,w=gmsh.model.mesh.getIntegrationPoints(9,'Gauss4');w=np.array(w);N=np.array(gmsh.model.mesh.getBasisFunctions(9,q,'Lagrange')[1]).reshape(-1,6);dN=np.array(gmsh.model.mesh.getBasisFunctions(9,q,'GradLagrange')[1]).reshape(-1,6,3)
    finally:gmsh.finalize()
    for opposite,face in enumerate(faces):
        xx=xyz[face];d=np.einsum('qij,ik->qjk',dN,xx);area=np.cross(d[:,0],d[:,1]);center=N@xx
        if np.mean(np.sum(area*(corners[opposite]-center),axis=1))>0:area=-area
        traction=area@sigma.T;F[face]+=np.einsum('qi,qj,q->ij',N,traction,w)
    np.savez_compressed(out/'input.npz',xyz=xyz,cells=np.arange(10).reshape(1,10),C_local=C,Q=Q,B=B,force_N=F,expected_U=U,expected_stress_xyz=stress)
    deck=['*HEADING','C02 rotated orthotropic 3D traction patch','*NODE,NSET=NALL']+[str(i+1)+','+','.join(f'{v:.16g}' for v in x) for i,x in enumerate(xyz)]
    deck+=['*ELEMENT,TYPE=C3D10,ELSET=EALL','1,1,2,3,4,5,6,7,8,10,9','*ORIENTATION,NAME=FIBER,SYSTEM=RECTANGULAR',','.join(f'{v:.16g}' for v in np.r_[Q[:,0],Q[:,1]]),'3,0','*MATERIAL,NAME=PLY','*ELASTIC,TYPE=ENGINEERING CONSTANTS',','.join(str(m[k]) for k in ['E1','E2','E3','nu12','nu13','nu23','G12','G13']),str(m['G23']),'*SOLID SECTION,ELSET=EALL,MATERIAL=PLY,ORIENTATION=FIBER','*STEP','*STATIC','1.,1.','*BOUNDARY']
    for i in range(4):
        for j in range(3):deck.append(f'{i+1},{j+1},{j+1},{U[i,j]:.16g}')
    deck+=['*CLOAD']
    for i in range(10):
        for j in range(3):
            if abs(F[i,j])>1e-10:deck.append(f'{i+1},{j+1},{F[i,j]:.12g}')
    deck+=['*NODE PRINT,NSET=NALL','U,RF','*EL PRINT,ELSET=EALL,GLOBAL=YES','S','*END STEP'];(out/'patch.inp').write_text('\n'.join(deck)+'\n',encoding='ascii');tick=time.perf_counter()
    with (out/'ccx.log').open('w') as f:r=subprocess.run([str(args.ccx.resolve()),'-i','patch'],cwd=out,stdout=f,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    if r.returncode:raise RuntimeError('CalculiX patch failed; inspect ccx.log')
    blocks=read_dat(out/'patch.dat');uu=next(b['rows'] for b in blocks if b['kind']=='U');ss=next(b['rows'] for b in blocks if b['kind']=='S');pred=uu[np.argsort(uu[:,0]),1:];stress_global=ss[:,2:][:,[0,1,2,5,4,3]];report={'U_relative_L2':float(np.linalg.norm(pred-U)/np.linalg.norm(U)),'stress_relative_L2':float(np.linalg.norm(stress_global-stress)/np.linalg.norm(np.broadcast_to(stress,stress_global.shape))),'seconds':time.perf_counter()-tick,'constant_rotated_orthotropic_stress_MPa':stress.tolist(),'load_force_balance_N':F.sum(0).tolist(),'load_moment_balance_Nmm':np.cross(xyz,F).sum(0).tolist()};report['passed']=report['U_relative_L2']<1e-5 and report['stress_relative_L2']<1e-5;(out/'calculix_review.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    if not report['passed']:raise RuntimeError('CalculiX analytical patch gate failed')
if __name__=='__main__':main()
