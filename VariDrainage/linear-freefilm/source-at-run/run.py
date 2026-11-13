"""Freeze, run and preserve falsifiable linear free-film controls."""
from pathlib import Path
import json,hashlib,datetime,subprocess,time,argparse
import numpy as np
from scipy.linalg import eigvals,eigvalsh
from reference import build,trajectory
R=Path(__file__).resolve().parents[3];O=R/'results/linear_freefilm_20261113';HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False),encoding='utf-8')
def prepare():
    assert not (O/'protocol.json').exists();s=json.loads((R/'results/global_space_20261113/status.json').read_text());rec=next(c for c in s['cases'] if c['name']=='u10_dt003125');assert rec['native_gate'] and rec['completed']
    for f in rec['files']:assert sha(R/f['path'])==f['sha256']
    initial=next(R/f['path'] for f in rec['files'] if Path(f['path']).name=='initial_states.csv');a=np.genfromtxt(initial,delimiter=',',names=True);i=int(np.argmin((a['x_m']-80e-6)**2+(a['y_m']-80e-6)**2));m=rec['summary'];theta=float(a['theta'][i]);E=m['RT_J_mol']*m['Gamma_inf_mol_m2']*theta/(1-theta)
    save(O/'protocol.json',dict(prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),sources={str(f.relative_to(R)):sha(f) for f in [Path(__file__),HERE/'reference.py',HERE/'PLAN.md']},initial_path=str(initial.relative_to(R)),initial_sha256=sha(initial),node_index=i,theta=theta,rho=m['rho_kg_m3'],mu=m['viscosity_Pa_s'],gamma=float(a['tension_N_m'][i]),E=E,h=6e-6,wavelengths_um=[40,80,160],orders=[8,16,32],times_us=[0,.5,1,2,5,10,20,50,100,200],g0=1e-6,budget_s=300,scope='Linear symmetric free-film reference, no curved edge/bath/native coupling or experiment validation'))
    print('Prepared linear free-film reference, no evolution run.')
def main():
    p=json.loads((O/'protocol.json').read_text());assert not (O/'summary.json').exists();assert all(sha(R/f)==h for f,h in p['sources'].items());assert sha(R/p['initial_path'])==p['initial_sha256'];start=time.monotonic();controls=[]
    for wavelength in [40,80]:
        k=2*np.pi/(wavelength*1e-6);exact=np.sqrt(p['gamma']/p['rho']*k**3*np.tanh(k*p['h']/2))
        for N in p['orders']:
            model=build(N,p['h'],p['rho'],0,p['gamma'],0,k);omega=max(abs(eigvals(model['A']).imag));controls.append(dict(kind='inviscid_capillary',wavelength_um=wavelength,N=N,relative_error=float(abs(omega/exact-1)),passed=bool(abs(omega/exact-1)<1e-8)))
    shear=build(32,p['h'],p['rho'],p['mu'],p['gamma'],0,0);rates=eigvalsh(shear['L']);exact=p['mu']/p['rho']*(2*np.pi*np.arange(1,4)/p['h'])**2;err=float(max(abs(rates[1:4]/exact-1)));controls.append(dict(kind='zero_k_shear',relative_error=err,zero_mode_normalized=float(abs(rates[0])/exact[0]),passed=err<1e-6 and abs(rates[0])/exact[0]<1e-10))
    allrows={};curves={};comparisons=[];times=np.array(p['times_us'])*1e-6
    for wavelength in p['wavelengths_um']:
        k=2*np.pi/(wavelength*1e-6)
        for N in p['orders']:
            model=build(N,p['h'],p['rho'],p['mu'],p['gamma'],p['E'],k);A=model['A'];identity=A+A.conj().T;identity[:N,:N]+=2*model['L'];identity_error=float(np.linalg.norm(identity)/np.linalg.norm(A));psd=float(min(eigvalsh(model['L']))/np.linalg.norm(model['L']));initial,states,rows=trajectory(model,times,p['g0']);key=f'lambda{wavelength}_N{N}';allrows[key]=rows
            np.savez_compressed(O/(key+'.npz'),times_s=times,energy_states=states,initial=initial,generator=A,mass=model['M'],dissipation=model['D'])
            tang=np.linalg.norm([r['tangential_residual'] for r in rows[1:]])/max(1e-60,np.linalg.norm([r['tangential_scale'] for r in rows[1:]]));normal=np.linalg.norm([r['normal_residual'] for r in rows[1:]])/max(1e-60,np.linalg.norm([r['normal_scale'] for r in rows[1:]]));E0=rows[0]['energy'];increase=max(0,max(np.diff([r['energy'] for r in rows])))/E0
            controls.append(dict(kind='coupled_freefilm',wavelength_um=wavelength,N=N,energy_identity_error=identity_error,normalized_minimum_damping_eigenvalue=psd,normalized_maximum_energy_increase=float(increase),tangential_stress_RMS_error=float(tang),normal_stress_RMS_error=float(normal),passed=identity_error<1e-10 and psd>=-1e-12 and increase<1e-8 and (N!=32 or max(tang,normal)<.01)))
            curves[key]=[dict(t_us=t,**{name:[float(r[name].real),float(r[name].imag)] for name in ['eta','g','u_surface','u_mean','flux']},energy_J_m2=r['energy'],dissipation_W_m2=r['dissipation']) for t,r in zip(p['times_us'],rows)]
        errors=[]
        for a,b in [(8,16),(16,32)]:
            x,y=allrows[f'lambda{wavelength}_N{a}'],allrows[f'lambda{wavelength}_N{b}'];delta={name:float(np.linalg.norm([r[name]-s[name] for r,s in zip(x,y)])/max(1e-60,np.linalg.norm([r[name] for r in y]))) for name in ['eta','g','u_surface','u_mean']};errors.append(max(delta.values()));comparisons.append(dict(wavelength_um=wavelength,N_candidate=a,N_reference=b,relative_L2=delta))
        controls.append(dict(kind='order_refinement',wavelength_um=wavelength,coarser_error=errors[0],finer_error=errors[1],passed=errors[1]<.01 and errors[1]<errors[0]))
    save(O/'controls.json',controls);save(O/'curves.json',curves);save(O/'comparisons.json',comparisons);save(O/'summary.json',dict(all_reference_controls_passed=all(c['passed'] for c in controls),controls=len(controls),failed_controls=[c for c in controls if not c['passed']],wall_s=time.monotonic()-start,protocol_sha256=sha(O/'protocol.json'),full_C1_complete=False,actual_edge_physical_flux_validated=False,scope=p['scope']));print(json.dumps(json.loads((O/'summary.json').read_text()),indent=2))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--prepare',action='store_true');a=parser.parse_args();prepare() if a.prepare else main()
