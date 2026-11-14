"""A new 5um pair using the frozen, qualified same-candidate reference runner."""
from pathlib import Path
import argparse, importlib.util, json, subprocess, tomllib

HERE=Path(__file__).resolve().parent;R=HERE.parents[2]
spec=importlib.util.spec_from_file_location('frozen_retained_runner',HERE.parent/'retained_reference/run.py')
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
O=R/'results/retained_five_20261113';RAW=R/'results/local-runs/retained_five_20261113';LOG=R/'results/logs/retained_five_20261113'
REF=R/'results/retained_reference_20261113'
base.O,base.RAW,base.LOG=O,RAW,LOG
sha,read,save,now,progress=base.sha,base.read,base.save,base.now,base.progress
THREAD='01a1c8ee-abdc-7c71-a299-31f9cb5ed265'
CONFIG=Path.home()/'.codex/automations/drift-5/automation.toml'

def references():
    p=read(REF/'protocol.json');s=read(REF/'status.json');a=read(REF/'audit_coverage.json');c=read(REF/'convergence_summary.json')
    assert s['state']=='completed' and not s['skipped'] and len(s['cases'])==4
    assert a['all4_native_cases_independently_audited'] and a['scalar_saved_steps']==a['mechanical_saved_steps']==36
    assert c['first_checkpoint_time_pairs_passed'] and not c['full_C1_complete']
    assert all(sha(R/f)==h for f,h in p['source_hashes'].items())
    seals={}
    for case in s['cases']:
        assert case['completed'] and case['native_gate'] and case['returncode']==0
        assert all(sha(R/f['path'])==f['sha256'] for f in case['files'])
        assert sha(R/case['log'])==case['log_sha256'] and sha(R/case['monitor_log'])==case['monitor_log_sha256']
        d=REF/'audits/cases'/case['name'];v=read(d/'verified.json');w=read(d/'wrapper_verified.json')
        assert w['passed'] and w['scalar_saved_steps']==w['mechanical_saved_steps']==9
        assert w['inherited_seal_sha256']==sha(d/'verified.json')
        assert all(sha(R/f)==h for f,h in w['sources'].items())
        assert all(sha(d/f)==h for f,h in v['outputs'].items())
        seals[case['name']]=sha(d/'verified.json')
    return p,dict(reference_protocol_sha256=sha(REF/'protocol.json'),reference_status_sha256=sha(REF/'status.json'),audit_coverage_sha256=sha(REF/'audit_coverage.json'),comparison_sha256=sha(REF/'comparisons.json'),convergence_summary_sha256=sha(REF/'convergence_summary.json'),audit_seals=seals)

def prepare():
    assert not O.exists() and not RAW.exists() and not LOG.exists()
    p,evidence=references();exe=Path(p['executable']);assert sha(exe)==p['binary_sha256']
    assert p['binary_sha256']=='fae0ce667a1456afd81850fe2c833fb95287b291c0319589d355df37c4e34c2c'
    assert all(sha(Path(f))==h for f,h in p['runtime_dlls'].items())
    files=[Path(__file__),HERE/'PLAN.md',HERE/'monitor.py']
    cases=[dict(name='u5_dt'+tag,mode='uniform',dx_um=5,N=2112,dt_us=dt,T_us=10,radius_um=500,amplitude=0,c_bulk_mM=6,preparation='variational',save_stride=stride,refine_start_um=200,observer=True,timeout_s=cap) for tag,dt,stride,cap in [('00625',.0625,20,18000),('003125',.03125,40,28000)]]
    O.mkdir();RAW.mkdir();LOG.mkdir()
    save(O/'protocol.json',dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),prepared_utc=now(),source_hashes={**p['source_hashes'],**{str(f.relative_to(R)):sha(f) for f in files}},binary_sha256=p['binary_sha256'],executable=str(exe),runtime_dlls=p['runtime_dlls'],cases=cases,whole_native_cap_s=46000,monitor_interval_s=1800,gates=p['gates'],reference_evidence=evidence,environment=dict(VF_THREADS='2',caller_VF='excluded',runtime_PATH=base.ENV['PATH']),required_monitor=dict(id='drift-5',status='ACTIVE',rrule='FREQ=MINUTELY;INTERVAL=30',target_thread_id=THREAD,campaign='retained_five_20261113'),scope='NEW same-candidate5um10us pair; original failed trajectory not restarted; fullC1 incomplete'))
    save(O/'status.json',dict(state='prepared',cases=[],skipped=[],full_C1_complete=False))
    print('Prepared NEW5um pair; no native process launched.')

def monitor_config():
    cfg=tomllib.loads(CONFIG.read_text(encoding='utf-8'));required=read(O/'protocol.json')['required_monitor']
    for k in ['id','status','rrule','target_thread_id']:assert cfg[k]==required[k],k
    assert required['campaign'] in cfg['prompt'] and 'retained_five' in cfg['prompt']
    return cfg

def authorize():
    assert read(O/'status.json')['state']=='prepared' and not (O/'monitor_authorization.json').exists()
    cfg=monitor_config()
    save(O/'monitor_authorization.json',dict(verified_utc=now(),protocol_sha256=sha(O/'protocol.json'),configuration_sha256=sha(CONFIG),configuration_path=str(CONFIG),id=cfg['id'],name=cfg['name'],status=cfg['status'],rrule=cfg['rrule'],target_thread_id=cfg['target_thread_id'],prompt=cfg['prompt']))
    print('Verified scheduled30min monitor ACTIVE for this new campaign; no native launch.')

def run():
    assert read(O/'status.json')['state']=='prepared'
    p=read(O/'protocol.json');seal=read(O/'monitor_authorization.json');cfg=monitor_config()
    assert seal['protocol_sha256']==sha(O/'protocol.json') and seal['configuration_sha256']==sha(CONFIG)
    assert seal['prompt']==cfg['prompt'] and references()[1]==p['reference_evidence']
    save(O/'launch_record.json',dict(utc=now(),execution_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),parent_pid=base.os.getpid(),protocol_sha256=sha(O/'protocol.json'),monitor_authorization_sha256=sha(O/'monitor_authorization.json'),binary_sha256=p['binary_sha256'],native_started=False,scope='Pre-Popen admission evidence; actual live PID must be checked separately'))
    base.run()

if __name__=='__main__':
    a=argparse.ArgumentParser();g=a.add_mutually_exclusive_group(required=True)
    for flag in ['prepare','authorize','run']:g.add_argument('--'+flag,action='store_true')
    args=a.parse_args();prepare() if args.prepare else authorize() if args.authorize else run()
