"""Original bounded chunk packager over terminal native, audit, log and source bytes."""
from pathlib import Path
import json,shutil,subprocess,sys
from run import R,O,RAW,LOG,read,save,sha
HERE=Path(__file__).resolve().parent
def main():
    eligible();state=read(O/'status.json');assert not (O/'download_manifest.json').exists()
    p=read(O/'protocol.json');assert all(sha(R/f)==h for f,h in p['source_hashes'].items()) and sha(Path(p['executable']))==p['binary_sha256']
    for rec in state['cases']:
        assert all(sha(R/f['path'])==f['sha256'] for f in rec['files'])
        for key in ['log','monitor_log']:assert sha(R/rec[key])==rec[key+'_sha256']
    snapshots={}
    for category in ['logs','audits','sources']:
        d=RAW/('_archive_'+category);assert not d.exists();d.mkdir();snapshots[category]=d
    for f in LOG.iterdir():shutil.copy2(f,snapshots['logs']/f.name)
    for f in (O/'audits').rglob('*'):
        if f.is_file():dest=snapshots['audits']/f.relative_to(O/'audits');dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,dest)
    source_map=dict(p['source_hashes'])
    for proto in [O/'postprocess_protocol.json',O/'comparison_protocol.json',O/'archive_report_protocol.json',R/'results/retained_checks_20261113/protocol.json']:source_map.update(read(proto)['sources'])
    pack=HERE.parent/'prepared_motion/package.py';restore=R/'results/edge_preparation_20261113/restore_raw.py'
    for f in [Path(__file__),HERE/'report.py',pack,restore]:source_map[str(f.relative_to(R))]=sha(f)
    for name,h in source_map.items():
        f=R/name;assert sha(f)==h;dest=snapshots['sources']/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,dest)
    shutil.copytree(snapshots['sources'],O/'source-at-run')
    roots=[('native','retained_five_20261113',O/'status.json')]
    for category,d in snapshots.items():
        inv=O/(category+'_inventory.json');save(inv,dict(files=[dict(path=str(f.relative_to(R)),bytes=f.stat().st_size,sha256=sha(f)) for f in sorted(d.rglob('*')) if f.is_file()]));roots.append((category,str(d.relative_to(R/'results/local-runs').as_posix()),inv))
    code=pack.read_text(encoding='utf-8');old="roots=[('native','prepared_motion_20261113',O/'status.json')]";assert code.count(old)==1
    expression='roots=['+','.join('('+repr(a)+','+repr(b)+',O/'+repr(c.name)+')' for a,b,c in roots)+']';code=code.replace(old,expression)
    old="shutil.copy2(Path(__file__).with_name('README.md'),O/'METHOD.md')";assert code.count(old)==1;code=code.replace(old,"shutil.copy2(Path(__file__).with_name('PLAN.md'),O/'METHOD.md')")
    ns=dict(__file__=str(Path(__file__)),__name__='retained_five_bounded_packager');exec(compile(code,str(pack),'exec'),ns);ns['main']()
    target=(R/'results/local-runs/retained_five_restore_20261113').resolve();assert target.is_relative_to((R/'results/local-runs').resolve()) and not target.exists()
    cmd=[sys.executable,str(O/'restore_raw.py'),str(O),str(target)];done=subprocess.run(cmd,cwd=R,capture_output=True,text=True,timeout=120);save(O/'restore_execution.json',dict(command=cmd,returncode=done.returncode,stdout=done.stdout,stderr=done.stderr));assert done.returncode==0
    manifest=read(O/'download_manifest.json');assert all(sha(target/f['path'])==f['sha256'] for f in manifest['files'])
    save(O/'archive_verified.json',dict(files=manifest['file_count'],raw_bytes=manifest['raw_bytes'],parts=len(manifest['parts']),all_restored_byte_exact=True,manifest_sha256=sha(O/'download_manifest.json'),restore_directory=str(target.relative_to(R)),full_C1_complete=False))
    print('Restored',manifest['file_count'],'native/log/audit/source files exactly.')
def eligible():
    state=read(O/'status.json');assert state['state']=='completed','Native matrix is live; no archive or live-frame reads'
    live=read(O/'live.json');assert live['state']=='terminal_batch'
    assert len(state['cases'])+len(state['skipped'])==2 and state['cases']
    p=read(O/'archive_report_protocol.json')
    assert p['native_protocol_sha256']==sha(O/'protocol.json') and p['sources']==sources()
    # A fully successful native pair must have completed independent audits and
    # original comparisons before the final archive. Failure evidence can be
    # preserved without manufacturing an audit/comparison for the failed case.
    if len(state['cases'])==2 and all(c['native_gate'] for c in state['cases']):
        import compare_terminal
        compare_terminal.eligible();assert (O/'comparisons.json').exists()
    return state

def sources():
    files=[Path(__file__),HERE/'report.py',HERE/'ARCHIVE_REPORT.md',HERE.parent/'prepared_motion/package.py',R/'results/edge_preparation_20261113/restore_raw.py']
    return {str(f.relative_to(R)):sha(f) for f in files}

def prepare():
    import datetime
    p=O/'archive_report_protocol.json';assert not p.exists() and not (O/'download_manifest.json').exists() and not (O/'index.html').exists()
    save(p,dict(prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),native_state_at_prepare=read(O/'status.json')['state'],sources=sources(),native_protocol_sha256=sha(O/'protocol.json'),scope='Terminal-only bounded evidence archive and actual-data report; no new native work',full_C1_complete=False))
    print('Prepared archive/report contract; no native frames read or result figures generated.')

if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser();g=a.add_mutually_exclusive_group(required=True)
    for name in ['prepare','readiness','run']:g.add_argument('--'+name,action='store_true')
    args=a.parse_args()
    if args.prepare:prepare()
    elif args.run:main()
    else:
        try:eligible()
        except (AssertionError,FileNotFoundError) as e:print(json.dumps(dict(ready=False,reason=str(e),archive_exists=(O/'download_manifest.json').exists(),new_native_executions=0)))
        else:print(json.dumps(dict(ready=True,new_native_executions=0)))

