"""Bounded deduplicated ZIP chunks, all native attempts including failures."""
from pathlib import Path
import json,hashlib,zipfile,shutil
from run import R,O,LOG,save,sha
def main():
    roots=[('native','prepared_motion_20261113',O/'status.json')]
    assert not (O/'download_manifest.json').exists()
    files=[];chunks={}
    for namespace,directory,statefile in roots:
        state=json.loads(statefile.read_text());records=state.get('cases',[state]);root=R/'results/local-runs'/directory
        for record in records:
            for f in record['files']:
                p=R/f['path'];assert sha(p)==f['sha256'];data=p.read_bytes();item=dict(path=namespace+'/'+p.relative_to(root).as_posix(),sha256=f['sha256'],bytes=len(data),chunks=[])
                for offset in range(0,len(data),2_000_000):
                    block=data[offset:offset+2_000_000];digest=hashlib.sha256(block).hexdigest();entry='chunks/'+digest
                    if digest not in chunks:chunks[digest]=dict(entry=entry,data=block,sha256=digest,bytes=len(block))
                    item['chunks'].append(dict(entry=entry,sha256=digest,bytes=len(block),offset=offset))
                files.append(item)
    dest=O/'downloads';dest.mkdir(exist_ok=True);parts=[];batch=[];size=0
    def flush():
        nonlocal batch,size
        if not batch:return
        p=dest/f'raw-part-{len(parts)+1:03d}.zip'
        with zipfile.ZipFile(p,'w',zipfile.ZIP_DEFLATED,compresslevel=9) as z:
            for c in batch:z.writestr(c['entry'],c['data']);c['archive']=p.name
        assert p.stat().st_size<5_000_000;parts.append(dict(path=p.relative_to(O).as_posix(),bytes=p.stat().st_size,sha256=sha(p),entries=[c['entry'] for c in batch]));batch=[];size=0
    for c in chunks.values():
        if size+c['bytes']>4_000_000:flush()
        batch.append(c);size+=c['bytes']
    flush()
    for f in files:
        digest=hashlib.sha256()
        for c in f['chunks']:
            c['archive']=chunks[c['sha256']]['archive']
            with zipfile.ZipFile(dest/c['archive']) as z:data=z.read(c['entry'])
            assert hashlib.sha256(data).hexdigest()==c['sha256'];digest.update(data)
        assert digest.hexdigest()==f['sha256']
    save(O/'download_manifest.json',dict(parts=parts,files=files,file_count=len(files),raw_bytes=sum(f['bytes'] for f in files),unique_chunk_bytes=sum(c['bytes'] for c in chunks.values()),deduplication='SHA256-identical chunks shared across original/replay files; restore preserves every path and byte'))
    shutil.copytree(LOG,O/'logs',dirs_exist_ok=True);shutil.copy2(Path(__file__).with_name('README.md'),O/'METHOD.md');shutil.copy2(R/'results/edge_preparation_20261113/restore_raw.py',O/'restore_raw.py')
    (O/'logs/index.html').write_text('<!doctype html><meta charset="utf-8"><title>Physical diagnostics logs</title><h1>Native logs</h1>'+''.join(f'<p><a href="{p.relative_to(O/"logs").as_posix()}">{p.relative_to(O/"logs").as_posix()}</a></p>' for p in sorted((O/'logs').rglob('*.log'))),encoding='utf-8')
    print('Verified files',len(files),'parts',len(parts),'rawMB',sum(f['bytes'] for f in files)/1e6,'uniqueMB',sum(c['bytes'] for c in chunks.values())/1e6)
if __name__=='__main__':main()
