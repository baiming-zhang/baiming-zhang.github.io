from pathlib import Path,PurePosixPath
import sys,json,hashlib,zipfile
source=Path(sys.argv[1]);target=Path(sys.argv[2]);m=json.loads((source/'download_manifest.json').read_text())
for p in m['parts']:assert hashlib.sha256((source/p['path']).read_bytes()).hexdigest()==p['sha256']
for f in m['files']:
 r=PurePosixPath(f['path']);assert not r.is_absolute() and '..' not in r.parts;p=target/Path(*r.parts);assert not p.exists()
 with zipfile.ZipFile(source/'downloads'/f['archive']) as z:data=z.read(f['path'])
 assert len(data)==f['bytes'] and hashlib.sha256(data).hexdigest()==f['sha256'];p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
print('Restored and verified',len(m['files']),'files')
