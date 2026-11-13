from pathlib import Path
import sys,json,zipfile,hashlib
root=Path(sys.argv[1]);out=Path(sys.argv[2]).resolve();assert not out.exists()
m=json.loads((root/'download_manifest.json').read_text());p=m['parts'][0];data=(root/p['path']).read_bytes();assert hashlib.sha256(data).hexdigest()==p['sha256']
with zipfile.ZipFile(root/p['path']) as z:
 for f in m['files']:
  target=(out/f['path']).resolve();assert target.is_relative_to(out);b=z.read(f['path']);assert len(b)==f['bytes'] and hashlib.sha256(b).hexdigest()==f['sha256']
 for f in m['files']:
  target=(out/f['path']).resolve();target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(z.read(f['path']))
print('Verified restored files',len(m['files']))
