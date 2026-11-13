"""Presentation-only structural patch, preserving the preregistered renderer."""
from pathlib import Path
import json, hashlib, shutil

R=Path(__file__).resolve().parents[3]
D=R/'results/uniform5_cut_20261113'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
p=D/'index.html'
assert not (D/'presentation_patch.json').exists()
source=p.read_text(encoding='utf-8')
assert source.startswith('<!doctype html><html lang="zh-CN"><meta')
assert source.endswith('</main></html>')
before=digest(p)
changed=source.replace('<html lang="zh-CN"><meta','<html lang="zh-CN"><head><meta',1).replace('</style><main>','</style></head><body><main>',1).replace('</main></html>','</main></body></html>')
assert source.count('</style><main>')==1
oldfig={f.name:digest(f) for f in D.glob('cut-flux.*')}
p.write_text(changed,encoding='utf-8')
assert oldfig=={f.name:digest(f) for f in D.glob('cut-flux.*')}
shutil.copy2(__file__,D/'cut_markup.py')
(D/'presentation_patch.json').write_text(json.dumps(dict(original_html_sha256=before,patched_html_sha256=digest(p),patch_source_sha256=digest(Path(__file__)),figures_unchanged=oldfig,scope='Explicit head/body prevents deployment resource router preceding doctype; no numerical rerun, original renderer/protocol retained'),indent=2)+'\n',encoding='utf-8')
