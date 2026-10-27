"""Stage the completed q8 alongside immutable q6; do not publish before QA."""
from pathlib import Path
import hashlib,json,shutil,re
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3';S=W/'showcase-public-site/GeoBC-show';O=B/'curved_viewer_q8_v1'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert read(B/'gpu_validation_queue.json')['status']=='completed_candidate_validation'
payload=O/'assets/copv-c04-repair-fields/copv_c04_0000-q8';export=read(payload/'export_review.json');manifest=read(payload/'manifest.json');assert export['actual_native_loads_checked'] and export['locked_components']==6 and export['pressure_arrows']==28 and manifest['summary']['assembly_quadrature_degree']==8
proof=O/'preview_prepared_v1.json';assert not proof.exists(),'Keep preview preparation versions'
for name in ['access.css','access.js','style.css','fe-glyphs-3d.js','copv-c04-curved-fields.js','copv-c04-curved-bc.js']:
 dest=O/'assets'/name;assert not dest.exists();shutil.copyfile(S/'assets'/name,dest)
shutil.copytree(S/'assets/copv-c04-repair-fields/copv_c04_0000-q6',O/'assets/copv-c04-repair-fields/copv_c04_0000-q6')
page=(S/'viewers/copv-c04-curved.html').read_text(encoding='utf-8');page,n=re.subn(r'<option value="8" disabled>[^<]*</option>','<option value="8">8阶</option>',page);assert n==1
(O/'viewers').mkdir();(O/'viewers/copv-c04-curved.html').write_text(page,encoding='utf-8');proof.write_text(json.dumps(dict(status='prepared_q8_pending_browser_QA',q6_original_site_hashes={p.name:sha(p) for p in (S/'assets/copv-c04-repair-fields/copv_c04_0000-q6').iterdir() if p.is_file()},q8_export_sha256=sha(payload/'export_review.json'),frontend_shared_sources_unchanged=True,production_eligible=False),indent=2)+'\n');print(json.dumps(dict(staged=str(O),published=False)))
