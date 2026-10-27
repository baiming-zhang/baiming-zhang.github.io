"""Publish verified q8 assets, keeping every immutable q6 payload untouched."""
from pathlib import Path
import hashlib,json,shutil
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3';P=B/'curved_viewer_q8_v1';S=W/'showcase-public-site/GeoBC-show'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
qa=read(P/'browser_qa_v1/review.json');assert qa['assembly_degree']==8 and qa['degree_selector_roundtrip'] and not qa['errors'] and 'SwiftShader' in qa['software_renderer'] and qa['after']['bcOn']=='true' and qa['after']['loadOn']=='true'
for name in ['access.css','access.js','style.css','fe-glyphs-3d.js','copv-c04-curved-fields.js','copv-c04-curved-bc.js']:assert (P/'assets'/name).read_bytes()==(S/'assets'/name).read_bytes(),'Shared frontend changed after QA'
proof=read(P/'preview_prepared_v1.json')
for name,digest in proof['q6_original_site_hashes'].items():assert sha(S/'assets/copv-c04-repair-fields/copv_c04_0000-q6'/name)==digest
accept=P/'browser_qa_v1/visual_acceptance_v1.json';assert not accept.exists();accept.write_text(json.dumps(dict(status='accepted_visual',assembly_degree=8,reviewed_screenshots=['01_field.png','02_BC_load.png','03_rotated.png','04_clipped_mesh.png'],scope='Actual q8 fields, preserved true3D rotating BC/load arrows without white boxes, six constraint components, curved cap forces, clipping/wireframe and q8-to-q6-to-q8 selector roundtrip checked. No affine cut fill. Original access guard retained.',browser_QA_sha256=sha(P/'browser_qa_v1/review.json')),indent=2)+'\n')
src=P/'assets/copv-c04-repair-fields/copv_c04_0000-q8';dest=S/'assets/copv-c04-repair-fields/copv_c04_0000-q8';assert not dest.exists();shutil.copytree(src,dest);shutil.copyfile(P/'viewers/copv-c04-curved.html',S/'viewers/copv-c04-curved.html')
print(json.dumps(dict(promoted_q8=True,q6_unchanged=True,viewer='viewers/copv-c04-curved.html?case=copv_c04_0000&q=8&bc=1&loads=1')))
