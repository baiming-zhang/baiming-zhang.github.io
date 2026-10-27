"""Promote only the verified q6 viewer; preserve original viewers/shared assets."""
from pathlib import Path
import hashlib,json,shutil
W=Path(__file__).resolve().parents[1];B=W/'data/copv/c04_repair_v3';P=B/'curved_viewer_v1';S=W/'showcase-public-site/GeoBC-show'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
qa=read(P/'browser_qa_v3/review.json');assert qa['status']=='browser_checks_passed_pending_visual_inspection' and not qa['errors'] and 'SwiftShader' in qa['software_renderer'];assert qa['after']['bcOn']=='true' and qa['after']['loadOn']=='true';assert qa['before']['basis']!=qa['after']['basis'];export=read(P/'assets/copv-c04-repair-fields/copv_c04_0000-q6/export_review.json');assert export['actual_native_loads_checked'] and export['locked_components']==6 and export['pressure_arrows']==28
visual=dict(status='accepted_visual',reviewed_screenshots=['01_field.png','02_BC_load.png','03_rotated.png','04_clipped_mesh.png'],scope='Field legend, cap and pressure solid3D arrows without white boxes, actualA/B/C translation supports, rotation-dependent world directions, clip/mesh controls visually inspected. Software-WebGL avoided physical GPU use. Original access guard retained in production files.',browser_QA_sha256=hashlib.sha256((P/'browser_qa_v3/review.json').read_bytes()).hexdigest());vp=P/'browser_qa_v3/visual_acceptance_v1.json';assert not vp.exists();vp.write_text(json.dumps(visual,indent=2)+'\n')
for name in ['access.css','access.js','style.css','fe-glyphs-3d.js']:assert (P/'assets'/name).read_bytes()==(S/'assets'/name).read_bytes(),'Shared asset changed after QA: '+name
for folder,name in [('viewers','copv-c04-curved.html'),('assets','copv-c04-curved-fields.js'),('assets','copv-c04-curved-bc.js')]:
 target=S/folder/name;assert not target.exists(),'Initial promotion must not overwrite existing viewer';shutil.copyfile(P/folder/name,target)
src=P/'assets/copv-c04-repair-fields/copv_c04_0000-q6';dest=S/'assets/copv-c04-repair-fields/copv_c04_0000-q6';dest.parent.mkdir(exist_ok=True);shutil.copytree(src,dest)
print(json.dumps(dict(viewer='viewers/copv-c04-curved.html?case=copv_c04_0000&q=6&bc=1&loads=1',status='QA_complete_pending_publication')))
