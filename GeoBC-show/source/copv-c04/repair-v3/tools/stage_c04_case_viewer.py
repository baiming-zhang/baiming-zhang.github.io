"""Stage/promote only a newly exported original case; existing payloads immutable."""
from pathlib import Path
import argparse,hashlib,json,shutil
W=Path(__file__).resolve().parents[1];S=W/'showcase-public-site/GeoBC-show';B=W/'data/copv/c04_repair_v3'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,v):
 assert not p.exists(),str(p)+' already exists'
 p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
ap=argparse.ArgumentParser();ap.add_argument('--case',required=True,choices=['copv_c04_0011','copv_c04_0002','copv_c04_0018']);ap.add_argument('--stage',type=Path,required=True);ap.add_argument('--promote',action='store_true');a=ap.parse_args();O=a.stage.resolve();case=a.case
names=['access.css','access.js','style.css','fe-glyphs-3d.js','copv-c04-curved-fields.js','copv-c04-curved-bc.js'];payloads=[O/'assets/copv-c04-repair-fields'/f'{case}-q{q}' for q in [6,8]]
for q,p in zip([6,8],payloads):
 m=read(p/'manifest.json');e=read(p/'export_review.json');assert m['summary']['case_id']==case and m['summary']['assembly_quadrature_degree']==q and not m['production_eligible'] and e['actual_native_loads_checked'] and e['locked_components']==6 and e['pressure_arrows']==28
if not a.promote:
 assert not (O/'preview_prepared_v1.json').exists()
 for name in names:
  dest=O/'assets'/name;assert not dest.exists();shutil.copyfile(S/'assets'/name,dest)
 (O/'viewers').mkdir(exist_ok=False);shutil.copyfile(S/'viewers/copv-c04-curved.html',O/'viewers/copv-c04-curved.html')
 old={str(p.relative_to(S)).replace('\\','/'):sha(p) for folder in (S/'assets/copv-c04-repair-fields').glob('copv_c04_*-q[68]') for p in folder.iterdir() if p.is_file()}
 put(O/'preview_prepared_v1.json',dict(case_id=case,status='staged_not_published',unchanged_old_payload_SHA256=old,export_review_sha256=[sha(p/'export_review.json') for p in payloads],production_eligible=False))
else:
 proof=read(O/'preview_prepared_v1.json');qa=read(O/'browser_qa_v1/review.json');v=read(O/'browser_qa_v1/visual_acceptance_v1.json')
 assert qa['case_id']==case and qa['assembly_degree']==8 and qa['degree_selector_roundtrip'] and 'SwiftShader' in qa['software_renderer'] and not qa['errors'] and v['status']=='accepted_visual' and v['browser_QA_sha256']==sha(O/'browser_qa_v1/review.json')
 for name in names:assert (O/'assets'/name).read_bytes()==(S/'assets'/name).read_bytes()
 assert (O/'viewers/copv-c04-curved.html').read_bytes()==(S/'viewers/copv-c04-curved.html').read_bytes()
 for p,h in proof['unchanged_old_payload_SHA256'].items():assert sha(S/p)==h,'Existing curved payload changed'
 targets=[S/'assets/copv-c04-repair-fields'/p.name for p in payloads];assert not any(p.exists() for p in targets)
 for p,t in zip(payloads,targets):shutil.copytree(p,t)
 for p,h in proof['unchanged_old_payload_SHA256'].items():assert sha(S/p)==h
 put(O/'promotion_v1.json',dict(case_id=case,new_payload_directories=[p.name for p in targets],old_payloads_unchanged=True,shared_frontend_unchanged=True,viewer=f'viewers/copv-c04-curved.html?case={case}&q=8&bc=1&loads=1',production_eligible=False))
print(json.dumps(dict(case_id=case,promoted=a.promote)))
