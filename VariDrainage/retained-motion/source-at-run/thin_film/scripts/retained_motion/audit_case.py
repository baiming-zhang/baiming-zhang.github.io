"""Reuse original independent SI/GEOS formulas; path/metadata adaptation only."""
from pathlib import Path
import sys
from run import O,read,save,sha
SOURCE=Path(__file__).parents[1]/'global_space/audit_ready.py'
def main(name):
    p=read(O/'protocol.json');assert p['sources'][str(Path(__file__).relative_to(Path(__file__).parents[3]))]==sha(Path(__file__))
    code=SOURCE.read_text(encoding='utf-8');old='global_space_audit_inputs_20261113';assert code.count(old)==1;code=code.replace(old,'retained_motion_audit_inputs_20261113')
    old='24previously executed unit checks verified via frozen hashes; not rerun here';assert code.count(old)==1;code=code.replace(old,'24candidate fixture checks verified via frozen hashes; not rerun by this auditor')
    ns=dict(__file__=str(Path(__file__)),__name__='retained_motion_independent_auditor');exec(compile(code,str(SOURCE),'exec'),ns)
    s=read(O/'status.json');rec=next(r for r in s['cases'] if r['name']==name);assert rec['passed'] and rec['native_gate'];d=ns['one'](rec)
    save(d/'wrapper_verified.json',dict(original_auditor_sha256=sha(SOURCE),protocol_sha256=sha(O/'protocol.json'),sealed_outputs_sha256=sha(d/'verified.json'),passed=True,scope='Original saved-step SI/GEOS scalar/mechanical formulas; no native replay or full dynamic qualification'))
if __name__=='__main__':main(sys.argv[1])
