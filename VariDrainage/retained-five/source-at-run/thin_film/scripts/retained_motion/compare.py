"""Original numerical-equivalence thresholds, explicitly scoped to saved step1."""
from pathlib import Path
SOURCE=Path(__file__).parents[1]/'local_assembly/equivalence.py'
def load():
    code=SOURCE.read_text(encoding='utf-8')
    old=';entry=dict(same_structure='
    assert code.count(old)==1
    code=code.replace(old,";a=a[a['step']==1] if filename in ['states.csv','faces.csv'] else a;b=b[b['step']==1] if filename in ['states.csv','faces.csv'] else b;entry=dict(same_structure=",1)
    old="for filename in ['summary.json','mesh.json','preparation.json']:";assert code.count(old)==1
    code=code.replace(old,"for filename in ['mesh.json','preparation.json']:")
    ns=dict(__file__=str(Path(__file__)),__name__='retained_first_step_equivalence');exec(compile(code,str(SOURCE),'exec'),ns)
    return ns['check']
