"""Stage2-only wrapper around verified validation Pareto evaluation."""
from pathlib import Path
import sys,os
sys.dont_write_bytecode=True
import hashlib,importlib.util,json,traceback
R=Path(__file__).resolve().parent;PARENT=R.parent;W=PARENT.parent
sys.path.insert(0,str(R))
from export_stage2 import setup_aliases,save,sha

def load_validator():
    base=PARENT/'validation_pareto.py';spec=importlib.util.spec_from_file_location('_stage2_validation_base',base)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    module.R=R;module.W=W;module.OUT=R/'validation_pareto';module.__file__=__file__
    module.PROTOCOL='resolved_stage2_validation_gpu_native_v1_base'+sha(base)[:12]
    save(R/'validation_wrapper_metadata.json',dict(stage='stage2_courant_objective',wrapper_path=str(Path(__file__).relative_to(W)),wrapper_sha256=sha(__file__),
        reused_base_script=str(base.relative_to(W)),base_script_sha256=sha(base),protocol=module.PROTOCOL,
        data_aliases='input_aliases.json',output_root=str((R/'validation_pareto').relative_to(W)),workspace=str(W),
        numerical_timing_accuracy_code_unchanged=True,validation_only=True,stage1_results_modified=False))
    return module
def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--formal',action='store_true');parser.add_argument('--repeats',type=int,default=31)
    parser.add_argument('--describe',action='store_true');parser.add_argument('--cpu-check',action='store_true');args=parser.parse_args();assert args.repeats>=31
    setup_aliases();module=load_validator();module.main(args)
if __name__=='__main__':
    try:main()
    except Exception:
        save(R/'validation_pareto/status.json',dict(state='failed',pid=os.getpid(),traceback=traceback.format_exc()));raise
