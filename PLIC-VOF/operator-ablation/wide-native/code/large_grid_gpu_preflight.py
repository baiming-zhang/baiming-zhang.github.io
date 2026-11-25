"""Functional CUDA transfer/chunk check under shared load; never speed evidence."""
from pathlib import Path
import os,sys,json,importlib.util,gc,time
for env in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[env]='1'
import numpy as np
import torch
from large_grid_operator_fixtures import ROOT,R,save
sys.path.insert(0,str(ROOT/'operator-width-native-baselines/02-native-bench'))
sys.path.insert(0,str(ROOT/'operator-ablation-20261122/shared'))
from nn_registry import registry,cpu32_session
from gpu_lease import gpu_lease

def main():
    destination=R/'gpu_functional_preflight.json'
    if destination.exists():print('Existing GPU functional checks retained');return
    torch.set_num_threads(1);torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    spec=importlib.util.spec_from_file_location('large_grid_functional_session',ROOT/'operator-ablation-20261122/02-inference-pareto/harness.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h);h.HERE=R/'functional_preflight';h.HERE.mkdir(exist_ok=True)
    candidates=[m for m in registry() if m.get('historical_reference') or (m.get('seed')==4101 and m.get('data_strategy')=='control' and m['stencil_size'] in [3,11])]
    assert len(candidates)==6
    full=np.load(R/'fixtures/many_drops_N4096.npz')['stencils'];indices=np.linspace(0,len(full)-1,37,dtype=int);results=[];begin=time.perf_counter()
    with gpu_lease('large-grid-functional-preflight',formal_timing=False):
        stream=torch.cuda.Stream()
        for m in candidates:
            k=int(m['stencil_size']);a=(11-k)//2;x=np.ascontiguousarray(full[indices,a:a+k,a:a+k]);c=np.full(len(x),.5,dtype=np.float32)
            cpu=cpu32_session(m);expected=cpu.run(None,{inp.name:c if len(inp.shape)==1 else x for inp in cpu.get_inputs()})[0].reshape(-1)
            gpu=h.session(ROOT/m['onnx_path'],stream);host=torch.from_numpy(x).pin_memory();hc=torch.from_numpy(c).pin_memory();output=torch.empty(len(x),pin_memory=True);chunks=[]
            with torch.cuda.stream(stream):
                for start in range(0,len(x),17):
                    end=min(start+17,len(x));xx=torch.empty((end-start,k,k),device='cuda');cc=torch.empty(end-start,device='cuda');yy=torch.empty(end-start,device='cuda');binding=gpu.io_binding()
                    for inp in gpu.get_inputs():
                        t=cc if len(inp.shape)==1 else xx;binding.bind_input(inp.name,'cuda',0,np.float32,tuple(t.shape),t.data_ptr())
                    binding.bind_output(gpu.get_outputs()[0].name,'cuda',0,np.float32,(end-start,),yy.data_ptr())
                    xx.copy_(host[start:end],non_blocking=True);cc.copy_(hc[start:end],non_blocking=True);chunks.append((start,end,xx,cc,yy,binding))
            stream.synchronize();values=[]
            for mode in ['cpu_gpu_cpu','resident']:
                with torch.cuda.stream(stream):
                    for start,end,xx,cc,yy,binding in chunks:
                        if mode=='cpu_gpu_cpu':xx.copy_(host[start:end],non_blocking=True);cc.copy_(hc[start:end],non_blocking=True)
                        gpu.run_with_iobinding(binding);output[start:end].copy_(yy,non_blocking=True)
                stream.synchronize();values.append(output.numpy().copy())
            diff=float(abs(values[0]-expected).max());mode_diff=float(abs(values[0]-values[1]).max());assert diff<1e-5 and mode_diff==0,(m['model_id'],diff,mode_diff)
            results.append(dict(model=m['model_id'],operator=m['operator'],stencil=k,queries=37,chunks=[17,17,3],cpu_gpu_max_abs=diff,resident_roundtrip_max_abs=mode_diff))
            del cpu,gpu,chunks,host,hc,output;gc.collect();torch.cuda.empty_cache()
    save(destination,dict(state='passed',purpose='Functional correctness under shared GPU load only;not formal speed measurements',models=len(results),checks=results,elapsed_seconds=time.perf_counter()-begin,formal_timing=False))
    print(json.dumps(dict(state='passed',models=len(results),maximum_difference=max(r['cpu_gpu_max_abs'] for r in results),elapsed_seconds=time.perf_counter()-begin)))

if __name__=='__main__':main()
