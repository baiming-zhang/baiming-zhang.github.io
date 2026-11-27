"""Frozen flux CUDA Graph API. No selection,data generation or output-file IO.

Caller holds the cooperative GPU lease. Inputs are CPU float32 already cropped
to the locked model stencil. Batch is captured exactly,no padding/chunk search.
This wraps the same Runner implementation already validated on all six stage3
models,including captured H2D/D2H and changed-input probes.
"""
from pathlib import Path
import hashlib,importlib.util,sys
import numpy as np
import torch

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
CORE_PATH=HERE/'runtime_cuda_graph.py'
FROZEN_CORE_SHA256='964e3b3edaf359d9b7427c4adac96d8fa838fb49f6cf2cec95626c46a16216f0'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def _array(value,ndim):
    if isinstance(value,torch.Tensor):
        assert value.device.type=='cpu' and value.dtype==torch.float32
        value=value.detach().numpy()
    assert isinstance(value,np.ndarray) and value.dtype==np.float32 and value.ndim==ndim
    assert np.isfinite(value).all() and ((value>=0)&(value<=1)).all()
    return np.ascontiguousarray(value)
def _implementation():
    assert sha(CORE_PATH)==FROZEN_CORE_SHA256,'Frozen validated graph implementation changed'
    name='_frozen_flux_graph_core'
    if name not in sys.modules:
        spec=importlib.util.spec_from_file_location(name,CORE_PATH);module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return sys.modules[name]
class FluxCUDAGraphAdapter:
    """One frozen artifact,one exact shape,one registered deployment boundary.

    set_inputs is input-buffer preparation outside operator timing. run includes
    graph launch and completion synchronization;cpu_gpu_cpu graph contains
    actual H2D,all unchanged neural arithmetic,and actual D2H. predict_numpy is
    an accuracy convenience and includes any requested CPU post-copy. It must
    not be substituted for run when timing the resident boundary.
    """
    def __init__(self,artifact_path,stencils_f32,courants_f32,mode='cpu_gpu_cpu',expected_artifact_sha256=None,warmup=12,memory_fraction=.10):
        x=_array(stencils_f32,3);c=_array(courants_f32,1)
        assert x.shape[0]==len(c) and x.shape[1]==x.shape[2]==5
        assert mode in ['cpu_gpu_cpu','resident']
        assert warmup==12,'Frozen validated capture uses12 preparation calls;no post-test warmup tuning'
        self.artifact_path=Path(artifact_path);self.artifact_sha256=sha(self.artifact_path)
        if expected_artifact_sha256 is not None:assert self.artifact_sha256==expected_artifact_sha256
        torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False;torch.set_float32_matmul_precision('highest')
        if memory_fraction is not None:torch.cuda.set_per_process_memory_fraction(memory_fraction)
        self.mode=mode;self.shape=tuple(x.shape);self.stream=torch.cuda.Stream();assert self.stream.cuda_stream!=0
        self.model=torch.jit.load(str(self.artifact_path),map_location='cuda').eval()
        assert all(t.device.type=='cuda' and (not t.is_floating_point() or t.dtype==torch.float32) for t in [*self.model.parameters(),*self.model.buffers()])
        self._runner=_implementation().Runner(self.model,x,c,mode,self.stream)
        self.metadata=dict(runtime='frozen_torchscript_CUDAgraph',adapter_path=str(Path(__file__)),adapter_sha256=sha(__file__),
            validated_core_path=str(CORE_PATH),validated_core_sha256=FROZEN_CORE_SHA256,artifact_path=str(self.artifact_path),artifact_sha256=self.artifact_sha256,
            input_shape=list(self.shape),courant_shape=[len(c)],padding=False,precision='float32;TF32off;noAMP',mode=mode,stream_pointer=int(self.stream.cuda_stream),
            capture_seconds=self._runner.capture_seconds,graph_includes_H2D_D2H=mode=='cpu_gpu_cpu',caller_owns_gpu_lease=True,
            boundary='pinnedCPU ready input -> H2D -> exact frozen model -> D2H -> completedCPU output' if mode=='cpu_gpu_cpu' else 'GPU ready input -> exact frozen model -> completedGPU output')
    def set_inputs(self,stencils_f32,courants_f32):
        x=_array(stencils_f32,3);c=_array(courants_f32,1)
        assert tuple(x.shape)==self.shape and len(c)==self.shape[0],'Capture shape is fixed;new registered shape requires a new adapter'
        self._runner.replace_inputs(x,c)
    def run(self):
        self._runner.invoke('cuda_graph')
        return self._runner.ho if self.mode=='cpu_gpu_cpu' else self._runner.graph_output
    def predict_numpy(self,stencils_f32=None,courants_f32=None):
        if stencils_f32 is not None or courants_f32 is not None:
            assert stencils_f32 is not None and courants_f32 is not None
            self.set_inputs(stencils_f32,courants_f32)
        output=self.run()
        return output.detach().cpu().numpy().copy()
    def close(self):
        self.stream.synchronize();self._runner=None;self.model=None

def cpu_contract_check():
    assert sha(CORE_PATH)==FROZEN_CORE_SHA256
    x=np.zeros((17,5,5),dtype=np.float32);c=np.linspace(0,1,17,dtype=np.float32)
    assert _array(x,3).shape==(17,5,5) and _array(c,1).shape==(17,)
    assert not torch.cuda.is_initialized(),'CPU adapter contract check must not initialize CUDA'
    return dict(passed=True,cuda_initialized=False,adapter_sha256=sha(__file__),validated_core_sha256=FROZEN_CORE_SHA256,
        GPU_parity_reused='stage3 runtime_summary accuracy:6 frozen models,full16384 rows+capacity and changed-input probe,graph/eager maxdiff0',API='FluxCUDAGraphAdapter(...).run()/predict_numpy()/set_inputs()')
if __name__=='__main__':
    import json
    print(json.dumps(cpu_contract_check(),indent=2))
