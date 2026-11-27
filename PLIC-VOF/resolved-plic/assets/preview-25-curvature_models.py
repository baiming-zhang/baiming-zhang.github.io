"""Fixed pure-network symmetry averages; no geometry inputs or fitted weights."""
import torch
from torch import nn

METHODS=('raw','d4','phase_odd','d4_phase_odd')

def d4_views(x):
    return [torch.rot90(x,r,(-2,-1)) for r in range(4)]+[torch.rot90(x.flip(-1),r,(-2,-1)) for r in range(4)]

class CurvatureProjection(nn.Module):
    def __init__(self,core,method):
        super().__init__()
        if method not in METHODS:raise ValueError(method)
        self.core=core;self.method=method
    def forward(self,stencils):
        if self.method=='raw':return self.core(stencils)
        views=d4_views(stencils) if self.method in ('d4','d4_phase_odd') else [stencils]
        spatial=torch.cat(views,0)
        if self.method in ('phase_odd','d4_phase_odd'):
            # Complement actual raw VOF, including its ordinary FP32 rounding.
            both=self.core(torch.cat((spatial,1-spatial),0)).reshape(2,len(views),-1)
            return ((both[0]-both[1])*.5).mean(0)
        return self.core(spatial).reshape(len(views),-1).mean(0)

def view_counts(method):
    calls={'raw':1,'d4':8,'phase_odd':2,'d4_phase_odd':16}[method]
    return dict(logical_source_predictions_per_sample=calls,source_internal_row_reflection_views=2,
                net_rows_per_sample=2*calls,distinct_spatial_phase_views={'raw':2,'d4':8,'phase_odd':4,'d4_phase_odd':16}[method],
                note='D4 contains duplicate internal row-reflection views; implemented faithfully as ordinary-source averaging, not deduplicated or timed.')
