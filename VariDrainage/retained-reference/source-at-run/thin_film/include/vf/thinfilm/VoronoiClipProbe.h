// Isolated geometry-cost candidate: not used by any production solver.
// Skip allocation only when the ORIGINAL floating-point plane predicates
// classify every vertex inside. Keep original site/cut/vertex order and
// intersection arithmetic. Established plane-intersection screening, not novelty.
#pragma once
#include "vf/thinfilm/OverlapTransfer.h"
namespace vf::thinfilm {
struct VoronoiClipCounts {size_t attempted=0,unchanged=0,changed=0;};
inline std::vector<TransferPolygon> Probe_Voronoi_Unchanged_Cuts(
    const std::vector<TransferPoint>& points,const TransferPolygon& patch,
    VoronoiClipCounts* counts=nullptr)
{
    std::vector<TransferPolygon> out;out.reserve(points.size());
    for(size_t i=0;i<points.size();++i){auto p=patch;
        if(!points[i].allFinite())throw std::invalid_argument("nonfinite transfer point");
        for(size_t j=0;j<points.size()&&!p.empty();++j)if(i!=j){
            const TransferPoint d=points[j]-points[i];
            if(!(d.squaredNorm()>0.))throw std::invalid_argument("duplicate transfer point");
            const TransferPoint origin=.5*(points[i]+points[j]);
            bool unchanged=true;
            for(const auto& x:p)if(!((x-origin).dot(d)<=0.)){unchanged=false;break;}
            if(counts){++counts->attempted;if(unchanged)++counts->unchanged;else ++counts->changed;}
            if(!unchanged)p=Transfer_Clip(p,origin,d);
        }
        if(!(Transfer_Area(p)>0.))throw std::invalid_argument("transfer point has no patch support");
        out.push_back(std::move(p));
    }
    return out;
}
}
