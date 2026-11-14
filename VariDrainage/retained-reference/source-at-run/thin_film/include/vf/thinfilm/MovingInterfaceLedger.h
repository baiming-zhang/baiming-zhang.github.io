// Regional bookkeeping on the ACTUAL moving Voronoi state. No force/flux law
// is changed. Fractional control-volume ownership is a piecewise-constant
// reconstruction, not an independently resolved material boundary velocity.
#pragma once
#include "vf/thinfilm/MovingVoronoiFilm.h"
namespace vf::thinfilm {
inline std::vector<real> Moving_Interface_Weights(const VoronoiDrainageGeometry& g,real cut){
    std::vector<real> w;for(size_t i=0;i<g.cells.size();++i)
        w.push_back(std::clamp(Transfer_Area(Transfer_Clip(g.cells[i],{cut,0.},{-1.,0.}))/g.graph.a[i],real(0),real(1)));
    return w; // right half-plane x >= cut
}
struct MovingInterfaceLedger {
    real volume_change=0.,material_volume=0.,relative_volume=0.,volume_residual=0.;
    real surface_change=0.,material_surface=0.,surface_residual=0.;
    real capillary_change=0.,liquid_work=0.,area_work=0.,shape_work=0.,metric_work=0.,ownership_work=0.,capillary_residual=0.;
    real surface_energy_change=0.,surface_area_work=0.,surface_ownership_work=0.,surface_energy_residual=0.;
};
// S remains attached to material cells in this insoluble split solver. The
// fractional regional transfer is nonzero even without centroid label flips.
// Energy reference is EXACTLY the fixed-domain uniform-monolayer reference
// used in Differentiate_Voronoi_Energy, including the factor 2 for two faces.
inline MovingInterfaceLedger Moving_Interface_Ledger(const FilmMatState& a,const FilmMatState& b,
    const VoronoiDrainageGeometry& ga,const VoronoiDrainageGeometry& gb,
    const std::vector<real>& wa,const std::vector<real>& wb,const FilmMechParams& p,
    const std::vector<real>* flux=nullptr,const std::vector<real>* reservoir=nullptr){
    const int n=a.Size();if(b.Size()!=n||wa.size()!=n||wb.size()!=n)throw std::invalid_argument("interface ledger size");
    MovingInterfaceLedger r;real totalS=0.,totalA=0.;for(int i=0;i<n;++i){totalS+=a.S[i];totalA+=ga.graph.a[i];
        if(a.S[i]!=b.S[i])throw std::invalid_argument("interface ledger requires insoluble material S");}
    const real mean=totalS/totalA;auto surface=[&](real area,real amount){const real gamma=amount/area;
        if(p.surf.Gamma_inf>0.){const real t=gamma/p.surf.Gamma_inf,m=mean/p.surf.Gamma_inf,d=t-m;
            return 2*area*p.surf.RT*p.surf.Gamma_inf*(t*std::log1p(d/m)+(1-t)*std::log1p(-d/(1-m)));}
        return 2*area*p.surf.RT*(gamma*std::log1p((gamma-mean)/mean)-(gamma-mean));};
    for(int i=0;i<n;++i){r.volume_change+=wb[i]*b.V[i]-wa[i]*a.V[i];r.material_volume+=(wb[i]-wa[i])*a.V[i];
        r.surface_change+=wb[i]*b.S[i]-wa[i]*a.S[i];r.material_surface+=(wb[i]-wa[i])*a.S[i];
        const real ea=surface(ga.graph.a[i],a.S[i]),eb=surface(gb.graph.a[i],b.S[i]);
        r.surface_energy_change+=wb[i]*eb-wa[i]*ea;r.surface_area_work+=.5*(wa[i]+wb[i])*(eb-ea);
        r.surface_ownership_work+=.5*(ea+eb)*(wb[i]-wa[i]);}
    if(flux){if(flux->size()!=gb.graph.pairs.Size())throw std::invalid_argument("interface flux size");
        for(int q=0;q<gb.graph.pairs.Size();++q)r.relative_volume+=(wb[gb.graph.pairs.pj[q]]-wb[gb.graph.pairs.pi[q]])*(*flux)[q];}
    // Union of material pair IDs also covers a face disappearing/appearing.
    std::map<std::pair<int,int>,std::array<real,2>> faces;
    auto fill=[&](const ThinFilmGraph& g,int stage){for(int q=0;q<g.pairs.Size();++q){int i=g.pairs.pi[q],j=g.pairs.pj[q];if(i>j)std::swap(i,j);
        if(reservoir&&(!std::isnan((*reservoir)[i])||!std::isnan((*reservoir)[j])))continue;faces[{i,j}][stage]=g.c[q];}};
    fill(ga.graph,0);fill(gb.graph,1);const real k=p.tf.sigma/8.;
    for(const auto& entry:faces){const auto [i,j]=entry.first;const real ca=entry.second[0],cb=entry.second[1];
        const real hai=a.V[i]/ga.graph.a[i],haj=a.V[j]/ga.graph.a[j],hbi=b.V[i]/gb.graph.a[i],hbj=b.V[j]/gb.graph.a[j];
        const real da=hai-haj,db=hbi-hbj,ea=k*ca*da*da,eb=k*cb*db*db;
        const real wi=.5*(wa[i]+wb[i]),wj=.5*(wa[j]+wb[j]),alpha=.5*(wi+wj);
        const real t=2*k*.5*(ca+cb)*.5*(da+db),g=k*.5*(da*da+db*db);
        const real dVi=b.V[i]-a.V[i],dVj=b.V[j]-a.V[j];
        const real bi=.5*(1/ga.graph.a[i]+1/gb.graph.a[i]),bj=.5*(1/ga.graph.a[j]+1/gb.graph.a[j]);
        r.capillary_change+=.5*(wb[i]+wb[j])*eb-.5*(wa[i]+wa[j])*ea;
        r.liquid_work+=t*(wi*bi*dVi-wj*bj*dVj);
        r.area_work+=t*(wi*.5*(a.V[i]+b.V[i])*(1/gb.graph.a[i]-1/ga.graph.a[i])-wj*.5*(a.V[j]+b.V[j])*(1/gb.graph.a[j]-1/ga.graph.a[j]));
        r.shape_work+=.5*t*(wj-wi)*((hbi-hai)+(hbj-haj));
        r.metric_work+=alpha*g*(cb-ca);
        r.ownership_work+=.5*(ea+eb)*.5*(wb[i]+wb[j]-wa[i]-wa[j]);
    }
    r.volume_residual=r.volume_change-r.material_volume-r.relative_volume;
    r.surface_residual=r.surface_change-r.material_surface;
    r.capillary_residual=r.capillary_change-r.liquid_work-r.area_work-r.shape_work-r.metric_work-r.ownership_work;
    r.surface_energy_residual=r.surface_energy_change-r.surface_area_work-r.surface_ownership_work;
    return r;
}
} // namespace vf::thinfilm
