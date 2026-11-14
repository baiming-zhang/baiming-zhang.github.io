// Orthogonal planar finite-volume geometry for the DRIFT drainage subproblem.
// Uses the inherited energy and Onsager flux step; no change to kernel mechanics.
// A single tessellation spans coarse film sites and dense local-edge sites.
#pragma once
#include <map>
#include <sstream>
#include "vf/thinfilm/OverlapTransfer.h"
#include "vf/thinfilm/ThinFilm.h"
namespace vf::thinfilm {
// A vanishing Voronoi face has an ill-conditioned ANGLE but well-conditioned
// absolute plane coordinates. Preserve the usual angle bound on resolved faces,
// and include a coordinate-roundoff bound on tiny faces; never move vertices.
inline real Voronoi_Orthogonality_Allowance(real length,real distance,real coordinate,real tolerance){
    return std::max(10*tolerance*length,64*std::numeric_limits<real>::epsilon()*coordinate*(1.+length/distance));
}
// Half-plane clipping can emit consecutive copies of a vertex separated only by
// coordinate rounding. Their nearly zero edges have meaningless turning angles.
// Remove ONLY consecutive coordinate duplicates at machine precision, never take
// a convex hull or relax a genuinely nonconvex corner. Area is audited by caller.
inline int Clean_Drainage_Polygon(TransferPolygon& p)
{
    real coordinate=0.;for(const auto& x:p)coordinate=std::max(coordinate,x.norm());
    const real eps=32*std::numeric_limits<real>::epsilon()*coordinate;
    TransferPolygon out;out.reserve(p.size());
    for(const auto& x:p)if(out.empty()||(x-out.back()).norm()>eps)out.push_back(x);
    if(out.size()>1&&(out.back()-out.front()).norm()<=eps)out.pop_back();
    const int removed=int(p.size()-out.size());p=std::move(out);return removed;
}
struct VoronoiDrainageGeometry {
    ThinFilmGraph graph;
    std::vector<TransferPolygon> cells;
    struct Face {int i,j;TransferPoint midpoint,normal;real length;}; // j<0: closed exterior face
    std::vector<Face> faces;
    real normal_closure=0.,orthogonality_error=0.,area_closure=0.;
    int duplicate_vertices_removed=0;
    real maximum_cleanup_area_change=0.;
    real orthogonality_certificate=0.,worst_angle_face_length=0.;
};
inline VoronoiDrainageGeometry Make_Voronoi_Drainage(const std::vector<TransferPoint>& sites,
    const TransferPolygon& domain,real tolerance=1e-10)
{
    if(sites.empty()||!(tolerance>0.&&tolerance<1e-3))throw std::invalid_argument("invalid Voronoi drainage geometry input");
    VoronoiDrainageGeometry r;r.cells=Transfer_Voronoi(sites,domain);
    for(auto& cell:r.cells){const real area=Transfer_Area(cell);r.duplicate_vertices_removed+=Clean_Drainage_Polygon(cell);
        r.maximum_cleanup_area_change=std::max(r.maximum_cleanup_area_change,std::abs(Transfer_Area(cell)/area-1.));}
    if(r.maximum_cleanup_area_change>1e-12)throw std::runtime_error("duplicate cleanup changed geometric area");
    Validate_Transfer_Partition(r.cells,domain,tolerance);
    const auto bounds=TransferBounds(domain);const real scale=(bounds.hi-bounds.lo).norm(),eps=tolerance*scale;
    real coordinate=scale;for(auto x:domain)coordinate=std::max(coordinate,x.norm());
    auto& G=r.graph;G.N=int(sites.size());G.a.resize(G.N);G.height.resize(G.N);
    struct Sides {real length[2]={0.,0.};TransferPoint midpoint=TransferPoint::Zero();int count[2]={0,0};};
    std::map<std::pair<int,int>,Sides> common;
    real sum_area=0.;
    for(int i=0;i<G.N;++i){
        G.x.emplace_back(sites[i].x(),sites[i].y(),0.);G.height[i]=sites[i].y();
        G.a[i]=Transfer_Area(r.cells[i]);sum_area+=G.a[i];
        const auto& cell=r.cells[i];
        for(size_t k=0;k<cell.size();++k){
            const TransferPoint a=cell[k],b=cell[(k+1)%cell.size()],e=b-a;const real length=e.norm();
            if(length<=eps)continue; // only machine-scale degenerate segments
            if(Transfer_Cross(e,sites[i]-a)<-eps*length)throw std::invalid_argument("Voronoi generator outside its cell");
            bool exterior=false;
            for(size_t d=0;d<domain.size();++d){const TransferPoint v=domain[(d+1)%domain.size()]-domain[d];const real vd=v.norm();
                if(vd>0.&&std::abs(Transfer_Cross(v,a-domain[d]))/vd<=eps&&std::abs(Transfer_Cross(v,b-domain[d]))/vd<=eps){exterior=true;break;}}
            if(exterior){r.faces.push_back({i,-1,.5*(a+b),TransferPoint(e.y(),-e.x())/length,length});continue;}
            int best=-1;real error=std::numeric_limits<real>::infinity();
            for(int j=0;j<G.N;++j)if(j!=i){const TransferPoint n=sites[j]-sites[i],m=.5*(sites[i]+sites[j]);
                const real err=std::max(std::abs((a-m).dot(n)),std::abs((b-m).dot(n)))/n.norm();
                if(err<error){error=err;best=j;}}
            if(best<0||error>eps)throw std::runtime_error("unmatched interior Voronoi face");
            const int lo=std::min(i,best),hi=std::max(i,best),side=i==lo?0:1;auto& f=common[{lo,hi}];
            f.length[side]+=length;f.midpoint+=.5*length*(a+b);++f.count[side];
            const real distance=(sites[best]-sites[i]).norm(),projection=std::abs(e.dot(sites[best]-sites[i]))/distance,angle=projection/length;
            if(angle>r.orthogonality_error){r.orthogonality_error=angle;r.worst_angle_face_length=length;}
            r.orthogonality_certificate=std::max(r.orthogonality_certificate,projection/Voronoi_Orthogonality_Allowance(length,distance,coordinate,tolerance));
        }
    }
    std::vector<TransferPoint> closure(G.N,TransferPoint::Zero());std::vector<real> perimeter(G.N,0.);
    for(const auto& entry:common){const int i=entry.first.first,j=entry.first.second;const auto& f=entry.second;
        if(!f.count[0]||!f.count[1]||std::abs(f.length[0]-f.length[1])>eps)
            throw std::runtime_error("Voronoi face sides do not agree");
        const real length=.5*(f.length[0]+f.length[1]),distance=(sites[j]-sites[i]).norm();
        const TransferPoint n=(sites[j]-sites[i])/distance;
        G.pairs.pi.push_back(i);G.pairs.pj.push_back(j);G.c.push_back(length/distance);
        r.faces.push_back({i,j,f.midpoint/(f.length[0]+f.length[1]),n,length});
    }
    G.pairs.Build_Incidence(G.N);
    for(const auto& f:r.faces){closure[f.i]+=f.length*f.normal;perimeter[f.i]+=f.length;
        if(f.j>=0){closure[f.j]-=f.length*f.normal;perimeter[f.j]+=f.length;}}
    for(int i=0;i<G.N;++i)r.normal_closure=std::max(r.normal_closure,closure[i].norm()/perimeter[i]);
    r.area_closure=std::abs(sum_area/Transfer_Area(domain)-1.);
    if(r.normal_closure>10*tolerance||r.orthogonality_certificate>1.){std::ostringstream diagnostic;
        diagnostic.precision(17);diagnostic<<"inconsistent Voronoi face geometry: normal_closure="<<r.normal_closure<<" angle_error="<<r.orthogonality_error<<" face_length="<<r.worst_angle_face_length<<" roundoff_certificate="<<r.orthogonality_certificate;
        throw std::runtime_error(diagnostic.str());}
    return r;
}
// Region a->b, positive if a loses liquid; the SAME pair flux enters both ledgers.
struct VoronoiInterfaceTransfer {real liquid=0.,dissolved=0.;int forward=0,reverse=0;};
inline VoronoiInterfaceTransfer Voronoi_Interface_Transfer(const ThinFilmGraph& G,
    const std::vector<int>& labels,const std::vector<real>& pair_volume,
    const std::vector<real>& final_V,const std::vector<real>& final_B,int a=0,int b=1)
{
    if(a==b||labels.size()!=size_t(G.N)||final_V.size()!=labels.size()||final_B.size()!=labels.size()||pair_volume.size()!=G.c.size())
        throw std::invalid_argument("invalid Voronoi interface ledger input");
    VoronoiInterfaceTransfer t;
    for(int q=0;q<G.pairs.Size();++q){const int i=G.pairs.pi[q],j=G.pairs.pj[q];
        const int sign=labels[i]==a&&labels[j]==b?1:labels[i]==b&&labels[j]==a?-1:0;
        if(!sign)continue;
        const real F=pair_volume[q],oriented=sign*F;const int donor=F>0.?i:j;
        if(!(final_V[donor]>0.))throw std::invalid_argument("nonpositive interface donor volume");
        t.liquid+=oriented;t.dissolved+=oriented*final_B[donor]/final_V[donor];
        if(oriented>0.)++t.forward;else if(oriented<0.)++t.reverse;
    }return t;
}
} // namespace vf::thinfilm
