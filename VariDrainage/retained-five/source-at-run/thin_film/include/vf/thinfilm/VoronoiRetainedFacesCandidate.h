// Opt-in diagnostic candidate; original geometry/defaults unchanged.
#pragma once
#include "vf/thinfilm/VoronoiDrainage.h"
namespace vf::thinfilm {
inline VoronoiDrainageGeometry Make_Voronoi_Drainage_RetainedFaces(const std::vector<TransferPoint>& sites,
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
            if(!(length>0.))continue; // candidate: retain every positive edge after original cleanup
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
}
