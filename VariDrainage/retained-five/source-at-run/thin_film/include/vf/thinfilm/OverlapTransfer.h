// Planar convex-control-volume transfer for the DRIFT local/particle interface.
// First-order geometric remapping is established ALE machinery, not a novelty claim.
// No drainage, mechanical time integration, bath or adsorption law lives here.
#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "vf/core/Types.h"

namespace vf::thinfilm {
using TransferPoint = Eigen::Matrix<real,2,1>;
using TransferPolygon = std::vector<TransferPoint>;

inline real Transfer_Cross(const TransferPoint& a,const TransferPoint& b){return a.x()*b.y()-a.y()*b.x();}
inline real Transfer_Area(const TransferPolygon& p)
{
    if(p.size()<3)return 0.;
    long double a=0.;
    for(size_t i=1;i+1<p.size();++i){const TransferPoint u=p[i]-p[0],v=p[i+1]-p[0];
        a+=(long double)u.x()*v.y()-(long double)u.y()*v.x();}
    return real(a*.5L); // signed: callers must supply convex CCW polygons
}
// Keep dot(x-origin,normal)<=0; using a nearby origin avoids small-cell cancellation.
inline TransferPolygon Transfer_Clip(const TransferPolygon& p,const TransferPoint& origin,const TransferPoint& normal)
{
    TransferPolygon q;if(p.empty())return q;
    TransferPoint a=p.back();real da=(a-origin).dot(normal);
    for(const TransferPoint& b:p){const real db=(b-origin).dot(normal);
        if((da<=0.)!=(db<=0.))q.push_back(a+(da/(da-db))*(b-a));
        if(db<=0.)q.push_back(b);
        a=b;da=db;
    }
    return q;
}
inline TransferPolygon Transfer_Intersection(TransferPolygon a,const TransferPolygon& b)
{
    for(size_t k=0;k<b.size()&&!a.empty();++k){const TransferPoint e=b[(k+1)%b.size()]-b[k];
        a=Transfer_Clip(a,b[k],TransferPoint(e.y(),-e.x()));}
    return a;
}
inline TransferPolygon Transfer_Rect(real x0,real y0,real x1,real y1)
{
    if(!(x1>x0&&y1>y0))throw std::invalid_argument("nonpositive transfer rectangle");
    return {TransferPoint(x0,y0),TransferPoint(x1,y0),TransferPoint(x1,y1),TransferPoint(x0,y1)};
}
struct TransferBounds {
    TransferPoint lo,hi;
    explicit TransferBounds(const TransferPolygon& p):lo(p.at(0)),hi(lo){for(const auto& x:p){lo=lo.cwiseMin(x);hi=hi.cwiseMax(x);}}
    bool Intersects(const TransferBounds& b) const{return (hi.array()>b.lo.array()).all()&&(b.hi.array()>lo.array()).all();}
};
// Voronoi cells clipped to a convex patch, with no periodic-image assumptions.
inline std::vector<TransferPolygon> Transfer_Voronoi(const std::vector<TransferPoint>& points,const TransferPolygon& patch)
{
    std::vector<TransferPolygon> out;out.reserve(points.size());
    for(size_t i=0;i<points.size();++i){auto p=patch;
        if(!points[i].allFinite())throw std::invalid_argument("nonfinite transfer point");
        for(size_t j=0;j<points.size()&&!p.empty();++j)if(i!=j){
            const TransferPoint d=points[j]-points[i];
            if(!(d.squaredNorm()>0.))throw std::invalid_argument("duplicate transfer point");
            p=Transfer_Clip(p,.5*(points[i]+points[j]),d);
        }
        if(!(Transfer_Area(p)>0.))throw std::invalid_argument("transfer point has no patch support");
        out.push_back(std::move(p));
    }
    return out;
}
inline void Validate_Transfer_Polygon(const TransferPolygon& p)
{
    if(p.size()<3||!(Transfer_Area(p)>0.))throw std::invalid_argument("empty or reversed transfer polygon");
    for(size_t k=0;k<p.size();++k){if(!p[k].allFinite())throw std::invalid_argument("nonfinite transfer polygon");
        const TransferPoint a=p[(k+1)%p.size()]-p[k],b=p[(k+2)%p.size()]-p[(k+1)%p.size()];
        if(Transfer_Cross(a,b)<-1e-12*a.norm()*b.norm())throw std::invalid_argument("nonconvex transfer polygon");}
}
inline void Validate_Transfer_Partition(const std::vector<TransferPolygon>& cells,const TransferPolygon& domain,real tol)
{
    Validate_Transfer_Polygon(domain);const real total=Transfer_Area(domain);real sum=0.;
    std::vector<TransferBounds> boxes;boxes.reserve(cells.size());
    for(const auto& p:cells){Validate_Transfer_Polygon(p);const real area=Transfer_Area(p);sum+=area;
        if(std::abs(Transfer_Area(Transfer_Intersection(p,domain))-area)>tol*area)
            throw std::invalid_argument("transfer cell outside domain");
        boxes.emplace_back(p);
    }
    if(std::abs(sum-total)>tol*total)throw std::invalid_argument("transfer partition area mismatch");
    for(size_t i=0;i<cells.size();++i)for(size_t j=0;j<i;++j)if(boxes[i].Intersects(boxes[j])){
        const real overlap=Transfer_Area(Transfer_Intersection(cells[i],cells[j]));
        if(overlap>tol*std::min(Transfer_Area(cells[i]),Transfer_Area(cells[j])))
            throw std::invalid_argument("overlapping transfer cells");
    }
}
// V is liquid volume; S is amount on ONE face; Q=2S+B. P is physical momentum.
// For constant liquid density, a shared nonnegative transfer of V and P preserves
// uniform velocity and cannot increase sum |P|^2/(2 rho V) (Jensen inequality).
struct TransferInventory {
    real V=0.,S=0.,B=0.;Vec3 P=Vec3::Zero();
    TransferInventory& operator+=(const TransferInventory& b){V+=b.V;S+=b.S;B+=b.B;P+=b.P;return *this;}
    TransferInventory Scaled(real w) const{TransferInventory r;r.V=w*V;r.S=w*S;r.B=w*B;r.P=w*P;return r;}
};
struct OverlapTransfer {
    struct Entry{int source,destination;real area,fraction;};
    std::vector<real> source_area,destination_area;
    std::vector<Entry> entries;
    real source_closure=0.,destination_closure=0.;
    std::vector<TransferInventory> Apply(const std::vector<TransferInventory>& input) const
    {
        if(input.size()!=source_area.size())throw std::invalid_argument("transfer inventory size");
        for(const auto& q:input)if(!(q.V>0.)||!(q.S>=0.)||!(q.B>=0.)||!std::isfinite(q.V)||!std::isfinite(q.S)||!std::isfinite(q.B)||!q.P.allFinite())
            throw std::invalid_argument("invalid transfer inventory");
        std::vector<TransferInventory> output(destination_area.size());
        for(const auto& e:entries)output[e.destination]+=input[e.source].Scaled(e.fraction);
        return output;
    }
};
inline OverlapTransfer Make_Overlap_Transfer(const std::vector<TransferPolygon>& source,
    const std::vector<TransferPolygon>& destination,const TransferPolygon& domain,real tolerance=1e-10)
{
    if(!(tolerance>0.&&tolerance<1.))throw std::invalid_argument("transfer tolerance");
    Validate_Transfer_Partition(source,domain,tolerance);Validate_Transfer_Partition(destination,domain,tolerance);
    OverlapTransfer r;std::vector<TransferBounds> sb,db;
    for(const auto& p:source){r.source_area.push_back(Transfer_Area(p));sb.emplace_back(p);}
    for(const auto& p:destination){r.destination_area.push_back(Transfer_Area(p));db.emplace_back(p);}
    std::vector<real> rows(source.size(),0.),cols(destination.size(),0.);
    for(size_t i=0;i<source.size();++i)for(size_t j=0;j<destination.size();++j)if(sb[i].Intersects(db[j])){
        const real a=Transfer_Area(Transfer_Intersection(source[i],destination[j]));
        if(a>0.){r.entries.push_back({int(i),int(j),a,a/r.source_area[i]});rows[i]+=a;cols[j]+=a;}
    }
    for(size_t i=0;i<source.size();++i)r.source_closure=std::max(r.source_closure,std::abs(rows[i]/r.source_area[i]-1.));
    for(size_t j=0;j<destination.size();++j)r.destination_closure=std::max(r.destination_closure,std::abs(cols[j]/r.destination_area[j]-1.));
    if(r.source_closure>tolerance||r.destination_closure>tolerance)throw std::runtime_error("unclosed geometric transfer; not renormalized");
    return r;
}
// Signed extensive inventory changing ownership a->b. This is the geometric/ALE
// component only, separate from relative drainage and material advection.
inline TransferInventory Transfer_Ownership(const OverlapTransfer& T,const std::vector<TransferInventory>& q,
    const std::vector<int>& before,const std::vector<int>& after,int a,int b)
{
    if(q.size()!=T.source_area.size()||before.size()!=q.size()||after.size()!=T.destination_area.size())
        throw std::invalid_argument("ownership transfer size");
    TransferInventory flux;
    for(const auto& e:T.entries){const int i=before[e.source],j=after[e.destination];
        if(i==a&&j==b)flux+=q[e.source].Scaled(e.fraction);
        else if(i==b&&j==a)flux+=q[e.source].Scaled(-e.fraction);}
    return flux;
}
} // namespace vf::thinfilm
