#pragma once
#include "vf/thinfilm/CurvedEdgeMotion.h"
#include "vf/thinfilm/CurvedEdgeRegional.h"
#include "vf/thinfilm/MovingInterfaceLedger.h"
namespace vf::thinfilm {
struct CurvedRegionMechanicalWork {
    real dE=0,dK=0,ownershipE=0,ownershipK=0,inertial_loss=0,viscous_dissipation=0;
    real inertial_work=0,viscous_node_work=0,capillary_node_work=0,capillary_local_work=0;
    real capillary_remainder=0,capillary_port=0,viscous_port=0,residual=0;
};
struct CurvedMechanicalObservation {
    std::vector<real> w0,w1,node_viscous;std::vector<TransferPoint> inertial,viscous,capillary,regional;
    CurvedEdgeFields f0,f1;std::array<CurvedRegionMechanicalWork,2> regions;
    real reconstructed_force=0,viscous_action_difference=0,coordinate_roundtrip=0,viscous_total=0,viscous_recovered_port=0,cut_length=0;
};
// Cut trace from the LEFT. The half-open rule avoids counting a shared triangle
// edge twice when the physical cut exactly passes through an old trace site.
inline std::pair<real,real> Mechanical_Triangle_Cut(const TransferPolygon& tri,real cut){
    real lo=std::numeric_limits<real>::infinity(),hi=-lo;bool left=false,notLeft=false;std::vector<real> y;
    for(auto x:tri){left|=x.x()<cut;notLeft|=x.x()>=cut;}
    if(!left||!notLeft)return {0,0};
    for(int k=0;k<3;++k){const auto a=tri[k],b=tri[(k+1)%3];if(a.x()==cut)y.push_back(a.y());
        if((a.x()<cut&&b.x()>cut)||(a.x()>cut&&b.x()<cut))y.push_back(a.y()+(b.y()-a.y())*(cut-a.x())/(b.x()-a.x()));}
    for(real x:y){lo=std::min(lo,x);hi=std::max(hi,x);}return hi>lo?std::make_pair(hi-lo,.5*(lo+hi)):std::make_pair(real(0),real(0));
}
// Independent scalar reconstruction of the inherited FULL face-pyramid metric.
// No matrix, force or solver state is modified. Dissipation is integrated over
// actual clipped old-geometry pyramids; recovered traction is a separate test.
inline void Observe_Curved_Viscosity(const FilmMatState& a,const FilmMatState& mid,const FilmMechParams& p,const VoronoiDrainageGeometry& geom,const CurvedEdgeBoundary& b,real dt,real cut,CurvedMechanicalObservation& out){
    const int n=a.Size();const auto& G=geom.graph;std::vector<Eigen::Matrix<real,2,2>> gradients(n,Eigen::Matrix<real,2,2>::Zero());
    if(!a.drag.empty()||!a.drag_n.empty())throw std::invalid_argument("observer requires inherited wall/full-gradient metric without extra drag");
    for(int i=0;i<n;++i)if(!b.Bath(i)){
        Eigen::Matrix<real,2,2> gram=Eigen::Matrix<real,2,2>::Zero(),rhs=gram;const TransferPoint xi=G.x[i].head<2>(),ui=mid.v[i].head<2>();
        for(const auto& face:geom.faces){int j=-2;if(face.i==i)j=face.j;else if(face.j==i)j=face.i;else continue;
            const bool dirichlet=j<0||a.fixed[j]||b.Bath(j);const TransferPoint d=dirichlet?TransferPoint(face.midpoint-xi):TransferPoint(G.x[j].head<2>()-xi),du=dirichlet?TransferPoint(-ui):TransferPoint(mid.v[j].head<2>()-ui);const real w=face.length/d.norm();gram.noalias()+=w*d*d.transpose();rhs.noalias()+=w*du*d.transpose();}
        gradients[i]=rhs*gram.inverse();
    }
    out.node_viscous.assign(n,0);real rightD=0;
    for(int i=0;i<n;++i)if(!b.Bath(i)){
        const TransferPoint xi=G.x[i].head<2>(),ui=mid.v[i].head<2>();const auto& g=gradients[i];
        for(const auto& face:geom.faces){int j=-2;TransferPoint normal=face.normal;if(face.i==i)j=face.j;else if(face.j==i){j=face.i;normal=-normal;}else continue;
            const TransferPoint dr=face.midpoint-xi;const real distance=dr.dot(normal);if(!(distance>0))throw std::runtime_error("observer nonpositive pyramid");TransferPoint residual;
            if(j>=0&&!a.fixed[j]&&!b.Bath(j)){const TransferPoint skew=face.midpoint-.5*(xi+TransferPoint(G.x[j].head<2>()));residual=.5*(TransferPoint(mid.v[j].head<2>())-ui)+.5*(g+gradients[j])*skew-g*dr;}
            else residual=-ui-g*dr;
            const Eigen::Matrix<real,2,2> H=g+residual*normal.transpose()/distance;const real div=H.trace(),beta=p.bulk_viscosity?1.:.5;
            const real strain=H(0,0)*H(0,0)+H(1,1)*H(1,1)+.5*std::pow(H(0,1)+H(1,0),2)+beta*div*div;
            const real area=.5*face.length*distance,density=2*p.Mu_S()*a.V[i]/G.a[i];const TransferPoint tangent(-normal.y(),normal.x());
            const TransferPolygon tri={xi,face.midpoint-.5*face.length*tangent,face.midpoint+.5*face.length*tangent};const real part=Transfer_Area(Transfer_Clip(tri,{cut,0},{-1.,0}));const real D=dt*density*area*strain;out.node_viscous[i]+=D;rightD+=dt*density*part*strain;
            const auto trace=Mechanical_Triangle_Cut(tri,cut);if(trace.first>0){const TransferPoint u=ui+H*(TransferPoint(cut,trace.second)-xi);const real tauXX=density*(H(0,0)+beta*div),tauYX=.5*density*(H(0,1)+H(1,0));out.viscous_recovered_port-=dt*trace.first*(u.x()*tauXX+u.y()*tauYX);out.cut_length+=trace.first;}
        }
    }
    for(real D:out.node_viscous)out.viscous_total+=D;out.regions[1].viscous_dissipation=rightD;out.regions[0].viscous_dissipation=out.viscous_total-rightD;
}
inline CurvedMechanicalObservation Observe_Curved_Mechanics(const FilmMatState& a,const FilmMatState& mid,const FilmMechParams& p,const TransferPolygon& domain,const CurvedEdgeBoundary& b,real dt,real cut){
    const int n=a.Size();if(a.V!=mid.V||a.S!=mid.S)throw std::invalid_argument("observer mechanical endpoint must precede drainage");
    CurvedMechanicalObservation out;const auto g0=Moving_Voronoi_Geometry(a,domain,p),g1=Moving_Voronoi_Geometry(mid,domain,p);out.w0=Moving_Interface_Weights(g0,cut);out.w1=Moving_Interface_Weights(g1,cut);
    out.f0=CurvedEdgeEnergy(g0.graph,p,a.S,b).Fields(a.V,false);out.f1=CurvedEdgeEnergy(g1.graph,p,mid.S,b).Fields(mid.V,false);
    const auto global=Differentiate_Curved_Edge(mid,p,domain,b),regional=Differentiate_Curved_Edge_Weighted(mid,p,domain,b,out.w0);
    MovingVoronoiProblem metric(a,p,domain,dt,&b.potential,{true,true});VecX q=VecX::Zero(metric.Dim()),mg;
    for(size_t k=0;k<metric.Free().size();++k){const int i=metric.Free()[k];q.segment<2>(2*k)=dt*mid.v[i].head<2>();out.coordinate_roundtrip=std::max(out.coordinate_roundtrip,(TransferPoint(mid.x[i].head<2>()-a.x[i].head<2>())-q.segment<2>(2*k)).norm());}
    metric.Metric_Gradient(q,mg);out.inertial.assign(n,TransferPoint::Zero());out.viscous=out.inertial;out.capillary=global.gradient;out.regional=regional.gradient;
    for(size_t k=0;k<metric.Free().size();++k){const int i=metric.Free()[k];out.inertial[i]=p.rho*a.V[i]*(TransferPoint(mid.v[i].head<2>()-a.v[i].head<2>()))/dt;out.viscous[i]=mg.segment<2>(2*k)-out.inertial[i];out.reconstructed_force=std::max(out.reconstructed_force,(out.inertial[i]+out.viscous[i]+out.capillary[i]).lpNorm<Eigen::Infinity>());}
    Observe_Curved_Viscosity(a,mid,p,g0,b,dt,cut,out);
    for(int region=0;region<2;++region){auto& r=out.regions[region];for(int i=0;i<n;++i){const real wa=region?out.w0[i]:1-out.w0[i],wb=region?out.w1[i]:1-out.w1[i];const TransferPoint v0=a.v[i].head<2>(),v1=mid.v[i].head<2>(),dx=dt*v1;
            const real k0=.5*p.rho*a.V[i]*v0.squaredNorm(),k1=.5*p.rho*a.V[i]*v1.squaredNorm();r.dK+=wb*k1-wa*k0;r.ownershipK+=(wb-wa)*k1;r.inertial_loss+=wa*.5*p.rho*a.V[i]*(v1-v0).squaredNorm();
            r.dE+=wb*out.f1.node_energy[i]-wa*out.f0.node_energy[i];r.ownershipE+=(wb-wa)*out.f1.node_energy[i];r.inertial_work+=wa*out.inertial[i].dot(dx);r.viscous_node_work+=wa*out.viscous[i].dot(dx);r.capillary_node_work+=wa*out.capillary[i].dot(dx);
            const TransferPoint local=region?out.regional[i]:TransferPoint(out.capillary[i]-out.regional[i]);r.capillary_local_work+=local.dot(dx);}
        r.capillary_remainder=r.capillary_local_work-(r.dE-r.ownershipE);r.capillary_port=r.capillary_local_work-r.capillary_node_work;r.viscous_port=r.viscous_dissipation-r.viscous_node_work;
        r.residual=r.dE+r.dK+r.inertial_loss+r.viscous_dissipation+r.capillary_remainder-r.ownershipE-r.ownershipK-r.capillary_port-r.viscous_port;
    }
    out.viscous_action_difference=out.viscous_total-out.regions[0].viscous_node_work-out.regions[1].viscous_node_work;return out;
}
} // namespace vf::thinfilm
