#include "vf/core/Check.h"
#include "vf/thinfilm/CurvedEdgePreparation.h"
#include "vf/thinfilm/CurvedEdgeMechanicalObserver.h"
using namespace vf;using namespace vf::thinfilm;
int main(){Set_Threads(2);FilmMatState a;const auto domain=Transfer_Rect(0,0,5,3);const real cut=2.75,dt=1e-4;
    for(int j=0;j<3;++j)for(real x:{.5,1.5,2.5,3.5,4.5,5.}){const int k=a.Size();a.x.emplace_back(x+(x==2.5||x==5?0:.015*std::sin(k)),j+.5+.015*std::cos(k),0);a.fixed.push_back(x==5);}
    a.v.assign(a.Size(),Vec3::Zero());FilmMechParams p;p.surf.Gamma_inf=1;p.surf.RT=.3;p.rho=.7;p.mu_s=.02;p.tf.mu=.2;
    const auto geom=Moving_Voronoi_Geometry(a,domain,p);for(int i=0;i<a.Size();++i){a.V.push_back(geom.graph.a[i]);a.S.push_back(.4*geom.graph.a[i]);}
    CurvedEdgeBoundary b;b.aspect=.2;b.clean_tension=1;b.height.assign(a.Size(),1.3);b.potential.assign(a.Size(),std::numeric_limits<real>::quiet_NaN());for(int i=0;i<a.Size();++i)if(a.fixed[i])b.potential[i]=-.01;
    Prepare_Curved_Edge(a,geom.graph,p,b,2.5,1.,.4);FilmMatState full=a,mid;const auto native=Curved_Edge_Coupled_Step(full,p,domain,dt,b,&mid);
    const auto x0=a.x,x1=mid.x;const auto o=Observe_Curved_Mechanics(a,mid,p,domain,b,dt,cut);
    const auto wa=Moving_Interface_Weights(geom,cut);auto complement=wa;for(auto& w:complement)w=1-w;
    const auto regional=Differentiate_Curved_Edge_Weighted(mid,p,domain,b,wa),left=Differentiate_Curved_Edge_Weighted(mid,p,domain,b,complement),global=Differentiate_Curved_Edge(mid,p,domain,b);
    real sumerror=std::abs(regional.energy+left.energy-global.energy);for(int i=0;i<a.Size();++i)sumerror=std::max(sumerror,(regional.gradient[i]+left.gradient[i]-global.gradient[i]).norm());
    VF_CHECK(sumerror<1e-12,"complementary weighted energy/gradient recovers common global energy");
    auto value=[&](const FilmMatState& s){const auto g=Moving_Voronoi_Geometry(s,domain,p);const auto f=CurvedEdgeEnergy(g.graph,p,s.S,b).Fields(s.V,false);real e=0;for(int i=0;i<s.Size();++i)e+=wa[i]*f.node_energy[i];return e;};
    real error=0,scale=0;const real eps=2e-6;for(int i=0;i<a.Size();++i)if(!a.fixed[i])for(int axis=0;axis<2;++axis){auto plus=mid,minus=mid;plus.x[i][axis]+=eps;minus.x[i][axis]-=eps;const real fd=(value(plus)-value(minus))/(2*eps);error=std::max(error,std::abs(fd-regional.gradient[i][axis]));scale=std::max(scale,std::abs(fd));}
    VF_CHECK(error<1e-7*scale,"regional energy derivative agrees with direct-energy finite differences: %.3e",error/scale);
    VF_CHECK(std::abs(o.viscous_total-native.mechanics.rayleigh_action)<1e-10*std::max(real(1e-20),o.viscous_total)&&std::abs(o.viscous_action_difference)<1e-10*std::max(real(1e-20),o.viscous_total),"scalar pyramid dissipation agrees with frozen native matrix action");
    VF_CHECK(o.regions[0].viscous_dissipation>=0&&o.regions[1].viscous_dissipation>=0&&std::abs(o.cut_length/3-1)<1e-9,"positive clipped regional dissipation and complete physical-cut length");
    real residual=0,workscale=0;for(int i=0;i<a.Size();++i){residual+=(o.inertial[i]+o.viscous[i]+o.capillary[i]).dot(dt*TransferPoint(mid.v[i].head<2>()));workscale+=std::abs(o.capillary[i].dot(dt*TransferPoint(mid.v[i].head<2>())));}
    VF_CHECK(std::abs(o.regions[0].residual+o.regions[1].residual-residual)<1e-9*workscale&&native.mechanics.newton.converged&&native.drainage.converged,"regional mechanical work closes against actual force residual");
    bool unchanged=a.x==x0&&mid.x==x1;VF_CHECK(unchanged&&a.V==mid.V&&a.S==mid.S,"observer leaves input states and inventories unchanged");
    return test::Report();}
