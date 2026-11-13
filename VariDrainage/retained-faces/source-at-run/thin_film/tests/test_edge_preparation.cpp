#include "vf/core/Check.h"
#include "vf/thinfilm/CurvedEdgePreparation.h"
#include "vf/thinfilm/CurvedEdgeMotion.h"
using namespace vf;using namespace vf::thinfilm;
int main(){Set_Threads(2);FilmMatState s;const auto domain=Transfer_Rect(0,0,5,2);
    for(int j=0;j<2;++j)for(real x:{.5,1.5,2.5,3.5,4.5,5.}){s.x.emplace_back(x,j+.5,0);s.fixed.push_back(x==5.);}
    s.v.assign(s.Size(),Vec3::Zero());FilmMechParams p;p.surf.Gamma_inf=1;p.surf.RT=.3;p.rho=.7;p.mu_s=.02;p.tf.mu=.2;
    const auto geom=Moving_Voronoi_Geometry(s,domain,p);for(int i=0;i<s.Size();++i){s.V.push_back(geom.graph.a[i]);s.S.push_back(.4*geom.graph.a[i]);}
    CurvedEdgeBoundary b;b.aspect=.2;b.clean_tension=1;b.height.assign(s.Size(),1.3);b.potential.assign(s.Size(),std::numeric_limits<real>::quiet_NaN());for(int i=0;i<s.Size();++i)if(s.fixed[i])b.potential[i]=-.01;
    const auto old=s;auto shifted=s;for(int i=0;i<s.Size();++i)if(b.Bath(i))shifted.V[i]-=1e6;
    const auto prep=Prepare_Curved_Edge(s,geom.graph,p,b,2.5,1.,.4);Prepare_Curved_Edge(shifted,geom.graph,p,b,2.5,1.,.4);
    VF_CHECK(prep.pressure_residual<1e-8,"prepared pressure matches actual variable-S energy: %.3e",prep.pressure_residual);
    const auto f=CurvedEdgeEnergy(geom.graph,p,s.S,b).Fields(s.V,false);real error=0,nonlinear=0;bool trace=true;
    for(int i=0;i<s.Size();++i){error=std::max(error,std::abs(s.S[i]/f.area[i]-.4));if(std::abs(s.x[i].x()-2.5)<1e-10)trace&=!s.fixed[i]&&s.V[i]==old.V[i];if(s.x[i].x()>2.5&&!b.Bath(i))nonlinear=std::max(nonlinear,std::abs(s.V[i]/geom.graph.a[i]-(1+.3*(s.x[i].x()-2.5)/2.5)));}
    VF_CHECK(error<1e-14&&nonlinear>1e-3&&trace,"actual-area coverage; solved interior; initial trace remains mobile");
    bool gauge=true;for(int i=0;i<s.Size();++i)gauge&=s.S[i]==shifted.S[i]&&(b.Bath(i)||s.V[i]==shifted.V[i]);VF_CHECK(gauge,"preparation independent of signed bath counter origin");
    auto moving=s;const auto step=Curved_Edge_Coupled_Step(moving,p,domain,1e-4,b);real dv=0,ds=0,move=0;
    for(int i=0;i<s.Size();++i){dv+=moving.V[i]-s.V[i];ds+=moving.S[i]-s.S[i];move=std::max(move,(moving.x[i]-s.x[i]).norm());}
    VF_CHECK(step.mechanics.newton.converged&&step.drainage.converged&&std::abs(dv)<1e-12&&ds==0&&move>0,"released prepared state evolves with conserved V/S");
    return test::Report();}
