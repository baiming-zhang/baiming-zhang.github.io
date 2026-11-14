#include "vf/core/Check.h"
#include "vf/thinfilm/CurvedEdgeMotion.h"
using namespace vf;using namespace vf::thinfilm;
int main(){Set_Threads(2);const auto domain=Transfer_Rect(0,0,5,4);FilmMatState s;
    for(int j=0;j<4;++j)for(int i=0;i<5;++i){s.x.emplace_back(i+.5+.025*std::sin(i+2*j),j+.5+.025*std::cos(2*i+j),0);s.fixed.push_back(i==0||i==4||j==0||j==3);}
    s.v.assign(s.Size(),Vec3::Zero());FilmMechParams p;p.surf.Gamma_inf=1;p.surf.RT=.3;p.rho=.7;p.mu_s=.02;p.tf.mu=.2;
    const auto geom=Moving_Voronoi_Geometry(s,domain,p);for(int i=0;i<s.Size();++i){s.V.push_back(geom.graph.a[i]*(1.5+.07*std::sin(s.x[i].x())));s.S.push_back(geom.graph.a[i]*(.4+.008*std::cos(s.x[i].x())*std::sin(s.x[i].y())));}
    CurvedEdgeBoundary b;b.aspect=.2;b.clean_tension=1;b.height.assign(s.Size(),2.);b.potential.assign(s.Size(),std::numeric_limits<real>::quiet_NaN());for(int i=0;i<s.Size();++i)if(s.x[i].x()>4)b.potential[i]=-.02;
    CurvedEdgeEnergy E(geom.graph,p,s.S,b);const auto f=E.Fields(s.V);real graderror=0,hesserror=0,gradscale=0,hessscale=0;const real eps=2e-6;
    for(int i=0;i<s.Size();++i){auto a=s.V,c=s.V;a[i]+=eps;c[i]-=eps;const real fd=(E.Value(a)-E.Value(c))/(2*eps);graderror=std::max(graderror,std::abs(fd-f.pressure[i]));gradscale=std::max(gradscale,std::abs(fd));const auto fa=E.Fields(a,false),fc=E.Fields(c,false);
        for(int j=0;j<s.Size();++j){const real h=(fa.pressure[j]-fc.pressure[j])/(2*eps);hesserror=std::max(hesserror,std::abs(h-f.hessian(j,i)));hessscale=std::max(hessscale,std::abs(h));}}
    VF_CHECK(graderror<1e-7*gradscale,"curved surface energy analytic pressure matches volume finite differences: %.2e",graderror/gradscale);
    VF_CHECK(hesserror<1e-7*hessscale&&(f.hessian-f.hessian.transpose()).norm()<1e-15,"full dilution/capillary Hessian agrees with finite differences and symmetry");
    bool bath=true;for(int i=0;i<s.Size();++i)if(b.Bath(i))bath&=f.pressure[i]==b.potential[i]&&f.hessian.col(i).norm()==0.;VF_CHECK(bath,"far height and pressure are independent boundary data; bath volume never sets h");
    real mindilution=1.;for(int i=0;i<s.Size();++i)mindilution=std::min(mindilution,geom.graph.a[i]/f.area[i]);VF_CHECK(mindilution<1-1e-8&&mindilution>0,"nonflat state dilutes monolayer using curved surface area: area ratio %.9f",mindilution);
    const auto ad=Differentiate_Curved_Edge(s,p,domain,b);real xe=0,xscale=0;
    for(int i=0;i<s.Size();++i)for(int axis=0;axis<2;++axis){auto a=s,c=s;a.x[i][axis]+=eps;c.x[i][axis]-=eps;const auto ga=Moving_Voronoi_Geometry(a,domain,p),gc=Moving_Voronoi_Geometry(c,domain,p);const real fd=(CurvedEdgeEnergy(ga.graph,p,a.S,b).Value(a.V)-CurvedEdgeEnergy(gc.graph,p,c.S,b).Value(c.V))/(2*eps);xe=std::max(xe,std::abs(fd-ad.gradient[i][axis]));xscale=std::max(xscale,std::abs(fd));}
    VF_CHECK(std::abs(ad.energy-ad.double_energy)<1e-13&&xe<1e-7*xscale,"reverse curved geometry force matches independently evaluated energy: %.2e",xe/xscale);
    auto flat=s;auto closed=b;closed.potential.clear();closed.height.clear();for(int i=0;i<s.Size();++i){flat.S[i]=.4*geom.graph.a[i];flat.V[i]=1.5*geom.graph.a[i];}
    const auto zero=Differentiate_Curved_Edge(flat,p,domain,closed);real force=0;for(auto x:zero.gradient)force=std::max(force,x.norm());VF_CHECK(std::abs(zero.energy)<1e-12&&force<1e-12,"flat uniform monolayer retains zero reference energy and force");
    for(int i=0;i<s.Size();++i)flat.V[i]*=1+1e-3*std::sin(s.x[i].x());p.tf.sigma=2*(1+.3*std::log(.6))*b.aspect*b.aspect;
    const real curved=CurvedEdgeEnergy(geom.graph,p,flat.S,closed).Value(flat.V),quadratic=ThinFilmEnergy(geom.graph,p.tf).Value(flat.V);
    VF_CHECK(std::abs(curved/quadratic-1)<1e-6,"small-slope constant-coverage limit matches inherited capillary energy: %.2e",curved/quadratic-1);
    auto moving=s;FilmMatState mid;const auto result=Curved_Edge_Coupled_Step(moving,p,domain,.01,b,&mid);
    VF_CHECK(result.mechanics.newton.converged&&result.drainage.converged,"actual curved-area mechanics and common-state drainage converge");
    real move=0,sv=0,ss=0;for(int i=0;i<s.Size();++i){move=std::max(move,(moving.x[i]-s.x[i]).norm());sv+=moving.V[i]-s.V[i];ss+=moving.S[i]-s.S[i];}
    VF_CHECK(move>1e-8&&std::abs(sv)<1e-12&&ss==0,"nonzero motion conserves liquid including bath accounting and material S");
    const real k0=Moving_Voronoi_Kinetic(s,p.rho);VF_CHECK(result.energy_after+result.kinetic_after+result.mechanics.rayleigh_action+result.half_drain_action+result.momentum_mixing_loss<=result.energy_before+k0+1e-11,"common free energy and inherited dissipative metric satisfy global inequality");
    const auto g1=Moving_Voronoi_Geometry(moving,domain,p);auto inferred=s.V;for(int q=0;q<g1.graph.pairs.Size();++q){inferred[g1.graph.pairs.pi[q]]-=result.flux[q];inferred[g1.graph.pairs.pj[q]]+=result.flux[q];}real node=0;for(int i=0;i<s.Size();++i)node=std::max(node,std::abs(inferred[i]-moving.V[i]));VF_CHECK(node<1e-12,"every shared face uses exactly one signed liquid flux");
    // A curved physical concentration can be admissible even when the amount
    // divided by PROJECTED area is >Gamma_inf. No cancelled flat-energy call
    // may silently rule out that valid curved state.
    auto steep=b;for(auto& h:steep.height)h=20;auto high=s;const auto af=CurvedEdgeEnergy(geom.graph,p,s.S,steep).Fields(s.V,false);int probe=0;for(int i=1;i<s.Size();++i)if(af.area[i]/geom.graph.a[i]>af.area[probe]/geom.graph.a[probe])probe=i;
    high.S[probe]=.8*af.area[probe];CurvedEdgeMotionProblem admissible(high,p,domain,.01,steep);VecX zeroq=VecX::Zero(admissible.Dim()),gradient;admissible.Gradient(zeroq,gradient);
    VF_CHECK(high.S[probe]/geom.graph.a[probe]>1&&std::isfinite(admissible.Value(zeroq))&&gradient.allFinite(),"physical metric does not impose an invalid projected-area Langmuir domain");
    auto signedState=s;for(int i=0;i<s.Size();++i)if(b.Bath(i))signedState.V[i]=-2.;const auto signedFields=CurvedEdgeEnergy(geom.graph,p,signedState.S,b).Fields(signedState.V,false);real shift=0;for(int i=0;i<s.Size();++i)if(b.Bath(i))shift+=b.potential[i]*(signedState.V[i]-s.V[i]);
    const auto originalFields=CurvedEdgeEnergy(geom.graph,p,s.S,b).Fields(s.V,false);real invariance=0;for(int i=0;i<s.Size();++i)invariance=std::max(invariance,std::abs(signedFields.pressure[i]-originalFields.pressure[i])+std::abs(signedFields.area[i]-originalFields.area[i]));
    VF_CHECK(invariance==0&&std::abs(signedFields.energy-originalFields.energy-shift)<1e-13,"signed bath counter changes only prescribed-potential accounting, not physical h/area/pressure");
    auto positiveRun=s,negativeRun=signedState;const auto pr=Curved_Edge_Coupled_Step(positiveRun,p,domain,.01,b),nr=Curved_Edge_Coupled_Step(negativeRun,p,domain,.01,b);real trajectory=0;
    for(int i=0;i<s.Size();++i){trajectory=std::max(trajectory,(positiveRun.x[i]-negativeRun.x[i]).norm());if(!b.Bath(i))trajectory=std::max(trajectory,std::abs(positiveRun.V[i]-negativeRun.V[i]));}
    VF_CHECK(pr.mechanics.newton.converged&&nr.mechanics.newton.converged&&pr.drainage.converged&&nr.drainage.converged&&trajectory<1e-10,"physical coupled trajectory is independent of bath bookkeeping origin: %.2e",trajectory);
    return vf::test::Report();
}
