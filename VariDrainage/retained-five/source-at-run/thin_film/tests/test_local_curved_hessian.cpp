#include "vf/core/Check.h"
#include "vf/thinfilm/CurvedEdgeLocalAssembly.h"
using namespace vf;using namespace vf::thinfilm;
int main(){Set_Threads(2);const auto domain=Transfer_Rect(0,0,5,4);FilmMatState state;for(int j=0;j<4;++j)for(int i=0;i<5;++i)state.x.emplace_back(i+.5+.025*std::sin(i+2*j),j+.5+.025*std::cos(2*i+j),0);
    FilmMechParams base;base.surf.Gamma_inf=1;base.surf.RT=.3;base.rho=.7;base.mu_s=.02;base.tf.mu=.2;const auto G=Moving_Voronoi_Geometry(state,domain,base).graph;
    for(int variant=0;variant<5;++variant){auto p=base;CurvedEdgeBoundary b;b.aspect=.2;b.clean_tension=1;b.height.assign(G.N,variant==4?20.:2.);b.potential.assign(G.N,std::numeric_limits<real>::quiet_NaN());std::vector<real> v(G.N),s(G.N);
        for(int i=0;i<G.N;++i){v[i]=G.a[i]*(1.5+.07*std::sin(G.x[i].x()));s[i]=G.a[i]*(.4+.008*std::cos(G.x[i].x())*std::sin(G.x[i].y()));if(G.x[i].x()>4)b.potential[i]=-.02;}
        if(variant==1){b.potential.clear();b.height.clear();for(int i=0;i<G.N;++i){v[i]=1.5*G.a[i];s[i]=.4*G.a[i];}}
        if(variant==2)for(int i=0;i<G.N;++i)if(b.Bath(i))v[i]=-2.;
        if(variant==3){p.tf.disj.eps=.007;p.tf.disj.hb=.7;p.tf.disj.c8=.003;p.tf.disj.c2=.002;p.tf.disj.beta=.01;p.tf.disj.lambda=.6;p.tf.rho_g=.04;}
        if(variant==4){const auto f=CurvedEdgeEnergy(G,p,s,b).Fields(v,false);int max=0;for(int i=1;i<G.N;++i)if(f.area[i]/G.a[i]>f.area[max]/G.a[max])max=i;s[max]=.8*f.area[max];}
        CurvedEdgeEnergy ref(G,p,s,b);CurvedEdgeLocalAssembly candidate(G,p,s,b);const auto a=ref.Fields(v),c=candidate.Fields(v);const real scale=std::max(real(1e-30),a.hessian.cwiseAbs().maxCoeff());
        VF_CHECK(a.energy==c.energy&&a.pressure==c.pressure&&a.area==c.area&&a.node_energy==c.node_energy&&a.tension==c.tension,"variant%d first fields exact",variant);
        VF_CHECK((a.hessian-c.hessian).cwiseAbs().maxCoeff()/scale<=2e-13,"variant%d all Hessian entries",variant);
        bool bath=true;for(int i=0;i<G.N;++i)if(b.Bath(i))bath&=c.hessian.col(i).norm()==0;VF_CHECK(bath,"variant%d signed bath columns zero",variant);
        real action=0,diag=0;std::vector<real> w(G.N),out;
        for(int k=0;k<3;++k){for(int i=0;i<G.N;++i)w[i]=std::sin((.17+.2*k)*(i+1));candidate.Hess_Apply(v,w,out);const Eigen::VectorXd expected=a.hessian*Eigen::Map<const Eigen::VectorXd>(w.data(),w.size());real error=0,norm=0;for(int i=0;i<G.N;++i){error=std::max(error,std::abs(out[i]-expected[i]));norm=std::max(norm,std::abs(expected[i]));}action=std::max(action,error/std::max(real(1e-30),norm));}
        VF_CHECK(action<=2e-13,"variant%d three Hessian actions",variant);for(int i=0;i<G.N;++i)diag=std::max(diag,std::abs(candidate.Hess_Diag(v,i)-a.hessian(i,i))/scale);VF_CHECK(diag<=2e-13,"variant%d all diagonal entries",variant);
        for(int i=0;i<G.N;++i)if(!b.Bath(i))v[i]*=1+.002*std::sin(.3*(i+1));const auto changed=ref.Fields(v);candidate.Hess_Apply(v,w,out);const Eigen::VectorXd expected=changed.hessian*Eigen::Map<const Eigen::VectorXd>(w.data(),w.size());real error=0,norm=0;for(int i=0;i<G.N;++i){error=std::max(error,std::abs(out[i]-expected[i]));norm=std::max(norm,std::abs(expected[i]));}VF_CHECK(error/std::max(real(1e-30),norm)<=2e-13,"variant%d volume cache invalidates",variant);
        if(variant==0){auto invalid=v;invalid[0]=0.;bool oldthrow=false,newthrow=false;try{ref.Fields(invalid);}catch(const std::exception&){oldthrow=true;}try{candidate.Fields(invalid);}catch(const std::exception&){newthrow=true;}VF_CHECK(oldthrow&&newthrow,"invalid physical volume rejected by both");auto bad=b;bad.aspect=0.;oldthrow=newthrow=false;try{CurvedEdgeEnergy x(G,p,s,bad);}catch(const std::exception&){oldthrow=true;}try{CurvedEdgeLocalAssembly x(G,p,s,bad);}catch(const std::exception&){newthrow=true;}VF_CHECK(oldthrow&&newthrow,"invalid physical units rejected by both");}
    }
    return test::Report();}
