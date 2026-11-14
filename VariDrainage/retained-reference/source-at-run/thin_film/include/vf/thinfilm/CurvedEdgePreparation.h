#pragma once
#include "vf/thinfilm/CurvedEdge.h"
namespace vf::thinfilm {
struct CurvedEdgePreparationStats {int iterations=0;real pressure_residual=0,free_volume=0,surface_amount=0;};
// Initial equilibrium-activity capillary preparation, not an evolving closure.
// Fixed physical trace h at x=cut and fixed far h. Only free edge heights are
// solved; initial trace constraints are not added to mechanics/drainage.
// Surface inventory is theta*Gamma_inf*ACTUAL area, explicitly a new preparation.
inline CurvedEdgePreparationStats Prepare_Curved_Edge(FilmMatState& s,const ThinFilmGraph& G,const FilmMechParams& p,const CurvedEdgeBoundary& b,real cut,real bulk_h,real coverage){
    if(!(bulk_h>0&&coverage>0&&coverage<1)||s.Size()!=G.N)throw std::invalid_argument("edge preparation inputs");
    const int n=G.N;const real tension=b.clean_tension+p.surf.RT*p.surf.Gamma_inf*std::log1p(-coverage);
    if(!(tension>0))throw std::invalid_argument("edge preparation tension");
    std::vector<int> free;std::vector<real> h(n,bulk_h);real mu=0;int baths=0,traces=0;
    const real tol=1e-10*std::max(real(1),std::abs(cut)); // coordinate comparison only
    for(int i=0;i<n;++i){if(b.Bath(i)){h[i]=b.height.at(i);if(baths&&std::abs(mu-b.potential[i])>1e-12)throw std::invalid_argument("nonuniform initial bath potential");mu=b.potential[i];++baths;}
        else if(G.x[i].x()>cut+tol)free.push_back(i);else if(std::abs(G.x[i].x()-cut)<=tol)++traces;}
    if(!baths||!traces||free.empty()||mu==0)throw std::invalid_argument("edge preparation requires exact physical traces");
    const int m=int(free.size());std::vector<int> map(n,-1);for(int k=0;k<m;++k)map[free[k]]=k;
    // Newton starting guess only. No interior profile is prescribed or retained.
    real far_x=cut,far_h=bulk_h;for(int i=0;i<n;++i)if(b.Bath(i)&&G.x[i].x()>far_x){far_x=G.x[i].x();far_h=h[i];}
    for(int i:free)h[i]=bulk_h+(far_h-bulk_h)*(G.x[i].x()-cut)/(far_x-cut);
    auto evaluate=[&](const std::vector<real>& v,Eigen::VectorXd* gradient,Eigen::MatrixXd* Hessian){
        real E=0;if(gradient)*gradient=Eigen::VectorXd::Zero(m);if(Hessian)*Hessian=Eigen::MatrixXd::Zero(m,m);
        for(int q=0;q<G.pairs.Size();++q){const int i=G.pairs.pi[q],j=G.pairs.pj[q];const real distance=(G.x[i]-G.x[j]).norm(),length=G.c[q]*distance,z=.5*b.aspect*(v[i]-v[j])/distance,f=std::hypot(1.,z);
            E+=2*tension*length*distance*z*z/(f+1);const real force=tension*length*b.aspect*z/f,k=tension*length*b.aspect*b.aspect/(2*distance*f*f*f);const int a=map[i],c=map[j];
            if(gradient){if(a>=0)(*gradient)[a]+=force;if(c>=0)(*gradient)[c]-=force;}
            if(Hessian){if(a>=0)(*Hessian)(a,a)+=k;if(c>=0)(*Hessian)(c,c)+=k;if(a>=0&&c>=0){(*Hessian)(a,c)-=k;(*Hessian)(c,a)-=k;}}}
        for(int k=0;k<m;++k){const int i=free[k];E-=mu*G.a[i]*v[i];if(gradient)(*gradient)[k]-=mu*G.a[i];}return E;};
    CurvedEdgePreparationStats r;
    for(int iteration=0;iteration<50;++iteration){Eigen::VectorXd g;Eigen::MatrixXd H;const real E=evaluate(h,&g,&H);real residual=0;for(int k=0;k<m;++k)residual=std::max(residual,std::abs(g[k]/G.a[free[k]]/mu));r.iterations=iteration;r.pressure_residual=residual;
        if(residual<1e-9)break;const Eigen::VectorXd step=H.ldlt().solve(-g);if(!step.allFinite())throw std::runtime_error("initial edge tangent solve");real alpha=1;bool accepted=false;
        for(int ls=0;ls<40;++ls){auto trial=h;for(int k=0;k<m;++k)trial[free[k]]+=alpha*step[k];if(*std::min_element(trial.begin(),trial.end())>0&&evaluate(trial,nullptr,nullptr)<=E+1e-4*alpha*g.dot(step)){h=std::move(trial);accepted=true;break;}alpha*=.5;}
        if(!accepted)throw std::runtime_error("initial edge line search");}
    if(r.pressure_residual>=1e-8)throw std::runtime_error("initial edge not stationary");
    std::vector<real> area=G.a;for(int q=0;q<G.pairs.Size();++q){const int i=G.pairs.pi[q],j=G.pairs.pj[q];const real d=(G.x[i]-G.x[j]).norm(),z=.5*b.aspect*(h[i]-h[j])/d,e=G.c[q]*d*d*z*z/(std::hypot(1.,z)+1);area[i]+=.5*e;area[j]+=.5*e;}
    for(int i=0;i<n;++i){if(!b.Bath(i))s.V[i]=h[i]*G.a[i];s.S[i]=coverage*p.surf.Gamma_inf*area[i];r.surface_amount+=s.S[i];if(map[i]>=0)r.free_volume+=s.V[i];}
    // Validate against the variable-S model actually used after preparation.
    const auto fields=CurvedEdgeEnergy(G,p,s.S,b).Fields(s.V,false);real residual=0;for(int i:free)residual=std::max(residual,std::abs((fields.pressure[i]-mu)/mu));
    if(residual>1e-8)throw std::runtime_error("prepared variable-S pressure mismatch");r.pressure_residual=std::max(r.pressure_residual,residual);return r;
}
} // namespace vf::thinfilm
