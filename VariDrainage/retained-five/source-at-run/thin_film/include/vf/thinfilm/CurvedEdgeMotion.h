#pragma once
#include "vf/thinfilm/CurvedEdge.h"
namespace vf::thinfilm {
// Reuse the existing physical inertial/viscous metric. Replace ONLY its
// variational free energy by the common curved-area/monolayer energy.
class CurvedEdgeMotionProblem:public Problem {
public:
    CurvedEdgeMotionProblem(const FilmMatState& s,const FilmMechParams& p,const TransferPolygon& domain,real dt,const CurvedEdgeBoundary& b):base_(s,p,domain,dt,b.potential.empty()?nullptr:&b.potential,{true,true}),p_(p),domain_(domain),b_(b){}
    int Dim() const override{return base_.Dim();}
    FilmMatState Trial(const VecX& q) const{return base_.Trial(q);}
    const std::vector<int>& Free() const{return base_.Free();}
    real Length() const{return base_.Length();}real Rayleigh(const VecX& q) const{return base_.Rayleigh(q);}real Kinetic(const FilmMatState& s) const{return base_.Kinetic(s);}
    real Value(const VecX& q) override{try{return base_.Metric_Value(q)+Differentiate_Curved_Edge(Trial(q),p_,domain_,b_).energy;}catch(const std::exception&){return std::numeric_limits<real>::infinity();}}
    void Gradient(const VecX& q,VecX& g) override{base_.Metric_Gradient(q,g);const auto b=Evaluate(Trial(q),"motion gradient");for(size_t k=0;k<Free().size();++k)g.segment<2>(2*k)+=b.gradient[Free()[k]];}
    void Prepare_Hessian(const VecX& q) override{prepared_=q;}
    void Hess_Vec(const VecX& v,VecX& h) override{
        base_.Metric_Hess_Vec(v,h);const real norm=v.lpNorm<Eigen::Infinity>();if(norm==0)return;const real eps=1e-5*Length()/norm;
        const auto plus=Trial(prepared_+eps*v),minus=Trial(prepared_-eps*v);
        const auto na=Evaluate(plus,"motion Hessian plus"),nb=Evaluate(minus,"motion Hessian minus");
        for(size_t k=0;k<Free().size();++k){const int i=Free()[k];h.segment<2>(2*k)+=(na.gradient[i]-nb.gradient[i])/(2*eps);}}
    void Hess_Diag(VecX& d) override{base_.Hess_Diag(d);}
private:
    VoronoiEnergyDerivative Evaluate(const FilmMatState& s,const char* stage){try{return Differentiate_Curved_Edge(s,p_,domain_,b_);}catch(const std::exception& e){throw std::runtime_error(std::string(stage)+": "+e.what());}}
    MovingVoronoiProblem base_;FilmMechParams p_;TransferPolygon domain_;CurvedEdgeBoundary b_;VecX prepared_;
};
inline MovingVoronoiStats Curved_Edge_Mechanics(FilmMatState& s,const FilmMechParams& p,const TransferPolygon& domain,real dt,const CurvedEdgeBoundary& b){
    CurvedEdgeMotionProblem problem(s,p,domain,dt,b);VecX q=VecX::Zero(problem.Dim()),g;problem.Gradient(q,g);MovingVoronoiStats r;
    r.initial_gradient_inf=g.lpNorm<Eigen::Infinity>();r.energy_before=Differentiate_Curved_Edge(s,p,domain,b).energy;r.kinetic_before=problem.Kinetic(s);const real j0=problem.Value(q);
    NewtonOptions options;options.max_iters=35;options.max_cg_iters=100;options.forcing_max=.05;options.grad_tol=std::max(real(1e-10),r.initial_gradient_inf*1e-8);options.step_tol=1e-10*problem.Length();options.max_total_step=.2*problem.Length();options.sparse_preconditioner=false;
    r.newton=Newton_CG_Minimize(problem,q,options);r.maximum_displacement=q.lpNorm<Eigen::Infinity>();r.incremental_change=problem.Value(q)-j0;problem.Gradient(q,g);r.final_gradient_inf=g.lpNorm<Eigen::Infinity>();
    auto final=problem.Trial(q);for(auto& v:final.v)v.setZero();for(size_t k=0;k<problem.Free().size();++k)final.v[problem.Free()[k]].head<2>()=q.segment<2>(2*k)/dt;
    r.energy_after=Differentiate_Curved_Edge(final,p,domain,b).energy;r.kinetic_after=problem.Kinetic(final);r.rayleigh_action=2*dt*problem.Rayleigh(q);s=std::move(final);s.has_prev=false;return r;
}
inline MovingVoronoiCoupledStats Curved_Edge_Coupled_Step(FilmMatState& s,const FilmMechParams& p,const TransferPolygon& domain,real dt,const CurvedEdgeBoundary& b,FilmMatState* mechanical_endpoint=nullptr){
    MovingVoronoiCoupledStats r;r.energy_before=Differentiate_Curved_Edge(s,p,domain,b).energy;const auto oldV=s.V;r.mechanics=Curved_Edge_Mechanics(s,p,domain,dt,b);
    if(!r.mechanics.newton.converged||r.mechanics.newton.hit_total_step)return r;if(mechanical_endpoint)*mechanical_endpoint=s;
    const auto geom=Moving_Voronoi_Geometry(s,domain,p);const auto& G=geom.graph;CurvedEdgeEnergy energy(G,p,s.S,b);std::vector<real> k(G.pairs.Size());
    for(int q=0;q<G.pairs.Size();++q){const int i=G.pairs.pi[q],j=G.pairs.pj[q];
        // Two prescribed-potential boundary nodes do not represent a resolved
        // bath flow link. Their counters must not define negative mobilities.
        if(b.Bath(i)&&b.Bath(j)){k[q]=0;continue;}
        real h=.5*(s.V[i]/G.a[i]+s.V[j]/G.a[j]);if(b.Bath(i))h=s.V[j]/G.a[j];if(b.Bath(j))h=s.V[i]/G.a[i];k[q]=p.tf.Mobility(h)*G.c[q];}
    NewtonOptions options;options.max_cg_iters=5000;r.drainage=H1_Flux_Step(energy,s.V,k,dt,options,1e-11,false,true,0.,nullptr,3,&r.flux);
    if(!r.drainage.converged)return r;for(int i=0;i<s.Size();++i)if(!std::isfinite(s.V[i])||(!b.Bath(i)&&!(s.V[i]>0)))return r;
    for(int q=0;q<G.pairs.Size();++q)if(k[q]>0)r.half_drain_action+=r.flux[q]*r.flux[q]/(2*dt*k[q]);
    real trueK=0;for(int i=0;i<s.Size();++i){trueK+=.5*p.rho*oldV[i]*s.v[i].squaredNorm();r.momentum_before_drain+=p.rho*oldV[i]*s.v[i];}
    r.pinned_impulse=Moving_Voronoi_Momentum_Transport(G,r.flux,oldV,s,p.rho,b.potential.empty()?nullptr:&b.potential);r.kinetic_after=Moving_Voronoi_Kinetic(s,p.rho);r.momentum_mixing_loss=trueK-r.kinetic_after;
    r.energy_after=energy.Value(s.V);const auto fields=energy.Fields(s.V,false);
    for(int i=0;i<s.Size();++i){r.momentum_after_drain+=p.rho*s.V[i]*s.v[i];if(b.Bath(i)){const real dv=s.V[i]-oldV[i];r.reservoir_volume_change+=dv;r.reservoir_work-=b.potential[i]*dv;}const real theta=s.S[i]/(p.surf.Gamma_inf*fields.area[i]);r.minimum_theta=std::min(r.minimum_theta,theta);r.maximum_theta=std::max(r.maximum_theta,theta);}
    r.film_energy_change=r.energy_after-r.energy_before+r.reservoir_work;return r;
}
} // namespace vf::thinfilm
