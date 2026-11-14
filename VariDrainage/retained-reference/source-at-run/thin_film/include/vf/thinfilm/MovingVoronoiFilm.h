// First-order inertial/viscous, insoluble planar particle-film prototype.
// Actual geometry derivatives, then conservative relative drainage/momentum.
// This is opt-in DRIFT development, not the default kernel mechanical solver.
#pragma once
#include "vf/thinfilm/VoronoiDifferentiation.h"
namespace vf::thinfilm {
struct MovingVoronoiDiscretization {
    bool full_velocity_gradient=false; // numerical choice, not a physical coefficient
    bool wall_velocity_samples=false; // Dirichlet samples at actual outer faces, not a mesh-wide pinned band
};
inline std::vector<TransferPoint> Moving_Voronoi_Sites(const FilmMatState& s){std::vector<TransferPoint> p;for(auto x:s.x)p.emplace_back(x.x(),x.y());return p;}
inline VoronoiDrainageGeometry Moving_Voronoi_Geometry(const FilmMatState& s,const TransferPolygon& domain,const FilmMechParams& p){
    auto g=Make_Voronoi_Drainage(Moving_Voronoi_Sites(s),domain);for(int i=0;i<s.Size();++i)g.graph.height[i]=p.up.dot(s.x[i]);return g;}
// Geometric area derivative at a shared face. Sum of these derivatives is zero
// for a fixed outer domain; wall-face positions themselves have zero derivative.
inline Eigen::MatrixXd Moving_Voronoi_Area_Jacobian(const VoronoiDrainageGeometry& g){
    const int n=g.graph.N;Eigen::MatrixXd J=Eigen::MatrixXd::Zero(n,2*n);
    for(const auto& f:g.faces)if(f.j>=0){const int i=f.i,j=f.j;const TransferPoint xi(g.graph.x[i].x(),g.graph.x[i].y()),xj(g.graph.x[j].x(),g.graph.x[j].y());
        const real ld=f.length/(xj-xi).norm();const TransferPoint di=ld*(f.midpoint-xi),dj=ld*(xj-f.midpoint);
        for(int k=0;k<2;++k){J(i,2*i+k)+=di[k];J(i,2*j+k)+=dj[k];J(j,2*i+k)-=di[k];J(j,2*j+k)-=dj[k];}}
    return J;
}
class MovingVoronoiProblem : public Problem {
public:
    MovingVoronoiProblem(const FilmMatState& initial,const FilmMechParams& p,const TransferPolygon& domain,real dt,const std::vector<real>* reservoir=nullptr,MovingVoronoiDiscretization discretization={})
        : old_(initial),prm_(p),domain_(domain),dt_(dt),res_(reservoir){
        const int n=old_.Size();if(!(dt>0.)||int(old_.v.size())!=n||int(old_.fixed.size())!=n||!(p.rho>0.))throw std::invalid_argument("moving prototype requires positive physical inertia and explicit fixed mask");
        id_.assign(n,-1);for(int i=0;i<n;++i)if(!old_.fixed[i]){id_[i]=int(free_.size());free_.push_back(i);}
        if(free_.empty())throw std::invalid_argument("moving prototype has no free particles");
        const auto geom=Moving_Voronoi_Geometry(old_,domain_,prm_);const auto& G=geom.graph;C_=Eigen::MatrixXd::Zero(Dim(),Dim());mass_.resize(Dim());velocity_.resize(Dim());
        auto bath=[&](int i){return res_&&!std::isnan((*res_)[i]);};
        auto row=[&](int i,TransferPoint v){VecX b=VecX::Zero(Dim());if(id_[i]>=0)b.segment<2>(2*id_[i])=v;return b;};
        if(!discretization.full_velocity_gradient){
        for(int q=0;q<G.pairs.Size();++q){const int i=G.pairs.pi[q],j=G.pairs.pj[q];const Vec3 d=G.x[j]-G.x[i];const TransferPoint e(d.x()/d.norm(),d.y()/d.norm());
            real hi=old_.V[i]/G.a[i],hj=old_.V[j]/G.a[j];if(bath(i)&&!bath(j))hi=hj;if(bath(j)&&!bath(i))hj=hi;
            const real cT=2.*p.Mu_S()*.5*(hi+hj)*G.c[q];const VecX b=row(j,e)-row(i,e);C_.noalias()+=2*cT*(b*b.transpose());}
        if(p.bulk_viscosity){const auto J=Moving_Voronoi_Area_Jacobian(geom);
            for(int i=0;i<n;++i)if(!bath(i)){VecX b(Dim());for(size_t k=0;k<free_.size();++k)b.segment<2>(2*k)=J.block<1,2>(i,2*free_[k]).transpose()/G.a[i];C_.noalias()+=p.Mu_S()*old_.V[i]*(b*b.transpose());}}
        }else{
            // Full symmetric gradient on face pyramids. The LS reconstruction
            // is affine exact; a shared face correction controls alternating
            // cell velocities that a centered cell gradient cannot see.
            // Integrate with the geometric pyramid area (no penalty coefficient).
            auto addOuter=[&](const VecX& b,real w){std::vector<int> nz;for(int k=0;k<b.size();++k)if(b[k]!=0.)nz.push_back(k);
                for(int j:nz)for(int i:nz)C_(i,j)+=w*b[i]*b[j];};
            using GradientRows=std::array<VecX,4>;
            std::vector<GradientRows> gradient(n);
            for(int i=0;i<n;++i){for(auto& b:gradient[i])b=VecX::Zero(Dim());if(bath(i))continue;
                struct Sample{int other;TransferPoint offset;real weight;};std::vector<Sample> samples;
                Eigen::Matrix<real,2,2> gram=Eigen::Matrix<real,2,2>::Zero();const TransferPoint xi=G.x[i].head<2>();
                for(const auto& face:geom.faces){int j=-2;if(face.i==i)j=face.j;else if(face.j==i)j=face.i;else continue;
                    if(j<0&&!discretization.wall_velocity_samples)continue;
                    const bool dirichlet=j<0||id_[j]<0||bath(j);
                    const TransferPoint d=dirichlet?TransferPoint(face.midpoint-xi):TransferPoint(G.x[j].head<2>()-xi);
                    const real w=face.length/d.norm();samples.push_back({dirichlet?-1:j,d,w});gram.noalias()+=w*(d*d.transpose());}
                if(!(gram.determinant()>128*std::numeric_limits<real>::epsilon()*gram.trace()*gram.trace()))throw std::runtime_error("rank-deficient velocity-gradient support");
                const auto inverse=gram.inverse().eval();auto& b=gradient[i];
                for(const auto& sample:samples){const TransferPoint coefficient=inverse*(sample.weight*sample.offset);
                    auto add=[&](int particle,real sign){if(particle<0||id_[particle]<0)return;const int k=2*id_[particle];b[0][k]+=sign*coefficient.x();b[1][k]+=sign*coefficient.y();b[2][k+1]+=sign*coefficient.x();b[3][k+1]+=sign*coefficient.y();};
                    add(sample.other,1.);add(i,-1.);}
            }
            for(int i=0;i<n;++i)if(!bath(i)){const TransferPoint xi=G.x[i].head<2>();const auto& g=gradient[i];
                for(const auto& face:geom.faces){int j=-2;TransferPoint normal=face.normal;if(face.i==i)j=face.j;else if(face.j==i){j=face.i;normal=-normal;}else continue;
                    const TransferPoint dr=face.midpoint-xi;const real distance=dr.dot(normal);
                    if(!(distance>0.))throw std::runtime_error("nonpositive Voronoi face pyramid");
                    GradientRows h=g;VecX residualX=VecX::Zero(Dim()),residualY=residualX;
                    if(j>=0&&id_[j]>=0&&!bath(j)){
                        const TransferPoint skew=face.midpoint-.5*(xi+TransferPoint(G.x[j].head<2>()));
                        residualX=.5*(row(j,{1.,0.})-row(i,{1.,0.}))+.5*skew.x()*(g[0]+gradient[j][0])+.5*skew.y()*(g[1]+gradient[j][1])-dr.x()*g[0]-dr.y()*g[1];
                        residualY=.5*(row(j,{0.,1.})-row(i,{0.,1.}))+.5*skew.x()*(g[2]+gradient[j][2])+.5*skew.y()*(g[3]+gradient[j][3])-dr.x()*g[2]-dr.y()*g[3];
                    }else if(j>=0||discretization.wall_velocity_samples){
                        residualX=-row(i,{1.,0.})-dr.x()*g[0]-dr.y()*g[1];residualY=-row(i,{0.,1.})-dr.x()*g[2]-dr.y()*g[3];}
                    h[0]+=normal.x()/distance*residualX;h[1]+=normal.y()/distance*residualX;
                    h[2]+=normal.x()/distance*residualY;h[3]+=normal.y()/distance*residualY;
                    const real pyramidArea=.5*face.length*distance,w=2*p.Mu_S()*old_.V[i]*pyramidArea/G.a[i];
                    addOuter(h[0],w);addOuter(h[3],w);addOuter(h[1]+h[2],.5*w);addOuter(h[0]+h[3],(p.bulk_viscosity?1.:.5)*w);
                }
            }
        }
        for(size_t k=0;k<free_.size();++k){const int i=free_[k];mass_.segment<2>(2*k).setConstant(p.rho*old_.V[i]);velocity_.segment<2>(2*k)=old_.v[i].head<2>();
            if(int(old_.drag.size())==n)C_.block<2,2>(2*k,2*k).diagonal().array()+=old_.drag[i];
            if(int(old_.drag_n.size())==n&&int(old_.drag_dir.size())==n){const TransferPoint e=old_.drag_dir[i].head<2>();C_.block<2,2>(2*k,2*k).noalias()+=old_.drag_n[i]*(e*e.transpose());}}
        diagonal_=mass_/(dt_*dt_)+C_.diagonal()/dt_;length_=std::numeric_limits<real>::infinity();
        for(int q=0;q<G.pairs.Size();++q)length_=std::min(length_,(G.x[G.pairs.pi[q]]-G.x[G.pairs.pj[q]]).norm());
    }
    int Dim() const override{return int(2*free_.size());}
    // Additive access to the frozen PHYSICAL metric, independent of any
    // particular free energy or monolayer concentration domain. Existing
    // mechanics methods below keep their original implementation/behavior.
    real Metric_Value(const VecX& q) const{const VecX d=q/dt_-velocity_;return .5*(mass_.array()*d.array().square()).sum()+.5*q.dot(C_*q)/dt_;}
    void Metric_Gradient(const VecX& q,VecX& g) const{g=mass_.array()*(q.array()/dt_-velocity_.array())/dt_;g.noalias()+=C_*q/dt_;}
    void Metric_Hess_Vec(const VecX& v,VecX& h) const{h=mass_.array()*v.array()/(dt_*dt_);h.noalias()+=C_*v/dt_;}
    FilmMatState Trial(const VecX& q) const{FilmMatState s=old_;for(size_t k=0;k<free_.size();++k){s.x[free_[k]].x()+=q[2*k];s.x[free_[k]].y()+=q[2*k+1];}return s;}
    real Value(const VecX& q) override{try{const auto E=Evaluate(q);const VecX d=q/dt_-velocity_;return E.energy+.5*(mass_.array()*d.array().square()).sum()+.5*q.dot(C_*q)/dt_;}catch(const std::exception&){return std::numeric_limits<real>::infinity();}}
    void Gradient(const VecX& q,VecX& g) override{const auto E=Evaluate(q);g.resize(Dim());for(size_t k=0;k<free_.size();++k)g.segment<2>(2*k)=E.gradient[free_[k]];
        g.array()+=mass_.array()*(q.array()/dt_-velocity_.array())/dt_;g.noalias()+=C_*q/dt_;}
    void Prepare_Hessian(const VecX& q) override{prepared_=q;}
    // Exact first derivative; Hessian-vector products are centered differences of
    // that derivative plus the analytic frozen inertial/viscous metric.
    void Hess_Vec(const VecX& v,VecX& h) override{const real vnorm=v.lpNorm<Eigen::Infinity>();if(vnorm==0.){h=VecX::Zero(Dim());return;}const real step=1e-5*length_/vnorm;
        const auto a=Evaluate(prepared_+step*v),b=Evaluate(prepared_-step*v);h.resize(Dim());
        for(size_t k=0;k<free_.size();++k)h.segment<2>(2*k)=(a.gradient[free_[k]]-b.gradient[free_[k]])/(2*step);
        h.array()+=mass_.array()*v.array()/(dt_*dt_);h.noalias()+=C_*v/dt_;}
    void Hess_Diag(VecX& d) override{d=diagonal_;}
    real Length() const{return length_;}
    real Rayleigh(const VecX& q) const{return .5*q.dot(C_*q)/(dt_*dt_);}
    real Kinetic(const FilmMatState& s) const{real k=0.;for(int i:free_)k+=.5*prm_.rho*s.V[i]*s.v[i].squaredNorm();return k;}
    const std::vector<int>& Free() const{return free_;}
private:
    VoronoiEnergyDerivative Evaluate(const VecX& q){if(cached_.size()==q.size()&&(cached_.array()==q.array()).all())return evaluation_;auto fresh=Differentiate_Voronoi_Energy(Trial(q),prm_,domain_,res_);cached_=q;evaluation_=std::move(fresh);return evaluation_;}
    FilmMatState old_;FilmMechParams prm_;TransferPolygon domain_;real dt_,length_;const std::vector<real>* res_;std::vector<int> id_,free_;
    Eigen::MatrixXd C_;VecX mass_,velocity_,diagonal_,prepared_,cached_;VoronoiEnergyDerivative evaluation_;
};
struct MovingVoronoiStats {NewtonStats newton;real energy_before=0.,energy_after=0.,kinetic_before=0.,kinetic_after=0.,rayleigh_action=0.,incremental_change=0.,maximum_displacement=0.,initial_gradient_inf=0.,final_gradient_inf=0.;};
inline MovingVoronoiStats Moving_Voronoi_Mechanics(FilmMatState& s,const FilmMechParams& prm,const TransferPolygon& domain,real dt,const std::vector<real>* reservoir=nullptr,MovingVoronoiDiscretization discretization={}){
    MovingVoronoiProblem problem(s,prm,domain,dt,reservoir,discretization);VecX q=VecX::Zero(problem.Dim()),g;problem.Gradient(q,g);MovingVoronoiStats r;
    r.initial_gradient_inf=g.lpNorm<Eigen::Infinity>();r.energy_before=Differentiate_Voronoi_Energy(s,prm,domain,reservoir).energy;r.kinetic_before=problem.Kinetic(s);const real j0=problem.Value(q);
    NewtonOptions options;options.max_iters=35;options.max_cg_iters=100;options.forcing_max=.05;options.grad_tol=std::max(real(1e-10),g.lpNorm<Eigen::Infinity>()*1e-8);options.step_tol=1e-10*problem.Length();options.max_total_step=.2*problem.Length();options.sparse_preconditioner=false;
    r.newton=Newton_CG_Minimize(problem,q,options);r.maximum_displacement=q.lpNorm<Eigen::Infinity>();r.incremental_change=problem.Value(q)-j0;problem.Gradient(q,g);r.final_gradient_inf=g.lpNorm<Eigen::Infinity>();
    auto final=problem.Trial(q);for(int i=0;i<s.Size();++i)final.v[i].setZero();for(size_t k=0;k<problem.Free().size();++k)final.v[problem.Free()[k]].head<2>()=q.segment<2>(2*k)/dt;
    r.energy_after=Differentiate_Voronoi_Energy(final,prm,domain,reservoir).energy;r.kinetic_after=problem.Kinetic(final);r.rayleigh_action=2*dt*problem.Rayleigh(q);
    s=std::move(final);s.has_prev=false;return r;
}
// Conservative momentum transport by the SAME implicit donor matrix as liquid.
// Unlike keeping each point's old momentum while changing its mass, this leaves
// a uniform velocity unchanged and is non-increasing in passive kinetic energy.
inline Vec3 Moving_Voronoi_Momentum_Transport(const ThinFilmGraph& G,const std::vector<real>& F,
    const std::vector<real>& oldV,FilmMatState& s,real rho,const std::vector<real>* reservoir=nullptr){
    Vec3 supplied=Vec3::Zero();
    for(int axis=0;axis<2;++axis){std::vector<real> P(G.N),bath(G.N,std::numeric_limits<real>::quiet_NaN());
        for(int i=0;i<G.N;++i){P[i]=rho*oldV[i]*s.v[i][axis];if((int(s.fixed.size())==G.N&&s.fixed[i])||(reservoir&&!std::isnan((*reservoir)[i])))bath[i]=0.;}
        supplied[axis]=Solute_Transport(G,F,s.V,P,&bath);for(int i=0;i<G.N;++i)s.v[i][axis]=std::isnan(bath[i])?P[i]/(rho*s.V[i]):0.;}
    return supplied; // free-particle momentum change from zero-velocity pinned supports
}
struct MovingVoronoiCoupledStats {
    MovingVoronoiStats mechanics;NewtonStats drainage;
    real energy_before=0.,energy_after=0.,kinetic_after=0.,half_drain_action=0.,momentum_mixing_loss=0.,reservoir_volume_change=0.,reservoir_work=0.,film_energy_change=0.,minimum_theta=1.,maximum_theta=0.;
    std::vector<real> flux;Vec3 pinned_impulse=Vec3::Zero(),momentum_before_drain=Vec3::Zero(),momentum_after_drain=Vec3::Zero();
};
inline real Moving_Voronoi_Kinetic(const FilmMatState& s,real rho){real k=0.;for(int i=0;i<s.Size();++i)k+=.5*rho*s.V[i]*s.v[i].squaredNorm();return k;}
inline MovingVoronoiCoupledStats Moving_Voronoi_Coupled_Step(FilmMatState& s,const FilmMechParams& p,const TransferPolygon& domain,real dt,const std::vector<real>* reservoir=nullptr,MovingVoronoiDiscretization discretization={}){
    MovingVoronoiCoupledStats r;r.energy_before=Differentiate_Voronoi_Energy(s,p,domain,reservoir).energy;
    const auto oldV=s.V;r.mechanics=Moving_Voronoi_Mechanics(s,p,domain,dt,reservoir,discretization);
    if(!r.mechanics.newton.converged||r.mechanics.newton.hit_total_step)return r;
    const auto geometry=Moving_Voronoi_Geometry(s,domain,p);const auto& G=geometry.graph;
    ThinFilmEnergy energy(G,p.tf);energy.Set_Reservoir(reservoir);std::vector<real> k(G.pairs.Size());
    auto bath=[&](int i){return reservoir&&!std::isnan((*reservoir)[i]);};
    // The bath's bookkeeping volume is NOT a film thickness. On a bath/film face
    // use the one-sided film thickness for physical lubrication mobility. Bath/
    // bath faces have equal prescribed potentials, hence zero physical flux.
    for(int q=0;q<G.pairs.Size();++q){const int i=G.pairs.pi[q],j=G.pairs.pj[q];real h=.5*(s.V[i]/G.a[i]+s.V[j]/G.a[j]);
        if(bath(i)&&!bath(j))h=s.V[j]/G.a[j];if(bath(j)&&!bath(i))h=s.V[i]/G.a[i];k[q]=p.tf.Mobility(h)*G.c[q];}
    NewtonOptions options;options.max_cg_iters=5000;r.drainage=H1_Flux_Step(energy,s.V,k,dt,options,1e-11,false,true,0.,nullptr,3,&r.flux);
    if(!r.drainage.converged||std::any_of(s.V.begin(),s.V.end(),[](real v){return !(v>0.);}))return r;
    for(int q=0;q<G.pairs.Size();++q)r.half_drain_action+=r.flux[q]*r.flux[q]/(2*dt*k[q]);
    real trueOldK=0.;for(int i=0;i<s.Size();++i){trueOldK+=.5*p.rho*oldV[i]*s.v[i].squaredNorm();r.momentum_before_drain+=p.rho*oldV[i]*s.v[i];}
    r.pinned_impulse=Moving_Voronoi_Momentum_Transport(G,r.flux,oldV,s,p.rho,reservoir);r.kinetic_after=Moving_Voronoi_Kinetic(s,p.rho);r.momentum_mixing_loss=trueOldK-r.kinetic_after;
    for(int i=0;i<s.Size();++i)r.momentum_after_drain+=p.rho*s.V[i]*s.v[i];
    r.energy_after=Differentiate_Voronoi_Energy(s,p,domain,reservoir).energy;
    for(int i=0;i<s.Size();++i){if(bath(i)){const real dv=s.V[i]-oldV[i];r.reservoir_volume_change+=dv;r.reservoir_work-=(*reservoir)[i]*dv;}
        const real theta=s.S[i]/G.a[i];r.minimum_theta=std::min(r.minimum_theta,theta);r.maximum_theta=std::max(r.maximum_theta,theta);}
    r.film_energy_change=r.energy_after-r.energy_before+r.reservoir_work;
    return r;
}
} // namespace vf::thinfilm
