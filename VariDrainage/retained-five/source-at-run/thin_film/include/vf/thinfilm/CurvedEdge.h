// Full-state curved edge on the moving Voronoi film. Experimental DRIFT model:
// directional dual-face area quadrature; long-wave inherited dissipation.
// No fitted law, imposed interior profile, reduced basis or kinetic coefficient.
#pragma once
#include <iomanip>
#include <sstream>
#include "vf/thinfilm/MovingVoronoiFilm.h"
namespace vf::thinfilm {
struct CurvedEdgeBoundary {
    real aspect=1.,clean_tension=1.; // H/L and measured gamma0/S from FilmScales/FilmMaterial
    std::vector<real> height,potential; // fixed far heights and liquid potentials; NaN = evolved cell
    // Bath V is a SIGNED exchange counter, never a finite inventory or h.
    // Physical non-bath V and every surface S must remain positive.
    bool Bath(int i) const{return !potential.empty()&&!std::isnan(potential.at(i));}
    real Height(int i,real v,real a) const{return Bath(i)?height.at(i):v/a;}
};
struct CurvedEdgeFields {real energy=0.;std::vector<real> area,pressure,node_energy,tension;Eigen::MatrixXd hessian;};
// Fixed-geometry variational drainage energy. Exact analytic gradient/Hessian
// includes surface dilution by curved area, including cross-face couplings.
class CurvedEdgeEnergy {
public:
    CurvedEdgeEnergy(const ThinFilmGraph& g,const FilmMechParams& p,const std::vector<real>& surface,const CurvedEdgeBoundary& b):g_(g),p_(p),s_(surface),b_(b){
        if(s_.size()!=g_.N||p_.sol.On()||!(b_.aspect>0&&b_.clean_tension>0)||p_.surf.Gamma_inf<=0)throw std::invalid_argument("curved edge requires positive units and insoluble Langmuir surface");
        real amount=0,area=0;for(int i=0;i<g_.N;++i){amount+=s_[i];area+=g_.a[i];}mean_=amount/area/p_.surf.Gamma_inf;
        if(!(mean_>0&&mean_<1))throw std::invalid_argument("curved edge reference coverage");reference_=b_.clean_tension+p_.surf.RT*p_.surf.Gamma_inf*std::log1p(-mean_);
    }
    const ThinFilmGraph& Graph() const{return g_;}
    CurvedEdgeFields Fields(const std::vector<real>& v,bool second=true) const{
        const int n=g_.N;CurvedEdgeFields r;r.area=g_.a;r.pressure.assign(n,0.);r.node_energy.assign(n,0.);r.hessian=Eigen::MatrixXd::Zero(n,n);
        if(v.size()!=n)throw std::invalid_argument("curved edge V size");
        std::vector<real> h(n),extra(n,0.),tension(n);std::vector<Eigen::VectorXd> da;
        for(int i=0;i<n;++i){if(!std::isfinite(v[i])||!(s_[i]>0)||(!b_.Bath(i)&&!(v[i]>0)))throw std::invalid_argument("nonpositive physical curved edge inventory");h[i]=b_.Height(i,v[i],g_.a[i]);da.push_back(Eigen::VectorXd::Zero(n));}
        struct Face{int i,j;real first,second;Eigen::VectorXd direction;};std::vector<Face> faces;
        for(int q=0;q<g_.pairs.Size();++q){const int i=g_.pairs.pi[q],j=g_.pairs.pj[q];const real distance=(g_.x[i]-g_.x[j]).norm(),length=g_.c[q]*distance;
            const real z=.5*b_.aspect*(h[i]-h[j])/distance,f=std::hypot(1.,z),e=length*distance*z*z/(f+1.);
            extra[i]+=.5*e;extra[j]+=.5*e;Eigen::VectorXd direction=Eigen::VectorXd::Zero(n);
            if(!b_.Bath(i))direction[i]=1/g_.a[i];if(!b_.Bath(j))direction[j]=-1/g_.a[j];
            const real first=.5*length*b_.aspect*z/f,second=length*b_.aspect*b_.aspect/(4*distance*f*f*f);
            da[i]+=.5*first*direction;da[j]+=.5*first*direction;faces.push_back({i,j,first,second,std::move(direction)});
        }
        for(int i=0;i<n;++i){r.area[i]+=extra[i];const real theta=s_[i]/(p_.surf.Gamma_inf*r.area[i]);if(!(theta>0&&theta<1))throw std::invalid_argument("curved monolayer outside Langmuir domain");
            const real d=theta-mean_,kl=p_.surf.RT*p_.surf.Gamma_inf*(theta*std::log1p(d/mean_)+(1-theta)*std::log1p(-d/(1-mean_)));
            r.node_energy[i]=2*(r.area[i]*kl+reference_*extra[i]);tension[i]=b_.clean_tension+p_.surf.RT*p_.surf.Gamma_inf*std::log1p(-theta);
            if(!(tension[i]>0))throw std::invalid_argument("negative curved surface tension");
            for(int j=0;j<n;++j)r.pressure[j]+=2*tension[i]*da[i][j];
            if(second)r.hessian.noalias()+=(2*p_.surf.RT*p_.surf.Gamma_inf*theta/(r.area[i]*(1-theta)))*(da[i]*da[i].transpose());
            if(b_.Bath(i)){r.node_energy[i]+=b_.potential[i]*v[i];r.pressure[i]+=b_.potential[i];}
            else{r.node_energy[i]+=g_.a[i]*p_.tf.disj.P(h[i])+p_.tf.rho_g*v[i]*g_.Height(i);r.pressure[i]+=p_.tf.disj.dP(h[i])+p_.tf.rho_g*g_.Height(i);r.hessian(i,i)+=p_.tf.disj.d2P(h[i])/g_.a[i];}
            r.energy+=r.node_energy[i];
        }
        if(second)for(const auto& f:faces)r.hessian.noalias()+=(tension[f.i]+tension[f.j])*f.second*(f.direction*f.direction.transpose());
        r.tension=std::move(tension);return r;
    }
    real Value(const std::vector<real>& v) const{try{return Fields(v,false).energy;}catch(const std::exception&){return std::numeric_limits<real>::infinity();}}
    void Pressure(const std::vector<real>& v,std::vector<real>& p) const{p=Fields(v,false).pressure;}
    void Hess_Apply(const std::vector<real>& v,const std::vector<real>& w,std::vector<real>& out) const{const auto& f=Tangent(v);const Eigen::VectorXd y=f.hessian*Eigen::Map<const Eigen::VectorXd>(w.data(),w.size());out.assign(y.data(),y.data()+y.size());}
    real Hess_Diag(const std::vector<real>& v,int i) const{return Tangent(v).hessian(i,i);}
    void Hess_Triplets(const std::vector<real>& v,std::vector<Eigen::Triplet<real>>& t,const std::vector<real>* =nullptr) const{
        const auto& f=Tangent(v);for(int j=0;j<g_.N;++j)for(int i=0;i<g_.N;++i)if(f.hessian(i,j)!=0)t.emplace_back(i,j,f.hessian(i,j));
        // Full exact tangent also used in the small-problem preconditioner;
        // sparsified pair scaling is deliberately not an accelerated variant.
    }
private:const CurvedEdgeFields& Tangent(const std::vector<real>& v) const{if(v!=cached_v_){auto f=Fields(v);cached_v_=v;cached_=std::move(f);}return cached_;}
    const ThinFilmGraph& g_;FilmMechParams p_;std::vector<real> s_;CurvedEdgeBoundary b_;real mean_,reference_;mutable std::vector<real> cached_v_;mutable CurvedEdgeFields cached_;
};
// Differentiate area geometry on its current Voronoi topology. Boundary heights
// remain fixed, interior h=V/A and material S is conserved during mechanics.
inline VoronoiEnergyDerivative Differentiate_Curved_Edge(const FilmMatState& st,const FilmMechParams& p,const TransferPolygon& domain,const CurvedEdgeBoundary& boundary){
    using namespace voronoi_ad;const int n=st.Size();const auto geom=Moving_Voronoi_Geometry(st,domain,p);const auto& G=geom.graph;const auto sites=Moving_Voronoi_Sites(st);
    if(p.sol.On()||!st.B.empty()||!st.Vs.empty()||!st.gscale.empty()||!st.fext.empty()||!st.bath_w.empty()||!st.sbath.empty())throw std::invalid_argument("curved edge unsupported soluble/closure state");
    Tape tape;std::vector<Point> x,box;std::vector<int> xi,yi;for(auto q:st.x){x.push_back({C(tape,q.x()),C(tape,q.y())});xi.push_back(x.back().x.id);yi.push_back(x.back().y.id);}
    for(auto q:domain)box.push_back({C(tape,q.x()),C(tape,q.y())});std::vector<std::vector<int>> neighbors(n);for(int q=0;q<G.pairs.Size();++q){neighbors[G.pairs.pi[q]].push_back(G.pairs.pj[q]);neighbors[G.pairs.pj[q]].push_back(G.pairs.pi[q]);}
    std::vector<std::vector<Point>> cells;std::vector<Scalar> area,extra,h;real totalS=0;for(real s:st.S)totalS+=s;const real mean=totalS/(Transfer_Area(domain)*p.surf.Gamma_inf),reference=boundary.clean_tension+p.surf.RT*p.surf.Gamma_inf*std::log1p(-mean);
    for(int i=0;i<n;++i){auto cell=box;std::sort(neighbors[i].begin(),neighbors[i].end());for(int j:neighbors[i])cell=Clip(cell,(x[i]+x[j])*.5,x[j]-x[i]);Clean(cell);const auto a=Area(cell,tape);
        if(!(a.Value()>0)||std::abs(a.Value()/G.a[i]-1)>1e-9)throw std::runtime_error("curved differentiated geometry disagrees");area.push_back(a);extra.push_back(C(tape,0));h.push_back(boundary.Bath(i)?C(tape,boundary.height.at(i)):st.V[i]/a);cells.push_back(std::move(cell));}
    const real eps=1e-10*(TransferBounds(domain).hi-TransferBounds(domain).lo).norm();
    for(int i=0;i<n;++i)for(size_t k=0;k<cells[i].size();++k){const auto a=cells[i][k],b=cells[i][(k+1)%cells[i].size()];const auto av=Value(a),bv=Value(b);if((bv-av).norm()<=eps)continue;int j=-1;real error=std::numeric_limits<real>::infinity();
        for(int other:neighbors[i]){const TransferPoint normal=sites[other]-sites[i],mid=.5*(sites[other]+sites[i]);const real e=std::max(std::abs((av-mid).dot(normal)),std::abs((bv-mid).dot(normal)))/normal.norm();if(e<error){error=e;j=other;}}
        if(error>eps||j<i)continue;const auto length=Sqrt(Dot(b-a,b-a)),distance=Sqrt(Dot(x[j]-x[i],x[j]-x[i]));const auto z=.5*boundary.aspect*(h[i]-h[j])/distance,u=z*z;
        const auto e=length*distance*u/(Sqrt(1+u)+1);extra[i]=extra[i]+.5*e;extra[j]=extra[j]+.5*e;
    }
    Scalar energy=C(tape,0);for(int i=0;i<n;++i){const auto surface=area[i]+extra[i],theta=st.S[i]/(p.surf.Gamma_inf*surface),delta=theta-mean;
        if(!(theta.Value()>0&&theta.Value()<1&&boundary.clean_tension+p.surf.RT*p.surf.Gamma_inf*std::log1p(-theta.Value())>0)){
            std::ostringstream message;message<<std::setprecision(17)<<"curved differentiated monolayer/tension domain node="<<i<<" theta="<<theta.Value()<<" projected_theta="<<st.S[i]/(p.surf.Gamma_inf*area[i].Value())<<" surface_area="<<surface.Value();throw std::invalid_argument(message.str());}
        const auto kl=p.surf.RT*p.surf.Gamma_inf*(theta*Log1p(delta/mean)+(1-theta)*Log1p((-1/(1-mean))*delta));energy=energy+2*(surface*kl+reference*extra[i]);
        if(boundary.Bath(i)){energy=energy+boundary.potential[i]*st.V[i];continue;}
        const auto& d=p.tf.disj;const auto ih=1/h[i],ih2=ih*ih,ih8=ih2*ih2*ih2*ih2;energy=energy+area[i]*(d.eps*(std::pow(d.hb,8)*ih8/3.-4*d.hb*d.hb*ih2/3.)+d.c8*ih8-d.c2*ih2+d.beta*Exp((-1/d.lambda)*h[i]));
        energy=energy+p.tf.rho_g*st.V[i]*(p.up.x()*x[i].x+p.up.y()*x[i].y);
    }
    VoronoiEnergyDerivative r;r.energy=energy.Value();r.double_energy=CurvedEdgeEnergy(G,p,st.S,boundary).Value(st.V);r.area=Transfer_Area(domain);r.tape_nodes=tape.nodes.size();r.gradient.resize(n);tape.Reverse(energy.id);for(int i=0;i<n;++i)r.gradient[i]={tape.nodes[xi[i]].adj,tape.nodes[yi[i]].adj};return r;
}
} // namespace vf::thinfilm
