// Reverse derivatives of the actual clipped planar Voronoi energy.
// This DRIFT-only tape is deliberately small; it does not modify core solvers.
// Generic Voronoi/power-particle mechanics is established prior art, not a novelty claim.
#pragma once
#include "vf/thinfilm/VoronoiDrainage.h"
#include "vf/thinfilm/ThinFilmMechanics.h"
namespace vf::thinfilm::voronoi_ad {
struct Tape {
    struct Node {real value,adj=0.;int a=-1,b=-1;real da=0.,db=0.;};
    std::vector<Node> nodes;
    int Constant(real v){nodes.push_back({v});return int(nodes.size())-1;}
    int Unary(real v,int a,real d){nodes.push_back({v,0.,a,-1,d,0.});return int(nodes.size())-1;}
    int Binary(real v,int a,int b,real da,real db){nodes.push_back({v,0.,a,b,da,db});return int(nodes.size())-1;}
    void Reverse(int root){nodes[root].adj=1.;for(int k=int(nodes.size())-1;k>=0;--k){const auto n=nodes[k];if(n.a>=0)nodes[n.a].adj+=n.adj*n.da;if(n.b>=0)nodes[n.b].adj+=n.adj*n.db;}}
};
struct Scalar {Tape* tape;int id;real Value() const{return tape->nodes[id].value;}};
inline Scalar C(Tape& t,real x){return {&t,t.Constant(x)};}
inline Scalar operator+(Scalar a,Scalar b){return {a.tape,a.tape->Binary(a.Value()+b.Value(),a.id,b.id,1.,1.)};}
inline Scalar operator-(Scalar a,Scalar b){return {a.tape,a.tape->Binary(a.Value()-b.Value(),a.id,b.id,1.,-1.)};}
inline Scalar operator*(Scalar a,Scalar b){return {a.tape,a.tape->Binary(a.Value()*b.Value(),a.id,b.id,b.Value(),a.Value())};}
inline Scalar operator/(Scalar a,Scalar b){const real ib=1./b.Value();return {a.tape,a.tape->Binary(a.Value()*ib,a.id,b.id,ib,-a.Value()*ib*ib)};}
inline Scalar operator*(Scalar a,real b){return {a.tape,a.tape->Unary(a.Value()*b,a.id,b)};}
inline Scalar operator*(real a,Scalar b){return b*a;}
inline Scalar operator+(Scalar a,real b){return {a.tape,a.tape->Unary(a.Value()+b,a.id,1.)};}
inline Scalar operator+(real a,Scalar b){return b+a;}
inline Scalar operator-(Scalar a,real b){return a+(-b);}
inline Scalar operator-(real a,Scalar b){return (-1.)*b+a;}
inline Scalar operator/(Scalar a,real b){return a*(1./b);}
inline Scalar operator/(real a,Scalar b){return C(*b.tape,a)/b;}
inline Scalar Log(Scalar x){return {x.tape,x.tape->Unary(std::log(x.Value()),x.id,1./x.Value())};}
inline Scalar Log1p(Scalar x){return {x.tape,x.tape->Unary(std::log1p(x.Value()),x.id,1./(1.+x.Value()))};}
inline Scalar Exp(Scalar x){const real e=std::exp(x.Value());return {x.tape,x.tape->Unary(e,x.id,e)};}
inline Scalar Sqrt(Scalar x){const real s=std::sqrt(x.Value());return {x.tape,x.tape->Unary(s,x.id,.5/s)};}
struct Point {Scalar x,y;};
inline Point operator+(Point a,Point b){return {a.x+b.x,a.y+b.y};}
inline Point operator-(Point a,Point b){return {a.x-b.x,a.y-b.y};}
inline Point operator*(Point a,Scalar b){return {a.x*b,a.y*b};}
inline Point operator*(Point a,real b){return {a.x*b,a.y*b};}
inline Scalar Dot(Point a,Point b){return a.x*b.x+a.y*b.y;}
inline Scalar Cross(Point a,Point b){return a.x*b.y-a.y*b.x;}
inline TransferPoint Value(Point a){return TransferPoint(a.x.Value(),a.y.Value());}
inline std::vector<Point> Clip(const std::vector<Point>& p,Point origin,Point normal){
    std::vector<Point> out;if(p.empty())return out;Point a=p.back();Scalar da=Dot(a-origin,normal);
    for(const auto& b:p){const Scalar db=Dot(b-origin,normal);
        if((da.Value()<=0.)!=(db.Value()<=0.))out.push_back(a+(b-a)*(da/(da-db)));
        if(db.Value()<=0.)out.push_back(b);a=b;da=db;
    }return out;
}
inline void Clean(std::vector<Point>& p){real coordinate=0.;for(auto x:p)coordinate=std::max(coordinate,Value(x).norm());const real eps=32*std::numeric_limits<real>::epsilon()*coordinate;
    std::vector<Point> out;for(auto x:p)if(out.empty()||(Value(x)-Value(out.back())).norm()>eps)out.push_back(x);
    if(out.size()>1&&(Value(out.back())-Value(out.front())).norm()<=eps)out.pop_back();p=std::move(out);}
inline Scalar Area(const std::vector<Point>& p,Tape& t){Scalar a=C(t,0.);for(size_t k=1;k+1<p.size();++k)a=a+.5*Cross(p[k]-p[0],p[k+1]-p[0]);return a;}
} // namespace vf::thinfilm::voronoi_ad

namespace vf::thinfilm {
struct VoronoiEnergyDerivative {real energy=0.,double_energy=0.,area=0.;std::vector<TransferPoint> gradient;size_t tape_nodes=0;};
inline VoronoiEnergyDerivative Differentiate_Voronoi_Energy(const FilmMatState& st,const FilmMechParams& prm,
    const TransferPolygon& domain,const std::vector<real>* reservoir=nullptr)
{
    using namespace voronoi_ad;const int n=st.Size();
    if(prm.sol.On()||!st.B.empty()||!st.Vs.empty())throw std::invalid_argument("moving Voronoi prototype requires an insoluble monolayer, no hidden bulk/column closure");
    if(!st.gscale.empty()||!st.fext.empty()||!st.bath_w.empty()||!st.sbath.empty())throw std::invalid_argument("moving prototype does not silently inherit kernel-frame closure fields");
    if(int(st.V.size())!=n||int(st.S.size())!=n||(reservoir&&int(reservoir->size())!=n))throw std::invalid_argument("moving Voronoi energy state size");
    std::vector<TransferPoint> sites;for(auto x:st.x){if(x.z()!=0.)throw std::invalid_argument("moving Voronoi energy is planar");sites.emplace_back(x.x(),x.y());}
    auto geom=Make_Voronoi_Drainage(sites,domain);auto& G=geom.graph;
    for(int i=0;i<n;++i){G.height[i]=prm.up.dot(st.x[i]);if(!(st.V[i]>0.&&st.S[i]>0.))throw std::invalid_argument("nonpositive moving inventory");}
    auto bath=[&](int i){return reservoir&&!std::isnan((*reservoir)[i]);};
    Tape tape;tape.nodes.reserve(size_t(n)*500);std::vector<Point> x,box;std::vector<int> xid,yid;
    for(auto p:sites){Point q{C(tape,p.x()),C(tape,p.y())};x.push_back(q);xid.push_back(q.x.id);yid.push_back(q.y.id);}
    for(auto p:domain)box.push_back({C(tape,p.x()),C(tape,p.y())});
    std::vector<std::vector<int>> neighbors(n);for(int q=0;q<G.pairs.Size();++q){neighbors[G.pairs.pi[q]].push_back(G.pairs.pj[q]);neighbors[G.pairs.pj[q]].push_back(G.pairs.pi[q]);}
    // Conserved S and fixed domain area permit a CONSTANT energy reference:
    // subtract the uniform-monolayer free energy. Stable relative entropy avoids
    // subtracting a huge constant from a tiny mechanical energy change. The
    // uniform bath-pressure area term is likewise constant on this domain.
    real totalS=0.;for(real v:st.S)totalS+=v;const real meanGamma=totalS/Transfer_Area(domain);
    std::vector<std::vector<Point>> cells;std::vector<Scalar> area,h;Scalar energy=C(tape,0.);VoronoiEnergyDerivative result;result.gradient.resize(n);
    for(int i=0;i<n;++i){auto p=box;std::sort(neighbors[i].begin(),neighbors[i].end());for(int j:neighbors[i])p=Clip(p,(x[i]+x[j])*.5,x[j]-x[i]);Clean(p);
        Scalar a=Area(p,tape);if(!(a.Value()>0.)||std::abs(a.Value()/G.a[i]-1.)>1e-9)throw std::runtime_error("differentiated Voronoi cell disagrees with primal geometry");
        area.push_back(a);h.push_back(st.V[i]/a);cells.push_back(std::move(p));result.area+=a.Value();
        const Scalar gamma=st.S[i]/a;Scalar g=C(tape,0.);
        if(prm.surf.Gamma_inf>0.){const Scalar theta=gamma/prm.surf.Gamma_inf;if(!(theta.Value()>0.&&theta.Value()<1.))throw std::invalid_argument("moving monolayer outside Langmuir domain");
            const real mean=meanGamma/prm.surf.Gamma_inf;const Scalar delta=theta-mean;
            g=prm.surf.RT*prm.surf.Gamma_inf*(theta*Log1p(delta/mean)+(1.-theta)*Log1p((-1./(1.-mean))*delta));
        }else g=prm.surf.RT*(gamma*Log1p((gamma-meanGamma)/meanGamma)-(gamma-meanGamma));
        energy=energy+2.*a*g;
        if(bath(i)){energy=energy+(*reservoir)[i]*st.V[i];continue;}
        const auto& d=prm.tf.disj;const Scalar ih=1./h.back(),ih2=ih*ih,ih8=ih2*ih2*ih2*ih2;
        energy=energy+a*(d.eps*(std::pow(d.hb,8)*ih8/3.-4.*d.hb*d.hb*ih2/3.)+d.c8*ih8-d.c2*ih2+d.beta*Exp((-1./d.lambda)*h.back()));
        energy=energy+prm.tf.rho_g*st.V[i]*(prm.up.x()*x[i].x+prm.up.y()*x[i].y);
    }
    const real eps=1e-10*(TransferBounds(domain).hi-TransferBounds(domain).lo).norm();
    for(int i=0;i<n;++i)if(!bath(i))for(size_t v=0;v<cells[i].size();++v){const Point a=cells[i][v],b=cells[i][(v+1)%cells[i].size()];const TransferPoint av=Value(a),bv=Value(b);if((bv-av).norm()<=eps)continue;
        int other=-1;real error=std::numeric_limits<real>::infinity();
        for(int j:neighbors[i]){const TransferPoint normal=sites[j]-sites[i],mid=.5*(sites[j]+sites[i]);const real e=std::max(std::abs((av-mid).dot(normal)),std::abs((bv-mid).dot(normal)))/normal.norm();if(e<error){error=e;other=j;}}
        if(error>eps||other<i||bath(other))continue;
        const Scalar length=Sqrt(Dot(b-a,b-a)),distance=Sqrt(Dot(x[other]-x[i],x[other]-x[i]));const Scalar dh=h[i]-h[other];
        energy=energy+(prm.tf.sigma/8.)*(length/distance)*dh*dh;
    }
    result.energy=energy.Value();result.tape_nodes=tape.nodes.size();if(!std::isfinite(result.energy))throw std::runtime_error("nonfinite moving energy");
    tape.Reverse(energy.id);for(int i=0;i<n;++i)result.gradient[i]=TransferPoint(tape.nodes[xid[i]].adj,tape.nodes[yid[i]].adj);
    ThinFilmEnergy ordinary(G,prm.tf);ordinary.Set_Reservoir(reservoir);result.double_energy=ordinary.Value(st.V);
    for(int i=0;i<n;++i){const real gamma=st.S[i]/G.a[i];real excess;
        if(prm.surf.Gamma_inf>0.){const real theta=gamma/prm.surf.Gamma_inf,mean=meanGamma/prm.surf.Gamma_inf,delta=theta-mean;
            excess=prm.surf.RT*prm.surf.Gamma_inf*(theta*std::log1p(delta/mean)+(1.-theta)*std::log1p(-delta/(1.-mean)));}
        else excess=prm.surf.RT*(gamma*std::log1p((gamma-meanGamma)/meanGamma)-(gamma-meanGamma));
        result.double_energy+=2.*G.a[i]*excess;}
    if(std::abs(result.energy-result.double_energy)>1e-9*std::max(real(1),std::abs(result.double_energy)))throw std::runtime_error("differentiated and inherited energies disagree");
    return result;
}
} // namespace vf::thinfilm
