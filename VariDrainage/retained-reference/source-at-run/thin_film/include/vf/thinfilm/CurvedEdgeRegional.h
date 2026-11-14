// Observer-only derivatives of a frozen weighted sum of the SAME node energies.
// Uses global reference coverage, never a separately fitted regional reference.
// Weights are held fixed during differentiation; moving ownership is separate.
#pragma once
#include "vf/thinfilm/CurvedEdge.h"
namespace vf::thinfilm {
inline VoronoiEnergyDerivative Differentiate_Curved_Edge_Weighted(const FilmMatState& st,const FilmMechParams& p,const TransferPolygon& domain,const CurvedEdgeBoundary& boundary,const std::vector<real>& weights){
    if(weights.size()!=st.Size())throw std::invalid_argument("regional energy weights");
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
        const auto kl=p.surf.RT*p.surf.Gamma_inf*(theta*Log1p(delta/mean)+(1-theta)*Log1p((-1/(1-mean))*delta));energy=energy+weights[i]*2*(surface*kl+reference*extra[i]);
        if(boundary.Bath(i)){energy=energy+weights[i]*boundary.potential[i]*st.V[i];continue;}
        const auto& d=p.tf.disj;const auto ih=1/h[i],ih2=ih*ih,ih8=ih2*ih2*ih2*ih2;energy=energy+weights[i]*area[i]*(d.eps*(std::pow(d.hb,8)*ih8/3.-4*d.hb*d.hb*ih2/3.)+d.c8*ih8-d.c2*ih2+d.beta*Exp((-1/d.lambda)*h[i]));
        energy=energy+weights[i]*p.tf.rho_g*st.V[i]*(p.up.x()*x[i].x+p.up.y()*x[i].y);
    }
    VoronoiEnergyDerivative r;r.energy=energy.Value();const auto fields=CurvedEdgeEnergy(G,p,st.S,boundary).Fields(st.V,false);r.double_energy=0;for(int i=0;i<n;++i)r.double_energy+=weights[i]*fields.node_energy[i];r.area=Transfer_Area(domain);r.tape_nodes=tape.nodes.size();r.gradient.resize(n);tape.Reverse(energy.id);for(int i=0;i<n;++i)r.gradient[i]={tape.nodes[xi[i]].adj,tape.nodes[yi[i]].adj};return r;
}
} // namespace vf::thinfilm
