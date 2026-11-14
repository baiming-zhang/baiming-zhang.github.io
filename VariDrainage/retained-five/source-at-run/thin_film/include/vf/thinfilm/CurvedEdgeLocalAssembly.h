// Opt-in DRIFT energy-assembly candidate, not a new model or solver backend.
// All first fields use the original energy. Only exact local Hessian support
// is assembled; dense output and existing H1 interfaces are preserved.
#pragma once
#include <map>
#include "vf/thinfilm/CurvedEdge.h"
namespace vf::thinfilm {
class CurvedEdgeLocalAssembly {
public:
    CurvedEdgeLocalAssembly(const ThinFilmGraph& g,const FilmMechParams& p,const std::vector<real>& s,const CurvedEdgeBoundary& b):original_(g,p,s,b),g_(g),p_(p),s_(s),b_(b){}
    const ThinFilmGraph& Graph() const{return g_;}
    CurvedEdgeFields Fields(const std::vector<real>& v,bool second=true) const{
        auto r=original_.Fields(v,false);if(!second)return r;const int n=g_.N;r.hessian.setZero(n,n);
        std::vector<real> h(n);for(int i=0;i<n;++i)h[i]=b_.Height(i,v[i],g_.a[i]);
        std::vector<std::map<int,real>> da(n);
        struct Face {int i,j;real second;std::map<int,real> direction;};std::vector<Face> faces;faces.reserve(g_.pairs.Size());
        for(int q=0;q<g_.pairs.Size();++q){const int i=g_.pairs.pi[q],j=g_.pairs.pj[q];const real distance=(g_.x[i]-g_.x[j]).norm(),length=g_.c[q]*distance;
            const real z=.5*b_.aspect*(h[i]-h[j])/distance,f=std::hypot(1.,z),first=.5*length*b_.aspect*z/f,second_value=length*b_.aspect*b_.aspect/(4*distance*f*f*f);
            std::map<int,real> direction;if(!b_.Bath(i))direction[i]=1/g_.a[i];if(!b_.Bath(j))direction[j]=-1/g_.a[j];
            for(const auto& item:direction){const real value=.5*first*item.second;da[i][item.first]+=value;da[j][item.first]+=value;}
            faces.push_back({i,j,second_value,std::move(direction)});
        }
        auto outer=[&](const std::map<int,real>& a,real coefficient){
            for(const auto& col:a)if(col.second!=0)for(const auto& row:a)if(row.second!=0)
                r.hessian(row.first,col.first)+=coefficient*(row.second*col.second);
        };
        for(int i=0;i<n;++i){const real theta=s_[i]/(p_.surf.Gamma_inf*r.area[i]);
            outer(da[i],2*p_.surf.RT*p_.surf.Gamma_inf*theta/(r.area[i]*(1-theta)));
            // Same diagonal insertion order as the original, even for nonzero DLVO.
            if(!b_.Bath(i))r.hessian(i,i)+=p_.tf.disj.d2P(h[i])/g_.a[i];
        }
        for(const auto& face:faces)outer(face.direction,(r.tension[face.i]+r.tension[face.j])*face.second);
        return r;
    }
    real Value(const std::vector<real>& v) const{return original_.Value(v);}
    void Pressure(const std::vector<real>& v,std::vector<real>& p) const{original_.Pressure(v,p);}
    void Hess_Apply(const std::vector<real>& v,const std::vector<real>& w,std::vector<real>& out) const{const auto& f=Tangent(v);const Eigen::VectorXd y=f.hessian*Eigen::Map<const Eigen::VectorXd>(w.data(),w.size());out.assign(y.data(),y.data()+y.size());}
    real Hess_Diag(const std::vector<real>& v,int i) const{return Tangent(v).hessian(i,i);}
    void Hess_Triplets(const std::vector<real>& v,std::vector<Eigen::Triplet<real>>& t,const std::vector<real>* =nullptr) const{const auto& f=Tangent(v);for(int j=0;j<g_.N;++j)for(int i=0;i<g_.N;++i)if(f.hessian(i,j)!=0)t.emplace_back(i,j,f.hessian(i,j));}
private:
    const CurvedEdgeFields& Tangent(const std::vector<real>& v) const{if(v!=cached_v_){auto f=Fields(v);cached_v_=v;cached_=std::move(f);}return cached_;}
    CurvedEdgeEnergy original_;const ThinFilmGraph& g_;FilmMechParams p_;std::vector<real> s_;CurvedEdgeBoundary b_;mutable std::vector<real> cached_v_;mutable CurvedEdgeFields cached_;
};
}
