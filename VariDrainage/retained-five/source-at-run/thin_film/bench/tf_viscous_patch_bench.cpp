#include "vf/thinfilm/CurvedEdgeMechanicalObserver.h"
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <iostream>
using namespace vf;using namespace vf::thinfilm;
// Manufactured SI vector field, NOT a drainage evolution or experimental case.
// The native metric is applied unchanged; no body force enters any solver.
int main(int argc,char** argv){
    if(argc!=4&&argc!=5)return 2;const bool solve=argc==5&&std::stoi(argv[4])!=0;const int n=std::stoi(argv[1]);const real jitter=std::stod(argv[2]);if(n<4||n>48||jitter<0||jitter>.2)return 2;
    Set_Threads(2);const real L=320e-6,W=160e-6,cut=.63*L,thickness=6e-6,U=1e-3,pi=std::acos(-1.);FilmMatState a;FilmMechParams p;p.rho=1045;p.mu_s=.0015;p.bulk_viscosity=true;
    for(int j=0;j<n;++j)for(int i=0;i<2*n;++i){const int k=int(a.x.size());const real x=(i+.5+jitter*std::sin(1.7*k))*L/(2*n),y=(j+.5+jitter*std::cos(2.3*k))*W/n;const real X=pi*x/L,Y=pi*y/W;
        a.x.emplace_back(x,y,0);a.v.emplace_back(U*std::pow(std::sin(X),2)*std::sin(2*Y),-U*W/L*std::sin(2*X)*std::pow(std::sin(Y),2),0);a.fixed.push_back(0);}
    const auto domain=Transfer_Rect(0,0,L,W);const auto g=Moving_Voronoi_Geometry(a,domain,p);a.S.assign(a.Size(),0);for(real area:g.graph.a)a.V.push_back(area*thickness);
    // dt=1 is only an algebraic action-to-power factor. No position advances.
    MovingVoronoiProblem metric(a,p,domain,1.,nullptr,{true,true});VecX q(metric.Dim()),force;for(int i=0;i<a.Size();++i)q.segment<2>(2*i)=a.v[i].head<2>();metric.Metric_Gradient(q,force);
    CurvedEdgeBoundary b;b.potential.assign(a.Size(),std::numeric_limits<real>::quiet_NaN());CurvedMechanicalObservation o;Observe_Curved_Viscosity(a,a,p,g,b,1.,cut,o);const auto weights=Moving_Interface_Weights(g,cut);real leftWork=0,action=0;
    const std::filesystem::path dest(argv[3]);if(std::filesystem::exists(dest))return 3;std::filesystem::create_directories(dest);
    if(solve){
        VecX mass(metric.Dim()),rhs(metric.Dim()),diagonal;metric.Hess_Diag(diagonal);
        for(int i=0;i<a.Size();++i){mass.segment<2>(2*i).setConstant(p.rho*a.V[i]);const real X=pi*a.x[i].x()/L,Y=pi*a.x[i].y()/W;
            const real lapu=U*(2*pi*pi/(L*L)*std::cos(2*X)-4*pi*pi/(W*W)*std::pow(std::sin(X),2))*std::sin(2*Y),lapv=-U*W/L*std::sin(2*X)*(-4*pi*pi/(L*L)*std::pow(std::sin(Y),2)+2*pi*pi/(W*W)*std::cos(2*Y));
            rhs.segment<2>(2*i)=-p.Mu_S()*thickness*g.graph.a[i]*TransferPoint(lapu,lapv);}
        diagonal-=mass;auto apply=[&](const VecX& v){VecX z;metric.Metric_Hess_Vec(v,z);z.array()-=mass.array()*v.array();return z;};
        VecX solution=VecX::Zero(metric.Dim()),r=rhs,z=(r.array()/diagonal.array()).matrix(),direction=z;real rz=r.dot(z);int it=0;
        while(r.norm()>1e-10*rhs.norm()&&it<2000){const VecX Ad=apply(direction);const real dAd=direction.dot(Ad);if(!(dAd>0))return 4;const real alpha=rz/dAd;solution+=alpha*direction;r-=alpha*Ad;z=(r.array()/diagonal.array()).matrix();const real next=r.dot(z);direction=z+(next/rz)*direction;rz=next;++it;}
        const real relative=(apply(solution)-rhs).norm()/rhs.norm();std::ofstream solved(dest/"solution.csv");solved<<std::setprecision(17)<<"id,ux_m_s,uy_m_s,rhs_x_N,rhs_y_N,res_x_N,res_y_N\n";const VecX residual=apply(solution)-rhs;
        for(int i=0;i<a.Size();++i)solved<<i<<','<<solution[2*i]<<','<<solution[2*i+1]<<','<<rhs[2*i]<<','<<rhs[2*i+1]<<','<<residual[2*i]<<','<<residual[2*i+1]<<'\n';
        std::ofstream solver(dest/"solver.json");solver<<std::setprecision(17)<<"{\"iterations\":"<<it<<",\"relative_residual\":"<<relative<<",\"converged\":"<<(relative<=1.01e-10?"true":"false")<<"}\n";if(relative>1.01e-10)return 5;
    }
    std::ofstream nodes(dest/"nodes.csv");nodes<<std::setprecision(17)<<"id,x_m,y_m,A_m2,ux_m_s,uy_m_s,fx_N,fy_N,left_weight\n";
    for(int i=0;i<a.Size();++i){const TransferPoint v=a.v[i].head<2>(),f=force.segment<2>(2*i);leftWork+=(1-weights[i])*v.dot(f);action+=v.dot(f);nodes<<i<<','<<a.x[i].x()<<','<<a.x[i].y()<<','<<g.graph.a[i]<<','<<v.x()<<','<<v.y()<<','<<f.x()<<','<<f.y()<<','<<1-weights[i]<<'\n';}
    const real DL=o.viscous_total-o.regions[1].viscous_dissipation;std::ofstream out(dest/"summary.json");out<<std::setprecision(17)<<"{\"n\":"<<n<<",\"sites\":"<<a.Size()<<",\"jitter\":"<<jitter<<",\"L_m\":"<<L<<",\"W_m\":"<<W<<",\"cut_m\":"<<cut<<",\"h_m\":"<<thickness<<",\"U_m_s\":"<<U<<",\"mu_Pa_s\":"<<p.Mu_S()<<",\"native_action_W\":"<<action<<",\"pyramid_action_W\":"<<o.viscous_total<<",\"left_dissipation_W\":"<<DL<<",\"left_node_work_W\":"<<leftWork<<",\"discrete_cut_power_W\":"<<DL-leftWork<<",\"recovered_cut_power_W\":"<<-o.viscous_recovered_port<<",\"cut_length_m\":"<<o.cut_length<<"}\n";
    std::cout<<"manufactured metric n="<<n<<" sites="<<a.Size()<<" jitter="<<jitter<<"\n";return 0;
}
