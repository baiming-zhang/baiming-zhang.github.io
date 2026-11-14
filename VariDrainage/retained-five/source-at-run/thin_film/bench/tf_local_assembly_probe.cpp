#include <chrono>
#include <iostream>
#include <iomanip>
#include "vf/thinfilm/FrozenCurvedProbe.h"
#include "vf/thinfilm/CurvedEdgeLocalAssembly.h"
using namespace vf;using namespace vf::thinfilm;
int main(int argc,char** argv){try{
    if(argc!=3)return 2;const std::filesystem::path input(argv[1]),output(argv[2]);if(std::filesystem::exists(output))return 3;Set_Threads(2);FrozenCurvedProbe c(input);const auto& s=c.state;const auto geom=Moving_Voronoi_Geometry(s,c.domain,c.params);
    CurvedEdgeEnergy ref(geom.graph,c.params,s.S,c.boundary);CurvedEdgeLocalAssembly candidate(geom.graph,c.params,s.S,c.boundary);auto a=ref.Fields(s.V),b=candidate.Fields(s.V);
    const bool fields=a.energy==b.energy&&a.area==b.area&&a.pressure==b.pressure&&a.tension==b.tension&&a.node_energy==b.node_energy;
    const real scale=std::max(real(1e-30),a.hessian.cwiseAbs().maxCoeff()),matrix=(a.hessian-b.hessian).cwiseAbs().maxCoeff()/scale,symmetry=(b.hessian-b.hessian.transpose()).cwiseAbs().maxCoeff()/scale;
    real diagonal=0,action=0;std::vector<real> w(s.Size()),old,newer;
    for(int k=0;k<3;++k){for(int i=0;i<s.Size();++i)w[i]=c.boundary.Bath(i)?0:.01*s.V[i]*std::sin((.17+.2*k)*(i+1));const Eigen::VectorXd x=a.hessian*Eigen::Map<const Eigen::VectorXd>(w.data(),w.size());candidate.Hess_Apply(s.V,w,newer);real error=0,norm=0;for(int i=0;i<s.Size();++i){error=std::max(error,std::abs(newer[i]-x[i]));norm=std::max(norm,std::abs(x[i]));}action=std::max(action,error/std::max(real(1e-30),norm));}
    for(int i=0;i<s.Size();++i)diagonal=std::max(diagonal,std::abs(candidate.Hess_Diag(s.V,i)-a.hessian(i,i))/scale);
    if(!fields||std::max({matrix,symmetry,diagonal,action})>2e-13)throw std::runtime_error("local Hessian equivalence guard failed");
    std::vector<double> rt,ct;real sink=0;
    auto measure=[&](std::vector<double>& times,auto&& fn){const auto start=std::chrono::steady_clock::now();fn();times.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());};
    for(int k=0;k<3;++k){auto oldcall=[&]{measure(rt,[&]{a=ref.Fields(s.V);sink+=a.hessian(0,0);});};auto newcall=[&]{measure(ct,[&]{b=candidate.Fields(s.V);sink+=b.hessian(0,0);});};if(k%2){newcall();oldcall();}else{oldcall();newcall();}}
    std::ofstream out(output);out<<std::setprecision(17)<<"{\"N\":"<<s.Size()<<",\"first_fields_exact\":"<<(fields?"true":"false")<<",\"matrix_relative_error\":"<<matrix<<",\"symmetry_relative_error\":"<<symmetry<<",\"diagonal_relative_error\":"<<diagonal<<",\"action_relative_error\":"<<action<<",\"sink\":"<<sink<<",\"original_wall_s\":[";for(size_t i=0;i<rt.size();++i){if(i)out<<',';out<<rt[i];}out<<"],\"candidate_wall_s\":[";for(size_t i=0;i<ct.size();++i){if(i)out<<',';out<<ct[i];}out<<"]}\n";
    std::cout<<"local Hessian equivalence N="<<s.Size()<<" matrix="<<matrix<<" action="<<action<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
