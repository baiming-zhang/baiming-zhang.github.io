#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <ctime>
#endif
#include <chrono>
#include <iostream>
#include <iomanip>
#include "vf/thinfilm/FrozenCurvedProbe.h"
using namespace vf;using namespace vf::thinfilm;
static double Cpu_Seconds(){
#ifdef _WIN32
    FILETIME c,e,k,u;if(!GetProcessTimes(GetCurrentProcess(),&c,&e,&k,&u))throw std::runtime_error("process timer failed");ULARGE_INTEGER a,b;a.LowPart=k.dwLowDateTime;a.HighPart=k.dwHighDateTime;b.LowPart=u.dwLowDateTime;b.HighPart=u.dwHighDateTime;return (a.QuadPart+b.QuadPart)*1e-7;
#else
    return double(std::clock())/CLOCKS_PER_SEC;
#endif
}
int main(int argc,char** argv){try{
    if(argc!=3)return 2;const std::filesystem::path input(argv[1]),output(argv[2]);auto pattern=output;pattern.replace_extension(".pattern.csv");if(std::filesystem::exists(output)||std::filesystem::exists(pattern))return 3;Set_Threads(2);
    FrozenCurvedProbe c(input);const auto& s=c.state;const auto geom=Moving_Voronoi_Geometry(s,c.domain,c.params);CurvedEdgeEnergy energy(geom.graph,c.params,s.S,c.boundary);CurvedEdgeFields full,plain;std::vector<real> pressure,out,w(s.Size());
    for(int i=0;i<s.Size();++i)w[i]=c.boundary.Bath(i)?0:.01*s.V[i]*std::sin(.37*(i+1));
    struct Sample {double wall,cpu;};std::map<std::string,std::vector<Sample>> times;real sink=0;
    auto measure=[&](const std::string& name,int n,auto&& fn,bool warm=true){if(warm)fn();for(int k=0;k<n;++k){const double cpu=Cpu_Seconds();const auto start=std::chrono::steady_clock::now();fn();const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();times[name].push_back({wall,Cpu_Seconds()-cpu});}};
    measure("value",3,[&]{sink+=energy.Value(s.V);});
    measure("pressure",3,[&]{energy.Pressure(s.V,pressure);sink+=pressure[0];});
    measure("fields_without_hessian",3,[&]{plain=energy.Fields(s.V,false);sink+=plain.energy;});
    measure("fields_with_hessian",3,[&]{full=energy.Fields(s.V,true);sink+=full.energy;});
    const bool equal=plain.energy==full.energy&&plain.area==full.area&&plain.pressure==full.pressure&&plain.tension==full.tension&&plain.node_energy==full.node_energy;if(!equal)throw std::runtime_error("second derivative assembly changed first fields");
    measure("cold_hessian_build_and_apply",1,[&]{energy.Hess_Apply(s.V,w,out);sink+=out[0];},false);
    measure("cached_hessian_apply",10,[&]{energy.Hess_Apply(s.V,w,out);sink+=out[0];},false);
    measure("cached_diagonal_scan",3,[&]{for(int i=0;i<s.Size();++i)sink+=energy.Hess_Diag(s.V,i);},false);
    const Eigen::VectorXd expected=full.hessian*Eigen::Map<const Eigen::VectorXd>(w.data(),w.size());real action_error=0,action_scale=0;for(int i=0;i<s.Size();++i){action_error=std::max(action_error,std::abs(out[i]-expected[i]));action_scale=std::max(action_scale,std::abs(expected[i]));}action_error/=std::max(real(1e-30),action_scale);
    const real maximum=full.hessian.cwiseAbs().maxCoeff(),symmetry=(full.hessian-full.hessian.transpose()).cwiseAbs().maxCoeff()/std::max(real(1e-30),maximum);if(action_error>1e-12||symmetry>1e-12||!std::isfinite(sink))throw std::runtime_error("native Hessian guard failed");
    size_t nnz=0;std::ofstream entries(pattern);entries<<std::setprecision(17)<<"row,col,H_code\n";for(int j=0;j<s.Size();++j)for(int i=0;i<s.Size();++i)if(full.hessian(i,j)!=0){++nnz;entries<<i<<','<<j<<','<<full.hessian(i,j)<<'\n';}
    std::ofstream result(output);result<<std::setprecision(17)<<"{\"N\":"<<s.Size()<<",\"area_import_error\":"<<c.area_import_error<<",\"energy_import_error\":"<<c.energy_import_error<<",\"first_fields_exact\":true,\"matrix_nonzeros\":"<<nnz<<",\"matrix_density\":"<<real(nnz)/(real(s.Size())*s.Size())<<",\"cached_action_relative_error\":"<<action_error<<",\"symmetry_relative_error\":"<<symmetry<<",\"sink\":"<<sink<<",\"timings\":{";bool first=true;for(const auto& row:times){if(!first)result<<',';first=false;result<<'"'<<row.first<<"\":[";for(size_t k=0;k<row.second.size();++k){if(k)result<<',';result<<"{\"wall_s\":"<<row.second[k].wall<<",\"cpu_s\":"<<row.second[k].cpu<<'}';}result<<']';}result<<"}}\n";std::cout<<"original drainage profile N="<<s.Size()<<" nnz="<<nnz<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
