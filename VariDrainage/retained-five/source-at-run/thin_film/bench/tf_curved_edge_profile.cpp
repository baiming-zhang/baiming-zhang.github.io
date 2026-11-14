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
#include <filesystem>
#include <fstream>
#include <sstream>
#include <map>
#include <memory>
#include <iomanip>
#include <iostream>
#include "vf/thinfilm/CurvedEdgeMotion.h"
#include "vf/thinfilm/CurvedEdgeMechanicalObserver.h"
#include "vf/thinfilm/Units.h"
using namespace vf;using namespace vf::thinfilm;
static double Cpu_Seconds(){
#ifdef _WIN32
    FILETIME created,exited,kernel,user;if(!GetProcessTimes(GetCurrentProcess(),&created,&exited,&kernel,&user))throw std::runtime_error("GetProcessTimes");ULARGE_INTEGER k,u;k.LowPart=kernel.dwLowDateTime;k.HighPart=kernel.dwHighDateTime;u.LowPart=user.dwLowDateTime;u.HighPart=user.dwHighDateTime;return (k.QuadPart+u.QuadPart)*1e-7;
#else
    return double(std::clock())/CLOCKS_PER_SEC;
#endif
}
struct ProfileSample {double wall,cpu;};
struct ProfileRow {std::string name;std::vector<ProfileSample> samples;};
int main(int argc,char** argv){try{
    if(argc!=5)return 2;const std::filesystem::path input(argv[1]),output(argv[4]);if(std::filesystem::exists(output))return 3;const real dtSI=std::stod(argv[2])*1e-6;const int repeats=std::stoi(argv[3]);if(!(dtSI>0&&repeats>=1&&repeats<=10))return 4;Set_Threads(2);
    auto M=SDS_Solution(6);M.soluble=false;const FilmScales sc;const auto cp=To_Code(M,sc);FilmMechParams p;cp.Apply(p);p.tf.disj={};p.tf.rho_g=0;const real vu=sc.H*sc.L*sc.L,su=M.Gamma_inf*sc.L*sc.L,eu=sc.S*sc.L*sc.L,dt=dtSI/sc.T;
    FilmMatState s;CurvedEdgeBoundary boundary;boundary.aspect=sc.H/sc.L;boundary.clean_tension=M.gamma0/sc.S;
    std::ifstream in(input);if(!in)return 5;std::string line;std::getline(in,line);std::stringstream head(line);std::string word;std::map<std::string,int> column;int index=0;while(std::getline(head,word,',')){if(!word.empty()&&word.back()=='\r')word.pop_back();column[word]=index++;}
    real observedE=0;std::vector<real> observedA;
    while(std::getline(in,line)){if(line.empty())continue;std::stringstream row(line);std::vector<real> v;while(std::getline(row,word,','))v.push_back(std::stod(word));auto get=[&](const std::string& key){return v.at(column.at(key));};if(int(get("id"))!=s.Size())throw std::runtime_error("nonsequential input ids");const bool bath=get("bath")!=0;
        s.x.emplace_back(get("x_m")/sc.L,get("y_m")/sc.L,0);s.v.emplace_back(get("ux_m_s")*sc.T/sc.L,get("uy_m_s")*sc.T/sc.L,0);s.fixed.push_back(bath);s.V.push_back(get("V_m3")/vu);s.S.push_back(get("S_mol")/su);boundary.height.push_back(bath?get("h_m")/sc.H:0);boundary.potential.push_back(bath?cp.p_code(get("pressure_Pa")):std::numeric_limits<real>::quiet_NaN());observedA.push_back(get("A_m2"));observedE+=get("node_energy_J");}
    if(s.Size()<3)throw std::runtime_error("too few input sites");const auto domain=Transfer_Rect(0,0,320e-6/sc.L,160e-6/sc.L);const auto geom=Moving_Voronoi_Geometry(s,domain,p);const auto energy=Differentiate_Curved_Edge(s,p,domain,boundary);real areaError=0;for(int i=0;i<s.Size();++i)areaError=std::max(areaError,std::abs(geom.graph.a[i]*sc.L*sc.L/observedA[i]-1));const real energyError=std::abs(energy.energy*eu-observedE)/std::max(real(1e-30),std::abs(observedE));
    if(areaError>1e-9||energyError>1e-9)throw std::runtime_error("SI snapshot import differs from native geometry/energy");
    real sink=0;std::vector<ProfileRow> rows;
    auto measure=[&](const std::string& name,int count,auto&& f,bool warm=true){if(warm)f();ProfileRow result;result.name=name;for(int i=0;i<count;++i){const double c=Cpu_Seconds();const auto t=std::chrono::steady_clock::now();f();const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();result.samples.push_back({wall,Cpu_Seconds()-c});}rows.push_back(std::move(result));};
    measure("geometry",repeats,[&]{auto g=Moving_Voronoi_Geometry(s,domain,p);sink+=g.graph.a[0];});
    measure("energy_gradient",repeats,[&]{auto e=Differentiate_Curved_Edge(s,p,domain,boundary);sink+=e.energy+e.gradient[0].squaredNorm();});
    std::unique_ptr<MovingVoronoiProblem> metric;VecX probe;
    measure("metric_build_plus_value",repeats,[&]{metric=std::make_unique<MovingVoronoiProblem>(s,p,domain,dt,&boundary.potential,MovingVoronoiDiscretization{true,true});probe=VecX::NullaryExpr(metric->Dim(),[](Eigen::Index k){return 1e-6*std::sin(.37*(k+1));});sink+=metric->Metric_Value(probe);});
    CurvedEdgeMotionProblem problem(s,p,domain,dt,boundary);const VecX zero=VecX::Zero(problem.Dim());VecX g,h;problem.Prepare_Hessian(zero);
    measure("problem_gradient",repeats,[&]{problem.Gradient(zero,g);sink+=g.squaredNorm();});
    measure("hessian_vector",repeats,[&]{problem.Hess_Vec(probe,h);sink+=h.squaredNorm();});
    FilmMatState mechanicsState;MovingVoronoiStats mechanics;
    measure("mechanics_step",1,[&]{mechanicsState=s;mechanics=Curved_Edge_Mechanics(mechanicsState,p,domain,dt,boundary);sink+=mechanics.energy_after;},false);
    FilmMatState final,mid;MovingVoronoiCoupledStats coupled;
    measure("coupled_step",1,[&]{final=s;coupled=Curved_Edge_Coupled_Step(final,p,domain,dt,boundary,&mid);sink+=coupled.energy_after;},false);
    if(!mechanics.newton.converged||mechanics.newton.hit_total_step||!coupled.mechanics.newton.converged||coupled.mechanics.newton.hit_total_step||!coupled.drainage.converged)throw std::runtime_error("profile single-step solve did not converge");
    real matching=0;for(int i=0;i<s.Size();++i)matching=std::max(matching,(mechanicsState.x[i]-mid.x[i]).norm());if(matching>1e-10*problem.Length())throw std::runtime_error("standalone/coupled mechanics endpoints differ");
    measure("saved_step_observer",1,[&]{const auto o=Observe_Curved_Mechanics(s,mid,p,domain,boundary,dt,240e-6/sc.L);sink+=o.regions[0].residual;},false);
    if(!std::isfinite(sink))throw std::runtime_error("nonfinite profile output");
    std::ofstream out(output);out<<std::setprecision(17)<<"{\"N\":"<<s.Size()<<",\"dt_s\":"<<dtSI<<",\"area_import_error\":"<<areaError<<",\"energy_import_error\":"<<energyError<<",\"tape_nodes\":"<<energy.tape_nodes<<",\"mechanics_endpoint_difference_m\":"<<matching*sc.L<<",\"mechanics_iters\":"<<coupled.mechanics.newton.iters<<",\"mechanics_cg_iters\":"<<coupled.mechanics.newton.total_cg_iters<<",\"force_target_ratio\":"<<coupled.mechanics.final_gradient_inf/std::max(real(1e-10),coupled.mechanics.initial_gradient_inf*1e-8)<<",\"sink\":"<<sink<<",\"timings\":[";
    for(size_t k=0;k<rows.size();++k){if(k)out<<',';out<<"{\"name\":\""<<rows[k].name<<"\",\"samples\":[";for(size_t j=0;j<rows[k].samples.size();++j){if(j)out<<',';out<<"{\"wall_s\":"<<rows[k].samples[j].wall<<",\"cpu_s\":"<<rows[k].samples[j].cpu<<'}';}out<<"]}";}out<<"]}\n";
    std::cout<<"profile N="<<s.Size()<<" import area="<<areaError<<" energy="<<energyError<<"\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
