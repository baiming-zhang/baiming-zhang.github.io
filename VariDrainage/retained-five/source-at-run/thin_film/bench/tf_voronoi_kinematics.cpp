#include "vf/thinfilm/MovingInterfaceLedger.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
using namespace vf;using namespace vf::thinfilm;
// Read-only geometry observer in SI. No film step, force or transport changes.
int main(int argc,char** argv){
    if(argc!=4)return 2;const std::filesystem::path input(argv[1]),output(argv[3]);if(std::filesystem::exists(output))return 3;const real cut=std::stod(argv[2]);
    std::ifstream in(input);if(!in)return 4;std::string line;std::getline(in,line);FilmMatState a;
    while(std::getline(in,line)){if(line.empty())continue;std::stringstream row(line);std::vector<real> values;std::string field;while(std::getline(row,field,','))values.push_back(std::stod(field));if(values.size()!=5)return 5;a.x.emplace_back(values[1],values[2],0);a.v.emplace_back(values[3],values[4],0);a.fixed.push_back(0);}
    if(a.Size()<3)return 6;FilmMechParams p;const auto g=Moving_Voronoi_Geometry(a,Transfer_Rect(0,0,320e-6,160e-6),p);const auto J=Moving_Voronoi_Area_Jacobian(g);VecX velocity(2*a.Size());for(int i=0;i<a.Size();++i)velocity.segment<2>(2*i)=a.v[i].head<2>();const VecX rate=J*velocity;const auto weights=Moving_Interface_Weights(g,cut);
    std::ofstream out(output);out<<std::setprecision(17)<<"id,A_m2,w,Jv_m2_s\n";for(int i=0;i<a.Size();++i)out<<i<<','<<g.graph.a[i]<<','<<weights[i]<<','<<rate[i]<<'\n';std::cout<<"native geometry observer N="<<a.Size()<<"\n";return 0;
}
