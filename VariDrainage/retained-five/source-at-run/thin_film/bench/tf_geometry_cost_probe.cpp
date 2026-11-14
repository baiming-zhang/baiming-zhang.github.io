#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <cstring>
#include "vf/thinfilm/VoronoiDrainage.h"
#include "vf/thinfilm/VoronoiClipProbe.h"
#include "vf/thinfilm/Units.h"
using namespace vf;using namespace vf::thinfilm;
int main(int argc,char** argv){try{
    if(argc!=4)return 2;std::filesystem::path input(argv[1]),output(argv[3]);if(std::filesystem::exists(output))return 3;
    const int repeats=std::stoi(argv[2]);if(repeats<3||repeats>15)return 4;Set_Threads(2);
    std::ifstream in(input);if(!in)return 5;std::string line,word;std::getline(in,line);std::stringstream head(line);std::map<std::string,int> columns;int k=0;
    while(std::getline(head,word,',')){if(!word.empty()&&word.back()=='\r')word.pop_back();columns[word]=k++;}
    FilmScales sc;std::vector<TransferPoint> sites;
    while(std::getline(in,line)){if(line.empty())continue;std::stringstream row(line);std::vector<real> values;while(std::getline(row,word,','))values.push_back(std::stod(word));sites.emplace_back(values.at(columns.at("x_m"))/sc.L,values.at(columns.at("y_m"))/sc.L);}
    const auto domain=Transfer_Rect(0,0,320e-6/sc.L,160e-6/sc.L);
    auto reference=Transfer_Voronoi(sites,domain);VoronoiClipCounts counts;auto candidate=Probe_Voronoi_Unchanged_Cuts(sites,domain,&counts);
    bool exact=reference.size()==candidate.size();size_t vertices=0;
    for(size_t i=0;i<reference.size()&&exact;++i){exact=reference[i].size()==candidate[i].size();
        for(size_t j=0;j<reference[i].size()&&exact;++j){exact=std::memcmp(reference[i][j].data(),candidate[i][j].data(),2*sizeof(real))==0;++vertices;}}
    if(!exact)throw std::runtime_error("candidate vertex bytes/order differ");
    Validate_Transfer_Partition(reference,domain,1e-10);Validate_Transfer_Partition(candidate,domain,1e-10);
    // Same post-cleanup input as full native validation; numeric cleanup is unchanged.
    auto cleaned=reference;for(auto& cell:cleaned)Clean_Drainage_Polygon(cell);
    std::map<std::string,std::vector<double>> times;real sink=0;
    auto time=[&](const std::string& name,auto&& fn){const auto t=std::chrono::steady_clock::now();fn();times[name].push_back(std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count());};
    // Warm both variants, then alternate their order to reduce monotonic load bias.
    Probe_Voronoi_Unchanged_Cuts(sites,domain);Transfer_Voronoi(sites,domain);Make_Voronoi_Drainage(sites,domain);
    for(int r=0;r<repeats;++r){
        auto old=[&]{time("reference_clip",[&]{auto p=Transfer_Voronoi(sites,domain);sink+=Transfer_Area(p[0]);});};
        auto next=[&]{time("candidate_clip",[&]{auto p=Probe_Voronoi_Unchanged_Cuts(sites,domain);sink+=Transfer_Area(p[0]);});};
        if(r%2){next();old();}else{old();next();}
        time("partition_validation",[&]{Validate_Transfer_Partition(cleaned,domain,1e-10);sink+=Transfer_Area(cleaned[0]);});
        time("full_reference_geometry",[&]{auto g=Make_Voronoi_Drainage(sites,domain);sink+=g.graph.a[0];});
    }
    std::ofstream out(output);out<<std::setprecision(17)<<"{\"N\":"<<sites.size()<<",\"vertex_bytes_equal\":true,\"vertices\":"<<vertices<<",\"attempted\":"<<counts.attempted<<",\"unchanged\":"<<counts.unchanged<<",\"changed\":"<<counts.changed<<",\"sink\":"<<sink<<",\"timings\":{";
    bool first=true;for(const auto& row:times){if(!first)out<<',';first=false;out<<'"'<<row.first<<"\":[";for(size_t i=0;i<row.second.size();++i){if(i)out<<',';out<<row.second[i];}out<<']';}out<<"}}\n";
    std::cout<<"exact geometry probe N="<<sites.size()<<" unchanged="<<counts.unchanged<<'/'<<counts.attempted<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
