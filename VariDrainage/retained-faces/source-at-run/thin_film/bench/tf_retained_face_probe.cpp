// Fixed original/candidate geometry observations, no physical time stepping.
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "vf/thinfilm/VoronoiRetainedFacesCandidate.h"
using namespace vf;using namespace vf::thinfilm;
int main(int argc,char** argv){try{
    if(argc!=3)return 2;std::ifstream input(argv[1]);
    real x0,y0,x1,y1;int n;input>>x0>>y0>>x1>>y1>>n;if(!input||n<1)throw std::runtime_error("invalid input header");
    const TransferPolygon domain={{x0,y0},{x1,y0},{x1,y1},{x0,y1}};
    std::vector<TransferPoint> sites;std::vector<real> areas;
    for(int i=0;i<n;++i){int id;real x,y,a;input>>id>>x>>y>>a;if(!input||id!=i)throw std::runtime_error("invalid site row");sites.emplace_back(x,y);areas.push_back(a);}
    const std::filesystem::path out(argv[2]);if(std::filesystem::exists(out))throw std::runtime_error("refuse prior output");std::filesystem::create_directories(out);
    std::ofstream csv(out/"observations.csv");csv<<std::setprecision(17)<<"variant,completed,normal_closure,area_closure,orthogonality_certificate,area_import_error,pairs,faces,minimum_face_length,message\n";
    std::vector<VoronoiDrainageGeometry> stored(2);bool ok[2]={false,false};
    for(int variant=0;variant<2;++variant){const std::string name=variant?"candidate":"original";
        try{auto g=variant?Make_Voronoi_Drainage_RetainedFaces(sites,domain):Make_Voronoi_Drainage(sites,domain);real error=0,minFace=std::numeric_limits<real>::infinity();
            for(int i=0;i<n;++i)if(areas[i]>0)error=std::max(error,std::abs(g.graph.a[i]/areas[i]-1));
            std::ofstream faces(out/(name+"_faces.csv"));faces<<std::setprecision(17)<<"i,j,x,y,nx,ny,length\n";
            for(const auto& f:g.faces){minFace=std::min(minFace,f.length);faces<<f.i<<','<<f.j<<','<<f.midpoint.x()<<','<<f.midpoint.y()<<','<<f.normal.x()<<','<<f.normal.y()<<','<<f.length<<'\n';}
            csv<<name<<",1,"<<g.normal_closure<<','<<g.area_closure<<','<<g.orthogonality_certificate<<','<<error<<','<<g.graph.pairs.Size()<<','<<g.faces.size()<<','<<minFace<<",\"\"\n";stored[variant]=std::move(g);ok[variant]=true;
        }catch(const std::exception& e){csv<<name<<",0,-1,-1,-1,-1,-1,-1,-1,"<<std::quoted(std::string(e.what()))<<'\n';}
    }
    std::ofstream comparison(out/"comparison.json");comparison<<std::setprecision(17);
    if(ok[0]&&ok[1]){const auto& a=stored[0];const auto& b=stored[1];bool cells=true,areasEqual=a.graph.a==b.graph.a;
        for(int i=0;i<n;++i){if(a.cells[i].size()!=b.cells[i].size()){cells=false;continue;}for(size_t k=0;k<a.cells[i].size();++k)cells&=(a.cells[i][k]-b.cells[i][k]).norm()==0.;}
        bool pairs=a.graph.pairs.pi==b.graph.pairs.pi&&a.graph.pairs.pj==b.graph.pairs.pj,conductance=a.graph.c==b.graph.c;
        comparison<<"{\"both_admitted\":true,\"cells_exact\":"<<(cells?"true":"false")<<",\"areas_exact\":"<<(areasEqual?"true":"false")<<",\"pairs_exact\":"<<(pairs?"true":"false")<<",\"conductances_exact\":"<<(conductance?"true":"false")<<",\"additional_pairs\":"<<int(b.graph.pairs.Size())-int(a.graph.pairs.Size())<<"}\n";
    }else comparison<<"{\"both_admitted\":false}\n";
    std::cout<<"Observed original and retained-face candidate on N="<<n<<"; no time steps.\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
