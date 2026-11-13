// Fixed manufactured observations, ORIGINAL geometry routine, no trajectories.
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "vf/thinfilm/VoronoiDrainage.h"
using namespace vf;using namespace vf::thinfilm;

int main(int argc,char** argv){try{
    if(argc!=2)return 2;
    const std::filesystem::path root(argv[1]);
    if(std::filesystem::exists(root))throw std::runtime_error("Refuse existing diagnostic output");
    std::filesystem::create_directories(root);
    const TransferPolygon domain={{-32,-16},{32,-16},{32,16},{-32,16}};
    constexpr real tolerance=1e-10;
    const real eps=tolerance*std::sqrt(real(64*64+32*32));
    std::ofstream out(root/"observations.csv"),vertices(root/"central_edges.csv"),points(root/"sites.csv");
    out<<std::setprecision(17)<<"case,delta,original_completed,original_normal_closure,central_full_closure,central_filtered_closure,central_perimeter,short_edge_count,short_edge_total,cutoff,message\n";
    vertices<<std::setprecision(17)<<"case,k,ax,ay,bx,by,length,omitted_by_original_cutoff\n";
    points<<std::setprecision(17)<<"case,id,x,y\n";
    for(const auto& entry:std::vector<std::pair<std::string,real>>{{"unperturbed",0},{"short_face",8e-9},{"resolved_face",8e-7}}){
        const auto& name=entry.first;const real delta=entry.second;
        std::vector<TransferPoint> sites;
        for(int y=-1;y<=1;++y)for(int x=-1;x<=1;++x)sites.emplace_back(x,y);
        sites[8].y()-=delta;
        for(size_t i=0;i<sites.size();++i)points<<name<<','<<i<<','<<sites[i].x()<<','<<sites[i].y()<<'\n';
        bool completed=false;real nativeClosure=-1;std::string message;
        try{const auto g=Make_Voronoi_Drainage(sites,domain,tolerance);completed=true;nativeClosure=g.normal_closure;}
        catch(const std::exception& e){message=e.what();const std::string label="normal_closure=";const auto pos=message.find(label);if(pos!=std::string::npos)nativeClosure=std::stod(message.substr(pos+label.size()));}
        auto cells=Transfer_Voronoi(sites,domain);
        for(auto& cell:cells)Clean_Drainage_Polygon(cell);
        const auto& cell=cells[4];TransferPoint full=TransferPoint::Zero(),filtered=TransferPoint::Zero();real perimeter=0,kept=0,shortTotal=0;int count=0;
        for(size_t k=0;k<cell.size();++k){
            const auto a=cell[k],b=cell[(k+1)%cell.size()];const TransferPoint edge=b-a,normalLength(edge.y(),-edge.x());const real length=edge.norm();
            full+=normalLength;perimeter+=length;
            if(length<=eps){++count;shortTotal+=length;}else{filtered+=normalLength;kept+=length;}
            vertices<<name<<','<<k<<','<<a.x()<<','<<a.y()<<','<<b.x()<<','<<b.y()<<','<<length<<','<<(length<=eps)<<'\n';
        }
        out<<name<<','<<delta<<','<<completed<<','<<nativeClosure<<','<<full.norm()/perimeter<<','<<filtered.norm()/kept<<','<<perimeter<<','<<count<<','<<shortTotal<<','<<eps<<','<<std::quoted(message)<<'\n';
    }
    std::cout<<"Executed three fixed nine-site observations; original failures retained. No native steps.\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
