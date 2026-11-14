// Export the rejected pre-validation geometry; no repair and no PDE solve.
#include <iostream>
#include <iomanip>
#include "vf/thinfilm/OverlapTransfer.h"
using namespace vf;using namespace vf::thinfilm;
int main(){std::vector<TransferPoint> sites;
    for(int j=0;j<8;++j)for(int i=0;i<10;++i)sites.emplace_back((i+.5)*96/10.,(j+.5)*80/8.);
    for(int j=0;j<16;++j)for(int i=0;i<5;++i)sites.emplace_back(96+(i+.5)*24/5.,(j+.5)*80/16.);
    const auto polys=Transfer_Voronoi(sites,Transfer_Rect(0,0,120,80));
    std::cout<<std::setprecision(17)<<"cell,vertex,x,y\n";
    for(size_t i=0;i<polys.size();++i)for(size_t j=0;j<polys[i].size();++j)std::cout<<i<<','<<j<<','<<polys[i][j].x()<<','<<polys[i][j].y()<<'\n';
}
