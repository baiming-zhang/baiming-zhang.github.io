// Passive common-transfer assay on genuine moving point-Voronoi support plus
// a refined rectangular edge patch. No force, drainage or physical time advance.
#include <filesystem>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "vf/thinfilm/OverlapTransfer.h"
#include "vf/thinfilm/Units.h"
using namespace vf;using namespace vf::thinfilm;
constexpr real W=.0016,Y=.0012,h0=2e-6;
struct Layout{std::vector<TransferPolygon> cells;std::vector<int> owner;std::vector<TransferPoint> sites;};
Layout layout(int n,int phase,bool moving){
    Layout r;const real cut=W*(moving?.62-.07*std::sin(phase*kPi/3):.62);
    const int ny=3*n/4,ne=n/2;std::vector<TransferPoint> seeds;
    for(int j=0;j<ny;++j)for(int i=0;i<n;++i){
        real x=(i+.5)*cut/n,y=(j+.5)*Y/ny;
        x+=.18*cut/n*std::sin(1.4*j+.7*i+.9*phase);
        y+=.18*Y/ny*std::cos(.8*j+1.1*i+.6*phase);
        seeds.emplace_back(x,y);
    }
    r.cells=Transfer_Voronoi(seeds,Transfer_Rect(0,0,cut,Y));r.sites=seeds;r.owner.assign(seeds.size(),0);
    // Nonmatching resolution along the shared interface; exact shared cut line.
    for(int j=0;j<2*ny;++j)for(int i=0;i<ne;++i){
        const real xa=cut+(W-cut)*i/ne,xb=cut+(W-cut)*(i+1)/ne;
        const real ya=Y*j/(2*ny),yb=Y*(j+1)/(2*ny);
        r.cells.push_back(Transfer_Rect(xa,ya,xb,yb));r.sites.emplace_back(.5*(xa+xb),.5*(ya+yb));r.owner.push_back(1);
    }
    return r;
}
TransferInventory field(const TransferPoint& x,int profile,bool shear,const FilmMaterial& M){
    const real a=x.x()/W,b=x.y()/Y;
    const real f=profile==0?0.:profile==1?.25*std::sin(2*kPi*a)*std::cos(2*kPi*b):
        -.4*std::exp(-std::pow((a-.62)/.06,2))*(1.+.1*std::cos(2*kPi*b));
    TransferInventory q;q.V=h0*(1.+f);q.S=.45*M.Gamma_inf*(1.+.2*f);q.B=M.c_bulk*(1.+.1*f)*q.V;
    const Vec3 u(.01+(shear?.004*std::sin(2*kPi*b):0.),shear?.003*std::cos(2*kPi*a):0.,0.);
    q.P=M.rho*q.V*u;return q;
}
// Degree-2 triangle quadrature for independent smooth prescribed-field integrals;
// these are diagnostic inputs, never a solver law. Refine geometry, not quadrature.
TransferInventory integrate(const TransferPolygon& p,int profile,bool shear,const FilmMaterial& M){
    TransferInventory q;
    for(size_t k=1;k+1<p.size();++k){const TransferPoint a=p[0],b=p[k],c=p[k+1];
        const real area=.5*Transfer_Cross(b-a,c-a);
        const std::array<TransferPoint,3> quadrature={(4*a+b+c)/6.,(a+4*b+c)/6.,(a+b+4*c)/6.};
        for(const auto& x:quadrature)q+=field(x,profile,shear,M).Scaled(area/3.);
    }return q;
}
TransferInventory total(const std::vector<TransferInventory>& q,const std::vector<int>& labels,int z=-1){
    TransferInventory t;for(size_t i=0;i<q.size();++i)if(z<0||labels[i]==z)t+=q[i];return t;
}
real kinetic(const std::vector<TransferInventory>& q,real rho){real k=0.;for(const auto& v:q)k+=v.P.squaredNorm()/(2*rho*v.V);return k;}
int main(int argc,char** argv){try{
    const std::filesystem::path folder=argc>1?argv[1]:"results/local-runs/overlap_transfer_20261111";
    std::filesystem::create_directories(folder);const auto M=SDS_Solution();const auto domain=Transfer_Rect(0,0,W,Y);
    std::cout<<std::setprecision(17)<<"name,n,profile,shear,moving,cycles,N0,N1,overlaps,geometric_source_residual,geometric_destination_residual,V_balance,Q_balance,P_balance,ledger_V,ledger_Q,ledger_P,h_L2,depletion_L2,theta_L2,c_L2,uniform_velocity_L2,boost_commutation,kinetic_change,feasible\n";
    for(int n:{8,16,32})for(bool moving:{false,true}){
        std::vector<Layout> mesh;for(int phase=0;phase<=4;++phase)mesh.push_back(layout(n,phase,moving));
        std::vector<OverlapTransfer> maps;real sr=0,dr=0;size_t overlaps=0;
        for(int k=0;k<4;++k){maps.push_back(Make_Overlap_Transfer(mesh[k].cells,mesh[k+1].cells,domain));
            sr=std::max(sr,maps.back().source_closure);dr=std::max(dr,maps.back().destination_closure);}
        for(int profile:{0,1,2})for(bool shear:{false,true})for(int cycles:{1,4}){
            const auto& first=mesh[0];const auto& last=mesh[cycles];std::vector<TransferInventory> q;
            for(const auto& p:first.cells)q.push_back(integrate(p,profile,shear,M));
            const auto initial=total(q,first.owner),edge0=total(q,first.owner,1);const real k0=kinetic(q,M.rho);
            TransferInventory flux;real boosterr=0.;overlaps=0;
            for(int k=0;k<cycles;++k){flux+=Transfer_Ownership(maps[k],q,mesh[k].owner,mesh[k+1].owner,0,1);
                auto boosted=q;const Vec3 boost(.017,-.013,.009);for(auto& v:boosted)v.P+=M.rho*v.V*boost;
                const auto shifted=maps[k].Apply(boosted);q=maps[k].Apply(q);overlaps+=maps[k].entries.size();
                for(size_t j=0;j<q.size();++j)boosterr=std::max(boosterr,(shifted[j].P-q[j].P-M.rho*q[j].V*boost).norm()/(M.rho*q[j].V*boost.norm()));
            }
            const auto final=total(q,last.owner),edge1=total(q,last.owner,1);real eh=0,dh=0,dep=0,et=0,dt=0,ec=0,dc=0,ev=0;bool feasible=true;
            const std::string name="n"+std::to_string(n)+"_p"+std::to_string(profile)+"_s"+std::to_string(shear)+"_m"+std::to_string(moving)+"_c"+std::to_string(cycles);
            std::ofstream cells;if(n==16&&profile==2&&!shear&&moving&&cycles==4){cells.open(folder/(name+".csv"));cells<<std::setprecision(17)<<"cell,owner,x_m,y_m,area_m2,V_m3,S_mol,B_mol,h_m,h_reference_m,theta,c_mM,ux_m_s,uy_m_s\n";}
            for(size_t j=0;j<q.size();++j){const real area=Transfer_Area(last.cells[j]);const auto exact=integrate(last.cells[j],profile,shear,M);
                const real h=q[j].V/area,ht=exact.V/area,theta=q[j].S/(area*M.Gamma_inf),tt=exact.S/(area*M.Gamma_inf),c=q[j].B/q[j].V,ct=exact.B/exact.V;
                const Vec3 vel=q[j].P/(M.rho*q[j].V);
                eh+=area*std::pow(h-ht,2);dh+=area*ht*ht;dep+=area*std::pow(ht-h0,2);et+=area*std::pow(theta-tt,2);dt+=area*tt*tt;
                ec+=area*std::pow(c-ct,2);dc+=area*ct*ct;ev+=q[j].V*(vel-Vec3(.01,0,0)).squaredNorm();
                feasible=feasible&&h>0.&&theta>0.&&theta<1.&&c>0.&&std::isfinite(h)&&vel.allFinite();
                if(cells.is_open())cells<<j<<','<<last.owner[j]<<','<<last.sites[j].x()<<','<<last.sites[j].y()<<','<<area<<','<<q[j].V<<','<<q[j].S<<','<<q[j].B<<','<<h<<','<<ht<<','<<theta<<','<<c<<','<<vel.x()<<','<<vel.y()<<'\n';
            }
            const real q0=2*initial.S+initial.B,pscale=M.rho*initial.V*.01;
            std::cout<<name<<','<<n<<','<<profile<<','<<shear<<','<<moving<<','<<cycles<<','<<first.cells.size()<<','<<last.cells.size()<<','<<overlaps<<','<<sr<<','<<dr<<','
                <<std::abs(final.V-initial.V)/initial.V<<','<<std::abs(2*final.S+final.B-q0)/q0<<','<<(final.P-initial.P).norm()/pscale<<','
                <<std::abs(edge1.V-edge0.V-flux.V)/initial.V<<','<<std::abs(2*(edge1.S-edge0.S-flux.S)+edge1.B-edge0.B-flux.B)/q0<<','<<(edge1.P-edge0.P-flux.P).norm()/pscale<<','
                <<std::sqrt(eh/dh)<<','<<(profile?std::sqrt(eh/dep):0.)<<','<<std::sqrt(et/dt)<<','<<std::sqrt(ec/dc)<<','<<(!shear?std::sqrt(ev/(final.V*.0001)):0.)<<','<<boosterr<<','<<kinetic(q,M.rho)/k0-1.<<','<<feasible<<'\n';
        }
        // Save actual support polygons once, for reproducible disc/geometry diagnostics.
        if(n==16&&moving){std::ofstream polygons(folder/"support_n16_moving.csv");polygons<<std::setprecision(17)<<"phase,cell,owner,vertex,x_m,y_m\n";
            for(int phase:{0,4})for(size_t j=0;j<mesh[phase].cells.size();++j)for(size_t v=0;v<mesh[phase].cells[j].size();++v)
                polygons<<phase<<','<<j<<','<<mesh[phase].owner[j]<<','<<v<<','<<mesh[phase].cells[j][v].x()<<','<<mesh[phase].cells[j][v].y()<<'\n';}
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
