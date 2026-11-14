// SI fixed-support capillary drainage plus passive dissolved-tracer transport.
// The local patch is an online fine part of ONE orthogonal Voronoi solve.
// No meniscus/free-surface motion, adsorption, DLVO or fitted drainage law.
#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "vf/thinfilm/VoronoiDrainage.h"
#include "vf/thinfilm/Units.h"
using namespace vf;using namespace vf::thinfilm;
constexpr real Width=.0012,Height=.0008,Hbase=6e-6;
using Pair=std::array<real,2>;
Pair prescribed(const TransferPoint& x,const FilmScales& sc,const CodeParams& cp,const FilmMaterial& M){
    const real a=x.x()*sc.L/Width,b=x.y()*sc.L/Height;
    const real h=Hbase*(1.+.08*std::cos(kPi*a)*std::cos(2*kPi*b)-.35*std::exp(-std::pow((a-.8)*Width/60e-6,2))*(1.+.12*std::cos(2*kPi*b)))/sc.H;
    const real c=M.c_bulk*(1.+.2*std::cos(kPi*a)*std::cos(kPi*b))/cp.sol.c_unit;
    return {h,h*c};
}
// Positive degree-five Dunavant triangle rule, diagnostic initialization only.
Pair integrate(const TransferPolygon& p,const FilmScales& sc,const CodeParams& cp,const FilmMaterial& M){
    Pair sum{};const std::array<std::array<real,4>,7> rule={{{1./3,1./3,1./3,.225},
        {.059715871789770,.470142064105115,.470142064105115,.132394152788506},
        {.470142064105115,.059715871789770,.470142064105115,.132394152788506},
        {.470142064105115,.470142064105115,.059715871789770,.132394152788506},
        {.797426985353087,.101286507323456,.101286507323456,.125939180544827},
        {.101286507323456,.797426985353087,.101286507323456,.125939180544827},
        {.101286507323456,.101286507323456,.797426985353087,.125939180544827}}};
    for(size_t k=1;k+1<p.size();++k){const real area=.5*Transfer_Cross(p[k]-p[0],p[k+1]-p[0]);
        for(const auto& v:rule){const auto q=prescribed(v[0]*p[0]+v[1]*p[k]+v[2]*p[k+1],sc,cp,M);sum[0]+=area*v[3]*q[0];sum[1]+=area*v[3]*q[1];}}
    return sum;
}
int main(int argc,char** argv){try{
    if(argc!=8){std::cerr<<"mode coarse_um fine_um patch_width_um dt_s T_s output_directory\n";return 2;}
    const bool local=std::string(argv[1])=="local";if(!local&&std::string(argv[1])!="uniform")return 2;
    const real coarse=std::stod(argv[2])*1e-6,fine=std::stod(argv[3])*1e-6,patch=std::stod(argv[4])*1e-6,dt=std::stod(argv[5]),T=std::stod(argv[6]);
    if(!(coarse>0&&fine>0&&patch>0&&patch<Width&&dt>0&&T>=dt))throw std::invalid_argument("invalid assay parameters");
    const int steps=int(std::llround(T/dt));if(std::abs(steps*dt-T)>1e-10*T)throw std::invalid_argument("unaligned physical horizon");
    const std::filesystem::path out=argv[7];std::filesystem::create_directories(out);Set_Threads(2);
    const auto M=SDS_Solution();const FilmScales sc;const auto cp=To_Code(M,sc);const real volume_unit=sc.H*sc.L*sc.L,amount_unit=M.Gamma_inf*sc.L*sc.L,energy_unit=sc.S*sc.L*sc.L;
    const real cut=(Width-patch)/sc.L;std::vector<TransferPoint> sites;
    auto append=[&](real x0,real x1,real spacing){const int nx=std::max(1,int(std::lround((x1-x0)/spacing))),ny=std::max(1,int(std::lround(Height/spacing)));
        for(int j=0;j<ny;++j)for(int i=0;i<nx;++i)sites.emplace_back((x0+(i+.5)*(x1-x0)/nx)/sc.L,(j+.5)*Height/ny/sc.L);};
    if(local){append(0.,Width-patch,coarse);append(Width-patch,Width,fine);}else append(0.,Width,coarse);
    const auto start=std::chrono::steady_clock::now();
    const auto geom=Make_Voronoi_Drainage(sites,Transfer_Rect(0,0,Width/sc.L,Height/sc.L));const auto& G=geom.graph;
    std::vector<real> V(G.N),B(G.N),S(G.N);std::vector<int> label(G.N);
    real V0=0.,B0=0.,S0=0.,edgeV0=0.,edgeB0=0.;
    for(int i=0;i<G.N;++i){const auto q=integrate(geom.cells[i],sc,cp,M);V[i]=q[0];B[i]=q[1];S[i]=M.theta0*G.a[i];label[i]=sites[i].x()>=cut;
        V0+=V[i];B0+=B[i];S0+=S[i];if(label[i]){edgeV0+=V[i];edgeB0+=B[i];}}
    const auto initialV=V,initialB=B;
    Drainage_Stepper solver;solver.prm=cp.tf;solver.prm.disj={};solver.prm.rho_g=0.; // horizontal capillary-only isolation
    solver.factored_precond=true;solver.rel_tol=1e-11;solver.newton.max_cg_iters=5000;solver.newton.max_iters=50;
    ThinFilmEnergy energy(G,solver.prm);const real E0=energy.Value(V),Q0=2*S0+B0;
    real maxV=0,maxB=0,maxQ=0,maxLV=0,maxLB=0,maxNodeV=0,maxNodeB=0,maxEnergy=0,maxStationary=0,minH=1e100,minC=1e100;
    real netV=0.,netB=0.,grossVforward=0.,grossVreverse=0.,totalDiss=0.;int failures=0,totalCG=0,forward=0,reverse=0,completed=0;
    std::vector<real> F,k,pressure;
    std::ofstream history(out/"history.csv");history<<std::setprecision(17)<<"step,t_s,V_m3,B_mol,Q_mol,edge_V_m3,edge_B_mol,E_J,half_action_J,interface_V_m3,interface_B_mol,V_balance,B_balance,edge_V_ledger,edge_B_ledger,node_V_residual,node_B_residual,energy_inequality,pressure_flux_residual,min_h_m,min_c_mM,converged,cg_iters\n";
    for(int step=0;step<steps;++step){const auto beforeV=V,beforeB=B;const real oldE=energy.Value(V);solver.Conductances(G,V,k);
        const auto ns=solver.Step(G,V,dt/sc.T,&F);totalCG+=ns.total_cg_iters;if(!ns.converged)++failures;
        Solute_Transport(G,F,V,B);energy.Pressure(V,pressure);
        const auto cross=Voronoi_Interface_Transfer(G,label,F,V,B);netV+=cross.liquid;netB+=cross.dissolved;forward+=cross.forward;reverse+=cross.reverse;
        auto predV=beforeV,predB=beforeB;real halfAction=0,stationarity=0;
        for(int q=0;q<G.pairs.Size();++q){const int i=G.pairs.pi[q],j=G.pairs.pj[q],donor=F[q]>0?i:j;const real solute=F[q]*B[donor]/V[donor];
            predV[i]-=F[q];predV[j]+=F[q];predB[i]-=solute;predB[j]+=solute;
            halfAction+=F[q]*F[q]/(2*(dt/sc.T)*k[q]);stationarity=std::max(stationarity,std::abs(F[q]-(dt/sc.T)*k[q]*(pressure[i]-pressure[j]))/V0);
            if(label[i]!=label[j]){const real f=(label[i]? -1.:1.)*F[q];if(f>0)grossVforward+=f;else grossVreverse-=f;}}
        real v=0,b=0,sv=0,sb=0,nodeV=0,nodeB=0,hmin=1e100,cmin=1e100;
        for(int i=0;i<G.N;++i){v+=V[i];b+=B[i];if(label[i]){sv+=V[i];sb+=B[i];}
            nodeV=std::max(nodeV,std::abs(V[i]-predV[i])/V0);nodeB=std::max(nodeB,std::abs(B[i]-predB[i])/B0);
            hmin=std::min(hmin,V[i]/G.a[i]*sc.H);cmin=std::min(cmin,B[i]/V[i]*cp.sol.c_unit);}
        const real en=energy.Value(V),vbal=std::abs(v-V0)/V0,bbal=std::abs(b-B0)/B0,qbal=std::abs(2*S0+b-Q0)/Q0;
        const real lv=std::abs(sv-edgeV0-netV)/V0,lb=std::abs(sb-edgeB0-netB)/B0,ineq=std::max(real(0),en-oldE+halfAction)/std::max(std::abs(E0),real(1e-300));
        maxV=std::max(maxV,vbal);maxB=std::max(maxB,bbal);maxQ=std::max(maxQ,qbal);maxLV=std::max(maxLV,lv);maxLB=std::max(maxLB,lb);maxNodeV=std::max(maxNodeV,nodeV);maxNodeB=std::max(maxNodeB,nodeB);
        maxEnergy=std::max(maxEnergy,ineq);maxStationary=std::max(maxStationary,stationarity);minH=std::min(minH,hmin);minC=std::min(minC,cmin);totalDiss+=halfAction;
        history<<step+1<<','<<(step+1)*dt<<','<<v*volume_unit<<','<<b*amount_unit<<','<<(2*S0+b)*amount_unit<<','<<sv*volume_unit<<','<<sb*amount_unit<<','<<en*energy_unit<<','<<halfAction*energy_unit<<','<<netV*volume_unit<<','<<netB*amount_unit<<','<<vbal<<','<<bbal<<','<<lv<<','<<lb<<','<<nodeV<<','<<nodeB<<','<<ineq<<','<<stationarity<<','<<hmin<<','<<cmin<<','<<ns.converged<<','<<ns.total_cg_iters<<'\n';
        completed=step+1;
        if(!ns.converged||!(hmin>0&&cmin>0)||!std::isfinite(en)){std::cerr<<"failed step "<<step+1<<'\n';break;}
    }
    std::ofstream cells(out/"cells.csv");cells<<std::setprecision(17)<<"cell,owner,x_m,y_m,area_m2,initial_V_m3,V_m3,initial_B_mol,B_mol,S_per_face_mol,initial_h_m,h_m,c_mM,pressure_Pa\n";
    for(int i=0;i<G.N;++i)cells<<i<<','<<label[i]<<','<<sites[i].x()*sc.L<<','<<sites[i].y()*sc.L<<','<<G.a[i]*sc.L*sc.L<<','<<initialV[i]*volume_unit<<','<<V[i]*volume_unit<<','<<initialB[i]*amount_unit<<','<<B[i]*amount_unit<<','<<S[i]*amount_unit<<','<<initialV[i]/G.a[i]*sc.H<<','<<V[i]/G.a[i]*sc.H<<','<<B[i]/V[i]*cp.sol.c_unit<<','<<cp.p_Pa(pressure[i])<<'\n';
    std::ofstream polygons(out/"polygons.csv");polygons<<std::setprecision(17)<<"cell,vertex,x_m,y_m\n";
    for(int i=0;i<G.N;++i)for(size_t v=0;v<geom.cells[i].size();++v)polygons<<i<<','<<v<<','<<geom.cells[i][v].x()*sc.L<<','<<geom.cells[i][v].y()*sc.L<<'\n';
    std::ofstream faces(out/"faces.csv");faces<<std::setprecision(17)<<"pair,i,j,x_m,y_m,nx,ny,length_m,interface,last_volume_m3\n";int q=0;
    for(const auto& f:geom.faces)if(f.j>=0){faces<<q<<','<<f.i<<','<<f.j<<','<<f.midpoint.x()*sc.L<<','<<f.midpoint.y()*sc.L<<','<<f.normal.x()<<','<<f.normal.y()<<','<<f.length*sc.L<<','<<(label[f.i]!=label[f.j])<<','<<F[q]*volume_unit<<'\n';++q;}
    std::ofstream summary(out/"summary.json");summary<<std::setprecision(17)<<"{\n\"N\":"<<G.N<<",\"pairs\":"<<G.pairs.Size()<<",\"steps_requested\":"<<steps<<",\"steps\":"<<completed<<",\"requested_T_s\":"<<T<<",\"T_s\":"<<completed*dt<<",\"dt_s\":"<<dt<<",\"local\":"<<(local?"true":"false")<<",\"coarse_um\":"<<coarse*1e6<<",\"fine_um\":"<<fine*1e6<<",\"width_um\":"<<patch*1e6
        <<",\"gamma_ref_N_m\":"<<M.gamma_ref()<<",\"mu_Pa_s\":"<<M.mu<<",\"rho_kg_m3\":"<<M.rho<<",\"Gamma_inf_mol_m2\":"<<M.Gamma_inf<<",\"theta0\":"<<M.theta0<<",\"bulk_scale_mM\":"<<M.c_bulk<<",\"S0_per_face_mol\":"<<S0*amount_unit<<",\"V0_m3\":"<<V0*volume_unit<<",\"B0_mol\":"<<B0*amount_unit<<",\"E0_J\":"<<E0*energy_unit<<",\"Efinal_J\":"<<energy.Value(V)*energy_unit<<",\"total_half_action_J\":"<<totalDiss*energy_unit
        <<",\"max_V_balance\":"<<maxV<<",\"max_B_balance\":"<<maxB<<",\"max_Q_balance\":"<<maxQ<<",\"max_edge_V_ledger\":"<<maxLV<<",\"max_edge_B_ledger\":"<<maxLB<<",\"max_node_V\":"<<maxNodeV<<",\"max_node_B\":"<<maxNodeB<<",\"max_energy_inequality\":"<<maxEnergy<<",\"max_pressure_flux\":"<<maxStationary
        <<",\"min_h_m\":"<<minH<<",\"min_c_mM\":"<<minC<<",\"net_interface_V_m3\":"<<netV*volume_unit<<",\"net_interface_B_mol\":"<<netB*amount_unit<<",\"forward_gross_V_m3\":"<<grossVforward*volume_unit<<",\"reverse_gross_V_m3\":"<<grossVreverse*volume_unit<<",\"forward_face_steps\":"<<forward<<",\"reverse_face_steps\":"<<reverse
        <<",\"Newton_failures\":"<<failures<<",\"cg_iterations\":"<<totalCG<<",\"area_closure\":"<<geom.area_closure<<",\"normal_closure\":"<<geom.normal_closure<<",\"orthogonality_error\":"<<geom.orthogonality_error<<",\"wall_s\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
    std::cout<<"N "<<G.N<<", "<<steps<<" physical steps, Newton failures "<<failures<<", maximum volume residual "<<maxV<<", capillary energy ratio "<<energy.Value(V)/E0<<'\n';
    return failures?1:0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
