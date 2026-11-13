#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include "vf/thinfilm/CurvedEdgeMotion.h"
#include "vf/thinfilm/CurvedEdgePreparation.h"
#include "vf/thinfilm/CurvedEdgeMechanicalObserver.h"
#include "vf/thinfilm/MovingInterfaceLedger.h"
#include "vf/thinfilm/Units.h"
using namespace vf;using namespace vf::thinfilm;
int main(int argc,char** argv){try{
    if(argc<8||argc>13){std::cerr<<"local/uniform dx_um dt_us T_us radius_um amplitude output [c_bulk_mM] [projected/balanced/variational] [save_stride] [refine_start_um] [mechanics_observer]\n";return 2;}
    const std::string mode=argv[1];if(mode!="local"&&mode!="uniform")return 2;const bool local=mode=="local";
    const real dx=std::stod(argv[2])*1e-6,dtSI=std::stod(argv[3])*1e-6,T=std::stod(argv[4])*1e-6,radius=std::stod(argv[5])*1e-6,amp=std::stod(argv[6]);
    constexpr real W=320e-6,Y=160e-6,Cut=240e-6,H0=6e-6;const int steps=int(std::llround(T/dtSI));
    if(!(dx>0&&dtSI>0&&steps>0&&std::abs(steps*dtSI-T)<1e-10*T&&radius>W-Cut&&std::abs(amp)<.01))throw std::invalid_argument("curved assay input");
    const std::filesystem::path out=argv[7];std::filesystem::create_directories(out);Set_Threads(2);
    const real concentration=argc>=9?std::stod(argv[8]):19.4;const std::string preparation=argc>=10?argv[9]:"projected";
    if(!(concentration>0)||(preparation!="projected"&&preparation!="balanced"&&preparation!="variational"))throw std::invalid_argument("material/preparation input");
    const int saveStride=argc>=11?std::stoi(argv[10]):1;if(saveStride<1)throw std::invalid_argument("save stride");
    real refineStart=argc>=12?std::stod(argv[11])*1e-6:Cut;
    if(std::abs(refineStart-Cut)<1e-15)refineStart=Cut;
    if(!(refineStart>0&&refineStart<=Cut))throw std::invalid_argument("refinement start outside film-side range");
    // Numerical mesh extent only: physical traces, bath, region cut and all
    // material/force/flux laws stay at their original physical locations.
    const bool observe=argc>=13&&std::stoi(argv[12])!=0;
    auto M=SDS_Solution(concentration);M.soluble=false;const FilmScales sc;const auto cp=To_Code(M,sc);FilmMechParams p;cp.Apply(p);p.tf.disj={};p.tf.rho_g=0;
    const real dt=dtSI/sc.T,vu=sc.H*sc.L*sc.L,su=M.Gamma_inf*sc.L*sc.L,eu=sc.S*sc.L*sc.L,pu=sc.S/sc.H;
    const auto domain=Transfer_Rect(0,0,W/sc.L,Y/sc.L);FilmMatState s;
    auto append=[&](real x0,real x1,real spacing){const int nx=int(std::lround((x1-x0)/spacing)),ny=int(std::lround(Y/spacing));for(int j=0;j<ny;++j)for(int i=0;i<nx;++i){const real hx=(x1-x0)/nx,hy=Y/ny;const real x=x0+(i+.5)*hx+.015*hx*std::sin(.7*i+1.3*j),y=(j+.5)*hy+.015*hy*std::cos(1.1*i+.8*j);s.x.emplace_back(x/sc.L,y/sc.L,0);s.fixed.push_back(0);}};
    if(local){append(0,refineStart,2*dx);append(refineStart,W,dx);}else append(0,W,dx);
    if(preparation=="variational"){if(amp!=0)throw std::invalid_argument("variational preparation uses uniform initial activity, amplitude must be zero");
        const int traceNy=int(std::lround(Y/dx));for(int j=0;j<traceNy;++j){s.x.emplace_back(Cut/sc.L,(j+.5)*(Y/traceNy)/sc.L,0);s.fixed.push_back(0);}}
    const int interior=s.Size(),ny=int(std::lround(Y/dx));for(int j=0;j<ny;++j){s.x.emplace_back(W/sc.L,(j+.5)*(Y/ny)/sc.L,0);s.fixed.push_back(1);}
    s.v.assign(s.Size(),Vec3::Zero());auto geometry=Moving_Voronoi_Geometry(s,domain,p);
    for(int i=0;i<s.Size();++i){const real x=s.x[i].x()*sc.L,y=s.x[i].y()*sc.L;s.V.push_back(H0/sc.H*geometry.graph.a[i]);s.S.push_back(M.theta0*(1+amp*std::cos(kPi*x/W)*std::cos(2*kPi*y/Y))*geometry.graph.a[i]);}
    CurvedEdgeBoundary boundary;boundary.aspect=sc.H/sc.L;boundary.clean_tension=M.gamma0/sc.S;boundary.height.assign(s.Size(),0);boundary.potential.assign(s.Size(),std::numeric_limits<real>::quiet_NaN());
    // Boundary geometry from specified radius at the physical x=W wall. Only
    // endpoint data: the interior starts flat and every interior h is evolved.
    const real length=W-Cut,farH=H0+2*length*length/(radius+std::sqrt(radius*radius-length*length)),mu=-M.gamma_ref()/radius;
    for(int i=interior;i<s.Size();++i){boundary.height[i]=farH/sc.H;boundary.potential[i]=cp.p_code(mu);}
    // Preparation control, not evolving adsorption: redistribute the SAME
    // initial total S over curved dual-face area, retaining the same small
    // spatial perturbation. Every amount remains positive; V/x unchanged.
    if(preparation=="balanced"){const auto areas=CurvedEdgeEnergy(geometry.graph,p,s.S,boundary).Fields(s.V,false).area;real oldAmount=0,newAmount=0;std::vector<real> weights(s.Size());
        for(int i=0;i<s.Size();++i){oldAmount+=s.S[i];weights[i]=s.S[i]*areas[i]/geometry.graph.a[i];newAmount+=weights[i];}
        for(int i=0;i<s.Size();++i)s.S[i]=oldAmount/newAmount*weights[i];}
    CurvedEdgePreparationStats prepared;
    if(preparation=="variational")prepared=Prepare_Curved_Edge(s,geometry.graph,p,boundary,Cut/sc.L,H0/sc.H,M.theta0);
    const auto initial=s;real v0=0,s0=0;for(int i=0;i<s.Size();++i){v0+=s.V[i];s0+=s.S[i];}const real e0=CurvedEdgeEnergy(geometry.graph,p,s.S,boundary).Value(s.V),k0=Moving_Voronoi_Kinetic(s,p.rho);
    std::ofstream fields(out/"states.csv"),faces(out/"faces.csv"),hist(out/"history.csv");fields<<std::setprecision(17)<<"step,stage,id,x_m,y_m,A_m2,As_m2,V_m3,S_mol,h_m,theta,pressure_Pa,tension_N_m,node_energy_J,w,bath,ux_m_s,uy_m_s\n";faces<<std::setprecision(17)<<"step,stage,i,j,c,F_m3\n";hist<<std::setprecision(17)<<"step,t_s,E_J,K_J,visc_J,halfDrain_J,mix_J,boundaryWork_J,global_budget,V_res,S_res,materialV_m3,relativeV_m3,regionV_res,materialS_mol,regionS_res,max_displacement_m\n";
    auto dump=[&](int step,int stage,const FilmMatState& state,const VoronoiDrainageGeometry& g,const std::vector<real>* F){const auto f=CurvedEdgeEnergy(g.graph,p,state.S,boundary).Fields(state.V,false);const auto w=Moving_Interface_Weights(g,Cut/sc.L);
        for(int i=0;i<s.Size();++i)fields<<step<<','<<stage<<','<<i<<','<<state.x[i].x()*sc.L<<','<<state.x[i].y()*sc.L<<','<<g.graph.a[i]*sc.L*sc.L<<','<<f.area[i]*sc.L*sc.L<<','<<state.V[i]*vu<<','<<state.S[i]*su<<','<<boundary.Height(i,state.V[i],g.graph.a[i])*sc.H<<','<<state.S[i]/f.area[i]<<','<<f.pressure[i]*pu<<','<<f.tension[i]*sc.S<<','<<f.node_energy[i]*eu<<','<<w[i]<<','<<boundary.Bath(i)<<','<<state.v[i].x()*sc.L/sc.T<<','<<state.v[i].y()*sc.L/sc.T<<'\n';
        for(int q=0;q<g.graph.pairs.Size();++q)faces<<step<<','<<stage<<','<<g.graph.pairs.pi[q]<<','<<g.graph.pairs.pj[q]<<','<<g.graph.c[q]<<','<<(F?(*F)[q]*vu:0)<<'\n';};
    if(preparation=="variational"){
        fields.flush();faces.flush();std::string sh,fh;{std::ifstream in(out/"states.csv");std::getline(in,sh);}{std::ifstream in(out/"faces.csv");std::getline(in,fh);}
        std::ofstream init(out/"initial_states.csv"),initFaces(out/"initial_faces.csv");init<<std::setprecision(17)<<sh<<'\n';initFaces<<std::setprecision(17)<<fh<<'\n';
        fields.swap(init);faces.swap(initFaces);dump(0,0,s,geometry,nullptr);fields.swap(init);faces.swap(initFaces);
        std::ofstream info(out/"preparation.json");info<<std::setprecision(17)<<"{\"iterations\":"<<prepared.iterations<<",\"pressure_residual\":"<<prepared.pressure_residual<<",\"free_volume_m3\":"<<prepared.free_volume*vu<<",\"S_mol\":"<<prepared.surface_amount*su<<",\"trace_x_m\":"<<Cut<<",\"trace_released\":true,\"constant_activity_initial_only\":true}\n";
    }
    std::ofstream mech,regional,nodes;
    if(observe){mech.open(out/"mechanics.csv");regional.open(out/"regional_mechanics.csv");nodes.open(out/"mechanical_nodes.csv");
        mech<<std::setprecision(17)<<"step,t_s,iters,cg_iters,negative_curvature,gradient_history_size,g0_N,g1_N,target_N,K0_J,K1_J,visc_J,inertia_loss_J,incremental_change_J,last_step_m,mixing_J\n";
        regional<<std::setprecision(17)<<"step,region,dE_J,dK_J,ownershipE_J,ownershipK_J,inertia_loss_J,viscous_D_J,inertial_work_J,viscous_node_J,capillary_node_J,capillary_local_J,capillary_remainder_J,capillary_port_J,viscous_port_J,residual_J,recovered_viscous_port_J,cut_length_m,reconstructed_force_N,viscous_action_difference_J,coordinate_roundtrip_m\n";
        nodes<<std::setprecision(17)<<"step,id,fixed,bath,x0_m,y0_m,x1_m,y1_m,A0_m2,A1_m2,V_m3,S_mol,w0,w1,u0x_m_s,u0y_m_s,u1x_m_s,u1y_m_s,e0_J,e1_J,gix_N,giy_N,gvx_N,gvy_N,gcx_N,gcy_N,grx_N,gry_N\n";}
    real visc=0,drain=0,mix=0,work=0,worstE=0,worstV=0,worstS=0,worstRV=0,worstRS=0,maxmove=0;int completed=0,failures=0;const auto start=std::chrono::steady_clock::now();
    for(int step=1;step<=steps;++step){const auto before=s;const auto previous=geometry;FilmMatState mid;MovingVoronoiCoupledStats r;
        try{r=Curved_Edge_Coupled_Step(s,p,domain,dt,boundary,&mid);}catch(const std::exception& e){++failures;std::ofstream failure(out/"failure.json");failure<<"{\"step\":"<<step<<",\"last_accepted_step\":"<<completed<<",\"message\":"<<std::quoted(std::string(e.what()))<<"}\n";std::cerr<<"step "<<step<<" after "<<completed<<" accepted: "<<e.what()<<'\n';break;}
        if(!r.mechanics.newton.converged||r.mechanics.newton.hit_total_step||!r.drainage.converged){++failures;
            // OBSERVER ONLY: accepted histories stay untouched. Save valid
            // mechanical input and rejected drainage candidate separately.
            if(r.mechanics.newton.converged&&!r.mechanics.newton.hit_total_step){
                fields.flush();faces.flush();std::string sh,fh;{std::ifstream in(out/"states.csv");std::getline(in,sh);}{std::ifstream in(out/"faces.csv");std::getline(in,fh);}
                std::ofstream rejected(out/"failed_states.csv"),rejectedFaces(out/"failed_faces.csv");rejected<<std::setprecision(17)<<sh<<'\n';rejectedFaces<<std::setprecision(17)<<fh<<'\n';
                fields.swap(rejected);faces.swap(rejectedFaces);
                try{const auto g=Moving_Voronoi_Geometry(mid,domain,p);dump(step,0,before,previous,nullptr);dump(step,1,mid,g,nullptr);dump(step,2,s,g,&r.flux);}catch(const std::exception& e){std::cerr<<"rejected snapshot: "<<e.what()<<'\n';}
                fields.swap(rejected);faces.swap(rejectedFaces);
            }
            std::ofstream failure(out/"failure.json");failure<<std::setprecision(17)<<"{\"step\":"<<step<<",\"last_accepted_step\":"<<completed<<",\"kind\":\"nonconverged\",\"mechanics_ok\":"<<(r.mechanics.newton.converged?"true":"false")<<",\"drainage_iterations\":"<<r.drainage.iters<<",\"gradient_start\":"<<(r.drainage.grad_inf_history.empty()?0:r.drainage.grad_inf_history.front())<<",\"gradient_last\":"<<(r.drainage.grad_inf_history.empty()?0:r.drainage.grad_inf_history.back())<<"}\n";
            std::cerr<<"incomplete step "<<step<<" mechanics "<<r.mechanics.newton.converged<<" drainage "<<r.drainage.converged<<'\n';break;}
        if(observe){const auto& m=r.mechanics;real iloss=0;for(int i=0;i<s.Size();++i)iloss+=.5*p.rho*before.V[i]*(mid.v[i]-before.v[i]).squaredNorm();
            mech<<step<<','<<step*dtSI<<','<<m.newton.iters<<','<<m.newton.total_cg_iters<<','<<m.newton.neg_curvature_hits<<','<<m.newton.grad_inf_history.size()<<','<<m.initial_gradient_inf*sc.S*sc.L<<','<<m.final_gradient_inf*sc.S*sc.L<<','<<std::max(real(1e-10),m.initial_gradient_inf*1e-8)*sc.S*sc.L<<','<<m.kinetic_before*eu<<','<<m.kinetic_after*eu<<','<<m.rayleigh_action*eu<<','<<iloss*eu<<','<<m.incremental_change*eu<<','<<(m.newton.step_inf_history.empty()?0:m.newton.step_inf_history.back()*sc.L)<<','<<r.momentum_mixing_loss*eu<<'\n';
            if(step==1||step%saveStride==0||step==steps){const auto o=Observe_Curved_Mechanics(before,mid,p,domain,boundary,dt,Cut/sc.L);
                for(int z=0;z<2;++z){const auto& w=o.regions[z];regional<<step<<','<<z<<','<<w.dE*eu<<','<<w.dK*eu<<','<<w.ownershipE*eu<<','<<w.ownershipK*eu<<','<<w.inertial_loss*eu<<','<<w.viscous_dissipation*eu<<','<<w.inertial_work*eu<<','<<w.viscous_node_work*eu<<','<<w.capillary_node_work*eu<<','<<w.capillary_local_work*eu<<','<<w.capillary_remainder*eu<<','<<w.capillary_port*eu<<','<<w.viscous_port*eu<<','<<w.residual*eu<<','<<(z?o.viscous_recovered_port:-o.viscous_recovered_port)*eu<<','<<o.cut_length*sc.L<<','<<o.reconstructed_force*sc.S*sc.L<<','<<o.viscous_action_difference*eu<<','<<o.coordinate_roundtrip*sc.L<<'\n';}
                const auto ga=Moving_Voronoi_Geometry(before,domain,p),gb=Moving_Voronoi_Geometry(mid,domain,p);
                for(int i=0;i<s.Size();++i)nodes<<step<<','<<i<<','<<before.fixed[i]<<','<<boundary.Bath(i)<<','<<before.x[i].x()*sc.L<<','<<before.x[i].y()*sc.L<<','<<mid.x[i].x()*sc.L<<','<<mid.x[i].y()*sc.L<<','<<ga.graph.a[i]*sc.L*sc.L<<','<<gb.graph.a[i]*sc.L*sc.L<<','<<before.V[i]*vu<<','<<before.S[i]*su<<','<<o.w0[i]<<','<<o.w1[i]<<','<<before.v[i].x()*sc.L/sc.T<<','<<before.v[i].y()*sc.L/sc.T<<','<<mid.v[i].x()*sc.L/sc.T<<','<<mid.v[i].y()*sc.L/sc.T<<','<<o.f0.node_energy[i]*eu<<','<<o.f1.node_energy[i]*eu<<','<<o.inertial[i].x()*sc.S*sc.L<<','<<o.inertial[i].y()*sc.S*sc.L<<','<<o.viscous[i].x()*sc.S*sc.L<<','<<o.viscous[i].y()*sc.S*sc.L<<','<<o.capillary[i].x()*sc.S*sc.L<<','<<o.capillary[i].y()*sc.S*sc.L<<','<<o.regional[i].x()*sc.S*sc.L<<','<<o.regional[i].y()*sc.S*sc.L<<'\n';}}
        geometry=Moving_Voronoi_Geometry(s,domain,p);const auto wa=Moving_Interface_Weights(previous,Cut/sc.L),wb=Moving_Interface_Weights(geometry,Cut/sc.L);
        real va=0,vb=0,sa=0,sb=0,mv=0,ms=0,rel=0,vt=0,st=0;for(int i=0;i<s.Size();++i){va+=wa[i]*before.V[i];vb+=wb[i]*s.V[i];sa+=wa[i]*before.S[i];sb+=wb[i]*s.S[i];mv+=(wb[i]-wa[i])*before.V[i];ms+=(wb[i]-wa[i])*before.S[i];vt+=s.V[i];st+=s.S[i];maxmove=std::max(maxmove,(s.x[i]-initial.x[i]).norm()*sc.L);}
        for(int q=0;q<geometry.graph.pairs.Size();++q)rel+=(wb[geometry.graph.pairs.pj[q]]-wb[geometry.graph.pairs.pi[q]])*r.flux[q];
        visc+=r.mechanics.rayleigh_action;drain+=r.half_drain_action;mix+=r.momentum_mixing_loss;work+=r.reservoir_work;
        const real scale=std::max(real(1e-12),std::abs(e0)+std::abs(work)+k0),budget=(r.energy_after-e0+r.kinetic_after-k0+visc+drain+mix)/scale;
        worstE=std::max(worstE,budget);worstV=std::max(worstV,std::abs(vt-v0)/v0);worstS=std::max(worstS,std::abs(st-s0)/s0);worstRV=std::max(worstRV,std::abs(vb-va-mv-rel)/v0);worstRS=std::max(worstRS,std::abs(sb-sa-ms)/s0);
        if(step==1||step%saveStride==0||step==steps){dump(step,0,before,previous,nullptr);dump(step,1,mid,geometry,nullptr);dump(step,2,s,geometry,&r.flux);}completed=step;
        hist<<step<<','<<step*dtSI<<','<<r.energy_after*eu<<','<<r.kinetic_after*eu<<','<<visc*eu<<','<<drain*eu<<','<<mix*eu<<','<<work*eu<<','<<budget<<','<<std::abs(vt-v0)/v0<<','<<std::abs(st-s0)/s0<<','<<mv*vu<<','<<rel*vu<<','<<std::abs(vb-va-mv-rel)/v0<<','<<ms*su<<','<<std::abs(sb-sa-ms)/s0<<','<<maxmove<<'\n';
    }
    std::ofstream summary(out/"summary.json");summary<<std::setprecision(17)<<"{\"refine_start_m\":"<<(local?refineStart:0)<<",\"save_stride\":"<<saveStride<<",\"signed_bath_counter\":true,\"c_bulk_mM\":"<<M.c_bulk<<",\"preparation\":"<<std::quoted(preparation)<<",\"N\":"<<s.Size()<<",\"interior\":"<<interior<<",\"steps\":"<<completed<<",\"requested_steps\":"<<steps<<",\"failures\":"<<failures<<",\"local\":"<<(local?"true":"false")<<",\"dx_m\":"<<dx<<",\"dt_s\":"<<dtSI<<",\"T_s\":"<<completed*dtSI<<",\"radius_m\":"<<radius<<",\"far_height_m\":"<<farH<<",\"bath_pressure_Pa\":"<<mu<<",\"amplitude\":"<<amp<<",\"Gamma_inf_mol_m2\":"<<M.Gamma_inf<<",\"theta0\":"<<M.theta0<<",\"clean_tension_N_m\":"<<M.gamma0<<",\"RT_J_mol\":"<<M.gibbs_n*si::Rgas*M.Temp<<",\"rho_kg_m3\":"<<M.rho<<",\"viscosity_Pa_s\":"<<M.mu<<",\"initial_energy_J\":"<<e0*eu<<",\"V0_m3\":"<<v0*vu<<",\"S0_mol\":"<<s0*su<<",\"max_budget_excess\":"<<worstE<<",\"max_V_res\":"<<worstV<<",\"max_S_res\":"<<worstS<<",\"max_regionV_res\":"<<worstRV<<",\"max_regionS_res\":"<<worstRS<<",\"max_displacement_m\":"<<maxmove<<",\"wall_s\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
    std::cout<<"N="<<s.Size()<<" accepted="<<completed<<'/'<<steps<<" failures="<<failures<<" energy excess="<<worstE<<'\n';return failures?1:0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
