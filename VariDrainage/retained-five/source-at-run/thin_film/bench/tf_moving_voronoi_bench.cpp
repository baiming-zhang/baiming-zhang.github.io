// Bounded SI startup dynamics: true mechanical motion + relative drainage.
// Insoluble limiting model, fixed-curvature reservoir (not resolved meniscus).
#include <fstream>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include "vf/thinfilm/MovingVoronoiFilm.h"
#include "vf/thinfilm/Meniscus.h"
#include "vf/thinfilm/MovingInterfaceLedger.h"
using namespace vf;using namespace vf::thinfilm;
constexpr real MW=.0012,MH=.0008,Patch=.0004,Hinitial=6e-6;
int main(int argc,char** argv){try{
    if(argc!=8&&argc!=9){std::cerr<<"uniform/local fine_um dt_us T_us amplitude radius_um output [ledger_cut_um]\n";return 2;}
    const std::string mode=argv[1];const bool wall=mode=="wall_local"||mode=="wall_uniform";
    const bool local=mode=="local"||mode=="wall_local";if(!local&&mode!="uniform"&&mode!="wall_uniform")return 2;
    const MovingVoronoiDiscretization discretization{wall,wall};
    const real fine=std::stod(argv[2])*1e-6,dtSI=std::stod(argv[3])*1e-6,horizon=std::stod(argv[4])*1e-6,amplitude=std::stod(argv[5]),radius=std::stod(argv[6])*1e-6;
    if(!(fine>0&&dtSI>0&&horizon>=dtSI&&std::abs(amplitude)<.02&&radius>=0))throw std::invalid_argument("invalid moving assay input");
    const int steps=int(std::llround(horizon/dtSI));if(std::abs(steps*dtSI-horizon)>1e-10*horizon)throw std::invalid_argument("unaligned startup horizon");
    const std::filesystem::path out=argv[7];std::filesystem::create_directories(out);Set_Threads(2);
    auto M=SDS_Solution();M.soluble=false;const FilmScales sc;const auto cp=To_Code(M,sc);FilmMechParams p;cp.Apply(p);p.tf.rho_g=0.;p.tf.disj={};
    const real dt=dtSI/sc.T,Vunit=sc.H*sc.L*sc.L,Sunit=M.Gamma_inf*sc.L*sc.L,Eunit=sc.S*sc.L*sc.L,velocityUnit=sc.L/sc.T;
    const auto domain=Transfer_Rect(0,0,MW/sc.L,MH/sc.L);const real cut=(MW-Patch)/sc.L;FilmMatState s;
    auto append=[&](real x0,real x1,real dx){const int nx=int(std::lround((x1-x0)/dx)),ny=int(std::lround(MH/dx));
        for(int j=0;j<ny;++j)for(int i=0;i<nx;++i){const real hx=(x1-x0)/nx,hy=MH/ny;
            real x=x0+(i+.5)*hx,y=(j+.5)*hy;
            // Deterministic geometric sampling, not physical forcing/noise.
            x+=.035*hx*std::sin(.7*i+1.3*j);y+=.035*hy*std::cos(1.1*i+.8*j);
            s.x.emplace_back(x/sc.L,y/sc.L,0.);s.fixed.push_back(j==0||j==ny-1||(x0==0.&&i==0)||(x1==MW&&i==nx-1));}}
    ;if(local){append(0,MW-Patch,2*fine);append(MW-Patch,MW,fine);}else append(0,MW,fine);
    if(wall)s.fixed.assign(s.Size(),0); // apply no-slip on the actual outer faces
    s.v.assign(s.Size(),Vec3::Zero());auto geometry=Moving_Voronoi_Geometry(s,domain,p);
    for(int i=0;i<s.Size();++i){const real x=s.x[i].x()*sc.L,y=s.x[i].y()*sc.L;s.V.push_back(Hinitial/sc.H*geometry.graph.a[i]);
        s.S.push_back(M.theta0*(1.+amplitude*std::cos(kPi*x/MW)*std::cos(2*kPi*y/MH))*geometry.graph.a[i]);}
    std::vector<real> bath(s.Size(),std::numeric_limits<real>::quiet_NaN());real mu=0.,arc=0.;MeniscusParams men;
    if(radius>0){men.r_ref=radius;mu=-M.gamma_ref()/radius;arc=Meniscus_Section::Width(men.fibre_radius,radius);
        for(int i=0;i<s.Size();++i)if(s.x[i].x()*sc.L>MW-.65*fine){bath[i]=cp.p_code(mu);s.fixed[i]=1;}}
    const std::vector<real>* res=radius>0?&bath:nullptr;
    auto region=[&](int i){return s.x[i].x()>=cut;};
    auto totals=[&](){std::array<real,5> z{};for(int i=0;i<s.Size();++i){z[0]+=s.V[i];z[1]+=s.S[i];if(region(i)){z[2]+=s.V[i];z[3]+=s.S[i];}if(res&&!std::isnan(bath[i]))z[4]+=bath[i]*s.V[i];}return z;};
    const auto initial=s;const auto initialTotals=totals();const real initialE=Differentiate_Voronoi_Energy(s,p,domain,res).energy,initialFilmE=initialE-initialTotals[4],initialK=Moving_Voronoi_Kinetic(s,p.rho);
    // Opt-in observer; no change to old modes, tolerances, state or fluxes.
    std::ofstream ledger,stages,faces;const real ledgerCut=argc==9?std::stod(argv[8])*1e-6/sc.L:cut;
    if(argc==9){ledger.open(out/"interface.csv");stages.open(out/"stages.csv");faces.open(out/"faces.csv");
        ledger<<std::setprecision(17)<<"step,stage,region,dV,materialV,relativeV,resV,dS,materialS,resS,dEcap,liquid,area,shape,metric,ownership,resCap,dEsurf,surfArea,surfOwnership,resSurf\n";
        stages<<std::setprecision(17)<<"step,stage,id,x,y,A,V,S,w\n";
        faces<<std::setprecision(17)<<"step,stage,i,j,c,F\n";}
    auto observe=[&](int step,const FilmMatState& before,const VoronoiDrainageGeometry& previous,const MovingVoronoiCoupledStats& result){
        if(argc!=9)return;auto mid=s;mid.V=before.V;const auto current=Moving_Voronoi_Geometry(s,domain,p);
        const auto w0=Moving_Interface_Weights(previous,ledgerCut),w1=Moving_Interface_Weights(current,ledgerCut);
        auto record=[&](int stage,const FilmMatState& a,const FilmMatState& b,const VoronoiDrainageGeometry& ga,const VoronoiDrainageGeometry& gb,const std::vector<real>& wa,const std::vector<real>& wb,const std::vector<real>* flux){
            for(int region=0;region<2;++region){auto u=wa,v=wb;if(region==0){for(auto& x:u)x=1-x;for(auto& x:v)x=1-x;}
                const auto r=Moving_Interface_Ledger(a,b,ga,gb,u,v,p,flux,res);
                ledger<<step<<','<<stage<<','<<region<<','<<r.volume_change*Vunit<<','<<r.material_volume*Vunit<<','<<r.relative_volume*Vunit<<','<<r.volume_residual*Vunit<<','<<r.surface_change*Sunit<<','<<r.material_surface*Sunit<<','<<r.surface_residual*Sunit<<','<<r.capillary_change*Eunit<<','<<r.liquid_work*Eunit<<','<<r.area_work*Eunit<<','<<r.shape_work*Eunit<<','<<r.metric_work*Eunit<<','<<r.ownership_work*Eunit<<','<<r.capillary_residual*Eunit<<','<<r.surface_energy_change*Eunit<<','<<r.surface_area_work*Eunit<<','<<r.surface_ownership_work*Eunit<<','<<r.surface_energy_residual*Eunit<<'\n';}}
        ;record(1,before,mid,previous,current,w0,w1,nullptr);record(2,mid,s,current,current,w1,w1,&result.flux);
        auto dump=[&](int stage,const FilmMatState& t,const VoronoiDrainageGeometry& g,const std::vector<real>& w){
            for(int i=0;i<t.Size();++i)stages<<step<<','<<stage<<','<<i<<','<<t.x[i].x()<<','<<t.x[i].y()<<','<<g.graph.a[i]<<','<<t.V[i]<<','<<t.S[i]<<','<<w[i]<<'\n';
            for(int q=0;q<g.graph.pairs.Size();++q){const int i=g.graph.pairs.pi[q],j=g.graph.pairs.pj[q];const bool isBath=res&&(!std::isnan(bath[i])||!std::isnan(bath[j]));
                faces<<step<<','<<stage<<','<<i<<','<<j<<','<<(isBath?0.:g.graph.c[q])<<','<<(stage==2?result.flux[q]:0.)<<'\n';}};
        dump(0,before,previous,w0);dump(1,mid,current,w1);dump(2,s,current,w1);
    };
    std::ofstream frames(out/"frames.csv"),history(out/"history.csv");frames<<std::setprecision(17)<<"step,t_s,id,x_m,y_m,area_m2,V_m3,S_per_face_mol,h_m,theta,ux_m_s,uy_m_s,zone,fixed,bath\n";
    history<<std::setprecision(17)<<"step,t_s,E_total_J,E_film_J,K_J,boundary_work_J,viscous_action_J,half_drain_action_J,momentum_mix_loss_J,energy_budget,V_balance,S_balance,region_V_ledger,region_S_ledger,node_V_residual,momentum_drain_residual,force_threshold_ratio,mechanics_ok,drainage_ok,max_displacement_m,min_h_m,min_theta,max_theta,material_crossings\n";
    auto snapshot=[&](int step,const VoronoiDrainageGeometry& g){for(int i=0;i<s.Size();++i)frames<<step<<','<<step*dtSI<<','<<i<<','<<s.x[i].x()*sc.L<<','<<s.x[i].y()*sc.L<<','<<g.graph.a[i]*sc.L*sc.L<<','<<s.V[i]*Vunit<<','<<s.S[i]*Sunit<<','<<s.V[i]/g.graph.a[i]*sc.H<<','<<s.S[i]/g.graph.a[i]<<','<<s.v[i].x()*velocityUnit<<','<<s.v[i].y()*velocityUnit<<','<<region(i)<<','<<int(s.fixed[i])<<','<<(res&&!std::isnan(bath[i]))<<'\n';};
    snapshot(0,geometry);const auto start=std::chrono::steady_clock::now();
    real work=0.,visc=0.,drain=0.,mix=0.,materialV=0.,materialS=0.,relativeV=0.,worstE=0.,worstV=0.,worstS=0.,worstRV=0.,worstRS=0.,worstNode=0.,worstMomentum=0.,worstForce=0.,maxMove=0.,minH=1e100,minTheta=1.,maxTheta=0.;int failures=0,completed=0,crossings=0,totalCG=0,totalNewton=0;
    for(int step=0;step<steps;++step){const auto before=s;std::vector<int> oldRegion(s.Size());for(int i=0;i<s.Size();++i)oldRegion[i]=region(i);
        const auto r=Moving_Voronoi_Coupled_Step(s,p,domain,dt,res,discretization);totalCG+=r.mechanics.newton.total_cg_iters+r.drainage.total_cg_iters;totalNewton+=r.mechanics.newton.iters+r.drainage.iters;
        if(!r.mechanics.newton.converged||r.mechanics.newton.hit_total_step||!r.drainage.converged){++failures;std::cerr<<"incomplete coupled step "<<step+1<<" mechanics "<<r.mechanics.newton.converged<<" limit "<<r.mechanics.newton.hit_total_step<<" drainage "<<r.drainage.converged<<'\n';break;}
        observe(step+1,before,geometry,r);
        geometry=Moving_Voronoi_Geometry(s,domain,p);const auto& G=geometry.graph;std::vector<int> labels(s.Size());std::vector<real> dummy(s.Size(),0.),predicted=before.V;
        for(int i=0;i<s.Size();++i){labels[i]=region(i);if(labels[i]!=oldRegion[i]){const real sign=labels[i]?1.:-1.;materialV+=sign*before.V[i];materialS+=sign*before.S[i];++crossings;}}
        const auto fluid=Voronoi_Interface_Transfer(G,labels,r.flux,s.V,dummy);relativeV+=fluid.liquid;
        for(int q=0;q<G.pairs.Size();++q){predicted[G.pairs.pi[q]]-=r.flux[q];predicted[G.pairs.pj[q]]+=r.flux[q];}
        real node=0.,hmin=1e100,tmin=1.,tmax=0.;for(int i=0;i<s.Size();++i){node=std::max(node,std::abs(s.V[i]-predicted[i])/initialTotals[0]);hmin=std::min(hmin,s.V[i]/G.a[i]*sc.H);const real t=s.S[i]/G.a[i];tmin=std::min(tmin,t);tmax=std::max(tmax,t);maxMove=std::max(maxMove,(s.x[i]-initial.x[i]).norm()*sc.L);}
        work+=r.reservoir_work;visc+=r.mechanics.rayleigh_action;drain+=r.half_drain_action;mix+=r.momentum_mixing_loss;
        const auto z=totals();const real filmE=r.energy_after-z[4],scale=std::max(real(1e-12),std::abs(initialFilmE)+std::abs(work)+initialK);
        const real budget=(filmE-initialFilmE+r.kinetic_after-initialK-work+visc+drain+mix)/scale;
        const real vb=std::abs(z[0]-initialTotals[0])/initialTotals[0],sb=std::abs(z[1]-initialTotals[1])/initialTotals[1],rv=std::abs(z[2]-initialTotals[2]-materialV-relativeV)/initialTotals[0],rs=std::abs(z[3]-initialTotals[3]-materialS)/initialTotals[1];
        const real forceRatio=r.mechanics.final_gradient_inf/std::max(real(1e-9),r.mechanics.initial_gradient_inf*1e-6);
        const real momentumResidual=(r.momentum_after_drain-r.momentum_before_drain-r.pinned_impulse).norm()/(p.rho*initialTotals[0]);
        worstE=std::max(worstE,budget);worstV=std::max(worstV,vb);worstS=std::max(worstS,sb);worstRV=std::max(worstRV,rv);worstRS=std::max(worstRS,rs);worstNode=std::max(worstNode,node);worstMomentum=std::max(worstMomentum,momentumResidual);worstForce=std::max(worstForce,forceRatio);minH=std::min(minH,hmin);minTheta=std::min(minTheta,tmin);maxTheta=std::max(maxTheta,tmax);
        completed=step+1;snapshot(completed,geometry);
        history<<completed<<','<<completed*dtSI<<','<<r.energy_after*Eunit<<','<<filmE*Eunit<<','<<r.kinetic_after*Eunit<<','<<work*Eunit<<','<<visc*Eunit<<','<<drain*Eunit<<','<<mix*Eunit<<','<<budget<<','<<vb<<','<<sb<<','<<rv<<','<<rs<<','<<node<<','<<momentumResidual<<','<<forceRatio<<','<<r.mechanics.newton.converged<<','<<r.drainage.converged<<','<<maxMove<<','<<hmin<<','<<tmin<<','<<tmax<<','<<crossings<<'\n';
    }
    std::ofstream result(out/"summary.json");result<<std::setprecision(17)<<"{\"N\":"<<s.Size()<<",\"local\":"<<(local?"true":"false")<<",\"fine_um\":"<<fine*1e6<<",\"dt_s\":"<<dtSI<<",\"requested_T_s\":"<<horizon<<",\"T_s\":"<<completed*dtSI<<",\"steps\":"<<completed<<",\"amplitude\":"<<amplitude<<",\"radius_m\":"<<radius<<",\"boundary_mu_Pa\":"<<mu<<",\"fibre_radius_m\":"<<men.fibre_radius<<",\"meniscus_arc_width_m\":"<<arc<<",\"gamma_ref_N_m\":"<<M.gamma_ref()<<",\"mu_Pa_s\":"<<M.mu<<",\"rho_kg_m3\":"<<M.rho<<",\"theta0\":"<<M.theta0
        <<",\"V0_m3\":"<<initialTotals[0]*Vunit<<",\"S0_per_face_mol\":"<<initialTotals[1]*Sunit<<",\"initial_film_excess_energy_J\":"<<initialFilmE*Eunit<<",\"boundary_work_J\":"<<work*Eunit<<",\"viscous_action_J\":"<<visc*Eunit<<",\"half_drain_action_J\":"<<drain*Eunit<<",\"momentum_mixing_loss_J\":"<<mix*Eunit<<",\"max_energy_budget_excess\":"<<worstE<<",\"max_V_balance\":"<<worstV<<",\"max_S_balance\":"<<worstS<<",\"max_region_V_ledger\":"<<worstRV<<",\"max_region_S_ledger\":"<<worstRS<<",\"max_node_V_residual\":"<<worstNode<<",\"max_momentum_drain_residual\":"<<worstMomentum<<",\"max_force_threshold_ratio\":"<<worstForce<<",\"max_displacement_m\":"<<maxMove<<",\"min_h_m\":"<<minH<<",\"min_theta\":"<<minTheta<<",\"max_theta\":"<<maxTheta<<",\"material_crossings\":"<<crossings<<",\"coupled_failures\":"<<failures<<",\"newton_iterations\":"<<totalNewton<<",\"cg_iterations\":"<<totalCG<<",\"wall_s\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
    std::cout<<"N="<<s.Size()<<" completed="<<completed<<"/"<<steps<<" failed="<<failures<<" movement(m)="<<maxMove<<" energy_excess="<<worstE<<'\n';return failures?1:0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
