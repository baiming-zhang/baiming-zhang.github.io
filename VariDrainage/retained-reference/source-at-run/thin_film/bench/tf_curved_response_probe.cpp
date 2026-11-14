// Read-only directional response of the ORIGINAL full curved-energy model.
// No native time stepping, model corrections, fitting or alternate geometry.
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include "vf/thinfilm/CurvedEdgeMotion.h"
#include "vf/thinfilm/Units.h"
using namespace vf;using namespace vf::thinfilm;
using PairSet=std::set<std::pair<int,int>>;
struct Observation {VoronoiEnergyDerivative e;std::vector<real> surface;PairSet pairs;};
static Observation Observe(const FilmMatState& s,const FilmMechParams& p,const TransferPolygon& box,const CurvedEdgeBoundary& b){
    const auto geom=Moving_Voronoi_Geometry(s,box,p);Observation r;r.e=Differentiate_Curved_Edge(s,p,box,b);r.surface=CurvedEdgeEnergy(geom.graph,p,s.S,b).Fields(s.V,false).area;
    for(int q=0;q<geom.graph.pairs.Size();++q)r.pairs.emplace(geom.graph.pairs.pi[q],geom.graph.pairs.pj[q]);return r;
}
int main(int argc,char** argv){try{
    if(argc!=4)return 2;const std::string input=argv[1];const real dxSI=std::stod(argv[2])*1e-6;const int nx=int(std::lround(320e-6/dxSI));if(nx!=16&&nx!=32)return 3;
    const std::filesystem::path output(argv[3]);auto meta=output;meta.replace_extension(".meta.json");auto nodes=output;nodes.replace_extension(".nodes.csv");if(std::filesystem::exists(output)||std::filesystem::exists(meta)||std::filesystem::exists(nodes))return 4;Set_Threads(2);
    const bool flat=input=="flat16"||input=="flat32";if(flat&&input!=(nx==16?"flat16":"flat32"))return 5;
    auto M=SDS_Solution(6);M.soluble=false;const FilmScales sc;const auto cp=To_Code(M,sc);FilmMechParams p;cp.Apply(p);p.tf.disj={};p.tf.rho_g=0;
    const real W=320e-6/sc.L,Y=160e-6/sc.L,dx=dxSI/sc.L,vu=sc.H*sc.L*sc.L,su=M.Gamma_inf*sc.L*sc.L,eu=sc.S*sc.L*sc.L;
    const auto domain=Transfer_Rect(0,0,W,Y);FilmMatState s;CurvedEdgeBoundary b;b.aspect=sc.H/sc.L;b.clean_tension=M.gamma0/sc.S;
    real observedE=0;std::vector<real> observedA;
    if(flat){
        for(int j=0;j<nx/2;++j)for(int i=0;i<nx;++i){s.x.emplace_back((i+.5)*dx,(j+.5)*dx,0);s.v.push_back(Vec3::Zero());s.fixed.push_back(0);}
        const auto G=Moving_Voronoi_Geometry(s,domain,p).graph;
        for(int i=0;i<s.Size();++i){s.V.push_back(6e-6/sc.H*G.a[i]);s.S.push_back(M.theta0*G.a[i]);b.height.push_back(0);b.potential.push_back(std::numeric_limits<real>::quiet_NaN());}
    }else{
        std::ifstream in{std::filesystem::path(input)};if(!in)return 6;std::string line,word;std::getline(in,line);std::stringstream head(line);std::map<std::string,int> columns;int id=0;
        while(std::getline(head,word,',')){if(!word.empty()&&word.back()=='\r')word.pop_back();columns[word]=id++;}
        while(std::getline(in,line)){if(line.empty())continue;std::stringstream row(line);std::vector<real> v;while(std::getline(row,word,','))v.push_back(std::stod(word));auto get=[&](const std::string& k){return v.at(columns.at(k));};if(int(get("id"))!=s.Size())throw std::runtime_error("nonsequential snapshot");const bool bath=get("bath")!=0;
            s.x.emplace_back(get("x_m")/sc.L,get("y_m")/sc.L,0);s.v.emplace_back(get("ux_m_s")*sc.T/sc.L,get("uy_m_s")*sc.T/sc.L,0);s.fixed.push_back(bath);s.V.push_back(get("V_m3")/vu);s.S.push_back(get("S_mol")/su);b.height.push_back(bath?get("h_m")/sc.H:0);b.potential.push_back(bath?cp.p_code(get("pressure_Pa")):std::numeric_limits<real>::quiet_NaN());observedA.push_back(get("A_m2"));observedE+=get("node_energy_J");}
    }
    const auto geom=Moving_Voronoi_Geometry(s,domain,p);const auto baseline=Observe(s,p,domain,b);const auto J=Moving_Voronoi_Area_Jacobian(geom);
    real areaError=0,energyError=0;if(!flat){for(int i=0;i<s.Size();++i)areaError=std::max(areaError,std::abs(geom.graph.a[i]*sc.L*sc.L/observedA[i]-1));energyError=std::abs(baseline.e.energy*eu-observedE)/std::max(real(1e-30),std::abs(observedE));if(areaError>1e-9||energyError>1e-9)throw std::runtime_error("snapshot import failed");}
    const real E0=p.surf.RT*p.surf.Gamma_inf*M.theta0/(1-M.theta0),c0=std::sqrt(2*E0*sc.S/(M.rho*6e-6));
    MovingVoronoiProblem metric(s,p,domain,1.,&b.potential,{true,true});
    std::ofstream out(output);out<<std::setprecision(17)<<"mode,epsilon,k_per_m,mass_kg,continuum_K_N_per_m,full_K_N_per_m,value_K_N_per_m,surface_K_N_per_m,projected_K_N_per_m,geometric_residual_K_N_per_m,damping_per_s,velocity_projection_m_per_s,kinetic_fraction,kinetic_available,force_projection_N,topology_plus_equal,topology_minus_equal,flat_expected_K_N_per_m,flat_J_error,base_energy_J,plus_energy_J,minus_energy_J\n";
    std::vector<int> modes=flat?std::vector<int>{1,8,nx}:nx==16?std::vector<int>{1,4,8,12,16}:std::vector<int>{1,4,8,12,16,24,32};
    for(int m:modes){
        const real k=m*kPi/W,kSI=k/sc.L;VecX q=VecX::Zero(2*s.Size());real mass=0,kinetic=0,velocityDot=0,force=0;
        for(int i=0;i<s.Size();++i)if(!b.Bath(i)){q[2*i]=std::sin(k*s.x[i].x());const real mi=M.rho*s.V[i]*vu;mass+=mi*q[2*i]*q[2*i];kinetic+=mi*s.v[i].squaredNorm()*sc.L*sc.L/(sc.T*sc.T);velocityDot+=mi*q[2*i]*s.v[i].x()*sc.L/sc.T;force+=baseline.e.gradient[i].x()*q[2*i]*eu/sc.L;}
        if(!(mass>0))throw std::runtime_error("empty mode mass");const real K0=mass*std::pow(c0*kSI,2),projection=velocityDot/mass,fraction=kinetic>0?velocityDot*velocityDot/(mass*kinetic):0;
        VecX freeQ(metric.Dim());for(size_t k0=0;k0<metric.Free().size();++k0)freeQ.segment<2>(2*k0)=q.segment<2>(2*metric.Free()[k0]);const real damping=2*metric.Rayleigh(freeQ)*sc.S*sc.T/mass;
        const VecX Jq=J*q;real projected=0,flatK=0,jError=0;std::vector<real> dh(s.Size());
        for(int i=0;i<s.Size();++i){const real theta=s.S[i]/(p.surf.Gamma_inf*baseline.surface[i]),Ei=p.surf.RT*p.surf.Gamma_inf*theta/(1-theta);projected+=2*Ei/baseline.surface[i]*Jq[i]*Jq[i];dh[i]=-(6e-6/sc.H)*Jq[i]/geom.graph.a[i];flatK+=2*E0/geom.graph.a[i]*Jq[i]*Jq[i];if(flat)jError=std::max(jError,std::abs(Jq[i]-dx*std::sin(k*dx)*std::cos(k*s.x[i].x()))/dx);}
        if(flat)for(int l=0;l<geom.graph.pairs.Size();++l){const int i=geom.graph.pairs.pi[l],j=geom.graph.pairs.pj[l];flatK+=.5*(M.gamma_ref()/sc.S)*b.aspect*b.aspect*geom.graph.c[l]*std::pow(dh[i]-dh[j],2);}
        for(real eps:{1e-3,1e-4,1e-5}){
            const real a=eps*dx;auto plus=s,minus=s;for(int i=0;i<s.Size();++i){plus.x[i].x()+=a*q[2*i];minus.x[i].x()-=a*q[2*i];}
            const auto pa=Observe(plus,p,domain,b),pb=Observe(minus,p,domain,b);real full=0,intrinsic=0;
            for(int i=0;i<s.Size();++i){full+=q[2*i]*(pa.e.gradient[i].x()-pb.e.gradient[i].x())/(2*a);const real ds=(pa.surface[i]-pb.surface[i])/(2*a),theta=s.S[i]/(p.surf.Gamma_inf*baseline.surface[i]),Ei=p.surf.RT*p.surf.Gamma_inf*theta/(1-theta);intrinsic+=2*Ei/baseline.surface[i]*ds*ds;}
            const real value=(pa.e.energy+pb.e.energy-2*baseline.e.energy)/(a*a);
            for(real finite:{full,value,intrinsic,projected,damping,projection,fraction})if(!std::isfinite(finite))throw std::runtime_error("nonfinite response");
            out<<m<<','<<eps<<','<<kSI<<','<<mass<<','<<K0<<','<<full*sc.S<<','<<value*sc.S<<','<<intrinsic*sc.S<<','<<projected*sc.S<<','<<(full-intrinsic)*sc.S<<','<<damping<<','<<projection<<','<<fraction<<','<<(kinetic>0)<<','<<force<<','<<(pa.pairs==baseline.pairs)<<','<<(pb.pairs==baseline.pairs)<<','<<(flat?flatK*sc.S:0)<<','<<jError<<','<<baseline.e.energy*eu<<','<<pa.e.energy*eu<<','<<pb.e.energy*eu<<'\n';
        }
    }
    std::ofstream nodeout(nodes);nodeout<<std::setprecision(17)<<"id,x_m,y_m,A_m2,As_m2,V_m3,S_mol,bath,ux_m_s,uy_m_s,gradient_x_N,gradient_y_N\n";
    for(int i=0;i<s.Size();++i)nodeout<<i<<','<<s.x[i].x()*sc.L<<','<<s.x[i].y()*sc.L<<','<<geom.graph.a[i]*sc.L*sc.L<<','<<baseline.surface[i]*sc.L*sc.L<<','<<s.V[i]*vu<<','<<s.S[i]*su<<','<<b.Bath(i)<<','<<s.v[i].x()*sc.L/sc.T<<','<<s.v[i].y()*sc.L/sc.T<<','<<baseline.e.gradient[i].x()*eu/sc.L<<','<<baseline.e.gradient[i].y()*eu/sc.L<<'\n';
    std::ofstream metadata(meta);metadata<<std::setprecision(17)<<"{\"N\":"<<s.Size()<<",\"flat\":"<<(flat?"true":"false")<<",\"dx_m\":"<<dxSI<<",\"area_import_error\":"<<areaError<<",\"energy_import_error\":"<<energyError<<",\"E0_N_per_m\":"<<E0*sc.S<<",\"gamma0_equilibrium_N_per_m\":"<<M.gamma_ref()<<",\"c0_m_per_s\":"<<c0<<",\"directions\":"<<modes.size()<<",\"epsilon_rows\":"<<3*modes.size()<<"}\n";
    std::cout<<"curved directional response N="<<s.Size()<<" modes="<<modes.size()<<" flat="<<flat<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
