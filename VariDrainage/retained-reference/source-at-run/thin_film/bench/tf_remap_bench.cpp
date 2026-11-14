// Passive remap accuracy assay. No mechanics, drainage, adsorption or time advance.
// Analytic fields are diagnostic inputs, never an empirical law in a solver.
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <filesystem>
#include "vf/thinfilm/Units.h"
using namespace vf;using namespace vf::thinfilm;
constexpr real W=.0016,H0=2e-6;
struct Totals {real V=0,S=0,B=0,A=0,E=0,mx=0,my=0;std::array<real,2> vr{},qr{};};
static int owner(const Vec3& x,real L){return x[0]<L/2?0:1;}
static real shape(real x,real y,int field){
	if(field==0)return 0.;
	if(field==1)return .25*std::sin(2*kPi*x)*std::cos(2*kPi*y);
	return -.4*std::exp(-std::pow((x-.5)/.05,2))*(1.+.1*std::cos(2*kPi*y));
}
static Totals measure(const FilmMatState& s,const ThinFilmGraph& g,const FilmMaterial& M,const CodeParams& cp,real Lx,real Ly){
	Totals t;ThinFilmEnergy e(g,cp.tf);t.E=e.Value(s.V);
	for(int i=0;i<s.Size();++i){
		const int z=owner(s.x[i],Lx);t.V+=s.V[i];t.S+=s.S[i];t.B+=s.B[i];t.A+=g.a[i];
		t.vr[z]+=s.V[i];t.qr[z]+=2*s.S[i]+s.B[i];
		t.mx+=s.V[i]*std::cos(2*kPi*s.x[i][0]/Lx);t.my+=s.V[i]*std::sin(2*kPi*s.x[i][1]/Ly);
		const real theta=s.S[i]/g.a[i];
		if(theta>0.&&theta<1.)t.E+=2*g.a[i]*(M.gamma0/cp.sc.S+cp.sol.g(theta))+s.V[i]*cp.sol.f(cp.sol.Solve(s.B[i]/s.V[i]));
		else t.E=std::numeric_limits<real>::quiet_NaN();
	}
	return t;
}
int main(int argc,char**argv){
	Set_Threads(2);const std::filesystem::path folder=argc>1?argv[1]:".";std::filesystem::create_directories(folder);
	const real boost_m_s=argc>2?std::stod(argv[2]):0.; // passive uniform translation, not a physical driving force
	std::cout<<std::setprecision(17)<<"name,nx,field,regional,cycles,N0,N1,splits,merges,volume_balance,surfactant_balance,regional_volume_change,regional_surfactant_change,h_L2,theta_L2,c_L2,depletion_L2,h_Linf,area_change,energy_change,moment_change,feasible,boost_m_s,momentum_balance,velocity_L2,kinetic_energy_change\n";
	for(int nx:{32,64,128})for(int field:{0,1,2})for(int regional:{0,1})for(int cycles:{1,4}){
		FilmScales sc;sc.L=W/nx;sc.H=1e-7;sc.T=1e-3;sc.S=.01;const FilmMaterial M=SDS_Solution();const auto cp=To_Code(M,sc);
		FilmBox box;box.Lx=nx;box.Ly=nx*std::sqrt(3.)/2.;Kernel2D kernel(4.);kernel.Set_Type(Kernel2D::WendlandC2);
		FilmMatState s;
		for(int j=0;j<nx;++j)for(int i=0;i<nx;++i){
			real x=i+.25+.5*(j&1),y=(j+.5)*std::sqrt(3.)/2.;x+=.85*box.Lx/(2*kPi)*std::sin(2*kPi*x/box.Lx);
			if(std::abs(x-.25*box.Lx)<1.05&&std::abs(y-.5*box.Ly)<2.)continue; // initial geometric hole, before assigning fields
			s.x.push_back(Vec3(x,y,0.));
		}
		auto graph=[&](const std::vector<Vec3>& x){return Make_Periodic_Patch(x,box.Lx,box.Ly,kernel);};
		auto g=graph(s.x);const real boost=boost_m_s*sc.T/sc.L;s.v.assign(s.Size(),Vec3(boost,0,0));s.fixed.assign(s.Size(),0);
		for(int i=0;i<s.Size();++i){real f=shape(s.x[i][0]/box.Lx,s.x[i][1]/box.Ly,field);
			s.V.push_back(H0/sc.H*(1.+f)*g.a[i]);s.S.push_back(.45*(1.+.2*f)*g.a[i]);s.B.push_back(cp.sol.c_bath*(1.+.1*f)*s.V[i]);}
		const int n0=s.Size();const auto before=measure(s,g,M,cp,box.Lx,box.Ly);
		RemeshParams rp;rp.absorb_fixed=false;rp.remap_fixed=false;
		if(regional)rp.remap_region=[&](const Vec3& x){return owner(x,box.Lx);};
		int splits=0,merges=0;
		for(int c=0;c<cycles;++c){const auto rs=Remesh_Film(s,g,rp,&box,nullptr,graph);splits+=rs.splits;merges+=rs.merges;g=graph(s.x);}
		const auto after=measure(s,g,M,cp,box.Lx,box.Ly);real eh=0,dh=0,et=0,dt=0,ec=0,dc=0,dep=0,linf=0;bool feasible=std::isfinite(after.E);
		const std::string name="n"+std::to_string(nx)+"_f"+std::to_string(field)+"_r"+std::to_string(regional)+"_c"+std::to_string(cycles);
		std::ofstream state;if(nx==64&&field==2&&cycles==4){state.open(folder/(name+".csv"));state<<std::setprecision(17)<<"x_m,y_m,h_m,h_exact_m,theta,theta_exact,c_mM,c_exact_mM,region\n";}
		for(int i=0;i<s.Size();++i){const real f=shape(s.x[i][0]/box.Lx,s.x[i][1]/box.Ly,field),a=g.a[i];
			const real h=s.V[i]/a*sc.H,ht=H0*(1.+f),theta=s.S[i]/a,tt=.45*(1.+.2*f),con=s.B[i]/s.V[i]*cp.sol.c_unit,ct=M.c_bulk*(1.+.1*f);
			eh+=a*(h-ht)*(h-ht);dh+=a*ht*ht;et+=a*(theta-tt)*(theta-tt);dt+=a*tt*tt;ec+=a*(con-ct)*(con-ct);dc+=a*ct*ct;dep+=a*(ht-H0)*(ht-H0);
			linf=std::max(linf,std::abs(h-ht)/H0);feasible=feasible&&h>0.&&theta>0.&&theta<1.&&con>0.;
			if(state.is_open())state<<s.x[i][0]*sc.L<<','<<s.x[i][1]*sc.L<<','<<h<<','<<ht<<','<<theta<<','<<tt<<','<<con<<','<<ct<<','<<owner(s.x[i],box.Lx)<<'\n';
		}
		real rv=0,rq=0;for(int z=0;z<2;++z){rv=std::max(rv,std::abs(after.vr[z]-before.vr[z])/before.V);rq=std::max(rq,std::abs(after.qr[z]-before.qr[z])/(2*before.S+before.B));}
		real momentum=0.,kinetic=0.,velocity_error=0.;
		for(int i=0;i<s.Size();++i){momentum+=s.V[i]*s.v[i][0];kinetic+=.5*s.V[i]*s.v[i].squaredNorm();velocity_error+=s.V[i]*(s.v[i]-Vec3(boost,0,0)).squaredNorm();}
		const real k0=.5*before.V*boost*boost;
		std::cout<<name<<','<<nx<<','<<field<<','<<regional<<','<<cycles<<','<<n0<<','<<s.Size()<<','<<splits<<','<<merges<<','
			<<std::abs(after.V-before.V)/before.V<<','<<std::abs(2*after.S+after.B-2*before.S-before.B)/(2*before.S+before.B)<<','<<rv<<','<<rq<<','
			<<std::sqrt(eh/dh)<<','<<std::sqrt(et/dt)<<','<<std::sqrt(ec/dc)<<','<<(dep>1e-30?std::sqrt(eh/dep):0.)<<','<<linf<<','
			<<after.A/before.A-1.<<','<<(after.E-before.E)/std::max(std::abs(before.E),1e-30)<<','<<std::hypot(after.mx-before.mx,after.my-before.my)/before.V<<','<<feasible<<','<<boost_m_s<<','
			<<(boost!=0.?std::abs(momentum-before.V*boost)/(before.V*std::abs(boost)):0.)<<','
			<<(boost!=0.?std::sqrt(velocity_error/(after.V*boost*boost)):0.)<<','<<(k0>0.?kinetic/k0-1.:0.)<<'\n';
	}
	return 0;
}
