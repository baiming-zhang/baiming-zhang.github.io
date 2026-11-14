// Optional diagnostics for the inherited moving film. No evolution or fitted closure.
// Material labels are used only across one mechanics substep, never across remeshing.
// Vsol is the effective surfactant solution-column volume, not an extra liquid inventory.
#pragma once
#include <array>
#include <string>
#include <cstdio>
#include <stdexcept>
#include "vf/thinfilm/ThinFilmMechanics.h"

namespace vf::thinfilm {
struct InterfaceAmounts {
	real V=0., S2=0., B=0., Vsol=0.; int count=0;
	real Q() const { return S2+B; }
};
using InterfaceZones=std::array<InterfaceAmounts,3>;
inline InterfaceZones Interface_Measure(const FilmMatState& st,const std::vector<int>& zones)
{
	InterfaceZones out{};
	for(int i=0;i<st.Size();++i) {
		auto& a=out[zones[i]];a.V+=st.V[i];a.S2+=2.*st.S[i];
		if((int)st.B.size()==st.Size())a.B+=st.B[i];
		a.Vsol+=st.Vsol(i);++a.count;
	}
	return out;
}
struct InterfaceLedger {
	struct Record { std::string kind,stage;int zone=-1,other=-1;InterfaceAmounts amount;
		real predictV=0.,predictQ=0.,energy=0.,work=0.; int inserted=0,merged=0; };
	std::vector<Record> records;
	std::vector<int> material_z;
	std::vector<real> material_V,material_S,material_B;
	void Clear(){records.clear();material_z.clear();material_V.clear();material_S.clear();material_B.clear();}
	void Capture(const char* stage,const FilmMatState& st,const std::vector<int>& z) {
		const auto sums=Interface_Measure(st,z);
		for(int i=0;i<3;++i){Record r;r.kind="state";r.stage=stage;r.zone=i;r.amount=sums[i];records.push_back(r);}
	}
	void Predict(const char* stage,const std::array<real,3>& v,const std::array<real,3>& q) {
		for(int i=0;i<3;++i){Record r;r.kind="prediction";r.stage=stage;r.zone=i;r.predictV=v[i];r.predictQ=q[i];records.push_back(r);}
	}
	void Scalar(const char* stage,real v,real q,real e=0.,real w=0.) {
		Record r;r.kind="scalar";r.stage=stage;r.predictV=v;r.predictQ=q;r.energy=e;r.work=w;records.push_back(r);
	}
	void Remesh_Event(int inserted,int merged){Record r;r.kind="event";r.stage="remesh";r.inserted=inserted;r.merged=merged;records.push_back(r);}
	void Begin_Material(const FilmMatState& st,const std::vector<int>& z) {
		material_z=z;material_V=st.V;material_S=st.S;material_B=st.B;
	}
	void End_Material(const std::vector<int>& z) {
		if(z.size()!=material_V.size())throw std::invalid_argument("Material IDs cannot cross a remeshing operation");
		std::array<real,3> v{},q{};std::array<InterfaceAmounts,9> cross{};
		for(size_t i=0;i<z.size();++i)if(material_z[i]!=z[i]){
			const int a=material_z[i],b=z[i];const real V=material_V[i],S2=2.*material_S[i],B=material_B.empty()?0.:material_B[i];
			v[a]-=V;v[b]+=V;q[a]-=S2+B;q[b]+=S2+B;
			auto& c=cross[3*a+b];c.V+=V;c.S2+=S2;c.B+=B;++c.count;
		}
		Predict("mechanics",v,q);
		for(int a=0;a<3;++a)for(int b=0;b<3;++b)if(cross[3*a+b].count){
			Record r;r.kind="material_crossing";r.stage="mechanics";r.zone=a;r.other=b;r.amount=cross[3*a+b];records.push_back(r);
		}
	}
	// These are INTEGRATED liquid pair volumes, positive pi -> pj.
	void Drainage(const ThinFilmGraph& G,const std::vector<real>& flux,const std::vector<int>& z) {
		std::array<real,3> v{},q{};
		for(int k=0;k<G.pairs.Size();++k){const int a=z[G.pairs.pi[k]],b=z[G.pairs.pj[k]];v[a]-=flux[k];v[b]+=flux[k];}
		Predict("drainage",v,q);
	}
	// A reservoir's B is not evolved by Solute_Transport. Its virtual counter-flux is
	// supplied/removed by the external bath. Use the final implicit donor concentration,
	// sampled BEFORE repartition changes B and S. Do not use an old concentration.
	void Transport(const ThinFilmGraph& G,const std::vector<real>& flux,const FilmMatState& st,
		const std::vector<int>& z,const std::vector<real>& cres) {
		std::array<real,3> v{},q{},pair{},bath{};
		for(int k=0;k<G.pairs.Size();++k){
			const real f=flux[k];if(f==0.)continue;
			const int d=f>0?G.pairs.pi[k]:G.pairs.pj[k],r=f>0?G.pairs.pj[k]:G.pairs.pi[k];
			const bool db=!std::isnan(cres[d]),rb=!std::isnan(cres[r]);
			const real c=db?cres[d]:st.B[d]/st.V[d],amount=std::abs(f)*c;
			pair[z[d]]-=amount;pair[z[r]]+=amount;
			if(db)bath[z[d]]+=amount;else q[z[d]]-=amount;
			if(rb)bath[z[r]]-=amount;else q[z[r]]+=amount;
		}
		Predict("transport",v,q);
		for(int i=0;i<3;++i){
			Record r;r.kind="solute_pair_and_bath";r.stage="transport";r.zone=i;
			r.predictQ=pair[i];r.amount.B=bath[i];records.push_back(r);
		}
	}
	void Write(FILE* f,int step,real t,real volume_unit,real amount_unit,real energy_unit)const {
		for(const auto& r:records)std::fprintf(f,"%d,%.17g,%s,%s,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d\n",
			step,t,r.kind.c_str(),r.stage.c_str(),r.zone,r.other,r.amount.count,
			r.amount.V*volume_unit,r.amount.S2*amount_unit,r.amount.B*amount_unit,r.amount.Vsol*volume_unit,
			r.predictV*volume_unit,r.predictQ*amount_unit,r.energy*energy_unit,r.work*energy_unit,r.inserted,r.merged);
	}
};
} // namespace vf::thinfilm
