// DRIFT-only frozen native snapshot importer for read-only numerical probes.
// Same SI/material/bath interpretation as the original native cost probe.
#pragma once
#include <filesystem>
#include <fstream>
#include <sstream>
#include <map>
#include "vf/thinfilm/CurvedEdge.h"
#include "vf/thinfilm/Units.h"
namespace vf::thinfilm {
struct FrozenCurvedProbe {
    FilmMaterial material=SDS_Solution(6);FilmScales scales;FilmMechParams params;
    FilmMatState state;CurvedEdgeBoundary boundary;TransferPolygon domain;
    real area_import_error=0,energy_import_error=0;
    explicit FrozenCurvedProbe(const std::filesystem::path& input){
        material.soluble=false;const auto cp=To_Code(material,scales);cp.Apply(params);params.tf.disj={};params.tf.rho_g=0;
        const real vu=scales.H*scales.L*scales.L,su=material.Gamma_inf*scales.L*scales.L,eu=scales.S*scales.L*scales.L;
        boundary.aspect=scales.H/scales.L;boundary.clean_tension=material.gamma0/scales.S;domain=Transfer_Rect(0,0,320e-6/scales.L,160e-6/scales.L);
        std::ifstream in(input);if(!in)throw std::runtime_error("cannot open snapshot");std::string line,word;std::getline(in,line);std::stringstream head(line);std::map<std::string,int> columns;int id=0;while(std::getline(head,word,',')){if(!word.empty()&&word.back()=='\r')word.pop_back();columns[word]=id++;}
        real observedE=0;std::vector<real> observedA;
        while(std::getline(in,line)){if(line.empty())continue;std::stringstream row(line);std::vector<real> v;while(std::getline(row,word,','))v.push_back(std::stod(word));auto get=[&](const std::string& k){return v.at(columns.at(k));};if(int(get("id"))!=state.Size())throw std::runtime_error("snapshot ids not sequential");const bool bath=get("bath")!=0;
            state.x.emplace_back(get("x_m")/scales.L,get("y_m")/scales.L,0);state.v.emplace_back(get("ux_m_s")*scales.T/scales.L,get("uy_m_s")*scales.T/scales.L,0);state.fixed.push_back(bath);state.V.push_back(get("V_m3")/vu);state.S.push_back(get("S_mol")/su);boundary.height.push_back(bath?get("h_m")/scales.H:0);boundary.potential.push_back(bath?cp.p_code(get("pressure_Pa")):std::numeric_limits<real>::quiet_NaN());observedA.push_back(get("A_m2"));observedE+=get("node_energy_J");}
        if(state.Size()<3)throw std::runtime_error("snapshot too small");const auto G=Moving_Voronoi_Geometry(state,domain,params).graph;const auto E=CurvedEdgeEnergy(G,params,state.S,boundary).Value(state.V);
        for(int i=0;i<state.Size();++i)area_import_error=std::max(area_import_error,std::abs(G.a[i]*scales.L*scales.L/observedA[i]-1));energy_import_error=std::abs(E*eu-observedE)/std::max(real(1e-30),std::abs(observedE));
        if(area_import_error>1e-9||energy_import_error>1e-9)throw std::runtime_error("snapshot native import differs");
    }
};
}
