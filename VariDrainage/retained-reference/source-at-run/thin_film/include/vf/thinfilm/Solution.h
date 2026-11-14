//////////////////////////////////////////////////////////////////////////
// Soluble ionic surfactant (SDS-like: a 1:1 salt, no added salt) in the film's liquid: monomers, counterions and
// micelles in mass-action equilibrium, and the Langmuir monolayer on the surfaces in equilibrium with them
// (thin_film/README.md, "Soluble surfactant").
//
// Thermodynamics: E. M. Blokhuis, Langmuir 41, 28521 (2025) (arXiv:2509.20110), eqs. 37-44; parameters from its
// Table 3, the fit of Elworthy & Mysels' surface tensions of pure SDS:
// - species, as mole fractions x_i = v0 c_i (v0 the solvent's molar volume): DS- monomers c1, free Na+ cNa, and the
//   surfactant in micelles cm (aggregation number m; a fraction r of the micelles' counterions bound). Balances
//   c = c1 + cm and c = cNa + r cm (all Na+ comes with the surfactant); micelles in mass-action equilibrium,
//   x_m = m (x1 xNa^r / x0)^m;
// - ideal dilute mixing: the chemical potential of a DS-Na pair is RT mu~ with mu~ = ln(K2 c1 cNa) (K2 also carries
//   the surface's standard potential), the osmotic pressure RT (c1 + cNa + cm / m) (a micelle counts once), and the
//   free energy density RT Phi(c), Phi = c mu~ - (c1 + cNa + cm / m): dPhi/dc = mu~ exactly (Gibbs-Duhem);
// - the surface: Langmuir free energy g(Gamma) = RT Gamma_inf [theta ln theta + (1 - theta) ln(1 - theta)] per area
//   (plus Gamma times the standard potential, in K2); in equilibrium with the bulk ln(theta / (1 - theta)) = mu~, i.e.
//   theta = K2 c1 cNa / (1 + K2 c1 cNa), and the surface pressure RT Gamma_inf ln(1 + K2 c1 cNa) (eq. 38). Below the
//   cmc c1 = cNa = c: the Gibbs factor 2 of an ionic surfactant without salt comes from the activity product.
// These are equilibrium properties of the solution, measured; no solution of the film's dynamics enters.
//
// Solubility (code units, used by FilmMechParams::sol): a particle carries the amount Q = 2 S + B, S per surface and
// B dissolved in its liquid V, in local equilibrium (one chemical potential): 2 a theta + V c = Q. The local
// equilibrium is a modelling choice: diffusion across a film of thickness h takes (h / 2)^2 / D ~ 1-3 ms at 2 um, and
// SDS micelles exchange monomers in ~ms near 20 mM (their slow relaxation time spans 1e-4 - 10 s with concentration,
// largest at 200 mM: Patist, Oh, Leung & Shah), against flows over seconds.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include "vf/core/Types.h"

namespace vf::thinfilm {

struct IonicSurfactantSolution {
	// SDS: Blokhuis (2025) Table 3, fit of Elworthy & Mysels' data (Tajima's: Gamma_inf 3.19e-6, 1 / K2 = 1.00 (mol/m^3)^2)
	real Gamma_inf = 3.70e-6;  // mol/m^2
	real K2 = 1. / 2.16;       // (m^3/mol)^2: 1 / K2 = 2.16e-6 (mol/L)^2
	real m = 65.;              // aggregation number
	real r = 0.75;             // bound counterions per micellar headgroup
	real x0 = 2.7e-7;          // micellization constant (mole-fraction units)
	real v0 = 1. / 55.3e3;     // m^3/mol, molar volume of water (55.3 mol/L)

	struct Species {
		real c1 = 0., cNa = 0., cm = 0.; // mol/m^3: monomers, free counterions, surfactant in micelles
		real mu = 0.;                    // mu~ = ln(K2 c1 cNa)
		real dmu = 0.;                   // dmu~/dc (m^3/mol)
		real s = 0.;                     // ln x1 (a warm start for the next solve)
	};
	// species at the total concentration c > 0 (mol/m^3); s_guess: ln x1 of a nearby solve
	Species Solve(real c, const real* s_guess = nullptr) const
	{
		Species sp;
		const real X = v0 * c, lx0 = std::log(x0), lm = std::log(m);
		// R(s) = ln m + m (s + r ln xNa - ln x0) - ln xm = 0, xm = X - x1, xNa = X - r xm: increasing in s = ln x1
		auto R = [&](real s, real& d) {
			const real x1 = std::exp(s), xm = X - x1;
			if (!(xm > 0.)) { d = 1.; return std::numeric_limits<real>::infinity(); }
			const real xNa = X - r * xm;
			d = m * (1. + r * r * x1 / xNa) + x1 / xm;
			return lm + m * (s + r * std::log(xNa) - lx0) - std::log(xm);
		};
		real hi = std::log(X), lo = hi - 1., d = 0.;
		while (R(lo, d) > 0.) { hi = lo; lo -= 1.; }
		real s = (s_guess && *s_guess > lo && *s_guess < hi) ? *s_guess : 0.5 * (lo + hi);
		for (int it = 0; it < 200; it++) {
			const real f = R(s, d);
			if (f > 0.) hi = s; else lo = s;
			real sn = s - f / d;
			if (!(sn > lo && sn < hi)) sn = 0.5 * (lo + hi);
			const bool done = std::abs(sn - s) <= 1e-15 * (1. + std::abs(s)) || hi - lo <= 1e-15 * (1. + std::abs(s));
			s = sn;
			if (done) break;
		}
		// micelles from their own equation where they are few (X - x1 would cancel)
		real x1 = std::exp(s), xm = X - x1;
		const real xm_eq = m * std::exp(m * (s + r * std::log(X - r * std::max<real>(xm, 0.)) - lx0));
		if (xm_eq < 0.5 * X) { xm = xm_eq; x1 = X - xm; }
		xm = std::max<real>(xm, 0.);
		const real xNa = X - r * xm;
		sp.s = std::log(x1);
		sp.c1 = x1 / v0; sp.cNa = xNa / v0; sp.cm = xm / v0;
		sp.mu = std::log(K2) + std::log(sp.c1) + std::log(sp.cNa);
		// implicit derivative: b = dxm/dX, a = dx1/dX = 1 - b, dmu~/dX = a / x1 + (1 - r b) / xNa (2 / X below the cmc)
		const real b = m * (1. / x1 + r / xNa) * xm / (1. + xm * (m * r * r / xNa + m / x1));
		sp.dmu = v0 * ((1. - b) / x1 + (1. - r * b) / xNa);
		return sp;
	}
	real Mu(real c) const { return Solve(c).mu; }
	real Phi(const Species& sp) const { return (sp.c1 + sp.cm) * sp.mu - (sp.c1 + sp.cNa + sp.cm / m); }
	real Phi(real c) const { return Phi(Solve(c)); }
	static real Theta(real mu) { return 1. / (1. + std::exp(-mu)); }
	static real Softplus(real mu) { return mu > 0. ? mu + std::log1p(std::exp(-mu)) : std::log1p(std::exp(mu)); } // ln(1 + e^mu)
	// equilibrium coverage and surface pressure (N/m, at temperature T) of a surface over the solution at c
	real Theta_Eq(real c) const { return Theta(Mu(c)); }
	real Pi_Eq(real c, real RT) const { return RT * Gamma_inf * Softplus(Mu(c)); }
};

// Solubility in code units: amounts in Gamma_inf L^2 (S per surface = theta a), volumes in H L^2, so concentrations
// in c_unit = Gamma_inf / H (mol/m^3), energies in S L^2; RT is Surfactant::RT = R T Gamma_inf / S (gibbs_n = 1:
// the factor 2 is in the bulk's activity product).
struct Solubility {
	IonicSurfactantSolution sol;
	real c_unit = 0.; // mol/m^3 per code concentration unit; 0: the surfactant is insoluble
	real RT = 1.;
	real c_bath = 0.; // code concentration of the menisci's liquid (their surfaces are a surfactant bath at it)
	bool On() const { return c_unit > 0.; }

	using Species = IonicSurfactantSolution::Species;
	Species Solve(real c, const real* s_guess = nullptr) const { return sol.Solve(c * c_unit, s_guess); }
	// free energy per volume of the solution at code concentration c (df/dc = RT mu~)
	real f(const Species& sp) const { return RT * sol.Phi(sp) / c_unit; }
	// surface free energy per area (Langmuir) and surface pressure at coverage theta
	real g(real th) const { return RT * (th * std::log(th) + (1. - th) * std::log1p(-th)); }
	real Pi(real th) const { return -RT * std::log1p(-th); }
	real Theta_Bath() const { return IonicSurfactantSolution::Theta(Solve(c_bath).mu); }
	real Pi_Bath() const { return RT * IonicSurfactantSolution::Softplus(Solve(c_bath).mu); }

	// Local equilibrium of the amount Q between the two surfaces (area a each) and the liquid V: 2 a theta(c) + V c = Q
	// with ln(theta / (1 - theta)) = mu~(c). guess: ln c of a nearby solve (updated).
	struct Part {
		real theta = 0., c = 0., dmu = 0.; // coverage, code concentration, dmu~/dc (code)
		Species sp;
	};
	Part Partition(real a, real V, real Q, real* guess = nullptr) const
	{
		Part p;
		real s_in = std::numeric_limits<real>::quiet_NaN();
		auto R = [&](real u, real& d) {
			const real c = std::exp(u);
			p.sp = Solve(c, std::isnan(s_in) ? nullptr : &s_in);
			s_in = p.sp.s;
			p.c = c; p.theta = IonicSurfactantSolution::Theta(p.sp.mu); p.dmu = p.sp.dmu * c_unit;
			d = c * (2. * a * p.theta * (1. - p.theta) * p.dmu + V);
			return 2. * a * p.theta + V * c - Q;
		};
		real hi = std::log(Q / V), lo = hi - 1., d = 0.;
		while (R(lo, d) > 0.) { hi = lo; lo -= 1.; }
		real u = (guess && *guess > lo && *guess < hi) ? *guess : 0.5 * (lo + hi);
		for (int it = 0; it < 200; it++) {
			const real F = R(u, d);
			if (F > 0.) hi = u; else lo = u;
			real un = u - F / d;
			if (!(un > lo && un < hi)) un = 0.5 * (lo + hi);
			const bool done = std::abs(un - u) <= 1e-14 * (1. + std::abs(u)) || hi - lo <= 1e-14 * (1. + std::abs(u));
			u = un;
			if (done) break;
		}
		R(u, d);
		if (guess) *guess = u;
		return p;
	}
	// The particle's free energy F(a) = 2 a g(theta) + V f(c) at fixed V and Q, as a function of its area:
	// dF/da = -2 Pi(theta) (envelope theorem: the partition is a minimum), d2F/da2 = 4 RT theta^2 mu~' / (2 a theta
	// (1 - theta) mu~' + V) (the insoluble limit V -> 0: 2 theta^2 g'' / a). Film elasticity a d2F/da2 = 2 E_G /
	// (1 + h / (2 dGamma/dc)), E_G = theta^2 g'' (Gibbs).
	void Energy(real a, real V, real Q, real& F, real& Fa, real& Faa, real* guess = nullptr) const
	{
		const Part p = Partition(a, V, Q, guess);
		const real th = p.theta;
		F = 2. * a * g(th) + V * f(p.sp);
		Fa = -2. * Pi(th);
		Faa = 4. * RT * th * th * p.dmu / (2. * a * th * (1. - th) * p.dmu + V);
	}
};

} // namespace vf::thinfilm
