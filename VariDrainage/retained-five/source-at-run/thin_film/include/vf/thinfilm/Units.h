//////////////////////////////////////////////////////////////////////////
// Physical units for the thin-film model.
//
// The model is the lubrication energy + dissipation of a free film (ThinFilm.h, ThinFilmMechanics.h): energy =
// surfactant free energy of both surfaces (insoluble, Langmuir) + disjoining energy P(h) + squeeze energy
// (sigma_f / 8) |grad h|^2 + weight; dissipation = Poiseuille drainage between the surfaces (h^3 / 12 mu) + in-plane
// (Trouton) viscosity 4 mu h of the surfaces' plug flow (+ inertia rho h, surface diffusion). Nothing else: no
// solution or empirical law of film drainage enters; those are diagnostics to compare against.
//
// The code's variables are dimensionless: thickness h / H, lateral position x / L, time t / T, energy / (S L^2),
// surfactant Gamma / Gamma_inf. Thickness and lateral lengths have their own units (a film is ~1e-3 to 1e-6 as thick
// as it is wide), so the code coefficients carry H / L where the physics does: the squeeze energy (H / L)^2, and
// the in-plane viscosity against the drainage's: mu_s,code / mu,code = (H / L)^2.
//
// FilmMaterial holds SI inputs; each says whether it is a measured constant of the liquid / surfactant / salt or a
// modelling choice.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <cmath>
#include "vf/thinfilm/ThinFilm.h"
#include "vf/thinfilm/ThinFilmMechanics.h"

namespace vf::thinfilm {

namespace si {
constexpr real kB = 1.380649e-23;        // J/K
constexpr real NA = 6.02214076e23;       // 1/mol
constexpr real Rgas = kB * NA;           // J/(mol K)
constexpr real e = 1.602176634e-19;      // C
constexpr real eps0 = 8.8541878128e-12;  // F/m
} // namespace si

struct FilmMaterial {
	// liquid (measured: water at 25 C)
	real mu = 0.89e-3;      // Pa s
	real rho = 997.;        // kg/m^3
	real g = 9.81;          // m/s^2
	real gamma0 = 0.072;    // N/m, clean surface tension
	real Temp = 298.15;     // K
	real eps_r = 78.4;      // relative permittivity
	// surfactant (insoluble Langmuir monolayer: Pi_s = n R T Gamma_inf ln(1 / (1 - theta)) per surface).
	// Gamma_inf and n: measured (Gibbs isotherm of the surfactant; n = 1 nonionic or in excess salt, 2 ionic without).
	// theta0: the reference coverage, from the bulk concentration through the isotherm. Insolubility is a modelling
	// choice: a soluble surfactant exchanging with the bulk has a smaller effective (Marangoni) elasticity.
	real Gamma_inf = 3.2e-6; // mol/m^2
	real gibbs_n = 1.;
	real theta0 = 0.9;
	real Ds = 0.;            // m^2/s, surface diffusivity
	// disjoining (DLVO): van der Waals with Hamaker constant A (measured / Lifshitz: ~3.7e-20 J for air-water-air),
	// electrostatic double layer for 1:1 salt of concentration c_salt (Debye length from it) and surface potential
	// psi0 (measured: electrokinetics or fitted to film-thickness isotherms), weak-overlap form. The short-range
	// repulsion c8 / h^8 that makes the Newton black film is phenomenological: its one parameter places the van der
	// Waals + repulsion minimum at h_nbf (measured NBF thickness, ~4-5 nm).
	real hamaker = 3.7e-20;  // J
	real c_salt = 1.;        // mol/m^3 (1 mM)
	real psi0 = 0.06;        // V (magnitude)
	real h_nbf = 4.5e-9;     // m
	// soluble surfactant (Solution.h): off unless `soluble`. Then the surfactant is the ionic one of `sol` (its Gamma_inf
	// must equal Gamma_inf above, gibbs_n = 1: the Gibbs factor is in the bulk's activity product), dissolved at the
	// total concentration c_bulk (mol/m^3) in the liquid the film is drawn from (the menisci's), and theta0 the
	// coverage in equilibrium with it.
	bool soluble = false;
	IonicSurfactantSolution sol;
	real c_bulk = 0.;
	// Modelling choice for the isolated symmetric surface-exchange diagnostic only (SurfaceExchange.h).
	// Effective reversible pair-activity desorption rate, NOT a measured SDS kinetic coefficient.
	// Zero keeps exchange disabled. The existing film stepper does not use this parameter yet.
	real surface_desorption_rate_per_s = 0.;

	real Pi_s(real theta) const { return gibbs_n * si::Rgas * Temp * Gamma_inf * -std::log1p(-theta); }
	real Gamma_ref() const { return theta0 * Gamma_inf; }
	real gamma_ref() const { return gamma0 - Pi_s(theta0); }                        // surface tension at theta0
	real Gibbs_E() const { return 2. * gibbs_n * si::Rgas * Temp * Gamma_inf * theta0 / (1. - theta0); } // film (both surfaces)
	real debye() const { return std::sqrt(eps_r * si::eps0 * si::kB * Temp / (2. * c_salt * si::NA * si::e * si::e)); }
	// disjoining energy per area (J/m^2): C8 / h^8 - A / (12 pi h^2) + B exp(-h / debye)
	real C8() const { return hamaker * std::pow(h_nbf, 6) / (48. * kPi); }
	real B() const
	{
		const real t = std::tanh(si::e * psi0 / (4. * si::kB * Temp));
		return 64. * c_salt * si::NA * si::kB * Temp * t * t * debye();
	}
	real P(real h) const { return C8() / std::pow(h, 8) - hamaker / (12. * kPi * h * h) + B() * std::exp(-h / debye()); }
	real Pi(real h) const { return 8. * C8() / std::pow(h, 9) - hamaker / (6. * kPi * h * h * h) + B() / debye() * std::exp(-h / debye()); }
};

// base units of the code's variables: thickness H (m), lateral L (m), time T (s), tension S (N/m)
struct FilmScales {
	real H = 1e-8, L = 1e-5, T = 1e-3, S = 1e-3;
};

// the code's coefficients for a material in given scales
struct CodeParams {
	ThinFilmParams tf;   // sigma (squeeze), mu (drainage), disjoining, rho_g
	real mu_s = 0.;      // in-plane viscosity (FilmMechParams::mu_s)
	real rho = 0.;       // inertia (FilmMechParams::rho)
	Surfactant surf;     // RT, Gamma_inf = 1 (Gamma in units of Gamma_inf)
	real G0 = 0.;        // reference concentration theta0
	real kT = 0.;        // thermal energy (thermal noise)
	real Ds = 0.;        // surface diffusivity
	Solubility sol;      // soluble surfactant (FilmMechParams::sol; off for an insoluble one)
	// Reference surface pressure Pi_s(theta0) per surface: the energy is taken grand-canonical with it
	// (FilmMechParams::bath_pressure), E + 2 Pi_ref sum_i a_i. For a closed or periodic film the total area is fixed and
	// the term is a constant in the continuum (as the clean tension gamma0, already left out); for an open film it is the
	// bath (meniscus) at the reference concentration. Its point is the particles: the kernel areas a_i are not exactly
	// conserved, and a film carrying its full pressure 2 Pi_s ~ 0.04 N/m on them is a particle crystal (or pairs: an
	// indefinite step Hessian) whose shear stiffness is 10-70x the physical one at ms steps; with the reference
	// subtracted the particles carry only the pressure differences (thin_film/README.md, "In-plane preconditioner").
	real Pi_ref = 0.;
	FilmScales sc;
	// conversions
	real h_m(real h) const { return h * sc.H; }
	real x_m(real x) const { return x * sc.L; }
	real t_s(real t) const { return t * sc.T; }
	real p_Pa(real p) const { return p * sc.S / sc.H; }           // liquid pressure dE/dV
	real p_code(real p_Pa) const { return p_Pa * sc.H / sc.S; }
	void Apply(FilmMechParams& m) const { m.tf = tf; m.mu_s = mu_s; m.rho = rho; m.surf = surf; m.bath_pressure = Pi_ref; m.sol = sol; }
};

inline CodeParams To_Code(const FilmMaterial& M, const FilmScales& s)
{
	CodeParams c;
	c.sc = s;
	const real H = s.H, L = s.L, T = s.T, S = s.S, eps2 = (H / L) * (H / L);
	c.tf.sigma = 2. * M.gamma_ref() * eps2 / S;           // squeeze energy (sigma_f / 8) |grad h|^2, sigma_f = 2 gamma
	c.tf.mu = M.mu * L * L / (T * S * H);                 // drainage mobility h^3 / (12 mu)
	c.mu_s = M.mu * H / (T * S);                          // in-plane viscosity 4 mu h (= mu * eps2 / ... of the drainage's)
	c.tf.rho_g = M.rho * M.g * H * L / S;                 // weight rho g h y
	c.rho = M.rho * H * L * L / (T * T * S);              // inertia rho h
	c.tf.disj = Disjoining();
	c.tf.disj.eps = 0.;
	c.tf.disj.c8 = M.C8() / (S * std::pow(H, 8));
	c.tf.disj.c2 = M.hamaker / (12. * kPi * S * H * H);
	c.tf.disj.beta = M.B() / S;
	c.tf.disj.lambda = M.debye() / H;
	c.surf.RT = M.gibbs_n * si::Rgas * M.Temp * M.Gamma_inf / S;
	c.surf.Gamma_inf = 1.;
	c.G0 = M.theta0;
	c.Pi_ref = c.surf.Pressure(c.G0); // = M.Pi_s(theta0) / S
	c.kT = si::kB * M.Temp / (S * L * L);
	c.Ds = M.Ds * T / (L * L);
	if (M.soluble) {
		c.sol.sol = M.sol;
		c.sol.c_unit = M.Gamma_inf / H;    // amounts in Gamma_inf L^2, volumes in H L^2
		c.sol.RT = c.surf.RT;
		c.sol.c_bath = M.c_bulk / c.sol.c_unit;
	}
	return c;
}

// The solution of Monier et al., Phys. Rev. Fluids 9, 124001 (2024): SDS 5.6 g/L (19.4 mM, ~2.4 cmc) in water-glycerol
// mixtures of viscosity 1-20 mPa s (their typical case 1.5 mPa s). Liquid: the viscosity mu as given and the density of
// the mixture with it (an input: ~1045 kg/m^3 at ~20 wt% glycerol, 1.5 mPa s); clean tension 72 mN/m, the water's with
// which the isotherm was fitted (glycerol's effect on the isotherm and the cmc is not modelled). Surfactant: Solution.h's
// SDS, theta0 in equilibrium with c_bulk; the double layer's ionic strength that of its free ions (micelles left out).
inline FilmMaterial SDS_Solution(real c_bulk = 19.4, real mu = 1.5e-3, real rho = 1045.)
{
	FilmMaterial M;
	M.mu = mu; M.rho = rho;
	M.soluble = true; M.c_bulk = c_bulk;
	M.Gamma_inf = M.sol.Gamma_inf; M.gibbs_n = 1.;
	const IonicSurfactantSolution::Species sp = M.sol.Solve(c_bulk);
	M.theta0 = IonicSurfactantSolution::Theta(sp.mu);
	M.c_salt = 0.5 * (sp.c1 + sp.cNa);
	return M;
}

} // namespace vf::thinfilm
