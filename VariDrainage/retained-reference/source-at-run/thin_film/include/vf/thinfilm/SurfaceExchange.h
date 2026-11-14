// Closed, homogeneous, symmetric two-surface SI diagnostic. No spatial film coupling.
// Inherits Solution.h's paired-salt equilibrium thermodynamics, including instantaneous
// bulk micellization. The effective kinetic law is a modelling choice, not a measured
// microscopic SDS adsorption mechanism: dtheta/dt = kd [ab (1-theta) - theta],
// ab = exp(mu_bulk). Its sign is consistent with chemical affinity at any coverage.
// With B = Q - 2 A Gamma_inf theta, total surfactant is conserved algebraically.
// At fixed A,V,Q: dF/dt = 2 A RT Gamma_inf kd (ab(1-theta)-theta)
//                                  * log(theta/(ab(1-theta))) <= 0.
// Backward Euler has a unique root in the feasible interval since activity increases
// with bulk concentration. Area changes do external work; compare energies only
// during the subsequent fixed-area relaxation, never across the imposed area jump.
#pragma once
#include <stdexcept>
#include "vf/thinfilm/Units.h"

namespace vf::thinfilm {

inline real Exchange_Activity(const FilmMaterial& M, real c)
{
	if (c <= 0.) return 0.;
	const auto sp = M.sol.Solve(c);
	return M.sol.K2 * sp.c1 * sp.cNa;
}

inline real Exchange_Rate(const FilmMaterial& M, real a, real V, real Q, real theta)
{
	const real c = (Q - 2. * a * M.Gamma_inf * theta) / V;
	return M.surface_desorption_rate_per_s * (Exchange_Activity(M, c) * (1. - theta) - theta);
}

inline real Exchange_Equilibrium(const FilmMaterial& M, real a, real V, real Q)
{
	real lo = 0., hi = std::min<real>(1., Q / (2. * a * M.Gamma_inf));
	for (int i = 0; i < 90; ++i) {
		const real th = .5 * (lo + hi), c = (Q - 2. * a * M.Gamma_inf * th) / V;
		if (Exchange_Activity(M, c) * (1. - th) > th) lo = th; else hi = th;
	}
	return .5 * (lo + hi);
}

inline real Exchange_Step(const FilmMaterial& M, real a, real V, real Q, real old, real dt)
{
	const real upper = std::min<real>(1., Q / (2. * a * M.Gamma_inf));
	if (!(a > 0. && V > 0. && Q > 0. && old > 0. && old < upper && dt >= 0.
		&& M.surface_desorption_rate_per_s >= 0.)) throw std::invalid_argument("infeasible surface exchange cell");
	if (dt == 0. || M.surface_desorption_rate_per_s == 0.) return old;
	real lo = 0., hi = upper, th = old;
	for (int i = 0; i < 100; ++i) {
		const real c = (Q - 2. * a * M.Gamma_inf * th) / V;
		const auto sp = M.sol.Solve(c);
		const real ab = M.sol.K2 * sp.c1 * sp.cNa, kd = M.surface_desorption_rate_per_s;
		const real f = th - old - dt * kd * (ab * (1. - th) - th);
		if (f > 0.) hi = th; else lo = th;
		const real df = 1. + dt * kd * (1. + ab + (1. - th) * ab * sp.dmu * 2. * a * M.Gamma_inf / V);
		real next = th - f / df;
		if (!(next > lo && next < hi)) next = .5 * (lo + hi);
		if (std::abs(next - th) < 2e-15) return next;
		th = next;
	}
	throw std::runtime_error("surface exchange root did not converge");
}

inline real Exchange_Free_Energy(const FilmMaterial& M, real a, real V, real Q, real theta)
{
	const real c = (Q - 2. * a * M.Gamma_inf * theta) / V;
	return si::Rgas * M.Temp * (2. * a * M.Gamma_inf *
		(theta * std::log(theta) + (1. - theta) * std::log1p(-theta)) + V * M.sol.Phi(c));
}

} // namespace vf::thinfilm
