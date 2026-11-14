#include <cstdio>
#include <algorithm>
#include "vf/thinfilm/SurfaceExchange.h"
using namespace vf;
using namespace vf::thinfilm;

int main()
{
	std::puts("h0_um,area_ratio,tau_bath_s,dt_s,t_s,theta,c_bulk_mM,gamma_N_m,equilibrium_gamma_N_m,insoluble_gamma_N_m,energy_J,balance_relative,energy_increase_relative");
	for (real h : {2e-6, 6e-6}) for (real area : {.99, .999, 1.001, 1.01}) for (real tau : {.0001, .001, .01, .1}) for (int half : {1, 2}) {
		FilmMaterial M = SDS_Solution();
		M.surface_desorption_rate_per_s = 1. / (tau * (1. + Exchange_Activity(M, M.c_bulk)));
		const real Q = 2.*M.Gamma_inf*M.theta0 + h*M.c_bulk, initial = M.theta0/area;
		const real equilibrium = Exchange_Equilibrium(M, area, h, Q);
		auto gamma = [&](real th) { return M.gamma0 - M.Pi_s(th); };
		const real dt = std::min(tau/20., 5e-5)/half, T = .02;
		const int n = (int)std::lround(T/dt), stride = std::max(1, n/400);
		real th = initial, energy = Exchange_Free_Energy(M, area, h, Q, th), up = 0.;
		for (int i = 0; i <= n; ++i) {
			const real bulk = Q - 2.*area*M.Gamma_inf*th, c = bulk/h;
			const real balance = std::abs(2.*area*M.Gamma_inf*th + bulk - Q)/Q;
			if (i % stride == 0 || i == n) std::printf("%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
				h*1e6,area,tau,dt,i*dt,th,c,gamma(th),gamma(equilibrium),gamma(initial),energy,balance,up);
			if (i == n) break;
			th = Exchange_Step(M, area, h, Q, th, dt);
			const real next = Exchange_Free_Energy(M, area, h, Q, th);
			up = std::max(up, (next-energy)/(si::Rgas*M.Temp*Q)); energy = next;
		}
	}
	return 0;
}
