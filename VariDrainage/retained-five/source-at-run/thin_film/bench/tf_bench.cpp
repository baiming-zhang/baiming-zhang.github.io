// Thin-film benchmarks: film thickness, drainage and black films (ThinFilm.h, thin_film/README.md). Immobile surfaces (stage 1).
//
// usage: tf_bench potential                 P(h), Pi(h), P''(h) of the disjoining potential (VF_EPS, VF_BETA, VF_LAMBDA)
//        tf_bench modes [nx ...]            growth / decay rate of thickness modes on a particle patch vs s(k) (dx = 64 / nx)
//        tf_bench front [dx ...]            a black film opening in a thick film: particle strip vs the 1D finite-difference
//                                           reference (dx = 1/8), same time steps
//        tf_bench spots [nx ...]            black spots nucleating and growing in a metastable thick film (2D patch)
//        tf_bench wave [dx ...]             longitudinal waves with mobile surfaces (Langmuir, Trouton) vs the exact dispersion (VF_MU)
//        tf_bench vmobile [dx]              vertical film with mobile surfaces: Marangoni support, drainage, surface diffusion
//                                           (VF_RT, VF_G0, VF_GINF, VF_DS, VF_RHO, VF_MU; VF_FRAMES + VF_EVERY: x y h Gamma frames)
//        tf_bench mbubble [dx]              spherical bubble with mobile surfaces (Langmuir, drainage, diffusion; VF_R, VF_RT, VF_DS)
//        tf_bench regenopen [dx]            marginal regeneration with an open meniscus (surface exchanged with a bath; remeshing)
//        tf_bench framed                    vertical film in a fibre frame with Plateau-border menisci, physical units (Meniscus.h;
//                                           VF_W_M, VF_H_M, VF_DX_M, VF_H0_M, VF_DT_S, VF_TMAX_S, VF_RREF_M, VF_FIBRE_M, VF_SIDES;
//                                           VF_SDS=1: soluble SDS, Monier et al.'s solution, VF_C_MM, VF_MU_PA, VF_RHO)
//        tf_bench regen [dx]                marginal regeneration at the bottom border of a vertical film (mobile surfaces;
//                                           border = reservoir VF_HM; VF_WIDTH, VF_HEIGHT, VF_NOISE, VF_RT, VF_MU, VF_REPORT)
//        tf_bench bubble [dx ...]           spherical bubble draining under gravity (immobile surfaces): apex law, axisymmetric
//                                           1D reference, black cap (VF_R, VF_CAP; VF_PREC=radius: sparsified factored preconditioner)
//        tf_bench vdrain [dx ...]           vertical film draining under gravity, immobile surfaces, closed ends: similarity
//                                           solution h = sqrt(4 mu y / (rho g t)) below the top, black film at the top (VF_RHOG, VF_H0)
// env:   VF_DT, VF_TMAX, VF_EVERY, VF_JITTER (lattice disorder in dx), VF_FRAMES=file (spots: thickness frames),
//        VF_SPARSE=1 (vdrain, vmobile: the exact sparse Hessian as preconditioner on the particle strip)
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include "vf/thinfilm/ThinFilm.h"
#include "vf/thinfilm/ThinFilmMechanics.h"
#include "vf/thinfilm/Meniscus.h"
#include "vf/thinfilm/Units.h"

using namespace vf;
using namespace vf::thinfilm;

static real Env(const char* k, real d) { const char* s = std::getenv(k); return s ? std::atof(s) : d; }

static std::vector<Vec3> Hex_Patch(int nx, int ny, real dx, real jitter, unsigned seed, real& Lx, real& Ly)
{
	std::mt19937 rng(seed);
	std::normal_distribution<real> N01;
	const real dy = dx * std::sqrt(3.) / 2.;
	Lx = nx * dx; Ly = ny * dy;
	std::vector<Vec3> x;
	for (int r = 0; r < ny; r++) for (int c = 0; c < nx; c++)
		x.push_back(Vec3((c + 0.5 * (r & 1)) * dx + jitter * dx * N01(rng), r * dy + jitter * dx * N01(rng), 0.));
	return x;
}

// Checkpoint of a film state plus a few counters (binary, native layout): to replay a step without rerunning to it.
template<class T> static void Put_Vec(FILE* f, const std::vector<T>& v) { const long long n = (long long)v.size(); std::fwrite(&n, sizeof n, 1, f); if (n) std::fwrite(v.data(), sizeof(T), (size_t)n, f); }
template<class T> static bool Get_Vec(FILE* f, std::vector<T>& v) { long long n = 0; if (std::fread(&n, sizeof n, 1, f) != 1) return false; v.resize((size_t)n); return !n || std::fread(v.data(), sizeof(T), (size_t)n, f) == (size_t)n; }
static void Save_State(const char* path, const FilmMatState& st, const std::vector<real>& counters)
{
	FILE* f = std::fopen(path, "wb");
	if (!f) { std::printf("  (cannot write %s)\n", path); return; }
	Put_Vec(f, st.x); Put_Vec(f, st.v); Put_Vec(f, st.x_prev); Put_Vec(f, st.v_prev); Put_Vec(f, st.V); Put_Vec(f, st.S); Put_Vec(f, st.fixed); Put_Vec(f, st.gscale);
	const real hp[2] = { st.dt_prev, st.has_prev ? 1. : 0. };
	std::fwrite(hp, sizeof(real), 2, f);
	Put_Vec(f, counters);
	std::fclose(f);
}
static bool Load_State(const char* path, FilmMatState& st, std::vector<real>& counters)
{
	FILE* f = std::fopen(path, "rb");
	if (!f) return false;
	real hp[2] = { 0., 0. };
	const bool ok = Get_Vec(f, st.x) && Get_Vec(f, st.v) && Get_Vec(f, st.x_prev) && Get_Vec(f, st.v_prev) && Get_Vec(f, st.V) && Get_Vec(f, st.S) && Get_Vec(f, st.fixed)
		&& Get_Vec(f, st.gscale) && std::fread(hp, sizeof(real), 2, f) == 2 && Get_Vec(f, counters);
	std::fclose(f);
	st.dt_prev = hp[0]; st.has_prev = hp[1] != 0.;
	return ok;
}

static Disjoining Disj_From_Env(real eps, real beta, real lambda)
{
	Disjoining d;
	d.eps = Env("VF_EPS", eps); d.hb = 1.; d.beta = Env("VF_BETA", beta); d.lambda = Env("VF_LAMBDA", lambda);
	return d;
}

static void Run_Potential()
{
	const Disjoining d = Disj_From_Env(1., 1., 1.5);
	std::printf("\n== Disjoining potential P = eps [(1/3)(hb/h)^8 - (4/3)(hb/h)^2] + beta exp(-h/lambda): eps %g, hb %g, beta %g, lambda %g ==\n", d.eps, d.hb, d.beta, d.lambda);
	std::printf("%8s %12s %12s %12s\n", "h", "P", "Pi", "P''");
	real prev = d.dP(0.5);
	for (real h = 0.6; h <= 20.; h *= 1.05) {
		const real dp = d.dP(h);
		const char* mark = (prev < 0. && dp >= 0.) ? "  <- minimum" : (prev > 0. && dp <= 0.) ? "  <- maximum (barrier)" : "";
		std::printf("%8.3f %12.5f %12.5f %12.5f%s\n", h, d.P(h), d.Pi(h), d.d2P(h), mark);
		prev = dp;
	}
}

// Thickness modes h0 + delta cos(k x) on a periodic particle patch of length 64; growth per backward-Euler step G,
// s_measured = (1 - 1/G) / dt (exact inversion for a linear mode), against s(k) = -M k^2 (sigma k^2 / 4 + P''(h0)).
static void Run_Modes(const std::vector<int>& nxs)
{
	const real L = 64., jitter = Env("VF_JITTER", 0.);
	std::printf("\n== Thickness modes, immobile surfaces (sigma = 1, mu = 1, eps = 1, hb = 1, beta = 0; patch length 64, kernel radius 4 dx; jitter %.2f dx) ==\n", jitter);
	std::printf("%6s %6s %8s %6s %8s %12s %12s %10s %8s\n", "nx", "dx", "N", "h0", "k", "s exact", "s measured", "rel err", "newton");
	for (int nx : nxs) {
		const real dx = L / nx;
		int ny = (int)std::ceil(2.2 * 4. * dx / (dx * std::sqrt(3.) / 2.));
		ny += ny & 1;
		real Lx, Ly;
		std::vector<Vec3> x = Hex_Patch(nx, ny, dx, jitter, 11, Lx, Ly);
		Quintic2D ker(4. * dx);
		ThinFilmGraph G = Make_Periodic_Patch(x, Lx, Ly, ker);
		Drainage_Stepper st;
		st.prm.sigma = 1.; st.prm.mu = 1.; st.prm.disj = Disj_From_Env(1., 0., 1.);
		st.rel_tol = 1e-8; // the perturbation is 1e-6: tighter is below round-off
		for (real h0 : { 1.1, 3. }) for (int m : { 1, 2, 4, 6 }) {
			const real k = 2. * kPi * m / Lx, delta = 1e-6;
			const real s_ex = -st.prm.Mobility(h0) * k * k * (st.prm.sigma / 4. * k * k + st.prm.disj.d2P(h0));
			// backward Euler stays a convex step only while s dt < 1 for every growing mode: limit dt by the fastest one too
			const real P2 = st.prm.disj.d2P(h0), s_max = P2 < 0. ? st.prm.Mobility(h0) * P2 * P2 / st.prm.sigma : 0.;
			const real dt = 0.1 / std::max(std::abs(s_ex), s_max);
			std::vector<real> V(G.N);
			for (int i = 0; i < G.N; i++) V[i] = (h0 + delta * std::cos(k * G.x[i][0])) * G.a[i];
			auto amp = [&]() {
				real num = 0., den = 0.;
				for (int i = 0; i < G.N; i++) { const real c = std::cos(k * G.x[i][0]); num += G.a[i] * (V[i] / G.a[i] - h0) * c; den += G.a[i] * c * c; }
				return num / den;
			};
			const real A0 = amp();
			const NewtonStats ns = st.Step(G, V, dt);
			const real A1 = amp(), s_m = (1. - A0 / A1) / dt;
			std::printf("%6d %6.3f %8d %6.2f %8.4f %12.4e %12.4e %+10.2e %8d\n", nx, dx, G.N, h0, k, s_ex, s_m, s_m / s_ex - 1., ns.iters);
		}
	}
}

// Black-film profile measures: black area (h < h_mid) per unit width and the thickest point
struct Profile { real black_halfwidth = 0., hmax = 0., hmin = 1e300, volume = 0., energy = 0.; };
static Profile Measure(const ThinFilmGraph& G, const std::vector<real>& V, real width, real hmid, const ThinFilmEnergy& E)
{
	Profile p;
	for (int i = 0; i < G.N; i++) {
		const real h = V[i] / G.a[i];
		if (h < hmid) p.black_halfwidth += G.a[i];
		p.hmax = std::max(p.hmax, h); p.hmin = std::min(p.hmin, h); p.volume += V[i];
	}
	p.black_halfwidth /= 2. * width;
	p.energy = E.Value(V);
	return p;
}

// A black film opening in a thick film (eps = 1, hb = 1, beta = 0, h0 = 3: the thick film is spinodally unstable, but
// the localized dimple ruptures first and the black film opens with a rim on each side). Domain 200, periodic;
// the dimple h0 - 1.8 exp(-(x - 100)^2 / 16). Particle strips at the given dx against the 1D grid at dx = 1/8.
static void Run_Front(const std::vector<real>& dxs)
{
	const real L = 200., h0 = 3., dt = Env("VF_DT", 0.5), tmax = Env("VF_TMAX", 400.), hmid = 2.;
	const int every = (int)Env("VF_EVERY", 50.);
	ThinFilmParams prm;
	prm.sigma = 1.; prm.mu = 1.; prm.disj = Disj_From_Env(1., 0., 1.);
	auto h_init = [&](real x) { return h0 - 1.8 * std::exp(-(x - 0.5 * L) * (x - 0.5 * L) / 16.); };
	std::printf("\n== Black film opening in a thick film (h0 = 3, eps = 1, hb = 1, sigma = 1, mu = 1; L = 200, dt = %g): half-width of the black film (h < 2) and rim height ==\n", dt);
	// reference: 1D grid
	struct Run { std::string name; ThinFilmGraph G; real width; std::vector<real> V; std::vector<Profile> prof; double sec = 0.; int newton = 0, cg = 0; bool mono = true; };
	std::vector<Run> runs;
	{
		Run r; r.name = "grid 1/8";
		const int n = (int)std::lround(L * 8.);
		r.G = Make_Grid_1D(n, L / n); r.width = 1.;
		r.V.resize(n);
		for (int i = 0; i < n; i++) r.V[i] = h_init(r.G.x[i][0]) * r.G.a[i];
		runs.push_back(std::move(r));
	}
	for (real dx : dxs) {
		Run r;
		char buf[64]; std::snprintf(buf, sizeof buf, "particles dx %g", dx); r.name = buf;
		int ny = (int)std::ceil(2.2 * 4. * dx / (dx * std::sqrt(3.) / 2.));
		ny += ny & 1;
		real Lx, Ly;
		std::vector<Vec3> x = Hex_Patch((int)std::lround(L / dx), ny, dx, Env("VF_JITTER", 0.), 3, Lx, Ly);
		Quintic2D ker(4. * dx);
		r.G = Make_Periodic_Patch(x, Lx, Ly, ker); r.width = Ly;
		r.V.resize(r.G.N);
		for (int i = 0; i < r.G.N; i++) r.V[i] = h_init(r.G.x[i][0]) * r.G.a[i];
		runs.push_back(std::move(r));
	}
	Drainage_Stepper st; st.prm = prm; st.newton.verbose = (int)Env("VF_VERBOSE", 0.);
	const int nsteps = (int)std::lround(tmax / dt);
	for (Run& r : runs) {
		ThinFilmEnergy E(r.G, prm);
		st.sparse_precond = r.name.rfind("grid", 0) == 0; // exact sparse Hessian on the 1D grid
		const auto t0 = std::chrono::steady_clock::now();
		real Eprev = E.Value(r.V);
		r.prof.push_back(Measure(r.G, r.V, r.width, hmid, E));
		for (int n = 1; n <= nsteps; n++) {
			const NewtonStats ns = st.Step(r.G, r.V, dt);
			r.newton += ns.iters; r.cg += ns.total_cg_iters;
			const real En = E.Value(r.V);
			if (En > Eprev + 1e-12 * std::abs(Eprev)) r.mono = false;
			Eprev = En;
			if (n % every == 0) r.prof.push_back(Measure(r.G, r.V, r.width, hmid, E));
		}
		r.sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
		std::printf("  (%s done: %.1f s, Newton %.1f / step, CG %.0f / step)\n", r.name.c_str(), r.sec, (double)r.newton / nsteps, (double)r.cg / nsteps);
	}
	std::printf("%8s", "t");
	for (const Run& r : runs) std::printf(" | %-24s", r.name.c_str());
	std::printf("\n%8s", "");
	for (size_t q = 0; q < runs.size(); q++) std::printf(" | %11s %12s", "half-width", "rim h");
	std::printf("\n");
	for (size_t s = 0; s < runs[0].prof.size(); s++) {
		std::printf("%8.1f", s * every * dt);
		for (const Run& r : runs) std::printf(" | %11.3f %12.4f", r.prof[s].black_halfwidth, r.prof[s].hmax);
		std::printf("\n");
	}
	for (const Run& r : runs) {
		const real dv = std::abs(r.prof.back().volume / r.prof.front().volume - 1.);
		std::printf("  %-18s N %6d  volume drift %.1e  energy %.6f -> %.6f (monotone %s)  min h %.4f  Newton %.1f / step, CG %.0f / step  %.1f s\n", r.name.c_str(), r.G.N, dv,
			r.prof.front().energy / r.width, r.prof.back().energy / r.width, r.mono ? "yes" : "NO", r.prof.back().hmin, (double)r.newton / nsteps, (double)r.cg / nsteps, r.sec);
	}
}

// connected black regions (h < hmid) over pairs closer than 1.5 dx
static int Count_Spots(const ThinFilmGraph& G, const std::vector<real>& V, real hmid, real dx, real Lx, real Ly)
{
	std::vector<int> par(G.N);
	std::iota(par.begin(), par.end(), 0);
	std::function<int(int)> find = [&](int i) { while (par[i] != i) { par[i] = par[par[i]]; i = par[i]; } return i; };
	auto wrap = [&](Vec3 d) { d[0] -= Lx * std::round(d[0] / Lx); d[1] -= Ly * std::round(d[1] / Ly); return d; };
	for (int q = 0; q < G.pairs.Size(); q++) {
		const int i = G.pairs.pi[q], j = G.pairs.pj[q];
		if (V[i] / G.a[i] < hmid && V[j] / G.a[j] < hmid && wrap(G.x[j] - G.x[i]).norm() < 1.5 * dx) par[find(i)] = find(j);
	}
	int n = 0;
	for (int i = 0; i < G.N; i++) if (V[i] / G.a[i] < hmid && find(i) == i) n++;
	return n;
}

// Black spots in a metastable thick film: eps = 1, hb = 1 with the double-layer barrier (beta, lambda); the thick
// film h0 sits in the local minimum of P. Nuclei: localized dimples deep enough to cross the barrier, at random
// positions and sizes; plus small noise. Square periodic patch of side nx dx (dx = 1).
static void Run_Spots(const std::vector<int>& nxs)
{
	const Disjoining d = Disj_From_Env(1., 1., 1.5);
	const real h0 = Env("VF_H0", 6.), dt = Env("VF_DT", 2.), tmax = Env("VF_TMAX", 2000.), hmid = Env("VF_HMID", 2.);
	const int every = (int)Env("VF_EVERY", 50.), nnuc = (int)Env("VF_NUCLEI", 6.);
	std::printf("\n== Black spots in a metastable thick film (h0 = %g; eps %g, beta %g, lambda %g; P(h0) = %.4f, P''(h0) = %.4f; sigma = 1, mu = 1; dt = %g) ==\n",
		h0, d.eps, d.beta, d.lambda, d.P(h0), d.d2P(h0), dt);
	for (int nx : nxs) {
		const real dx = 1.;
		real Lx, Ly;
		std::vector<Vec3> x = Hex_Patch(nx, (int)std::lround(nx * 2. / std::sqrt(3.)) & ~1, dx, Env("VF_JITTER", 0.1), 5, Lx, Ly);
		Quintic2D ker(4. * dx);
		ThinFilmGraph G = Make_Periodic_Patch(x, Lx, Ly, ker);
		ThinFilmParams prm; prm.sigma = 1.; prm.mu = 1.; prm.disj = d;
		std::mt19937 rng(17);
		std::uniform_real_distribution<real> U01(0., 1.);
		std::normal_distribution<real> N01;
		std::vector<Vec3> nuc;
		std::vector<real> rad;
		for (int q = 0; q < nnuc; q++) { nuc.push_back(Vec3(U01(rng) * Lx, U01(rng) * Ly, 0.)); rad.push_back(3. + 4. * U01(rng)); }
		auto wrap = [&](Vec3 v) { v[0] -= Lx * std::round(v[0] / Lx); v[1] -= Ly * std::round(v[1] / Ly); return v; };
		std::vector<real> V(G.N);
		for (int i = 0; i < G.N; i++) {
			real h = h0 * (1. + 0.01 * N01(rng));
			for (int q = 0; q < nnuc; q++) { const real r2 = wrap(G.x[i] - nuc[q]).squaredNorm(); h -= (h0 - 1.2) * std::exp(-r2 / (rad[q] * rad[q])); }
			V[i] = std::max(h, 1.05) * G.a[i];
		}
		ThinFilmEnergy E(G, prm);
		Drainage_Stepper st; st.prm = prm;
		real V0 = 0.;
		for (real v : V) V0 += v;
		real Eprev = E.Value(V);
		bool mono = true;
		FILE* fr = std::getenv("VF_FRAMES") ? std::fopen(std::getenv("VF_FRAMES"), "w") : nullptr;
		auto frame = [&](real t) {
			if (!fr) return;
			std::fprintf(fr, "t %.2f\n", t);
			for (int i = 0; i < G.N; i++) std::fprintf(fr, "%.3f %.3f %.4f\n", G.x[i][0], G.x[i][1], V[i] / G.a[i]);
		};
		std::printf("\n  patch %.0f x %.1f, %d particles, %d nuclei\n", Lx, Ly, G.N, nnuc);
		std::printf("  %8s %8s %10s %10s %10s %14s %10s %8s %8s\n", "t", "spots", "black frac", "min h", "max h", "energy", "vol drift", "newton", "sec");
		const auto t0 = std::chrono::steady_clock::now();
		auto report = [&](real t, int it) {
			real black = 0., area = 0., hmin = 1e300, hmax = 0., Vt = 0.;
			for (int i = 0; i < G.N; i++) { const real h = V[i] / G.a[i]; area += G.a[i]; if (h < hmid) black += G.a[i]; hmin = std::min(hmin, h); hmax = std::max(hmax, h); Vt += V[i]; }
			std::printf("  %8.1f %8d %10.4f %10.4f %10.4f %14.6f %10.1e %8d %8.1f\n", t, Count_Spots(G, V, hmid, dx, Lx, Ly), black / area, hmin, hmax, E.Value(V),
				std::abs(Vt / V0 - 1.), it, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
		};
		report(0., 0);
		frame(0.);
		const int nsteps = (int)std::lround(tmax / dt);
		for (int n = 1; n <= nsteps; n++) {
			const NewtonStats ns = st.Step(G, V, dt);
			const real En = E.Value(V);
			if (En > Eprev + 1e-12 * std::abs(Eprev)) mono = false;
			Eprev = En;
			if (n % every == 0 || n == nsteps) { report(n * dt, ns.iters); frame(n * dt); }
		}
		if (fr) std::fclose(fr);
		std::printf("  energy monotone: %s\n", mono ? "yes" : "NO");
	}
}

// A vertical film draining under gravity with immobile surfaces, closed at the top (y = H) and the bottom (y = 0).
// Away from the ends the top part follows the similarity solution of h_t = -(rho g / 4 mu) h^2 h_y' (y' = H - y, depth
// below the top): h = sqrt(4 mu y' / (rho g t)), valid while h < h0, i.e. for y' < rho g h0^2 t / (4 mu); the drained
// liquid collects at the bottom. With the barrier potential the top reaches the black film once it is thin enough.
// 1D grid (dx = 1/4) and a particle strip (dx = 1, periodic in x, closed in y).
static void Run_VDrain(const std::vector<real>& dxs)
{
	const real H = Env("VF_HEIGHT", 200.), h0 = Env("VF_H0", 6.), rho_g = Env("VF_RHOG", 0.01), dt = Env("VF_DT", 2.), tmax = Env("VF_TMAX", 2000.);
	ThinFilmParams prm;
	prm.sigma = 1.; prm.mu = 1.; prm.rho_g = rho_g; prm.disj = Disj_From_Env(1., 1., 1.5);
	std::printf("\n== Vertical film draining with immobile surfaces (height %g, h0 %g, rho g %g, mu 1, sigma 1; eps %g, beta %g, lambda %g; dt %g) ==\n",
		H, h0, rho_g, prm.disj.eps, prm.disj.beta, prm.disj.lambda, dt);
	std::printf("   similarity solution h = sqrt(4 mu y' / (rho g t)) below the top (y' = depth); error: rms of h / h_sim - 1 over y' in [0.15, 0.6] rho g h0^2 t / (4 mu)\n");
	struct Run { std::string name; ThinFilmGraph G; real width; std::vector<real> V; };
	std::vector<Run> runs;
	{
		Run r; r.name = "grid 1/4";
		const int n = (int)std::lround(H * 4.);
		r.G = Make_Grid_1D(n, H / n, false); r.width = 1.;
		r.G.height.resize(n);
		for (int k = 0; k < n; k++) r.G.height[k] = r.G.x[k][0];
		r.V.resize(n);
		for (int k = 0; k < n; k++) r.V[k] = h0 * r.G.a[k];
		runs.push_back(std::move(r));
	}
	for (real dx : dxs) {
		Run r;
		char buf[64]; std::snprintf(buf, sizeof buf, "particles dx %g", dx); r.name = buf;
		const real dy = dx * std::sqrt(3.) / 2.;
		const int nx = (int)std::ceil(2.2 * 4.), ny = (int)std::lround(H / dy);
		std::vector<Vec3> x;
		for (int rr = 0; rr < ny; rr++) for (int c = 0; c < nx; c++) x.push_back(Vec3((c + 0.5 * (rr & 1)) * dx, (rr + 0.5) * dy, 0.));
		Quintic2D ker(4. * dx);
		r.G = Make_Periodic_Patch(x, nx * dx, ny * dy, ker, false); r.width = nx * dx;
		r.G.height.resize(r.G.N);
		for (int i = 0; i < r.G.N; i++) r.G.height[i] = r.G.x[i][1];
		r.V.resize(r.G.N);
		for (int i = 0; i < r.G.N; i++) r.V[i] = h0 * r.G.a[i];
		runs.push_back(std::move(r));
	}
	const int nsteps = (int)std::lround(tmax / dt);
	const std::vector<real> report_t = { 250., 500., 1000., 2000., 4000., 8000. };
	for (Run& r : runs) {
		Drainage_Stepper st; st.prm = prm;
		st.sparse_precond = r.name.rfind("grid", 0) == 0 || Env("VF_SPARSE", 0.) > 0.; // exact sparse Hessian (banded on a strip)
		ThinFilmEnergy E(r.G, prm);
		real V0 = 0.;
		for (real v : r.V) V0 += v;
		real Eprev = E.Value(r.V);
		bool mono = true;
		const auto t0 = std::chrono::steady_clock::now();
		std::printf("\n  %s (%d unknowns)\n  %8s %12s %14s %12s %12s %10s %8s\n", r.name.c_str(), r.G.N, "t", "sim error", "black at top", "h(y'=20)", "h_sim(20)", "vol drift", "sec");
		size_t next = 0;
		for (int n = 1; n <= nsteps; n++) {
			st.Step(r.G, r.V, dt);
			const real En = E.Value(r.V);
			if (En > Eprev + 1e-12 * std::abs(Eprev)) mono = false;
			Eprev = En;
			const real t = n * dt;
			if (next < report_t.size() && std::abs(t - report_t[next]) < 0.5 * dt) {
				next++;
				// thickness profile vs depth: bin by depth (bins of 1)
				const int nb = (int)std::ceil(H);
				std::vector<real> sh(nb, 0.), sa(nb, 0.);
				real black = 0., Vt = 0.;
				for (int i = 0; i < r.G.N; i++) {
					const real yd = H - r.G.Height(i), hi = r.V[i] / r.G.a[i];
					const int b = std::clamp((int)yd, 0, nb - 1);
					sh[b] += r.V[i]; sa[b] += r.G.a[i];
					if (hi < 2.) black += r.G.a[i];
					Vt += r.V[i];
				}
				const real yf = rho_g * h0 * h0 * t / 4.;
				real e2 = 0.; int ne = 0;
				for (int b = 0; b < nb; b++) {
					const real yd = b + 0.5;
					if (yd < 0.15 * yf || yd > 0.6 * yf || sa[b] <= 0.) continue;
					const real hs = std::sqrt(4. * yd / (rho_g * t)), hm = sh[b] / sa[b];
					e2 += (hm / hs - 1.) * (hm / hs - 1.); ne++;
				}
				const int b20 = 20;
				std::printf("  %8.0f %12.4f %14.2f %12.4f %12.4f %10.1e %8.1f\n", t, ne ? std::sqrt(e2 / ne) : 0., black / r.width,
					sa[b20] > 0. ? sh[b20] / sa[b20] : 0., std::sqrt(4. * 20.5 / (rho_g * t)), std::abs(Vt / V0 - 1.),
					std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
			}
		}
		std::printf("  energy monotone: %s\n", mono ? "yes" : "NO");
	}
}

// Cost of the factored preconditioner on a sphere: build, factorizations, solves (diagnostics).
static void Run_Precond_Timing(real R)
{
	using clk = std::chrono::steady_clock;
	auto sec = [](clk::time_point a) { return std::chrono::duration<double>(clk::now() - a).count(); };
	const real dx = 1.;
	const int n = (int)std::lround(4. * kPi * R * R / (std::sqrt(3.) / 2. * dx * dx));
	Quintic2D ker(4. * dx);
	auto t0 = clk::now();
	ThinFilmGraph G = Make_Surface_Graph(Fibonacci_Sphere(n, R), ker);
	std::printf("sphere R %g: %d particles, %d pairs (%.1f per particle), graph %.2f s\n", R, n, G.pairs.Size(), 2. * G.pairs.Size() / n, sec(t0));
	ThinFilmParams prm; prm.rho_g = 0.01; prm.disj.eps = 1.; prm.disj.beta = 1.; prm.disj.lambda = 1.5;
	std::vector<real> V(n);
	for (int i = 0; i < n; i++) V[i] = 6. * G.a[i];
	Drainage_Stepper st; st.prm = prm;
	std::vector<real> k;
	st.Conductances(G, V, k);
	ThinFilmEnergy E(G, prm);
	// a film with a thickness bump, so that the step does real work
	std::vector<real> Vb(n);
	for (int i = 0; i < n; i++) Vb[i] = (6. + 1.5 * std::exp(-((G.x[i] - Vec3(R, 0., 0.)).squaredNorm()) / 25.)) * G.a[i];
	std::printf("%8s %12s %12s %10s %10s %10s\n", "radius", "1st prep s", "LU s", "apply s", "step CG", "step s");
	for (real rad : { 0., 2.5, 1.8, 1.5 }) {
		if (rad == 0. && n > 4000) continue;
		Drainage_Problem P(E, V, k, 2.);
		P.factored_precond = true; P.precond_radius = rad;
		VecX phi = VecX::Zero(n), r = VecX::Random(n), z;
		P.Prepare_Hessian(phi);
		t0 = clk::now(); P.Prepare_Preconditioner(); const double t1 = sec(t0);
		t0 = clk::now(); P.Prepare_Preconditioner(); const double t2 = sec(t0);
		t0 = clk::now(); P.Apply_Preconditioner(r, z); const double t3 = sec(t0);
		Drainage_Stepper sf = st; sf.factored_precond = true; sf.precond_radius = rad;
		std::vector<real> Vs = Vb;
		t0 = clk::now(); const NewtonStats ns = sf.Step(G, Vs, 2.); const double t4 = sec(t0);
		std::printf("%8.2f %12.2f %12.2f %10.3f %10d %10.2f%s\n", rad, t1, t2, t3, ns.total_cg_iters, t4, ns.converged ? "" : " (not converged)");
	}
	{
		std::vector<real> Vs = Vb;
		t0 = clk::now(); const NewtonStats ns = st.Step(G, Vs, 2.);
		std::printf("%8s %12s %12s %10s %10d %10.2f%s\n", "Jacobi", "-", "-", "-", ns.total_cg_iters, sec(t0), ns.converged ? "" : " (not converged)");
	}
}

// Two-level (fine/coarse particle) preconditioner against Jacobi on spheres of increasing fine resolution: one drainage
// step of a film with a bump and the barrier potential. Fine particles at dx (kernel 4 dx), coarse at VF_COARSE
// (prolongation radius 2 dx_coarse).
static void Run_Two_Level_Timing(const std::vector<real>& dxs)
{
	using clk = std::chrono::steady_clock;
	auto sec = [](clk::time_point a) { return std::chrono::duration<double>(clk::now() - a).count(); };
	const real R = Env("VF_R", 24.), dxc = Env("VF_COARSE", 2.), dt = Env("VF_DT", 2.);
	ThinFilmParams prm; prm.rho_g = 0.01; prm.disj.eps = 1.; prm.disj.beta = 1.; prm.disj.lambda = 1.5;
	const int nc = (int)std::lround(4. * kPi * R * R / (std::sqrt(3.) / 2. * dxc * dxc));
	const std::vector<Vec3> xc = Fibonacci_Sphere(nc, R);
	std::printf("\n== Two-level preconditioner on a sphere R %g: coarse dx %g (%d particles), Chebyshev degree %d; one drainage step dt %g ==\n",
		R, dxc, nc, (int)Env("VF_CHEB", 3.), dt);
	std::printf("%8s %8s %14s %10s %14s %10s %10s %12s\n", "fine dx", "N", "Jacobi CG", "Jacobi s", "2-level CG", "V-cycles", "2-level s", "dV / change");
	for (real dx : dxs) {
		const int n = (int)std::lround(4. * kPi * R * R / (std::sqrt(3.) / 2. * dx * dx));
		Quintic2D ker(4. * dx);
		const ThinFilmGraph G = Make_Surface_Graph(Fibonacci_Sphere(n, R), ker);
		std::vector<real> V0(n);
		for (int i = 0; i < n; i++) V0[i] = (6. + 1.5 * std::exp(-((G.x[i] - Vec3(R, 0., 0.)).squaredNorm()) / 25.)) * G.a[i];
		const SpMat P = Make_Prolongation(G.x, xc, 2. * dxc);
		Drainage_Stepper sj; sj.prm = prm; sj.rel_tol = 1e-9; sj.newton.max_cg_iters = 50000;
		std::vector<real> Vj = V0, Vt = V0;
		auto t0 = clk::now();
		const NewtonStats nj = sj.Step(G, Vj, dt);
		const double tj = sec(t0);
		Drainage_Stepper st = sj; st.prolong = &P; st.cheb_degree = (int)Env("VF_CHEB", 3.);
		t0 = clk::now();
		const NewtonStats ntw = st.Step(G, Vt, dt);
		const double tt = sec(t0);
		real dV = 0., ch = 0.;
		for (int i = 0; i < n; i++) { dV = std::max(dV, std::abs(Vj[i] - Vt[i])); ch = std::max(ch, std::abs(Vj[i] - V0[i])); }
		std::printf("%8.2f %8d %10d %3s %10.2f %10d %3s %10d %10.2f %12.1e\n", dx, n, nj.total_cg_iters, nj.converged ? "" : "nc", tj,
			ntw.total_cg_iters, ntw.converged ? "" : "nc", ntw.precond_factorizations, tt, dV / ch);
		// setup cost of one two-level preconditioner at the start of the step, and one V-cycle
		{
			std::vector<real> k;
			sj.Conductances(G, V0, k);
			ThinFilmEnergy E(G, prm);
			Drainage_Problem Pb(E, V0, k, dt);
			Pb.prolong = &P; Pb.cheb_degree = st.cheb_degree;
			VecX phi = VecX::Zero(n), r(n), z;
			for (int i = 0; i < n; i++) r[i] = std::sin(0.37 * i);
			r.array() -= r.mean();
			Pb.Prepare_Hessian(phi);
			Pb.Prepare_Preconditioner();
			t0 = clk::now(); Pb.Apply_Preconditioner(r, z); const double tv = sec(t0);
			std::printf("           setup: Galerkin %.2f s, coarse factorization %.2f s, power iteration %.2f s; one V-cycle %.3f s; nnz(P) %d\n",
				Pb.t_galerkin, Pb.t_factor, Pb.t_power, tv, (int)P.nonZeros());
		}
	}
}

// A spherical bubble (radius R) whose film drains under gravity with immobile surfaces. At the apex symmetry makes the
// thinning uniform: dh/dt = -(rho g / (6 mu R)) h^3 (squeeze and disjoining aside), so h_apex = (h0^-2 + rho g t / (3 mu
// R))^-1/2; the liquid collects at the bottom; with the barrier potential a black cap opens at the top. Particles on a
// Fibonacci sphere (factored preconditioner) against the axisymmetric 1D reference (zones of arc 1/4).
static void Run_Bubble(const std::vector<real>& dxs)
{
	const real R = Env("VF_R", 24.), h0 = Env("VF_H0", 6.), rho_g = Env("VF_RHOG", 0.01), dt = Env("VF_DT", 2.), tmax = Env("VF_TMAX", 2000.);
	ThinFilmParams prm;
	prm.sigma = 1.; prm.mu = 1.; prm.rho_g = rho_g; prm.disj = Disj_From_Env(1., 1., 1.5);
	const real cap = Env("VF_CAP", 3.); // apex readout: area-weighted mean thickness within this arc of the top
	std::printf("\n== Bubble draining with immobile surfaces (R %g, h0 %g, rho g %g, mu 1, sigma 1; eps %g, beta %g, lambda %g; dt %g) ==\n",
		R, h0, rho_g, prm.disj.eps, prm.disj.beta, prm.disj.lambda, dt);
	std::printf("   apex law (no squeeze, no disjoining) h = (h0^-2 + rho g t / (3 mu R))^-1/2; apex = mean within arc %g of the top;\n"
		"   black cap: arc from the top of the area with h < 2; bottom: mean within arc %g of the bottom\n", cap, cap);
	struct Run { std::string name; ThinFilmGraph G; std::vector<real> V; bool factored; };
	std::vector<Run> runs;
	{
		Run r; r.name = "axisymmetric 1D, arc 1/4"; r.factored = false;
		r.G = Make_Sphere_Axisym((int)std::lround(kPi * R * 4.), R);
		r.V.resize(r.G.N);
		for (int k = 0; k < r.G.N; k++) r.V[k] = h0 * r.G.a[k];
		runs.push_back(std::move(r));
	}
	for (real dx : dxs) {
		Run r; r.factored = true;
		const int n = (int)std::lround(4. * kPi * R * R / (std::sqrt(3.) / 2. * dx * dx));
		char buf[96]; std::snprintf(buf, sizeof buf, "particles dx %g (%d)", dx, n); r.name = buf;
		Quintic2D ker(4. * dx);
		r.G = Make_Surface_Graph(Fibonacci_Sphere(n, R), ker);
		r.V.resize(r.G.N);
		// VF_NOISE: relative random thickness perturbation (breaks the axial symmetry of the drainage; without it the
		// black film opens in rings)
		std::mt19937 rng(5);
		std::normal_distribution<real> N01;
		const real noise = Env("VF_NOISE", 0.);
		for (int i = 0; i < r.G.N; i++) r.V[i] = h0 * (1. + noise * N01(rng)) * r.G.a[i];
		runs.push_back(std::move(r));
	}
	const int nsteps = (int)std::lround(tmax / dt);
	std::vector<real> report_t = { 100., 250., 500., 1000., 2000., 4000., 8000. };
	if (const char* rt = std::getenv("VF_TIMES")) { // comma-separated report times
		report_t.clear();
		for (const char* p = rt; *p;) { char* e; const real v = std::strtod(p, &e); if (e == p) break; report_t.push_back(v); p = (*e == ',') ? e + 1 : e; }
	}
	const bool profile = Env("VF_PROFILE", 0.) > 0.; // the 1D run: print h along the arc from the top
	const bool only_1d = Env("VF_ONLY1D", 0.) > 0.;
	// particle runs: frames ("t T", then "x y z h" for every particle) from VF_FSTART every VF_EVERY
	const char* frames = std::getenv("VF_FRAMES");
	const real every = Env("VF_EVERY", 0.), fstart = Env("VF_FSTART", 0.);
	FILE* ff = frames ? std::fopen(frames, "w") : nullptr;
	for (Run& r : runs) {
		if (only_1d && r.factored) continue;
		Drainage_Stepper st; st.prm = prm;
		st.sparse_precond = !r.factored; // the 1D zones: exact banded Hessian
		st.factored_precond = r.factored && Env("VF_PREC", 0.) > 0.; // particles: Jacobi CG (cheap parallel Hessian products), or the sparsified factored preconditioner of this radius
		st.precond_radius = Env("VF_PREC", 0.);
		st.newton.max_cg_iters = 20000;
		st.newton.verbose = r.factored ? (int)Env("VF_VERBOSE", 0.) : 0;
		// VF_COARSE = coarse particle spacing: the two-level (MELP-like fine/coarse) preconditioner
		SpMat Pr, Wd; // prolongation (fine x coarse) and the area-weighted distribution of coarse volume changes to fine
		ThinFilmGraph Gc; // coarse graph (the thermal noise lives on it)
		std::vector<real> colw;
		if (r.factored && Env("VF_COARSE", 0.) > 0.) {
			const real dxc = Env("VF_COARSE", 0.);
			const std::vector<Vec3> xc = Fibonacci_Sphere((int)std::lround(4. * kPi * R * R / (std::sqrt(3.) / 2. * dxc * dxc)), R);
			Pr = Make_Prolongation(r.G.x, xc, 2. * dxc);
			st.prolong = &Pr;
			Gc = Make_Surface_Graph(xc, Quintic2D(4. * dxc));
			// W_lc = P_lc a_l / sum_l P_lc a_l: each column sums to 1 (a coarse volume change is shared exactly)
			colw.assign(xc.size(), 0.);
			for (int k = 0; k < Pr.outerSize(); k++) for (SpMat::InnerIterator it(Pr, k); it; ++it) colw[it.col()] += it.value() * r.G.a[it.row()];
			std::vector<Eigen::Triplet<real>> tw;
			for (int k = 0; k < Pr.outerSize(); k++) for (SpMat::InnerIterator it(Pr, k); it; ++it) tw.emplace_back((int)it.row(), (int)it.col(), it.value() * r.G.a[it.row()] / colw[it.col()]);
			Wd.resize(Pr.rows(), Pr.cols());
			Wd.setFromTriplets(tw.begin(), tw.end());
			std::printf("  (two-level preconditioner: %d coarse particles at spacing %g; thermal noise on the coarse particles)\n", (int)xc.size(), dxc);
		}
		ThinFilmEnergy E(r.G, prm);
		real V0 = 0.;
		for (real v : r.V) V0 += v;
		real Eprev = E.Value(r.V);
		bool mono = true;
		int cg = 0, nt = 0;
		const auto t0 = std::chrono::steady_clock::now();
		std::printf("\n  %s\n  %8s %10s %10s %10s %10s %10s %10s %8s\n", r.name.c_str(), "t", "h apex", "apex law", "black cap", "h bottom", "vol drift", "CG/step", "sec");
		size_t next = 0;
		// VF_KT: thermal fluctuations of the film (particle runs): after each drainage step a random liquid flux along every
		// pair with variance 2 kT k_q dt (k_q the pair conductance): the noise consistent with the dissipation of the
		// step (fluctuation-dissipation for the gradient flow in the metric of L; the stochastic thin-film equation,
		// Grun, Mecke & Rauscher 2006). Volume is conserved exactly; the energy then fluctuates.
		const real kT = r.factored ? Env("VF_KT", 0.) : 0.;
		std::mt19937 trng(17);
		std::normal_distribution<real> TN;
		std::vector<real> kq;
		for (int s = 1; s <= nsteps; s++) {
			// noise first, then the implicit step: the step relaxes the stiff (particle-scale) part of the kick in the same
			// step, the slow modes keep their thermal amplitude (noise added after an implicit step leaves the stiff modes
			// with the full kick, which grows as dx^-2 and can empty particles); a pair flux never takes more than half a
			// particle's liquid
			if (kT > 0. && Gc.N > 0) {
				// on the coarse particles (fluctuations coarse-grained at the coarse spacing; the particle-scale kick of a
				// step of this size would push the fine thickness into the disjoining wall): coarse thickness, conductances
				// and pair fluxes, the coarse volume changes shared among the fine particles by W (exactly conservative)
				VecX Vf = Eigen::Map<const VecX>(r.V.data(), r.G.N);
				const VecX Vc = Pr.transpose() * Vf;
				std::vector<real> hc(Gc.N), kc(Gc.pairs.Size());
				for (int c = 0; c < Gc.N; c++) hc[c] = Vc[c] / colw[c];
				VecX dVc = VecX::Zero(Gc.N);
				for (int q = 0; q < Gc.pairs.Size(); q++) {
					const int i = Gc.pairs.pi[q], j = Gc.pairs.pj[q];
					const real kcq = prm.Mobility(0.5 * (hc[i] + hc[j])) * Gc.c[q];
					const real F = std::sqrt(2. * kT * kcq * dt) * TN(trng);
					dVc[i] -= F; dVc[j] += F;
				}
				const VecX dVf = Wd * dVc;
				for (int i = 0; i < r.G.N; i++) r.V[i] = std::max(r.V[i] + dVf[i], 0.05 * r.V[i]);
			}
			else if (kT > 0.) {
				st.Conductances(r.G, r.V, kq);
				for (int q = 0; q < r.G.pairs.Size(); q++) {
					const int i = r.G.pairs.pi[q], j = r.G.pairs.pj[q];
					real F = std::sqrt(2. * kT * kq[q] * dt) * TN(trng);
					F = std::clamp(F, -0.5 * r.V[j], 0.5 * r.V[i]);
					r.V[i] -= F; r.V[j] += F;
				}
			}
			const NewtonStats ns = st.Step(r.G, r.V, dt);
			cg += ns.total_cg_iters; nt++;
			const real En = E.Value(r.V);
			if (kT == 0. && En > Eprev + 1e-12 * std::abs(Eprev)) mono = false;
			Eprev = En;
			const real t = s * dt;
			if (ff && r.factored && every > 0. && t >= fstart - 0.5 * dt && std::abs(t / every - std::round(t / every)) < 0.5 * dt / every) {
				std::fprintf(ff, "t %g\n", t);
				for (int i = 0; i < r.G.N; i++) std::fprintf(ff, "%.4f %.4f %.4f %.4f\n", r.G.x[i][0], r.G.x[i][1], r.G.x[i][2], r.V[i] / r.G.a[i]);
			}
			if (next < report_t.size() && std::abs(t - report_t[next]) < 0.5 * dt) {
				next++;
				real va = 0., aa = 0., vb = 0., ab = 0., black = 0., Vt = 0.;
				for (int i = 0; i < r.G.N; i++) {
					const real z = r.G.Height(i), th = std::acos(std::clamp(z / R, -1., 1.)), hi = r.V[i] / r.G.a[i];
					if (R * th < cap) { va += r.V[i]; aa += r.G.a[i]; }
					if (R * (kPi - th) < cap) { vb += r.V[i]; ab += r.G.a[i]; }
					if (hi < 2.) black += r.G.a[i];
					Vt += r.V[i];
				}
				const real law = 1. / std::sqrt(1. / (h0 * h0) + rho_g * t / (3. * prm.mu * R));
				const real arc = R * std::acos(std::clamp(1. - black / (2. * kPi * R * R), -1., 1.));
				std::printf("  %8.0f %10.4f %10.4f %10.2f %10.3f %10.1e %10.1f %8.1f\n", t, aa > 0. ? va / aa : 0., law, arc, ab > 0. ? vb / ab : 0.,
					std::abs(Vt / V0 - 1.), (double)cg / nt, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
				cg = nt = 0;
				if (profile && !r.factored) {
					std::printf("           h along the arc from the top:");
					for (int k = 0; k < r.G.N; k += 8) {
						const real arc_k = R * std::acos(std::clamp(r.G.Height(k) / R, -1., 1.));
						if (arc_k > 70.) break;
						std::printf(" %.0f:%.2f", arc_k, r.V[k] / r.G.a[k]);
					}
					std::printf("\n");
				}
			}
		}
		std::printf("  energy monotone: %s\n", mono ? "yes" : "NO");
	}
	if (ff) std::fclose(ff);
}

// A vertical film with mobile surfaces (Langmuir surfactant), drainage and surface diffusion (FilmCoupled_Stepper):
// a strip periodic in x, closed at y = 0 and y = H by frame rows of particles that do not move (thickness: one kernel
// radius). Initially uniform h0 and Gamma0 at rest. The film's weight pulls the surfaces down until the Marangoni stress
// carries it (2 dPi_s/dy = -rho g h, the top stretched, the bottom compressed): with enough elasticity the surfaces then
// stay put and the liquid drains between them as with immobile surfaces; surface diffusion relaxes the supporting
// surfactant gradient, so the surfaces keep sliding down at u ~ -Ds rho g h / E_G and the top thins by stretching too.
// Inertia rho is small (elastic waves are far faster than drainage; the surfaces follow quasi-statically).
static void Run_VMobile(real dx)
{
	const real H = Env("VF_HEIGHT", 200.), h0 = Env("VF_H0", 6.), rho_g = Env("VF_RHOG", 0.01), dt = Env("VF_DT", 2.), tmax = Env("VF_TMAX", 2000.);
	const real G0 = Env("VF_G0", 0.5), Ds = Env("VF_DS", 0.), rdx = Env("VF_RDX", 4.);
	FilmCoupled_Stepper cs;
	FilmMechParams& prm = cs.mech.prm;
	prm.tf.sigma = 1.; prm.tf.mu = Env("VF_MU", 1.); prm.tf.rho_g = rho_g; prm.tf.disj = Disj_From_Env(1., 1., 1.5);
	prm.surf.RT = Env("VF_RT", 20.); prm.surf.Gamma_inf = Env("VF_GINF", 1.);
	prm.rho = Env("VF_RHO", 1e-3);
	prm.up = Vec3(0., 1., 0.);
	cs.diff.Ds = Ds;
	cs.drain.sparse_precond = cs.diff.sparse_precond = Env("VF_SPARSE", 0.) > 0.; // exact sparse Hessians (banded on a strip)
	const real R = rdx * dx;
	cs.mech.kernel.Set_Radius(R);
	const real dy = dx * std::sqrt(3.) / 2.;
	const int nx = (int)std::ceil(2.2 * rdx), ny = (int)std::lround(H / dy);
	cs.mech.box.Lx = nx * dx; cs.mech.box.Ly = ny * dy; cs.mech.box.periodic_y = false;
	const real Ht = ny * dy, width = nx * dx;
	FilmMatState st;
	for (int rr = 0; rr < ny; rr++) for (int c = 0; c < nx; c++) st.x.push_back(Vec3((c + 0.5 * (rr & 1)) * dx, (rr + 0.5) * dy, 0.));
	const int n = st.Size();
	st.v.assign(n, Vec3::Zero()); st.fixed.assign(n, 0);
	for (int i = 0; i < n; i++) if (st.x[i][1] < R || st.x[i][1] > Ht - R) st.fixed[i] = 1;
	{
		const ThinFilmGraph G = cs.Graph(st);
		st.V.resize(n); st.S.resize(n);
		for (int i = 0; i < n; i++) { st.V[i] = h0 * G.a[i]; st.S[i] = G0 * G.a[i]; }
	}
	const real EG = 2. * G0 * G0 * prm.surf.d2g(G0);
	std::printf("\n== Vertical film, mobile surfaces (height %g, width %g, dx %g, frame rows within %g of the ends; h0 %g, rho g %g, mu %g, rho %g;\n"
		"   Langmuir RT %g, Gamma_inf %g, Gamma0 %g: E_G %.3f, Pi_s %.3f; Ds %g; eps 1, beta %g, lambda %g; dt %g) ==\n",
		Ht, width, dx, R, h0, rho_g, prm.tf.mu, prm.rho, prm.surf.RT, prm.surf.Gamma_inf, G0, EG, prm.surf.Pressure(G0), Ds, prm.tf.disj.beta, prm.tf.disj.lambda, dt);
	std::printf("   weight to carry: rho g h0 H / 2 = %.3f per surface (surface pressure drop top to bottom); sliding u ~ -Ds rho g h0 / E_G = %.2e\n",
		rho_g * h0 * Ht / 2., -Ds * rho_g * h0 / EG);
	std::printf("   immobile-surface similarity below the top: h = sqrt(4 mu y' / (rho g t)); error over y' in [0.15, 0.6] rho g h0^2 t / (4 mu)\n");
	std::printf("  %7s %9s %8s %9s %9s %9s %9s %9s %10s %9s %9s %6s %6s %7s\n", "t", "sim err", "black", "h(y'20)", "h_sim20", "G top", "G mid", "G bot",
		"v_y mid", "vol drft", "surf drft", "N mech", "N drn", "sec");
	real V0 = 0., S0 = 0.;
	for (int i = 0; i < st.Size(); i++) { V0 += st.V[i]; S0 += st.S[i]; }
	const char* frames = std::getenv("VF_FRAMES");
	FILE* ff = frames ? std::fopen(frames, "w") : nullptr;
	const real every = Env("VF_EVERY", 0.);
	const int nsteps = (int)std::lround(tmax / dt);
	const std::vector<real> report_t = { 10., 50., 250., 500., 1000., 2000., 4000., 8000. };
	size_t next = 0;
	int it_mech = 0, it_drain = 0, nrep = 0;
	// VF_REMESH: particle insertion and merging after every step (Remesh_Film; target area the initial lattice cell)
	const bool remesh = Env("VF_REMESH", 0.) > 0.;
	RemeshParams rp; rp.dx0 = dx; rp.kernel_radius = R;
	const real cell = dx * dy;
	int splits = 0, merges = 0;
	const auto t0 = std::chrono::steady_clock::now();
	for (int s = 1; s <= nsteps; s++) {
		const FilmCoupled_Stepper::Stats ss = cs.Step(st, dt);
		it_mech += ss.mech.iters; it_drain += ss.drain.iters; nrep++;
		if (!ss.mech.converged || !ss.drain.converged || !ss.diff.converged) {
			std::printf("  (step %d: Newton not converged: mech %d drain %d diff %d", s, (int)ss.mech.converged, (int)ss.drain.converged, (int)ss.diff.converged);
			// diagnostics: extremes of h and Gamma over the free particles, the closest pair, the last gradient norms
			const ThinFilmGraph Gd = cs.Graph(st);
			real hmn = 1e300, hmx = 0., gmn = 1e300, gmx = 0., rmin = 1e300; int ih = -1, ig = -1;
			for (int i = 0; i < st.Size(); i++) {
				if (st.fixed[i]) continue;
				const real h = st.V[i] / Gd.a[i], g = st.S[i] / Gd.a[i];
				if (h < hmn) { hmn = h; ih = i; }
				hmx = std::max(hmx, h);
				if (g > gmx) { gmx = g; ig = i; }
				gmn = std::min(gmn, g);
			}
			for (int q = 0; q < Gd.pairs.Size(); q++) rmin = std::min(rmin, cs.mech.box.Wrap(st.x[Gd.pairs.pj[q]] - st.x[Gd.pairs.pi[q]]).norm());
			std::printf("; h %.3f-%.3f (min at y %.1f), Gamma %.4f-%.4f (max at y %.1f), closest pair %.3f; mech iters %d |g| %.2e -> %.2e)\n", hmn, hmx,
				ih >= 0 ? st.x[ih][1] : -1., gmn, gmx, ig >= 0 ? st.x[ig][1] : -1., rmin, ss.mech.iters,
				ss.mech.grad_inf_history.empty() ? 0. : ss.mech.grad_inf_history.front(), ss.mech.grad_inf_history.empty() ? 0. : ss.mech.grad_inf_history.back());
		}
		if (remesh) {
			const ThinFilmGraph Gr = cs.Graph(st);
			const RemeshStats rs = Remesh_Film(st, Gr, rp, &cs.mech.box, nullptr, [&](const std::vector<Vec3>& xs) { FilmMatState tmp; tmp.x = xs; return cs.Graph(tmp); });
			splits += rs.splits; merges += rs.merges;
		}
		const int n = st.Size();
		const real t = s * dt;
		if (ff && every > 0. && std::abs(t / every - std::round(t / every)) < 0.5 * dt / every) {
			const ThinFilmGraph G = cs.Graph(st);
			std::fprintf(ff, "t %g\n", t);
			for (int i = 0; i < n; i++) std::fprintf(ff, "%.4f %.4f %.5f %.5f\n", st.x[i][0], st.x[i][1], st.V[i] / G.a[i], st.S[i] / G.a[i]);
		}
		if (next < report_t.size() && std::abs(t - report_t[next]) < 0.5 * dt) {
			next++;
			const ThinFilmGraph G = cs.Graph(st);
			if (remesh) {
				real amin = 1e300, amax = 0.;
				for (int i = 0; i < n; i++) if (!st.fixed[i]) { amin = std::min(amin, G.a[i]); amax = std::max(amax, G.a[i]); }
				std::printf("  (remeshing: %d particles; %d splits, %d merges so far; free particle areas %.2f-%.2f of the cell)\n", n, splits, merges, amin / cell, amax / cell);
			}
			const int nb = (int)std::ceil(Ht);
			std::vector<real> sh(nb, 0.), sa(nb, 0.);
			real black = 0., Vt = 0., St = 0., gt = 0., at = 0., gm = 0., am = 0., gbm = 0., ab = 0., vm = 0.;
			int nm = 0;
			for (int i = 0; i < n; i++) {
				const real y = st.x[i][1], yd = Ht - y, hi = st.V[i] / G.a[i], Gi = st.S[i] / G.a[i];
				const int b = std::clamp((int)yd, 0, nb - 1);
				sh[b] += st.V[i]; sa[b] += G.a[i];
				if (hi < 2.) black += G.a[i];
				Vt += st.V[i]; St += st.S[i];
				if (yd > R && yd < R + 10.) { gt += st.S[i]; at += G.a[i]; }
				if (std::abs(y - Ht / 2.) < 5.) { gm += st.S[i]; am += G.a[i]; vm += st.v[i][1]; nm++; }
				if (y > R && y < R + 10.) { gbm += st.S[i]; ab += G.a[i]; }
			}
			const real yf = rho_g * h0 * h0 * t / (4. * prm.tf.mu);
			real e2 = 0.; int ne = 0;
			for (int b = 0; b < nb; b++) {
				const real yd = b + 0.5;
				if (yd < 0.15 * yf || yd > 0.6 * yf || sa[b] <= 0.) continue;
				const real hs = std::sqrt(4. * prm.tf.mu * yd / (rho_g * t)), hm = sh[b] / sa[b];
				e2 += (hm / hs - 1.) * (hm / hs - 1.); ne++;
			}
			std::printf("  %7.0f %9.4f %8.2f %9.4f %9.4f %9.4f %9.4f %9.4f %+10.2e %9.1e %9.1e %6.1f %6.1f %7.1f\n", t, ne ? std::sqrt(e2 / ne) : 0., black / width,
				sa[20] > 0. ? sh[20] / sa[20] : 0., std::sqrt(4. * prm.tf.mu * 20.5 / (rho_g * t)), at > 0. ? gt / at : 0., am > 0. ? gm / am : 0., ab > 0. ? gbm / ab : 0.,
				nm ? vm / nm : 0., std::abs(Vt / V0 - 1.), std::abs(St / S0 - 1.), (double)it_mech / nrep, (double)it_drain / nrep,
				std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
			it_mech = it_drain = nrep = 0;
		}
	}
	if (ff) std::fclose(ff);
}

// A spherical bubble with mobile surfaces (Langmuir surfactant), drainage and surface diffusion: the film's weight
// pulls the surfaces down until the Marangoni stress carries it (top diluted, bottom compressed), then the liquid
// drains between them to the bottom and a black cap opens at the top. Against the immobile-surface bubble (bubble).
static void Run_MBubble(real dx)
{
	const real R = Env("VF_R", 16.), h0 = Env("VF_H0", 6.), rho_g = Env("VF_RHOG", 0.01), dt = Env("VF_DT", 2.), tmax = Env("VF_TMAX", 2000.);
	const real G0 = Env("VF_G0", 0.5), cap = Env("VF_CAP", 3.);
	FilmCoupled_Stepper cs;
	FilmMechParams& prm = cs.mech.prm;
	prm.tf.sigma = 1.; prm.tf.mu = Env("VF_MU", 1.); prm.tf.rho_g = rho_g; prm.tf.disj = Disj_From_Env(1., 1., 1.5);
	prm.surf.RT = Env("VF_RT", 20.); prm.surf.Gamma_inf = 1.;
	prm.rho = Env("VF_RHO", 1e-3);
	prm.up = Vec3(0., 0., 1.);
	cs.diff.Ds = Env("VF_DS", 0.);
	cs.drain.newton.max_cg_iters = cs.diff.newton.max_cg_iters = 20000;
	cs.mech.kernel.Set_Radius(4. * dx);
	cs.mech.on_sphere = true; cs.mech.sphere.R = R;
	const int n = (int)std::lround(4. * kPi * R * R / (std::sqrt(3.) / 2. * dx * dx));
	FilmMatState st;
	st.x = Fibonacci_Sphere(n, R);
	st.v.assign(n, Vec3::Zero());
	{
		const ThinFilmGraph G = cs.Graph(st);
		st.V.resize(n); st.S.resize(n);
		for (int i = 0; i < n; i++) { st.V[i] = h0 * G.a[i]; st.S[i] = G0 * G.a[i]; }
	}
	std::printf("\n== Bubble with mobile surfaces (R %g, %d particles; h0 %g, rho g %g, mu %g, rho %g; Langmuir RT %g, Gamma0 %g: E_G %.1f; Ds %g; dt %g) ==\n",
		R, n, h0, rho_g, prm.tf.mu, prm.rho, prm.surf.RT, G0, 2. * G0 * G0 * prm.surf.d2g(G0), cs.diff.Ds, dt);
	std::printf("   weight per surface rho g h0 R = %.3f against Pi_s0 = %.2f; apex law (immobile) h = (h0^-2 + rho g t / (3 mu R))^-1/2\n", rho_g * h0 * R, prm.surf.Pressure(G0));
	std::printf("  %7s %9s %9s %9s %9s %9s %9s %10s %9s %9s %6s %7s\n", "t", "h apex", "apex law", "black", "h bottom", "G apex", "G bottom", "max|v|", "vol drft", "surf drft", "N drn", "sec");
	real V0 = 0., S0 = 0.;
	for (int i = 0; i < st.Size(); i++) { V0 += st.V[i]; S0 += st.S[i]; }
	const int nsteps = (int)std::lround(tmax / dt);
	const std::vector<real> report_t = { 10., 50., 100., 250., 500., 1000., 2000., 4000. };
	size_t next = 0;
	int it_d = 0, nrep = 0;
	const auto t0 = std::chrono::steady_clock::now();
	for (int s = 1; s <= nsteps; s++) {
		const FilmCoupled_Stepper::Stats ss = cs.Step(st, dt);
		it_d += ss.drain.iters; nrep++;
		if (!ss.mech.converged || !ss.drain.converged || !ss.diff.converged)
			std::printf("  (step %d: Newton not converged: mech %d drain %d diff %d)\n", s, (int)ss.mech.converged, (int)ss.drain.converged, (int)ss.diff.converged);
		const real t = s * dt;
		if (next < report_t.size() && std::abs(t - report_t[next]) < 0.5 * dt) {
			next++;
			const ThinFilmGraph G = cs.Graph(st);
			real va = 0., sa = 0., aa = 0., vb = 0., sb = 0., ab = 0., black = 0., Vt = 0., St = 0., vmax = 0.;
			for (int i = 0; i < n; i++) {
				const real th = std::acos(std::clamp(st.x[i][2] / R, -1., 1.));
				if (R * th < cap) { va += st.V[i]; sa += st.S[i]; aa += G.a[i]; }
				if (R * (kPi - th) < cap) { vb += st.V[i]; sb += st.S[i]; ab += G.a[i]; }
				if (st.V[i] / G.a[i] < 2.) black += G.a[i];
				Vt += st.V[i]; St += st.S[i];
				vmax = std::max(vmax, st.v[i].norm());
			}
			const real law = 1. / std::sqrt(1. / (h0 * h0) + rho_g * t / (3. * prm.tf.mu * R));
			const real arc = R * std::acos(std::clamp(1. - black / (2. * kPi * R * R), -1., 1.));
			std::printf("  %7.0f %9.4f %9.4f %9.2f %9.3f %9.4f %9.4f %10.2e %9.1e %9.1e %6.1f %7.1f\n", t, va / aa, law, arc, vb / ab, sa / aa, sb / ab, vmax,
				std::abs(Vt / V0 - 1.), std::abs(St / S0 - 1.), (double)it_d / nrep, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
			it_d = nrep = 0;
		}
	}
}

// Marginal regeneration at the bottom border of a vertical film (mobile surfaces). The border (a Plateau border /
// meniscus) is a reservoir: the frame rows within one kernel radius of the bottom keep the thickness h_m and the
// concentration Gamma0 (liquid and surfactant exchanged with the bath, counted). Where the film meets the thick border
// the capillary (squeeze) term lowers the pressure, sucks liquid and pinches the film: thinner, lighter film under
// thicker film, which with mobile surfaces is unstable to 2D buoyancy (thin film rises in plumes, thick film sinks).
// Periodic in x, closed at the top (fixed rows). Reports the pinch (thinnest band above the border), the across-x
// variation of h in bands (plumes), the largest horizontal surface speed, the black area, the exchange.
static void Run_Regen(real dx)
{
	const real W = Env("VF_WIDTH", 64.), H = Env("VF_HEIGHT", 64.), h0 = Env("VF_H0", 6.), hm = Env("VF_HM", 16.), rho_g = Env("VF_RHOG", 0.01);
	const real dt = Env("VF_DT", 2.), tmax = Env("VF_TMAX", 1000.), G0 = Env("VF_G0", 0.5), rdx = Env("VF_RDX", 4.), noise = Env("VF_NOISE", 0.01);
	FilmCoupled_Stepper cs;
	FilmMechParams& prm = cs.mech.prm;
	prm.tf.sigma = Env("VF_SIGMA", 1.); prm.tf.mu = Env("VF_MU", 1.); prm.tf.rho_g = rho_g; prm.tf.disj = Disj_From_Env(1., 1., 1.5);
	prm.surf.RT = Env("VF_RT", 20.); prm.surf.Gamma_inf = 1.;
	prm.rho = Env("VF_RHO", 1e-3);
	prm.up = Vec3(0., 1., 0.);
	cs.diff.Ds = Env("VF_DS", 0.);
	cs.drain.newton.max_cg_iters = cs.diff.newton.max_cg_iters = 20000;
	const real R = rdx * dx, dy = dx * std::sqrt(3.) / 2.;
	cs.mech.kernel.Set_Radius(R);
	const int nx = (int)std::lround(W / dx), ny = (int)std::lround(H / dy);
	cs.mech.box.Lx = nx * dx; cs.mech.box.Ly = ny * dy; cs.mech.box.periodic_y = false;
	const real Lx = nx * dx, Ly = ny * dy;
	FilmMatState st;
	for (int rr = 0; rr < ny; rr++) for (int c = 0; c < nx; c++) st.x.push_back(Vec3((c + 0.5 * (rr & 1)) * dx, (rr + 0.5) * dy, 0.));
	const int n = st.Size();
	st.v.assign(n, Vec3::Zero()); st.fixed.assign(n, 0);
	std::vector<char> res(n, 0);
	for (int i = 0; i < n; i++) {
		if (st.x[i][1] < R) { st.fixed[i] = 1; res[i] = 1; }
		if (st.x[i][1] > Ly - R) st.fixed[i] = 1;
	}
	// the border's thickness profile: a parabola h = h_e + kappa (R - y)^2 / 2 rising from the film to h_m at the bottom:
	// its curvature kappa = 1 / r_m gives the suction -sigma kappa / 4 of the meniscus (marginal pinching)
	const real he = Env("VF_HE", h0), kappa = 2. * (hm - he) / (R * R);
	auto hres = [&](real y) { return he + 0.5 * kappa * (R - y) * (R - y); };
	std::mt19937 rng(11);
	std::normal_distribution<real> N01;
	{
		const ThinFilmGraph G = cs.Graph(st);
		st.V.resize(n); st.S.resize(n);
		for (int i = 0; i < n; i++) { st.V[i] = (res[i] ? hres(st.x[i][1]) : h0 * (1. + noise * N01(rng))) * G.a[i]; st.S[i] = G0 * G.a[i]; }
	}
	std::printf("\n== Marginal regeneration: vertical film %g x %g (dx %g, %d particles), bottom border = reservoir rising to h_m %g within %g (suction %.3f);\n"
		"   h0 %g (noise %g), rho g %g, mu %g, rho %g, sigma %g; Langmuir RT %g, Gamma0 %g (E_G %.2f); Ds %g; dt %g ==\n",
		Lx, Ly, dx, n, hm, R, -prm.tf.sigma * kappa / 4., h0, noise, rho_g, prm.tf.mu, prm.rho, prm.tf.sigma, prm.surf.RT, G0, 2. * G0 * G0 * prm.surf.d2g(G0), cs.diff.Ds, dt);
	std::printf("  %7s %8s %8s %10s %10s %10s %10s %8s %10s %6s %6s %7s\n", "t", "pinch h", "at y", "x-var pin", "x-var 3R", "x-var mid", "max|v_x|", "black%", "bath dV",
		"N mech", "N drn", "sec");
	const char* frames = std::getenv("VF_FRAMES");
	FILE* ff = frames ? std::fopen(frames, "w") : nullptr;
	const real every = Env("VF_EVERY", 0.), report = Env("VF_REPORT", 50.);
	const int nsteps = (int)std::lround(tmax / dt);
	real bath = 0.;
	int it_m = 0, it_d = 0, nrep = 0;
	// VF_KT: thermal fluctuations of the thickness after each step (random liquid fluxes along the pairs, variance
	// 2 kT k_q dt: the noise of the drainage step's dissipation; see bubble)
	const real kT = Env("VF_KT", 0.);
	std::mt19937 trng(23);
	std::normal_distribution<real> TN;
	std::vector<real> kq;
	if (kT > 0.) std::printf("   thermal fluctuations kT %g\n", kT);
	const auto t0 = std::chrono::steady_clock::now();
	for (int s = 1; s <= nsteps; s++) {
		const FilmCoupled_Stepper::Stats ss = cs.Step(st, dt);
		it_m += ss.mech.iters; it_d += ss.drain.iters; nrep++;
		if (!ss.mech.converged || !ss.drain.converged || !ss.diff.converged)
			std::printf("  (step %d: Newton not converged: mech %d drain %d diff %d)\n", s, (int)ss.mech.converged, (int)ss.drain.converged, (int)ss.diff.converged);
		const ThinFilmGraph G = cs.Graph(st);
		// thermal kick for the next step (relaxed by its implicit drainage, see bubble)
		if (kT > 0.) {
			cs.drain.Conductances(G, st.V, kq);
			for (int q = 0; q < G.pairs.Size(); q++) {
				const int i = G.pairs.pi[q], j = G.pairs.pj[q];
				const real F = std::clamp(std::sqrt(2. * kT * kq[q] * dt) * TN(trng), -0.5 * st.V[j], 0.5 * st.V[i]);
				st.V[i] -= F; st.V[j] += F;
			}
		}
		for (int i = 0; i < n; i++) if (res[i]) { const real hr = hres(st.x[i][1]); bath += st.V[i] - hr * G.a[i]; st.V[i] = hr * G.a[i]; st.S[i] = G0 * G.a[i]; }
		const real t = s * dt;
		if (ff && every > 0. && std::abs(t / every - std::round(t / every)) < 0.5 * dt / every) {
			std::fprintf(ff, "t %g\n", t);
			for (int i = 0; i < n; i++) std::fprintf(ff, "%.4f %.4f %.5f %.5f\n", st.x[i][0], st.x[i][1], st.V[i] / G.a[i], st.S[i] / G.a[i]);
		}
		if (std::abs(t / report - std::round(t / report)) < 0.5 * dt / report) {
			// bands of height 2 above the border: mean h, and the across-x variation (std over 8 x-bins / mean)
			const int nb = (int)std::ceil(Ly / 2.), nxb = 8;
			std::vector<real> sv(nb * nxb, 0.), sa(nb * nxb, 0.);
			real black = 0., area = 0., vx = 0.;
			for (int i = 0; i < n; i++) {
				if (st.fixed[i]) continue;
				const int b = std::clamp((int)(st.x[i][1] / 2.), 0, nb - 1), c = std::clamp((int)(st.x[i][0] / Lx * nxb), 0, nxb - 1);
				sv[b * nxb + c] += st.V[i]; sa[b * nxb + c] += G.a[i];
				if (st.V[i] / G.a[i] < 2.) black += G.a[i];
				area += G.a[i];
				vx = std::max(vx, std::abs(st.v[i][0]));
			}
			auto band = [&](int b, real& mean, real& var) {
				real m = 0., m2 = 0.; int k = 0;
				for (int c = 0; c < nxb; c++) if (sa[b * nxb + c] > 0.) { const real h = sv[b * nxb + c] / sa[b * nxb + c]; m += h; m2 += h * h; k++; }
				if (!k) { mean = 0.; var = 0.; return; }
				m /= k; mean = m; var = std::sqrt(std::max<real>(0., m2 / k - m * m)) / m;
			};
			int bp = -1; real hp = 1e300, vp = 0.;
			for (int b = (int)std::ceil(R / 2.); b < nb / 2; b++) { real m, v; band(b, m, v); if (m > 0. && m < hp) { hp = m; bp = b; vp = v; } }
			real m3, v3, mm, vm;
			band(std::min(nb - 1, (int)(3. * R / 2.)), m3, v3);
			band(nb / 2, mm, vm);
			std::printf("  %7.0f %8.3f %8.1f %10.2e %10.2e %10.2e %10.2e %8.2f %+10.3e %6.1f %6.1f %7.1f\n", t, hp, 2. * bp + 1., vp, v3, vm, vx, 100. * black / area, bath,
				(double)it_m / nrep, (double)it_d / nrep, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
			it_m = it_d = nrep = 0;
		}
	}
	if (ff) std::fclose(ff);
}

// Marginal regeneration with an open meniscus: the bottom border exchanges surface with a surfactant bath. The film's
// energy carries the bath term 2 sum_i Pi_b(y_i) a_i (Pi_b = Pi_m = Pi_s(Gamma0) at the open edge: a film at the bath
// concentration is in equilibrium with it; lower by delta in the film, see below; the stretched top carries the
// weight). The border band (free particles with y < R) is weightless (the bath carries the meniscus) and is reset every
// step to the bath concentration; for the mechanics it has the thickness of the film above it, for the drainage the
// meniscus profile (parabola rising to h_m at y = 0: the suction -sigma kappa / 4). Particles that move below y = 0 are
// absorbed by the bath; where the film pulls away from the bath, new particles are injected at the bottom (one per
// column of width dx whose lowest particle is above 0.9 dy); the film is remeshed (insertion, merging, conservative
// remap; nothing inserted below the bath surface, nothing absorbed by the frame). The top frame is 2R thick and a
// mechanical continuation of the film (ghost_frame = 2). Liquid and surfactant exchanged with the bath are counted.
// Why each of these choices: thin_film/README.md, "The open meniscus".
static void Run_Regen_Open(real dx)
{
	const real W = Env("VF_WIDTH", 64.), H = Env("VF_HEIGHT", 64.), h0 = Env("VF_H0", 6.), hm = Env("VF_HM", 16.), rho_g = Env("VF_RHOG", 0.01);
	const real dt = Env("VF_DT", 2.), tmax = Env("VF_TMAX", 1000.), G0 = Env("VF_G0", 0.5), rdx = Env("VF_RDX", 4.), noise = Env("VF_NOISE", 0.01);
	FilmCoupled_Stepper cs;
	FilmMechParams& prm = cs.mech.prm;
	prm.tf.sigma = Env("VF_SIGMA", 1.); prm.tf.mu = Env("VF_MU", 1.); prm.tf.rho_g = rho_g; prm.tf.disj = Disj_From_Env(1., 1., 1.5);
	// VF_MU: the liquid viscosity in the drainage (Poiseuille mobility h^3 / (12 mu)); VF_MU_S (default VF_MU): the one in
	// the in-plane (Trouton) viscosity 4 mu h of the surfaces' motion. In lubrication units both are one mu only if the
	// thickness and lateral units are the same; for a real film (h / L ~ 1e-3) the in-plane flow is faster than the
	// drainage by ~ (L / h)^2, and VF_MU_S << VF_MU restores part of that separation (marginal regeneration lives in it)
	const real mu_drain = prm.tf.mu;
	prm.mu_s = Env("VF_MU_S", mu_drain); // FilmMechParams::mu_s (Units.h derives it from H / L)
	prm.surf.RT = Env("VF_RT", 20.); prm.surf.Gamma_inf = 1.;
	prm.rho = Env("VF_RHO", 1e-3);
	prm.up = Vec3(0., 1., 0.);
	// the bath balances the film's whole in-plane stress at the reference state (h0, Gamma0): the surface pressure of
	// both surfaces and the disjoining part h P' - P (else the edge drifts into the bath), as the energy
	// 2 sum_i Pi_b(y_i) a_i with Pi_b = Pi_m in the band and Pi_m - delta in the film (VF_BATH_DELTA), joined by a
	// smoothstep over VF_BATH_W (frozen factor in bath_w, body force -2 Pi_b'(y_i) a_i in fext). In the continuum Pi_b(y)
	// drops out of the bulk exactly (its stress gradient cancels its body force): only Pi_b = Pi_m at the open edge
	// counts (zero stress there at Gamma0). On the particles it sets the stress they carry: 2 (Pi_s - Pi_b).
	// - delta = 0: the excess stress 2 (Pi_s - Pi_m), a tension wherever the film carries its weight (up to ~4.5 under
	//   the frame). Under tension the lattice there churns (insertions and merges ~35 per step each), the
	//   concentration scatters by +-0.05 between neighbours (a tension noise as large as the tension itself), and when
	//   the drained top reaches the disjoining barrier (h ~ 3.4) the black spots that nucleate tear it (t ~ 330);
	// - delta = Pi_m (the closed films' positive pressure 2 Pi_s ~ 28 everywhere above the band): the stress gradient
	//   and the body force cancel only to the kernel's accuracy, and a few percent of 2 Pi_m / w ~ 5 per area is more
	//   than the film's weight (0.06): the first step moves the film by ~10;
	// - delta ~ 3 (default) over w = 4R: the particles keep a positive pressure 2 delta - T > 0 everywhere, and the body
	//   force is 2 delta / w ~ 0.4 per area.
	const real Pi_m = prm.surf.Pressure(G0) + 0.5 * (h0 * prm.tf.disj.dP(h0) - prm.tf.disj.P(h0));
	prm.bath_pressure = Pi_m;
	const real bath_delta = Env("VF_BATH_DELTA", 3.);
	cs.diff.Ds = Env("VF_DS", 0.);
	cs.drain.newton.max_cg_iters = cs.diff.newton.max_cg_iters = 20000;
	cs.mech.newton.verbose = (int)Env("VF_VERBOSE", 0.);
	// the top frame as a wire (VF_GHOST, see FilmMech_Stepper::ghost_frame): without it the frame particles' own state
	// drifts (they drain to black, their concentration falls to 0.3 against the film's 0.45) and the frame pulls the
	// film into it; with the film's thickness and concentration (1) the frozen frame lattice, denser than the stretched
	// film below, pulls less than the film and the top rows come off; 2: the film's stress coefficient
	cs.mech.ghost_frame = (int)Env("VF_GHOST", 2.);
	cs.mech.max_halvings = (int)Env("VF_HALVINGS", 3.); // fast plumes: redo a step in halves rather than lose pairs
	const bool band_ghost = Env("VF_BAND_GHOST", 1.) != 0.;
	const real R = rdx * dx, dy = dx * std::sqrt(3.) / 2.;
	const real bath_w = Env("VF_BATH_W", 4. * R);
	cs.mech.kernel.Set_Radius(R);
	// VF_FRAME: the top frame's thickness in kernel radii. With 1 the frame particles next to the film are also next to
	// the frame's outer edge: their kernel areas are too large, and a pair's force (a_i^2 T_i + a_j^2 T_j) |W'| (T: the
	// film tension) pulls the film into the frame even at equal tension. With 2 those particles have full neighbourhoods;
	// the extra band is added above (the film keeps its height)
	const real frame = Env("VF_FRAME", 2.);
	const int nx = (int)std::lround(W / dx), ny = (int)std::lround((H + (frame - 1.) * R) / dy);
	cs.mech.box.Lx = nx * dx; cs.mech.box.Ly = ny * dy; cs.mech.box.periodic_y = false;
	const real Lx = nx * dx, Ly = ny * dy, Lmid = 0.5 * (Ly - (frame - 1.) * R);
	FilmMatState st;
	for (int rr = 0; rr < ny; rr++) for (int c = 0; c < nx; c++) st.x.push_back(Vec3((c + 0.5 * (rr & 1)) * dx, (rr + 0.5) * dy, 0.));
	int n = st.Size();
	st.v.assign(n, Vec3::Zero()); st.fixed.assign(n, 0); st.gscale.assign(n, 1.);
	for (int i = 0; i < n; i++) if (st.x[i][1] > Ly - frame * R) st.fixed[i] = 1; // the top frame
	const real he = Env("VF_HE", h0), kappa = 2. * (hm - he) / (R * R);
	auto hres = [&](real y) { const real d = std::max<real>(0., R - y); return he + 0.5 * kappa * d * d; };
	auto in_band = [&](int i) { return !st.fixed[i] && st.x[i][1] < R; };
	const bool reservoir_bath = Env("VF_RESERVOIR", 1.) != 0.;
	const real mu_res = -prm.tf.sigma * kappa / 4. + rho_g * R;
	std::mt19937 rng(11);
	std::normal_distribution<real> N01;
	{
		const ThinFilmGraph G = cs.Graph(st);
		st.V.resize(n); st.S.resize(n);
		// VF_INIT_EQ=1: the concentration graded so that the film starts in Marangoni equilibrium (the tension at height y
		// carries the weight below it: Pi_s(Gamma(y)) = Pi_s(Gamma0) - rho_g h0 (y - R) / 2); 0: uniform Gamma0, the film
		// unsupported at first (with a small in-plane viscosity it falls freely in the first step and crushes)
		const bool init_eq = Env("VF_INIT_EQ", 0.) != 0.;
		for (int i = 0; i < n; i++) {
			st.V[i] = (in_band(i) ? hres(st.x[i][1]) : h0 * (1. + noise * N01(rng))) * G.a[i];
			const real yy = std::max<real>(0., std::min(st.x[i][1], Ly - frame * R) - R);
			st.S[i] = (init_eq ? prm.surf.Concentration(prm.surf.Pressure(G0) - 0.5 * rho_g * h0 * yy) : G0) * G.a[i];
		}
	}
	std::printf("\n== Marginal regeneration, open meniscus: vertical film %g x %g (dx %g, %d particles); bath Pi_m %.3f (Gamma0 %g), meniscus rising to h_m %g\n"
		"   within %g (suction %.3f); h0 %g (noise %g), rho g %g, mu %g (in-plane %g), rho %g; Langmuir RT %g (E_G %.2f); Ds %g; dt %g ==\n",
		Lx, Ly, dx, n, Pi_m, G0, hm, R, -prm.tf.sigma * kappa / 4., h0, noise, rho_g, mu_drain, prm.Mu_S(), prm.rho, prm.surf.RT, 2. * G0 * G0 * prm.surf.d2g(G0),
		cs.diff.Ds, dt);
	const real kT = Env("VF_KT", 0.);
	std::mt19937 trng(23);
	std::normal_distribution<real> TN;
	std::vector<real> kq;
	if (kT > 0.) std::printf("   thermal fluctuations kT %g\n", kT);
	// x-var: the thickness's variation along x (std over 8 x-bins / mean) in a band of height 2: at the pinch, at
	// y = 3R (where plumes rising from the pinch show), and at mid-height
	std::printf("  %7s %7s %8s %8s %10s %10s %10s %10s %10s %8s %10s %10s %7s %7s %7s\n", "t", "N", "pinch h", "at y", "x-var pin", "x-var 3R", "x-var mid", "max|v_x|",
		"v_y mid", "black%", "bath dV", "bath dS", "in/out", "ins/mrg", "sec");
	const char* frames = std::getenv("VF_FRAMES");
	FILE* ff = frames ? std::fopen(frames, "w") : nullptr;
	const real every = Env("VF_EVERY", 0.), report = Env("VF_REPORT", 50.);
	const int nsteps = (int)std::lround(tmax / dt);
	real bathV = 0., bathS = 0.;
	int injected = 0, absorbed = 0, inserted = 0, merged = 0;
	RemeshParams rp; rp.dx0 = dx; rp.kernel_radius = R; rp.absorb_fixed = Env("VF_ABSORB_FRAME", 0.) != 0.;
	rp.insert_ok = [&](const Vec3& p) { return p[1] > 0.5 * dy; }; // nothing below the bath surface (it would be absorbed next step)
	auto graph_of = [&](const std::vector<Vec3>& xs) { FilmMatState tmp; tmp.x = xs; return cs.Graph(tmp); };
	auto remove_if = [&](const std::function<bool(int)>& drop) {
		FilmMatState out; out.dt_prev = st.dt_prev; out.has_prev = st.has_prev;
		const bool hist = st.has_prev && (int)st.x_prev.size() == st.Size();
		for (int i = 0; i < st.Size(); i++) {
			if (drop(i)) continue;
			out.x.push_back(st.x[i]); out.v.push_back(st.v[i]); out.V.push_back(st.V[i]); out.S.push_back(st.S[i]);
			out.fixed.push_back(st.fixed[i]); out.gscale.push_back(st.gscale[i]);
			if (hist) { out.x_prev.push_back(st.x_prev[i]); out.v_prev.push_back(st.v_prev[i]); }
		}
		st = std::move(out);
	};
	// diagnostics: VF_SAVE=file + VF_SAVE_AT=step (checkpoint after that step), VF_LOAD=file (resume from a checkpoint),
	// VF_STEPLOG=1 (one line per step), VF_DUMP=file + VF_DUMP_STEP=step (per-particle state before and after its mechanics)
	const char* save_file = std::getenv("VF_SAVE");
	const char* load_file = std::getenv("VF_LOAD");
	const char* dump_file = std::getenv("VF_DUMP");
	const int save_at = (int)Env("VF_SAVE_AT", 0.), dump_step = (int)Env("VF_DUMP_STEP", 0.), steplog = (int)Env("VF_STEPLOG", 0.);
	int s_first = 1;
	if (load_file) {
		std::vector<real> c;
		if (!Load_State(load_file, st, c) || c.size() < 7) { std::printf("  (cannot load %s)\n", load_file); return; }
		s_first = (int)c[0] + 1; bathV = c[1]; bathS = c[2]; injected = (int)c[3]; absorbed = (int)c[4]; inserted = (int)c[5]; merged = (int)c[6];
		std::printf("  (resumed from %s after step %d, %d particles)\n", load_file, s_first - 1, st.Size());
	}
	auto dump = [&](FILE* f, const char* tag) {
		const ThinFilmGraph Gd = cs.Graph(st);
		std::fprintf(f, "%s %d\n", tag, st.Size());
		for (int i = 0; i < st.Size(); i++)
			std::fprintf(f, "%.6f %.6f %.6e %.6e %.6e %.6e %.6e %d %g\n", st.x[i][0], st.x[i][1], st.v[i][0], st.v[i][1], st.V[i], st.S[i], Gd.a[i], (int)st.fixed[i],
				st.gscale.empty() ? 1. : st.gscale[i]);
	};
	const auto t0 = std::chrono::steady_clock::now();
	for (int s = s_first; s <= nsteps; s++) {
		n = st.Size();
		const std::vector<Vec3> x_step0 = st.x;
		FILE* fd = (dump_file && s == dump_step) ? std::fopen(dump_file, "w") : nullptr;
		// the band for the mechanics: weightless, bath concentration, and the thickness of the film just above it (mean over
		// x-bins of width 2 of the film in R <= y < R + 2; h0 where a bin is empty). The meniscus profile belongs to the bath,
		// not to the particles sliding over it: with h0 in the band, the squeeze energy (sigma / 8) c (h_i - h_j)^2 across
		// the band's edge (film there ~8, thicker than h0 from the liquid gathering at the meniscus) compressed the band and
		// stretched the film's lowest rows; the reset erased the band's compression every step, and the film sank into it
		{
			const ThinFilmGraph G = cs.Graph(st);
			const int nbx = std::max(1, (int)std::lround(Lx / 2.));
			std::vector<real> sv(nbx, 0.), sa(nbx, 0.);
			for (int i = 0; i < n; i++) {
				if (st.fixed[i] || st.x[i][1] < R || st.x[i][1] >= R + 2.) continue;
				const int b = std::clamp((int)(st.x[i][0] / Lx * nbx), 0, nbx - 1);
				sv[b] += st.V[i]; sa[b] += G.a[i];
			}
			for (int i = 0; i < n; i++) {
				st.gscale[i] = in_band(i) ? 0. : 1.;
				if (!in_band(i)) continue;
				const int b = std::clamp((int)(st.x[i][0] / Lx * nbx), 0, nbx - 1);
				const real hb = (band_ghost && sa[b] > 0.) ? sv[b] / sa[b] : h0;
				st.V[i] = hb * G.a[i]; st.S[i] = G0 * G.a[i];
			}
			// the bath coefficient Pi_b = Pi_m - delta (1 - chi(y)), chi = 1 up to y = R, smoothstep down to 0 at R + bath_w,
			// and its body force -2 Pi_b' a_i = -2 delta chi' a_i
			st.fext.assign(n, Vec3::Zero());
			st.bath_w.assign(n, 1.);
			if (bath_delta > 0. && bath_w > 0.) for (int i = 0; i < n; i++) {
				const real s = std::clamp((st.x[i][1] - R) / bath_w, 0., 1.);
				const real chi = 1. - s * s * (3. - 2. * s), dchi = (s > 0. && s < 1.) ? -6. * s * (1. - s) / bath_w : 0.;
				st.bath_w[i] = 1. - bath_delta / Pi_m * (1. - chi);
				st.fext[i] = Vec3(0., -2. * bath_delta * dchi * G.a[i], 0.);
			}
		}
		if (fd) dump(fd, "before");
		const NewtonStats nm = cs.mech.Step(st, dt);
		if (fd) { dump(fd, "after"); std::fclose(fd); }
		if (steplog) {
			real dmax = 0., vmax = 0., ke = 0., vy = 0., wy = 0.; int idm = -1, ivm = -1;
			for (int i = 0; i < n; i++) {
				if (st.fixed[i]) continue;
				const real d = cs.mech.box.Wrap(st.x[i] - x_step0[i]).norm(), v = st.v[i].norm();
				if (d > dmax) { dmax = d; idm = i; }
				if (v > vmax) { vmax = v; ivm = i; }
				ke += 0.5 * prm.rho * st.V[i] * v * v;
				if (!in_band(i)) { vy += st.V[i] * st.v[i][1]; wy += st.V[i]; }
			}
			std::printf("  step %4d t %6.0f N %5d  mech %2d it %s |g| %.1e  max move %.3f at (%.1f,%.1f)  max|v| %.3f at (%.1f,%.1f)  KE %.3e  <v_y> %+.4f\n", s, s * dt, n,
				nm.iters, nm.converged ? "ok" : "NC", nm.grad_inf_history.empty() ? 0. : nm.grad_inf_history.back(), dmax, idm >= 0 ? st.x[idm][0] : 0., idm >= 0 ? st.x[idm][1] : 0.,
				vmax, ivm >= 0 ? st.x[ivm][0] : 0., ivm >= 0 ? st.x[ivm][1] : 0., ke, wy > 0. ? vy / wy : 0.);
		}
		// the band for the drainage: the meniscus profile; the exchange with the bath is what the drainage moves into it
		const ThinFilmGraph G = cs.Graph(st);
		for (int i = 0; i < n; i++) if (in_band(i)) { st.V[i] = hres(st.x[i][1]) * G.a[i]; st.S[i] = G0 * G.a[i]; }
		real bandV0 = 0., bandS0 = 0.;
		for (int i = 0; i < n; i++) if (in_band(i)) { bandV0 += st.V[i]; bandS0 += st.S[i]; }
		cs.drain.prm = prm.tf;
		// VF_RESERVOIR=1 (default): the band is a bath at the meniscus's chemical potential -sigma kappa / 4 + rho_g R (its
		// capillary suction at the junction height); its thickness profile only sets the pair mobilities. Without it the
		// band's pressure comes from its own thickness profile, and its lowest rows (neighbours only above, all thinner)
		// get a large positive squeeze pressure instead of the suction: with the band reset every step, a liquid source
		// that thickens the film above the meniscus (h ~ 17 at y = R for h_m = 40)
		std::vector<real> resv;
		if (reservoir_bath) {
			resv.assign(n, std::numeric_limits<real>::quiet_NaN());
			for (int i = 0; i < n; i++) if (in_band(i)) resv[i] = mu_res;
			cs.drain.reservoir = &resv;
		}
		NewtonStats nd;
		{
			// a step that does not converge (the film next to the meniscus pulled into the black film's wall) is redone in
			// 4, then 16 substeps (fresh conductances each)
			const std::vector<real> V_start = st.V;
			nd = cs.drain.Step(G, st.V, dt);
			for (int sub = 4; !nd.converged && sub <= 16; sub *= 4) {
				st.V = V_start;
				bool ok = true;
				for (int k = 0; k < sub && ok; k++) { nd = cs.drain.Step(G, st.V, dt / sub); ok = nd.converged; }
			}
		}
		cs.drain.reservoir = nullptr;
		cs.diff.surf = prm.surf;
		const NewtonStats ns = cs.diff.Step(G, st.S, dt);
		real bandV1 = 0., bandS1 = 0.;
		for (int i = 0; i < n; i++) if (in_band(i)) { bandV1 += st.V[i]; bandS1 += st.S[i]; }
		bathV += bandV1 - bandV0; bathS += bandS1 - bandS0;
		if (!nm.converged || !nd.converged || !ns.converged) {
			real hmn = 1e300, hmx = 0., gmn = 1e300, gmx = 0., rmin = 1e300; real yh = -1., yg = -1.;
			for (int i = 0; i < n; i++) {
				if (st.fixed[i]) continue;
				const real h = st.V[i] / G.a[i], g = st.S[i] / G.a[i];
				if (h < hmn) { hmn = h; yh = st.x[i][1]; }
				hmx = std::max(hmx, h);
				if (g > gmx) { gmx = g; yg = st.x[i][1]; }
				gmn = std::min(gmn, g);
			}
			for (int q = 0; q < G.pairs.Size(); q++) rmin = std::min(rmin, cs.mech.box.Wrap(st.x[G.pairs.pj[q]] - st.x[G.pairs.pi[q]]).norm());
			std::printf("  (step %d: Newton not converged: mech %d drain %d diff %d; h %.3f-%.3f (min at y %.1f), Gamma %.4f-%.4f (max at y %.1f), closest pair %.3f;"
				" mech iters %d, |g| %.2e -> %.2e, line-search fails %s)\n", s, (int)nm.converged, (int)nd.converged, (int)ns.converged, hmn, hmx, yh, gmn, gmx, yg, rmin,
				nm.iters, nm.grad_inf_history.empty() ? 0. : nm.grad_inf_history.front(), nm.grad_inf_history.empty() ? 0. : nm.grad_inf_history.back(),
				nm.iters < cs.mech.newton.max_iters ? "yes" : "no (max iters)");
		}
		// absorb what went below the bath surface
		int na = 0;
		for (int i = 0; i < n; i++) if (!st.fixed[i] && st.x[i][1] < 0.) { bathV += st.V[i]; bathS += st.S[i]; na++; }
		if (na) { remove_if([&](int i) { return !st.fixed[i] && st.x[i][1] < 0.; }); absorbed += na; }
		// VF_KT: thermal fluctuations of the thickness (random liquid fluxes along the pairs, variance 2 kT k_q dt; as in
		// regen), relaxed by the next step's implicit drainage. Only between film particles: the band is reset to the
		// meniscus profile every step, so what a kick moves out of it stays in the film and what it moves in is erased,
		// and its conductances (k ~ h^3, h up to h_m) are huge: kicks across the band pump the bath's liquid into the film
		if (kT > 0.) {
			cs.drain.Conductances(G, st.V, kq);
			for (int q = 0; q < G.pairs.Size(); q++) {
				const int i = G.pairs.pi[q], j = G.pairs.pj[q];
				if (in_band(i) || in_band(j)) continue;
				const real F = std::clamp(std::sqrt(2. * kT * kq[q] * dt) * TN(trng), -0.5 * st.V[j], 0.5 * st.V[i]);
				st.V[i] -= F; st.V[j] += F;
			}
		}
		// inject where the film pulled away from the bath
		{
			std::vector<real> ymin(nx, 1e300);
			std::vector<int> imin(nx, -1);
			for (int i = 0; i < st.Size(); i++) {
				const int c = std::clamp((int)std::floor(st.x[i][0] / dx), 0, nx - 1);
				if (st.x[i][1] < ymin[c]) { ymin[c] = st.x[i][1]; imin[c] = i; }
			}
			const bool hist = st.has_prev && (int)st.x_prev.size() == st.Size();
			for (int c = 0; c < nx; c++) {
				if (ymin[c] <= 0.9 * dy) continue;
				const Vec3 x((c + 0.5) * dx, 0.45 * dy, 0.);
				const Vec3 v = imin[c] >= 0 ? st.v[imin[c]] : Vec3::Zero();
				const real cell = dx * dy;
				st.x.push_back(x); st.v.push_back(v); st.V.push_back(hres(x[1]) * cell); st.S.push_back(G0 * cell);
				st.fixed.push_back(0); st.gscale.push_back(0.);
				if (hist) { st.x_prev.push_back(x - dt * v); st.v_prev.push_back(v); }
				bathV -= hres(x[1]) * cell; bathS -= G0 * cell;
				injected++;
			}
		}
		// remesh the film (the remap keeps h and Gamma; the band is reset next step anyway)
		{
			const ThinFilmGraph Gr = cs.Graph(st);
			const RemeshStats rs = Remesh_Film(st, Gr, rp, &cs.mech.box, nullptr, graph_of);
			inserted += rs.splits; merged += rs.merges;
		}
		n = st.Size();
		if (save_file && s == save_at) {
			Save_State(save_file, st, { (real)s, bathV, bathS, (real)injected, (real)absorbed, (real)inserted, (real)merged });
			std::printf("  (checkpoint after step %d: %s)\n", s, save_file);
		}
		const real t = s * dt;
		if (ff && every > 0. && std::abs(t / every - std::round(t / every)) < 0.5 * dt / every) {
			const ThinFilmGraph Gf = cs.Graph(st);
			std::fprintf(ff, "t %g\n", t);
			for (int i = 0; i < n; i++) std::fprintf(ff, "%.4f %.4f %.5f %.5f\n", st.x[i][0], st.x[i][1], st.V[i] / Gf.a[i], st.S[i] / Gf.a[i]);
		}
		if (std::abs(t / report - std::round(t / report)) < 0.5 * dt / report) {
			const ThinFilmGraph Gp = cs.Graph(st);
			const int nb = (int)std::ceil(Lmid), nxb = 8;
			std::vector<real> sv(nb * nxb, 0.), sa(nb * nxb, 0.);
			real black = 0., area = 0., vx = 0., vym = 0.;
			int nym = 0;
			for (int i = 0; i < n; i++) {
				if (st.fixed[i] || in_band(i)) continue;
				const int b = std::clamp((int)(st.x[i][1] / 2.), 0, nb - 1), c = std::clamp((int)(st.x[i][0] / Lx * nxb), 0, nxb - 1);
				sv[b * nxb + c] += st.V[i]; sa[b * nxb + c] += Gp.a[i];
				if (st.V[i] / Gp.a[i] < 2.) black += Gp.a[i];
				area += Gp.a[i];
				vx = std::max(vx, std::abs(st.v[i][0]));
				if (std::abs(st.x[i][1] - Lmid) < 4.) { vym += st.v[i][1]; nym++; }
			}
			auto band = [&](int b, real& mean, real& var) {
				real m = 0., m2 = 0.; int k = 0;
				for (int c = 0; c < nxb; c++) if (sa[b * nxb + c] > 0.) { const real h = sv[b * nxb + c] / sa[b * nxb + c]; m += h; m2 += h * h; k++; }
				if (!k) { mean = 0.; var = 0.; return; }
				m /= k; mean = m; var = std::sqrt(std::max<real>(0., m2 / k - m * m)) / m;
			};
			int bp = -1; real hp = 1e300, vp = 0.;
			for (int b = (int)std::ceil(R / 2.); b < nb / 2; b++) { real m, v; band(b, m, v); if (m > 0. && m < hp) { hp = m; bp = b; vp = v; } }
			real mm, vm, m3, v3;
			band(nb / 2, mm, vm);
			band(std::min(nb - 1, (int)(3. * R / 2.)), m3, v3);
			std::printf("  %7.0f %7d %8.3f %8.1f %10.2e %10.2e %10.2e %10.2e %+10.2e %8.2f %+10.3e %+10.3e %3d/%-3d %3d/%-3d %7.1f\n", t, n, hp, 2. * bp + 1., vp, v3, vm, vx,
				nym ? vym / nym : 0., 100. * black / std::max<real>(area, 1e-300), bathV, bathS, injected, absorbed, inserted, merged,
				std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
		}
	}
	if (ff) std::fclose(ff);
}

// Longitudinal (elastic) waves of a film with mobile surfaces: displacement A sin(k x) on a periodic strip, released from
// rest. Linear theory (plug flow, Langmuir surfactant, disjoining, squeeze energy, Trouton viscosity):
//   rho h s^2 + 4 mu h k^2 s + (E_G + h^2 P''(h)) k^2 + (sigma h^2 / 4) k^4 = 0,   E_G = 2 RT Gamma / (1 - theta)
// i.e. a damped oscillation with rate b = 2 mu k^2 / rho and frequency omega_d = sqrt(omega_0^2 - b^2). Measured from the
// zero crossings (spacing pi / omega_d) and the ratio of successive extrema (exp(-b pi / omega_d)).
// Physical units: the material's code coefficients in given scales, and the same physical film (a flat periodic
// patch, h0 500 nm, thickness and surfactant modulations, gravity in the plane) stepped in two sets of scales, with
// the solver statistics and the difference of the physical results after each step (they must agree).
static void Run_Units()
{
	FilmMaterial M;
	M.theta0 = Env("VF_THETA0", M.theta0); M.c_salt = Env("VF_SALT", M.c_salt); M.psi0 = Env("VF_PSI0", M.psi0);
	std::printf("\n== Physical units: water + insoluble Langmuir surfactant (Gamma_inf %.2g mol/m^2, theta0 %g) + DLVO (A %.2g J, %g mM, psi0 %g mV, NBF %g nm) ==\n",
		M.Gamma_inf, M.theta0, M.hamaker, M.c_salt, 1e3 * M.psi0, 1e9 * M.h_nbf);
	std::printf("   gamma %.4f N/m, Pi_s %.4f N/m, Gibbs elasticity (film) %.4f N/m, Debye length %.2f nm, double-layer B %.3g J/m^2\n",
		M.gamma_ref(), M.Pi_s(M.theta0), M.Gibbs_E(), 1e9 * M.debye(), M.B());
	std::printf("   disjoining pressure Pi(h):");
	for (real hn : { 3., 4.5, 6., 10., 20., 30., 50., 100., 200. }) std::printf("  %g nm %.3g Pa", hn, M.Pi(hn * 1e-9));
	std::printf("\n");
	const real dx_m = Env("VF_DX_M", 20e-6), h0_m = Env("VF_H0_M", 500e-9), dt_s = Env("VF_DT_S", 1e-3);
	const int nx = (int)Env("VF_NX", 16.), nsteps = (int)Env("VF_STEPS", 5.);
	struct Out { std::vector<Vec3> d; std::vector<real> h; };
	auto run = [&](const FilmScales& sc, std::vector<Out>& outs) {
		const CodeParams cp = To_Code(M, sc);
		std::printf("   scales H %.3g m, L %.3g m, T %.3g s, S %.3g N/m: sigma %.4g, mu (drainage) %.4g, mu_s %.4g (ratio (H/L)^2 = %.3g), rho %.3g, rho_g %.4g,"
			" RT %.4g, kT %.3g; disjoining c8 %.3g c2 %.3g beta %.3g lambda %.3g\n", sc.H, sc.L, sc.T, sc.S, cp.tf.sigma, cp.tf.mu, cp.mu_s, cp.mu_s / cp.tf.mu,
			cp.rho, cp.tf.rho_g, cp.surf.RT, cp.kT, cp.tf.disj.c8, cp.tf.disj.c2, cp.tf.disj.beta, cp.tf.disj.lambda);
		FilmCoupled_Stepper cs;
		cp.Apply(cs.mech.prm);
		cs.mech.prm.up = Vec3(0., 1., 0.);
		cs.diff.Ds = cp.Ds;
		cs.mech.max_halvings = (int)Env("VF_HALVINGS", 0.);
		cs.mech.factored_precond = Env("VF_FPC", 1.) > 0.;
		cs.mech.precond_cholmod = Env("VF_CHOLMOD", 1.) > 0.;
		if (Env("VF_BATH", 1.) <= 0.) cs.mech.prm.bath_pressure = 0.; // VF_BATH=0: the full surface pressure on the particles (no reference)
		cs.mech.kernel.Set_Type((Kernel2D::Type)(int)Env("VF_KERNEL", 0.));
		const real dx = dx_m / sc.L;
		real Lx, Ly;
		std::vector<Vec3> x0 = Hex_Patch(nx, nx, dx, 0.05, 5, Lx, Ly);
		cs.mech.kernel.Set_Radius(Env("VF_RDX", 4.) * dx);
		cs.mech.box.Lx = Lx; cs.mech.box.Ly = Ly; cs.mech.box.periodic_y = true;
		cs.mech.rel_tol = cs.drain.rel_tol = cs.diff.rel_tol = Env("VF_RELTOL", 1e-12);
		FilmMatState st;
		st.x = x0;
		const int n = st.Size();
		st.v.assign(n, Vec3::Zero()); st.fixed.assign(n, 0);
		const ThinFilmGraph G = cs.Graph(st);
		st.V.resize(n); st.S.resize(n);
		for (int i = 0; i < n; i++) {
			st.V[i] = h0_m / sc.H * (1. + 0.05 * Env("VF_MOD", 1.) * std::cos(2. * kPi * x0[i][0] / Lx)) * G.a[i];
			st.S[i] = cp.G0 * (1. + 0.02 * Env("VF_MOD", 1.) * std::sin(2. * kPi * x0[i][1] / Ly)) * G.a[i];
		}
		for (int s = 0; s < nsteps; s++) {
			const FilmCoupled_Stepper::Stats ss = cs.Step(st, dt_s / sc.T);
			const ThinFilmGraph G1 = cs.Graph(st);
			Out o; o.d.resize(n); o.h.resize(n);
			for (int i = 0; i < n; i++) { o.d[i] = cs.mech.box.Wrap(st.x[i] - x0[i]) * sc.L; o.h[i] = st.V[i] / G1.a[i] * sc.H; }
			outs.push_back(o);
			std::printf("     step %d: mech %d it %s (|g| %.1e, cg %d, precond %.2f s%s), drain %d it %s, diff %d it %s\n", s + 1, ss.mech.iters, ss.mech.converged ? "ok" : "NC",
				ss.mech.grad_inf_history.empty() ? 0. : ss.mech.grad_inf_history.back(), ss.mech.total_cg_iters, ss.mech.precond_seconds, cs.mech.precond_stress_projected ? ", stress projected" : "",
				ss.drain.iters, ss.drain.converged ? "ok" : "NC", ss.diff.iters, ss.diff.converged ? "ok" : "NC");
		}
	};
	FilmScales A, B;
	A.H = 1e-8; A.L = 2e-5; A.T = 1e-3; A.S = 1e-2;
	B.H = 2.5e-8; B.L = 1e-5; B.T = 3e-4; B.S = 3e-2;
	std::vector<Out> oa, ob;
	run(A, oa);
	run(B, ob);
	for (size_t s = 0; s < oa.size(); s++) {
		real dmax = 0., ddiff = 0., hmax = 0., hdiff = 0.;
		for (size_t i = 0; i < oa[s].d.size(); i++) {
			dmax = std::max(dmax, oa[s].d[i].norm()); ddiff = std::max(ddiff, (oa[s].d[i] - ob[s].d[i]).norm());
			hmax = std::max(hmax, std::abs(oa[s].h[i] - h0_m)); hdiff = std::max(hdiff, std::abs(oa[s].h[i] - ob[s].h[i]));
		}
		std::printf("   after step %zu: max displacement %.3e m (scales differ by %.2e of it), max |h - h0| %.3e m (%.2e)\n", s + 1, dmax, ddiff / dmax, hmax,
			hdiff / hmax);
	}
}

// A vertical film in a frame of wetted fibres with Plateau-border menisci (Meniscus.h), in physical units. Released
// at rest with a uniform thickness and concentration; nothing prescribed at the borders but the menisci's geometry
// (fibre radius, meniscus radius at one height) and the friction closure. Reports the film's thickness, the
// vertical velocities of the meniscus strips, of the film next to them and of the film's centre (a circulation:
// meniscus surface rising, centre descending), the pinch along the lateral junctions, and the conservation of the
// monolayer. VF_SIDES: bit mask of the borders with a meniscus (1 left, 2 right, 4 bottom, 8 top; 15 all).
// VF_SDS=1: Monier et al.'s solution (Units.h SDS_Solution: soluble SDS at VF_C_MM mM, default 19.4, viscosity VF_MU_PA,
// default 1.5e-3, density VF_RHO, default 1045), released at rest in chemical equilibrium with the menisci's liquid.
static void Run_Framed()
{
	const bool sds = Env("VF_SDS", 0.) > 0.;
	FilmMaterial M = sds ? SDS_Solution(Env("VF_C_MM", 19.4), Env("VF_MU_PA", 1.5e-3), Env("VF_RHO", 1045.)) : FilmMaterial();
	if (!sds) { M.mu = Env("VF_MU_PA", M.mu); M.theta0 = Env("VF_THETA0", M.theta0); }
	const real W_m = Env("VF_W_M", 4e-3), H_m = Env("VF_H_M", 5e-3), dx_m = Env("VF_DX_M", 100e-6), h0_m = Env("VF_H0_M", 2e-6);
	const real dt_s = Env("VF_DT_S", 5e-3), tmax_s = Env("VF_TMAX_S", 1.), report_s = Env("VF_REPORT_S", 0.05), every_s = Env("VF_EVERY_S", 0.);
	MeniscusParams mp;
	mp.fibre_radius = Env("VF_FIBRE_M", mp.fibre_radius); mp.r_ref = Env("VF_RREF_M", mp.r_ref); mp.y_ref = Env("VF_YREF_M", mp.y_ref);
	mp.c_drag = Env("VF_CDRAG", mp.c_drag);
	mp.wedge = Env("VF_WEDGE", 1.) > 0.; mp.column = Env("VF_COLUMN", 1.) > 0.; // the wedge friction; SDS: meniscus surfaces over their liquid column (0: a bath)
	const int sides = (int)Env("VF_SIDES", 15.);
	mp.left = sides & 1; mp.right = sides & 2; mp.bottom = sides & 4; mp.top = sides & 8;
	FilmScales sc;
	sc.H = 1e-7; sc.L = dx_m; sc.T = 1e-3; sc.S = 1e-2;
	const CodeParams cp = To_Code(M, sc);
	Framed_Film_Stepper fs;
	FilmCoupled_Stepper& cs = fs.cs;
	fs.remesh=Env("VF_REMESH",1.)>0.;
	fs.equilibrate_region_remap=Env("VF_REGION_REMAP_EQUILIBRIUM",1.)>0.;
	cp.Apply(cs.mech.prm);
	cs.mech.prm.up = Vec3(0., 1., 0.);
	cs.diff.Ds = cp.Ds;
	cs.mech.factored_precond = Env("VF_FPC", 1.) > 0.;
	cs.mech.precond_cholmod = Env("VF_CHOLMOD", 1.) > 0.;
	// the reference pressure (Units.h, CodeParams::Pi_ref): for a closed frame any constant is the same continuum
	// model. VF_PREF = k sets it k times the film's weight per width (rho g h0 H, both faces) below Pi_s(theta0), so
	// that the particles carry no tension where the film hangs from its top (on particles a tension is the kernel's
	// tensile instability: pairing, negative curvature in the step). With a positive pressure the quintic kernel
	// pairs at 4 dx (Dehnen & Aly 2012); the Wendland kernels do not (VF_KERNEL, default 1: Wendland C2)
	const real pref_k = Env("VF_PREF", 1.);
	cs.mech.prm.bath_pressure = cp.Pi_ref - pref_k * 0.5 * M.rho * M.g * h0_m * H_m / sc.S;
	cs.mech.max_halvings = (int)Env("VF_HALVINGS", 3.);
	// pair-list margin (kernel radii): a step may move a particle by a quarter of it (FilmMech_Stepper::Step_Adaptive);
	// thin film rises at 1-2 cm/s, so 1 (100 um per step at dx 100 um) rather than the default 0.5
	cs.mech.list_margin = Env("VF_MARGIN", 1.);
	cs.mech.rel_tol = cs.drain.rel_tol = cs.diff.rel_tol = Env("VF_RELTOL", 1e-9);
	cs.drain.factored_precond = Env("VF_DRAIN_FPC", 1.) > 0.;
	cs.drain.newton.max_cg_iters = 20000;
	cs.mech.newton.verbose = (int)Env("VF_VERBOSE", 0.);
	// refactorize the in-plane preconditioner every k Newton iterations (0: once per step); with fast flow the Newton
	// iterates move far from the inertial prediction where it was factorized, and a stale one costs ~100 CG iterations each
	cs.mech.newton.precond_rebuild_every = (int)Env("VF_REBUILD", 0.);
	// a step whose Newton solve needs more than this is redone in halves (FilmMech_Stepper::max_halvings): as the flow
	// speeds up the film's stress varies by as much as the reference offset, the step Hessian turns indefinite at 5 ms
	// (lattice shear modes: the stress term against inertia rho h / dt^2) and Newton crawls; shorter steps are convex
	cs.mech.newton.max_iters = (int)Env("VF_MAXIT", 15.);
	const real dx = 1., R = Env("VF_RDX", 4.) * dx, dt = dt_s / sc.T;
	cs.mech.kernel.Set_Type((Kernel2D::Type)(int)Env("VF_KERNEL", 1.));
	cs.mech.kernel.Set_Radius(R);
	const Framed_Lattice lat(W_m / sc.L, H_m / sc.L, dx, Env("VF_FRAME", 2.) * R);
	cs.mech.box = lat.box;
	const Framed_Menisci men(M, cp, mp, lat.x0, lat.x1, lat.y0, lat.y1);
	if(Env("VF_REGION_REMAP",0.)>0.)fs.rp.remap_region=[&men](const Vec3& p){return men.Zone(p,false);};
	if (!men.Valid()) { std::printf("  (meniscus above the atmosphere's pressure at the bottom: r_ref too large for y_ref)\n"); return; }
	FilmMatState st;
	st.x = lat.x; st.fixed = lat.fixed;
	int n = st.Size();
	st.v.assign(n, Vec3::Zero());
	{
		const ThinFilmGraph G = cs.Graph(st);
		st.V.resize(n); st.S.resize(n);
		for (int i = 0; i < n; i++) { st.V[i] = h0_m / sc.H * G.a[i]; st.S[i] = cp.G0 * G.a[i]; }
		if (sds) { st.B.resize(n); for (int i = 0; i < n; i++) st.B[i] = cp.sol.c_bath * st.V[i]; }
	}
	auto free_surfactant = [&]() { real q = 0.; for (int i = 0; i < st.Size(); i++) if (!st.fixed[i]) q += sds ? 2. * st.S[i] + st.B[i] : st.S[i]; return q; };
	const real ymid = 0.5 * (lat.y0 + lat.y1), xmid = 0.5 * (lat.x0 + lat.x1);
	const real mu_Pa = cp.p_Pa(men.Mu());
	std::printf("\n== Framed vertical film with menisci (Meniscus.h): %g x %g mm inside fibres of radius %g um, dx %g um (R %g dx), %d particles (%d free) ==\n",
		1e3 * W_m, 1e3 * H_m, 1e6 * mp.fibre_radius, 1e6 * dx_m, R, n, (int)std::count(st.fixed.begin(), st.fixed.end(), 0));
	std::printf("   material: mu %.3g Pa s, gamma %.4f N/m, film Gibbs elasticity %.4f N/m (theta0 %g); h0 %g um; dt %g ms; borders with menisci: %s%s%s%s\n",
		M.mu, M.gamma_ref(), M.Gibbs_E(), M.theta0, 1e6 * h0_m, 1e3 * dt_s, mp.left ? "left " : "", mp.right ? "right " : "", mp.bottom ? "bottom " : "", mp.top ? "top" : "");
	if (sds) {
		// the soluble film's elasticity at h0 (local equilibrium with its liquid: Solution.h), against its weight
		const IonicSurfactantSolution::Species sp = M.sol.Solve(M.c_bulk);
		real F, Fa, Faa;
		cp.sol.Energy(1., h0_m / sc.H, 2. * cp.G0 + cp.sol.c_bath * h0_m / sc.H, F, Fa, Faa);
		std::printf("   SDS %.2f mM (%.0f%% in micelles; monomers %.2f mM, Na+ %.2f mM), rho %.0f kg/m^3: soluble film elasticity %.3g N/m at h0 against the weight rho g h0 H = %.3g N/m\n",
			M.c_bulk, 100. * sp.cm / M.c_bulk, sp.c1, sp.cNa, M.rho, Faa * sc.S, M.rho * M.g * h0_m * H_m);
	}
	std::printf("   meniscus: r %.3f mm at the bottom, %.3f mm at the top (r_ref %.3f mm at y %.2f mm); chemical potential %.2f Pa; strip width %.3f -> %.3f mm;"
		" friction 2 zeta %.3g -> %.3g Pa s/m (c_drag %g)\n", 1e3 * men.Radius(lat.y0) * sc.L, 1e3 * men.Radius(lat.y1) * sc.L, 1e3 * mp.r_ref, 1e3 * mp.y_ref, mu_Pa,
		1e3 * men.Width(lat.y0) * sc.L, 1e3 * men.Width(lat.y1) * sc.L, men.Zeta(lat.y0) * sc.S * sc.T / (sc.L * sc.L), men.Zeta(lat.y1) * sc.S * sc.T / (sc.L * sc.L), mp.c_drag);
	std::printf("   scale (not used by the solver): a strip carrying the film's weight gradient rho g h0 / 2 per face against its friction moves at ~ %.2g mm/s\n",
		1e3 * M.rho * M.g * h0_m / 2. * men.Radius(ymid) * sc.L / (mp.c_drag * M.mu));
	{
		// the menisci's liquid (both faces, all four borders) against the film's: mu_m is held fixed, so the film's
		// drainage must not change the menisci much
		const real a_m = mp.fibre_radius, nseg = 200.;
		real Vm = 0.;
		for (int k = 0; k < (int)nseg; k++) {
			const real y = lat.y0 + (k + 0.5) / nseg * (lat.y1 - lat.y0);
			Vm += 2. * 2. * Meniscus_Section::Liquid_Area(a_m, men.Radius(y) * sc.L) * H_m / nseg; // left and right
		}
		Vm += 2. * W_m * (Meniscus_Section::Liquid_Area(a_m, men.Radius(lat.y0) * sc.L) + Meniscus_Section::Liquid_Area(a_m, men.Radius(lat.y1) * sc.L));
		std::printf("   liquid in the menisci %.3g mm^3 against %.3g mm^3 in the film\n", 1e9 * Vm, 1e9 * h0_m * W_m * H_m);
	}
	{
		// the meniscus strips' friction as set up at the start: along the border (closure) and across the junction (wedge)
		FilmMatState tmp = st;
		const ThinFilmGraph G = cs.Graph(tmp);
		men.Prepare(tmp, G, cs.mech.kernel, cs.mech.box, &cs.mech.prm.sol);
		std::vector<int> zz;
		men.Zones(tmp, zz);
		int ns = 0, nw = 0;
		real wmax = 0., wsum = 0., dsum = 0.;
		const real to_si = sc.S * sc.T / (sc.L * sc.L);
		for (int i = 0; i < n; i++) {
			if (zz[i] != 1) continue;
			ns++; dsum += tmp.drag[i] / G.a[i] * to_si;
			const real e = (int)tmp.drag_n.size() == n ? tmp.drag_n[i] : 0.;
			if (e > 0.) { nw++; wsum += (e + tmp.drag[i]) / G.a[i] * to_si; wmax = std::max(wmax, (e + tmp.drag[i]) / G.a[i] * to_si); }
		}
		std::printf("   strip friction (both faces): %d strip particles, along the border %.3g Pa s/m on average; across the junction (wedge) on %d of them, %.3g on average, %.3g at most\n",
			ns, ns ? dsum / ns : 0., nw, nw ? wsum / nw : 0., wmax);
	}
	const char* frames = std::getenv("VF_FRAMES");
	FILE* ff = frames ? std::fopen(frames, "w") : nullptr;
	if (ff) std::fprintf(ff, "frame %g %g %g %g %g\n", lat.x0 * sc.L * 1e3, lat.x1 * sc.L * 1e3, lat.y0 * sc.L * 1e3, lat.y1 * sc.L * 1e3, 1e3 * sc.L);
	fs.rp.dx0 = dx;
	const char* audit_path=std::getenv("VF_INTERFACE_CSV");
	FILE* audit_file=audit_path?std::fopen(audit_path,"w"):nullptr;
	if(audit_path&&!audit_file)throw std::runtime_error("Cannot open VF_INTERFACE_CSV");
	if(audit_file)std::fprintf(audit_file,"step,time_s,kind,stage,zone,other_zone,count,V_m3,S2_mol,B_mol,Vsolution_m3,predict_dV_m3,predict_dQ_mol,energy_J,work_J,inserted,merged\n");
	InterfaceLedger audit;
	const real S_free0 = free_surfactant(); // the monolayer (insoluble) or all the surfactant (soluble: 2 S + B)
	real Q_film0 = 0.; // SDS: the film zone's surfactant at the start (the strips' columns are the menisci's liquid)
	{
		std::vector<int> z0;
		men.Zones(st, z0);
		if (sds) for (int i = 0; i < n; i++) if (z0[i] == 0) Q_film0 += 2. * st.S[i] + st.B[i];
	}
	int nc_mech = 0, nc_drain = 0;
	const int nsteps = (int)std::lround(tmax_s / dt_s);
	// Diagnostic control: decouple remesh-call frequency from integration dt.
	// Zero preserves the inherited every-step policy. Positive periods must align
	// with the outer time grid; this is not a physical material parameter.
	const real remesh_period_s=Env("VF_REMESH_PERIOD_S",0.);
	const int remesh_every=remesh_period_s>0.?std::max(1,(int)std::lround(remesh_period_s/dt_s)):1;
	if(remesh_period_s>0.&&std::abs(remesh_every*dt_s-remesh_period_s)>1e-10*remesh_period_s)
		throw std::invalid_argument("VF_REMESH_PERIOD_S must be an integer multiple of dt");
	const bool remesh_enabled=fs.remesh;
	std::printf("  %7s %6s %8s %8s %8s %8s %8s %9s %9s %9s %9s %8s %10s %6s %9s %7s\n", "t (s)", "N", "film h", "h centre", "h side", "pinch", "film %",
		"vy strip", "vy side", "vy centre", "max|v|", "dS/S", "to men.", "mech", "ins/mrg", "sec");
	std::printf("  %7s %6s %8s %8s %8s %8s %8s %9s %9s %9s %9s %8s %10s %6s %9s %7s\n", "", "", "(um)", "(um)", "(um)", "(um)", "", "(mm/s)", "(mm/s)", "(mm/s)",
		"(mm/s)", "", "(um^3/um2)", "it/cg", "", "");
	const auto t0 = std::chrono::steady_clock::now();
	std::vector<int> z;
	const bool steplog = Env("VF_STEPLOG", 0.) != 0.;
	for (int s = 1; s <= nsteps; s++) {
		const auto ts0 = std::chrono::steady_clock::now();
		fs.remesh=remesh_enabled&&(s%remesh_every==0);
		const Framed_Film_Stepper::Stats ss = fs.Step(st, men, dt,audit_file?&audit:nullptr);
		if(audit_file)audit.Write(audit_file,s,s*dt_s,sc.H*sc.L*sc.L,M.Gamma_inf*sc.L*sc.L,sc.S*sc.L*sc.L);
		const NewtonStats& nm = ss.mech;
		if (!nm.converged) nc_mech++;
		if (!ss.drain.converged) nc_drain++;
		if (steplog) std::printf("    step %d: mech %d it %s, cg %d, %d factorizations (%.2f s); drainage %d it %s; %.2f s in all\n", s, nm.iters, nm.converged ? "ok" : "NC",
			nm.total_cg_iters, nm.precond_factorizations, nm.precond_seconds, ss.drain.iters, ss.drain.converged ? "ok" : "NC",
			std::chrono::duration<double>(std::chrono::steady_clock::now() - ts0).count());
		n = st.Size();
		const real t_s = s * dt_s;
		const bool rep = std::abs(t_s / report_s - std::round(t_s / report_s)) < 0.5 * dt_s / report_s;
		const bool frm = ff && every_s > 0. && std::abs(t_s / every_s - std::round(t_s / every_s)) < 0.5 * dt_s / every_s;
		if (!rep && !frm && s != nsteps) continue;
		const ThinFilmGraph Gp = cs.Graph(st);
		men.Zones(st, z);
		if (frm) {
			std::fprintf(ff, "t %g\n", t_s);
			for (int i = 0; i < n; i++) {
				std::fprintf(ff, "%.5f %.5f %.5g %.5f %d %.4g %.4g", st.x[i][0] * sc.L * 1e3, st.x[i][1] * sc.L * 1e3, st.V[i] / Gp.a[i] * sc.H * 1e6, st.S[i] / Gp.a[i],
					z[i], st.v[i][0] * sc.L / sc.T * 1e3, st.v[i][1] * sc.L / sc.T * 1e3);
				if (sds) std::fprintf(ff, " %.4g", st.B[i] / st.Vsol(i) * cp.sol.c_unit); // SDS: + the dissolved concentration (mM)
				std::fprintf(ff, "\n");
			}
		}
		if (!rep && s != nsteps) continue;
		// film: mean thickness; centre column (|x - xmid| < W / 8) and side bands (film within 3 dx of a lateral strip) at
		// mid-height (|y - ymid| < H / 8); the pinch: the thinnest 3 dx band along the lateral junctions (bins of 2 R in y)
		real Vf = 0., Af = 0., hc = 0., ac = 0., hs = 0., as = 0., vys = 0., ws = 0., vyc = 0., wc = 0., vym = 0., wm = 0., vmax = 0., S_free = 0.;
		real cmin = 1e300, cmax = 0., thmin = 1., thmax = 0., ctop = 0., vtop = 0., cbot = 0., vbot = 0.; real Q_film = 0.; // SDS
		const int nb = std::max(1, (int)((lat.y1 - lat.y0) / (2. * R)));
		std::vector<real> pv(nb, 0.), pa(nb, 0.);
		// the film's width-averaged thickness in bins of 2 dx in y, for the fringes (isothickness lines, as the
		// experiments track them): the height where it first reaches h* coming down from the top
		const int nfb = std::max(1, (int)((lat.y1 - lat.y0) / (2. * dx)));
		std::vector<real> fv(nfb, 0.), fa(nfb, 0.);
		for (int i = 0; i < n; i++) {
			if (st.fixed[i]) continue;
			S_free += sds ? 2. * st.S[i] + st.B[i] : st.S[i];
			if (sds && z[i] == 0) Q_film += 2. * st.S[i] + st.B[i];
			vmax = std::max(vmax, st.v[i].norm());
			const real x = st.x[i][0], y = st.x[i][1], a = Gp.a[i];
			const bool mid = std::abs(y - ymid) < 0.125 * (lat.y1 - lat.y0);
			if (z[i] == 1) {
				const bool lateral = (mp.left && x - lat.x0 < men.Width(y)) || (mp.right && lat.x1 - x < men.Width(y));
				if (lateral && mid) { vym += a * st.v[i][1]; wm += a; }
				continue;
			}
			Vf += st.V[i]; Af += a;
			if (sds) {
				const real c = st.B[i] / st.Vsol(i) * cp.sol.c_unit, th = st.S[i] / a, yr = (y - lat.y0) / (lat.y1 - lat.y0);
				cmin = std::min(cmin, c); cmax = std::max(cmax, c); thmin = std::min(thmin, th); thmax = std::max(thmax, th);
				if (yr > 0.75) { ctop += st.B[i]; vtop += st.Vsol(i); }
				if (yr < 0.25) { cbot += st.B[i]; vbot += st.Vsol(i); }
			}
			{ const int b = std::clamp((int)((y - lat.y0) / (2. * dx)), 0, nfb - 1); fv[b] += st.V[i]; fa[b] += a; }
			const real dside = std::min(mp.left ? x - lat.x0 - men.Width(y) : 1e300, mp.right ? lat.x1 - x - men.Width(y) : 1e300);
			if (dside < 3. * dx) {
				const int b = std::clamp((int)((y - lat.y0) / (2. * R)), 0, nb - 1);
				pv[b] += st.V[i]; pa[b] += a;
				if (mid) { hs += st.V[i]; as += a; vys += a * st.v[i][1]; ws += a; }
			}
			if (std::abs(x - xmid) < 0.125 * (lat.x1 - lat.x0) && mid) { hc += st.V[i]; ac += a; vyc += a * st.v[i][1]; wc += a; }
		}
		real pinch = 1e300;
		for (int b = 1; b + 1 < nb; b++) if (pa[b] > 0.) pinch = std::min(pinch, pv[b] / pa[b]);
		const real um = sc.H * 1e6, mms = sc.L / sc.T * 1e3;
		std::printf("  %7.3f %6d %8.4f %8.4f %8.4f %8.4f %8.1f %+9.3f %+9.3f %+9.3f %9.3f %+8.1e %10.3e %2d/%-3d %4d/%-4d %7.1f%s\n", t_s, n, Vf / Af * um, ac > 0. ? hc / ac * um : 0.,
			as > 0. ? hs / as * um : 0., pinch < 1e300 ? pinch * um : 0., 100. * Af / ((lat.x1 - lat.x0) * (lat.y1 - lat.y0)), wm > 0. ? vym / wm * mms : 0., ws > 0. ? vys / ws * mms : 0.,
			wc > 0. ? vyc / wc * mms : 0., vmax * mms, sds ? Q_film / Q_film0 - 1. : S_free / S_free0 - 1., fs.drained * sc.H * sc.L * sc.L * 1e18 / (W_m * H_m * 1e12), nm.iters, nm.total_cg_iters, fs.inserted, fs.merged,
			std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), (nc_mech || nc_drain) ? " (Newton not converged in some steps)" : "");
		std::printf("          fringes z/H:");
		for (real lv : { 0.98, 0.95, 0.9, 0.8, 0.7, 0.5 }) {
			real zf = 0.;
			for (int b = nfb - 1; b >= 0; b--) if (fa[b] > 0. && fv[b] / fa[b] >= lv * h0_m / sc.H) { zf = std::min<real>(1., (b + 0.5) * 2. * dx / (lat.y1 - lat.y0)); break; }
			std::printf("  %.2f h0 %.3f", lv, zf);
		}
		std::printf("\n");
		if (sds) {
			auto gam = [&](real th) { return 1e3 * (M.gamma0 - cp.sol.Pi(th) * sc.S); };
			std::printf("          SDS: film c %.2f - %.2f mM (top quarter %.2f, bottom quarter %.2f), tension %.3f - %.3f mN/m; film zone's surfactant %+.3e; menisci gave %+.3e of all the free particles' (balance %.1e)\n",
				cmin, cmax, vtop > 0. ? ctop / vtop * cp.sol.c_unit : 0., vbot > 0. ? cbot / vbot * cp.sol.c_unit : 0., gam(thmax), gam(thmin), Q_film / Q_film0 - 1., fs.exchanged / S_free0,
				(S_free - S_free0 - fs.exchanged) / S_free);
		}
	}
	if (ff) std::fclose(ff);
	if(audit_file)std::fclose(audit_file);
	std::printf("   Newton failures: mechanics %d, drainage %d; %.0f s\n", nc_mech, nc_drain, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
}

// The factored in-plane preconditioner against the Hessian it inverts: a physical film (the units bench's state), the
// first Newton iterate of one BDF1 step. Residual |H P^-1 r - r| / |r| (0 if exact: no squeeze energy, VF_NOSQ=1),
// negative pivots vs auxiliary rows (equal: positive definite), and CG on H with each preconditioner.
static void Run_Precond_Check()
{
	FilmMaterial M;
	const real dx_m = Env("VF_DX_M", 20e-6), h0_m = Env("VF_H0_M", 500e-9), dt_s = Env("VF_DT_S", 1e-3);
	const int nx = (int)Env("VF_NX", 16.);
	FilmScales sc;
	sc.H = 1e-8; sc.L = 2e-5; sc.T = 1e-3; sc.S = 1e-2;
	const CodeParams cp = To_Code(M, sc);
	FilmMech_Stepper ms;
	cp.Apply(ms.prm);
	if (Env("VF_NOSQ", 0.) > 0.) ms.prm.tf.sigma = 0.;
	ms.prm.bulk_viscosity = Env("VF_BULK", 1.) > 0.;
	// VF_BATH=0: the full surface pressure on the particles (To_Code subtracts the reference: CodeParams::Pi_ref)
	if (Env("VF_BATH", 1.) <= 0.) ms.prm.bath_pressure = 0.;
	const real dx = dx_m / sc.L;
	real Lx, Ly;
	std::vector<Vec3> x0 = Hex_Patch(nx, nx, dx, Env("VF_JIT", 0.05), 5, Lx, Ly);
	ms.kernel.Set_Type((Kernel2D::Type)(int)Env("VF_KERNEL", 0.));
	ms.kernel.Set_Radius(Env("VF_RDX", 4.) * dx);
	ms.box.Lx = Lx; ms.box.Ly = Ly;
	const ThinFilmGraph G = Make_Periodic_Patch(x0, Lx, Ly, ms.kernel);
	FilmMatState st;
	st.x = x0;
	const int n = st.Size();
	st.v.assign(n, Vec3::Zero()); st.V.resize(n); st.S.resize(n);
	for (int i = 0; i < n; i++) {
		st.V[i] = h0_m / sc.H * (1. + 0.05 * Env("VF_MOD", 1.) * std::cos(2. * kPi * x0[i][0] / Lx)) * G.a[i];
		st.S[i] = cp.G0 * (1. + 0.02 * Env("VF_MOD", 1.) * std::sin(2. * kPi * x0[i][1] / Ly)) * G.a[i];
	}
	ShiftedPairs sp;
	Build_Shifted_Pairs(st.x, ms.box, 1.5 * ms.kernel.Radius(), sp);
	const real hh = dt_s / sc.T;
	FilmMech_Problem P(ms.kernel, ms.prm, sp, st, st.x, st.x, hh);
	P.precond_cholmod = Env("VF_CHOLMOD", 1.) > 0.;
	const int m = P.Dim();
	std::printf("\n== Factored in-plane preconditioner: %d particles, dx %.3g m, dt %.3g s, squeeze %s, bulk viscosity %s ==\n", n, dx_m, dt_s,
		ms.prm.tf.sigma > 0. ? "on" : "off", ms.prm.bulk_viscosity ? "on" : "off");
	VecX q = VecX::Zero(m), r(m);
	for (int k = 0; k < m; k++) r[k] = std::sin(1.3 * k + 0.7);
	{
		// effective shear modulus of the particle set: Rayleigh quotient of the energy Hessian (no inertia, no viscosity)
		// on shear waves v_x = sin(k y): v^T H v = mu_eff k^2 sum_i a_i v_i^2 (a fluid: 0), against the viscous 4 mu h / dt
		FilmMechParams p0 = ms.prm;
		p0.rho = 0.; p0.mu_s = 0.;
		FilmMech_Problem P0(ms.kernel, p0, sp, st, st.x, st.x, 1e30);
		P0.Prepare_Hessian(q);
		std::printf("   shear waves (energy Hessian only): ");
		for (int mk : { 1, 2, 4 }) {
			const real kk = 2. * kPi * mk / Ly;
			VecX v = VecX::Zero(m), Hv;
			real av = 0.;
			for (int i = 0; i < n; i++) { v[2 * i] = std::sin(kk * x0[i][1]); av += G.a[i] * v[2 * i] * v[2 * i]; }
			P0.Hess_Vec(v, Hv);
			const real mu_eff = v.dot(Hv) / (kk * kk * av); // code tension units S
			std::printf(" lambda = Ly/%d: mu_eff %+.2e N/m;", mk, mu_eff * sc.S);
		}
		std::printf("  viscous 4 mu h / dt = %.2e N/m\n", 4. * M.mu * h0_m / dt_s);
	}
	P.Prepare_Hessian(q);
	if (P.Factor_Precond(false)) {
		// relaxed: the stiffness of one implicit step (inertia, viscosity, energy) against shear waves when every other
		// particle motion adjusts: min x^T H x over x with (v, x)_M = (v, v)_M is (v, v)_M^2 / (v, H^-1 M v)_M
		std::printf("   shear waves, one step of %.3g s, non-affine relaxation allowed: ", dt_s);
		for (int mk : { 1, 2, 4 }) {
			const real kk = 2. * kPi * mk / Ly;
			VecX Mv = VecX::Zero(m), z;
			real vMv = 0.;
			for (int i = 0; i < n; i++) { const real v = std::sin(kk * x0[i][1]); Mv[2 * i] = G.a[i] * v; vMv += G.a[i] * v * v; }
			P.Apply_Preconditioner(Mv, z);
			std::printf(" Ly/%d: %.3e N/m;", mk, vMv / Mv.dot(z) / (kk * kk) * sc.S);
		}
		std::printf("\n");
	}
	auto pcg = [&](int mode, int& negc) { // CG on H z = r; mode 0 Jacobi, 1 factored
		VecX D, z = VecX::Zero(m), res = r, y, d, Hd;
		P.Hess_Diag(D);
		auto prec = [&](const VecX& a, VecX& b) { if (mode) P.Apply_Preconditioner(a, b); else b = a.cwiseQuotient(D); };
		prec(res, y); d = y;
		real ry = res.dot(y);
		negc = 0;
		int k = 0;
		for (; k < 5000 && res.norm() > 1e-10 * r.norm(); k++) {
			P.Hess_Vec(d, Hd);
			const real dHd = d.dot(Hd);
			if (!(dHd > 0.)) { negc++; break; }
			const real al = ry / dHd;
			z += al * d; res -= al * Hd;
			prec(res, y);
			const real ry2 = res.dot(y);
			d = y + (ry2 / ry) * d; ry = ry2;
		}
		return k;
	};
	int negc = 0;
	const int kj = pcg(0, negc);
	std::printf("   Jacobi CG: %d iterations to 1e-10%s\n", kj, negc ? " (negative curvature: H indefinite)" : "");
	for (int proj = 0; proj < 2; proj++) {
		const auto t0 = std::chrono::steady_clock::now();
		const bool ok = P.Factor_Precond(proj);
		const double ts = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
		VecX z, Hz;
		P.Apply_Preconditioner(r, z);
		P.Hess_Vec(z, Hz);
		const int kf = pcg(1, negc);
		std::printf("   [setFromTriplets %.2f s, ordering %.2f s, factorization %.2f s (cumulative); matrix lower %lld nonzeros, factor %lld]\n", P.precond_t_assemble,
			P.precond_t_analyze, P.precond_t_factor, P.precond_nnz_A, P.precond_nnz_L);
		std::printf("   stress %s: factorized %.3f s, %s (%d negative pivots, %d auxiliary rows), |H P^-1 r - r| / |r| = %.2e, r.P^-1 r %s 0; CG %d iterations%s\n",
			proj ? "projected" : "exact    ", ts, ok ? "positive definite" : "NOT positive definite", P.precond_negatives, P.precond_aux, (Hz - r).norm() / r.norm(),
			r.dot(z) > 0. ? ">" : "<=", kf, negc ? " (negative curvature)" : "");
	}
}

static void Run_Wave(const std::vector<real>& dxs)
{
	const real L = Env("VF_L", 64.), h0 = 6., G0 = Env("VF_G0", 0.5), A = Env("VF_AMP", 0.05) * L / 64., mu = Env("VF_MU", 0.1), periods = Env("VF_PERIODS", 3.), ppt = Env("VF_PPT", 80.);
	FilmMechParams prm;
	prm.tf.sigma = 1.; prm.tf.mu = mu; prm.tf.disj = Disj_From_Env(1., 1., 1.5);
	prm.surf.RT = 0.5; prm.surf.Gamma_inf = Env("VF_GINF", 1.);
	prm.rho = 1.;
	const real k = 2. * kPi / L;
	const real EG = 2. * G0 * G0 * prm.surf.d2g(G0), Keff = EG + h0 * h0 * prm.tf.disj.d2P(h0) + prm.tf.sigma * h0 * h0 / 4. * k * k;
	const real w0 = k * std::sqrt(Keff / (prm.rho * h0)), b = 2. * mu * k * k / prm.rho, wd = std::sqrt(w0 * w0 - b * b);
	const real T = 2. * kPi / wd, dt = T / ppt;
	std::printf("\n== Longitudinal waves, mobile surfaces (h0 %g, Gamma0 %g, RT 0.5, Gamma_inf %g: E_G %.4f; h^2 P'' %.4f; squeeze %.4f; mu %g; k %.4f; BDF2 dt = T/%g) ==\n",
		h0, G0, prm.surf.Gamma_inf, EG, h0 * h0 * prm.tf.disj.d2P(h0), prm.tf.sigma * h0 * h0 / 4. * k * k, mu, k, ppt);
	std::printf("   exact: omega_d %.6f, decay rate b %.6f\n", wd, b);
	std::printf("%6s %8s %12s %10s %12s %10s %8s %10s\n", "dx", "N", "omega_d", "rel err", "decay b", "rel err", "newton", "Hess RQ");
	for (real dx : dxs) {
		const int nx = (int)std::lround(L / dx);
		const real rdx = Env("VF_RDX", 4.);
		int ny = (int)std::ceil(2.2 * rdx * dx / (dx * std::sqrt(3.) / 2.));
		ny += ny & 1;
		real Lx, Ly;
		std::vector<Vec3> X = Hex_Patch(nx, ny, dx, 0., 1, Lx, Ly);
		FilmMech_Stepper ms;
		ms.kernel.Set_Radius(rdx * dx);
		ms.prm = prm;
		ms.box.Lx = Lx; ms.box.Ly = Ly;
		ThinFilmGraph G = Make_Periodic_Patch(X, Lx, Ly, ms.kernel);
		FilmMatState st;
		const int n = (int)X.size();
		st.x = X; st.v.assign(n, Vec3::Zero()); st.V.resize(n); st.S.resize(n);
		for (int i = 0; i < n; i++) { st.V[i] = h0 * G.a[i]; st.S[i] = G0 * G.a[i]; st.x[i][0] += A * std::sin(k * X[i][0]); }
		auto amp = [&]() { real s = 0.; for (int i = 0; i < n; i++) s += 2. / n * ms.box.Wrap(st.x[i] - X[i])[0] * std::sin(k * X[i][0]); return s; };
		// discrete stiffness: Rayleigh quotient of the energy Hessian (inertia, viscosity off) for the mode at rest
		real w_rq = 0.;
		{
			FilmMatState s0 = st;
			s0.x = X;
			ShiftedPairs sp;
			Build_Shifted_Pairs(s0.x, ms.box, 1.5 * ms.kernel.Radius(), sp);
			FilmMechParams p0 = prm; p0.tf.mu = 0.;
			FilmMech_Problem P(ms.kernel, p0, sp, s0, s0.x, s0.x, 1e12);
			VecX q = VecX::Zero(P.Dim()), v(P.Dim()), Hv;
			for (int i = 0; i < n; i++) { v[2 * i] = std::sin(k * X[i][0]); v[2 * i + 1] = 0.; }
			P.Prepare_Hessian(q);
			P.Hess_Vec(v, Hv);
			real m = 0.;
			for (int i = 0; i < n; i++) m += prm.rho * s0.V[i] * v[2 * i] * v[2 * i];
			w_rq = std::sqrt(std::max<real>(0., v.dot(Hv) / m - b * b));
		}
		std::vector<real> ts, us;
		ts.push_back(0.); us.push_back(amp());
		const int nsteps = (int)std::lround(periods * ppt);
		int newton = 0;
		for (int s = 1; s <= nsteps; s++) { newton += ms.Step(st, dt).iters; ts.push_back(s * dt); us.push_back(amp()); }
		// zero crossings (linear interpolation) and extrema between them
		std::vector<real> zc, ext;
		for (size_t s = 1; s < us.size(); s++) if ((us[s - 1] > 0.) != (us[s] > 0.)) zc.push_back(ts[s - 1] + (ts[s] - ts[s - 1]) * us[s - 1] / (us[s - 1] - us[s]));
		for (size_t s = 1; s + 1 < us.size(); s++) if (std::abs(us[s]) >= std::abs(us[s - 1]) && std::abs(us[s]) > std::abs(us[s + 1])) ext.push_back(std::abs(us[s]));
		real wm = 0., bm = 0.;
		if (zc.size() >= 2) wm = kPi * (zc.size() - 1) / (zc.back() - zc.front());
		if (ext.size() >= 2 && wm > 0.) bm = -std::log(ext.back() / ext.front()) / ((ext.size() - 1) * kPi / wm);
		std::printf("%6.2f %8d %12.6f %+10.2e %12.6f %+10.2e %8.1f %+10.2e\n", dx, n, wm, wm / wd - 1., bm, b > 0. ? bm / b - 1. : bm, (double)newton / nsteps, w_rq / wd - 1.);
	}
}

int main(int argc, char** argv)
{
	std::setvbuf(stdout, nullptr, _IONBF, 0);
	vf::Set_Threads(-1);
	const std::string what = argc > 1 ? argv[1] : "modes";
	std::vector<real> args;
	for (int a = 2; a < argc; a++) args.push_back(std::atof(argv[a]));
	if (what == "potential") Run_Potential();
	if (what == "modes") { std::vector<int> n; for (real v : args) n.push_back((int)v); if (n.empty()) n = { 32, 64, 128 }; Run_Modes(n); }
	if (what == "front") { if (args.empty()) args = { 1., 0.5 }; Run_Front(args); }
	if (what == "wave") { if (args.empty()) args = { 1., 0.5 }; Run_Wave(args); }
	if (what == "vdrain") { if (args.empty()) args = { 1. }; Run_VDrain(args); }
	if (what == "vmobile") Run_VMobile(args.empty() ? 1. : args[0]);
	if (what == "regen") Run_Regen(args.empty() ? 1. : args[0]);
	if (what == "regenopen") Run_Regen_Open(args.empty() ? 1. : args[0]);
	if (what == "units") Run_Units();
	if (what == "pcheck") Run_Precond_Check();
	if (what == "framed") Run_Framed();
	if (what == "mbubble") Run_MBubble(args.empty() ? 1. : args[0]);
	if (what == "bubble") { if (args.empty()) args = { 1. }; Run_Bubble(args); }
	if (what == "ptime") Run_Precond_Timing(args.empty() ? 12. : args[0]);
	if (what == "ml") { if (args.empty()) args = { 1., 0.5 }; Run_Two_Level_Timing(args); }
	if (what == "spots") { std::vector<int> n; for (real v : args) n.push_back((int)v); if (n.empty()) n = { 96 }; Run_Spots(n); }
	return 0;
}
