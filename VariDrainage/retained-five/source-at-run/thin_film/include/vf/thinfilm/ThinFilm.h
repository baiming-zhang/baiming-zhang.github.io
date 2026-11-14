//////////////////////////////////////////////////////////////////////////
// Thin film: film thickness, drainage and black films (thin_film/README.md).
//
// Particles (here at fixed positions: the surfaces are immobile) carry liquid volumes V_i; the film thickness is
// h_i = V_i / a_i with the particle area a_i. Energy (per film, sigma_f the film tension of both interfaces):
//
//   E(V) = sum_i a_i P(h_i) + (sigma_f / 8) sum_pairs c_p (h_i - h_j)^2 + rho_g sum_i V_i z_i
//
// - P(h): disjoining (interaction) energy per area, Pi = -P';
// - the squeeze energy: each interface displaced by h / 2 gains area |grad h|^2 / 8; sum_p c_p (f_i - f_j)^2 approximates
//   int |grad f|^2 (c_p: Brookshaw weights of the kernel, or 1 / dx on a uniform 1D grid);
// - gravity: the weight of the liquid at height z.
// The liquid drains along the pairs (Poiseuille flow between the immobile surfaces, mobility M = h^3 / (12 mu)): a flux
// F_p from i to j costs F_p^2 / (2 k_p) per unit time with k_p = M(h_p) c_p, so that dV_i / dt = sum_j k_ij (p_j - p_i)
// ~ a_i div(M grad p), p_i = dE / dV_i (the liquid pressure). One backward-Euler step is the minimizing movement
//   min_F  E(V0 - D F) + sum_p F_p^2 / (2 k_p dt);
// its minimizer is a potential flow F = K D^T phi, so with L = D K D^T (weighted graph Laplacian) the step minimizes
//   J(phi) = E(V0 - L phi) + phi^T L phi / (2 dt)
// over one unknown per particle (Newton-CG). Volume is conserved exactly, the energy cannot increase.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>
#include "vf/core/NewtonCG.h"
#include "vf/core/Parallel.h"
#include "vf/thinfilm/Particles.h"

namespace vf::thinfilm {

// Disjoining potential per unit area:
//   P(h) = eps [ (1/3) (hb/h)^8 - (4/3) (hb/h)^2 ] + beta exp(-h / lambda)
// eps term: minimum -eps at hb (black film), van der Waals tail -(4/3) eps hb^2 / h^2, steep wall for h -> 0;
// beta term: double-layer repulsion of range lambda (makes thick films metastable behind a barrier).
// P(h) = eps [(1/3)(hb/h)^8 - (4/3)(hb/h)^2] + c8 / h^8 - c2 / h^2 + beta exp(-h / lambda). The eps form is the
// dimensionless benchmark potential; c8, c2, beta, lambda the DLVO form in physical units (Units.h): short-range
// repulsion, van der Waals (c2 = A / 12 pi), electrostatic double layer (Debye length lambda).
struct Disjoining {
	real eps = 0., hb = 1., beta = 0., lambda = 1.;
	real c8 = 0., c2 = 0.;
	real P(real h) const
	{
		const real x = hb / h, x2 = x * x, x8 = x2 * x2 * x2 * x2, y2 = 1. / (h * h), y8 = y2 * y2 * y2 * y2;
		return eps * (x8 / 3. - 4. / 3. * x2) + c8 * y8 - c2 * y2 + (beta != 0. ? beta * std::exp(-h / lambda) : 0.);
	}
	real dP(real h) const
	{
		const real x = hb / h, x2 = x * x, x8 = x2 * x2 * x2 * x2, y2 = 1. / (h * h), y8 = y2 * y2 * y2 * y2;
		return eps * (-8. / 3. * x8 + 8. / 3. * x2) / h + (-8. * c8 * y8 + 2. * c2 * y2) / h - (beta != 0. ? beta / lambda * std::exp(-h / lambda) : 0.);
	}
	real d2P(real h) const
	{
		const real x = hb / h, x2 = x * x, x8 = x2 * x2 * x2 * x2, y2 = 1. / (h * h), y8 = y2 * y2 * y2 * y2;
		return eps * (24. * x8 - 8. * x2) / (h * h) + (72. * c8 * y8 - 6. * c2 * y2) / (h * h) + (beta != 0. ? beta / (lambda * lambda) * std::exp(-h / lambda) : 0.);
	}
	real Pi(real h) const { return -dP(h); }
};

// Particles of one film with fixed positions: areas, pairs and the pair weights of the gradient energy.
struct ThinFilmGraph {
	int N = 0;
	std::vector<real> a;     // particle areas
	PairList pairs;          // with incidence
	std::vector<real> c;     // sum_p c_p (f_i - f_j)^2 ~ int |grad f|^2
	std::vector<Vec3> x;     // positions (diagnostics)
	std::vector<real> height; // height of each particle for gravity (empty: x[i][2])
	real Height(int i) const { return height.empty() ? x[i][2] : height[i]; }
};

// Uniform 1D grid of n cells of width dx (the finite-difference reference): a_k = dx, c = 1 / dx. Periodic, or closed
// at both ends (no flux through them).
inline ThinFilmGraph Make_Grid_1D(int n, real dx, bool periodic = true)
{
	ThinFilmGraph G;
	G.N = n;
	G.a.assign(n, dx);
	for (int k = 0; k < n; k++) G.x.push_back(Vec3((k + 0.5) * dx, 0., 0.));
	for (int k = 0; k < (periodic ? n : n - 1); k++) { G.pairs.pi.push_back(k); G.pairs.pj.push_back((k + 1) % n); G.c.push_back(1. / dx); }
	G.pairs.Build_Incidence(n);
	return G;
}

// Flat patch [0, Lx) x [0, Ly) (z = 0) of the given particles, periodic in x and (optionally) in y: kernel densities
// with minimum-image distances, a_i = 1 / N_i, Brookshaw weights c_p = 2 |W'(r)| / (r N_i N_j). Without periodicity in y
// the patch is closed at y = 0 and y = Ly (no pairs across: no flux; the kernel sums there are one-sided).
// Requires Lx (and a periodic Ly) > 2 kernel radii.
inline ThinFilmGraph Make_Periodic_Patch(const std::vector<Vec3>& x, real Lx, real Ly, const Quintic2D& ker, bool periodic_y = true)
{
	ThinFilmGraph G;
	const int n = (int)x.size();
	G.N = n; G.x = x;
	const real h = ker.Radius();
	auto wrap = [&](Vec3 d) { d[0] -= Lx * std::round(d[0] / Lx); if (periodic_y) d[1] -= Ly * std::round(d[1] / Ly); return d; };
	const int ncx = std::max(1, (int)std::floor(Lx / h)), ncy = std::max(1, (int)std::floor(Ly / h));
	auto cell = [&](const Vec3& p) {
		int cx = (int)std::floor(p[0] / Lx * ncx), cy = (int)std::floor(p[1] / Ly * ncy);
		cx = ((cx % ncx) + ncx) % ncx; cy = periodic_y ? ((cy % ncy) + ncy) % ncy : std::clamp(cy, 0, ncy - 1);
		return cy * ncx + cx;
	};
	std::vector<std::vector<int>> cells(ncx * ncy);
	for (int i = 0; i < n; i++) cells[cell(x[i])].push_back(i);
	std::vector<std::vector<int>> nb(n);
	Parallel_For(n, [&](int i) {
		const int c0 = cell(x[i]), cx = c0 % ncx, cy = c0 / ncx;
		std::vector<int> seen;
		for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
			if (!periodic_y && (cy + dy < 0 || cy + dy >= ncy)) continue;
			const int c = (((cy + dy) % ncy + ncy) % ncy) * ncx + (((cx + dx) % ncx + ncx) % ncx);
			if (std::find(seen.begin(), seen.end(), c) != seen.end()) continue;
			seen.push_back(c);
			for (int j : cells[c]) if (j > i && wrap(x[j] - x[i]).norm() < h) nb[i].push_back(j);
		}
		std::sort(nb[i].begin(), nb[i].end());
	}, 64);
	for (int i = 0; i < n; i++) for (int j : nb[i]) { G.pairs.pi.push_back(i); G.pairs.pj.push_back(j); }
	G.pairs.Build_Incidence(n);
	const int np = G.pairs.Size();
	std::vector<real> r(np), nd(n, ker.W(0.));
	for (int p = 0; p < np; p++) r[p] = wrap(x[G.pairs.pj[p]] - x[G.pairs.pi[p]]).norm();
	for (int p = 0; p < np; p++) { const real w = ker.W(r[p]); nd[G.pairs.pi[p]] += w; nd[G.pairs.pj[p]] += w; }
	G.a.resize(n);
	for (int i = 0; i < n; i++) G.a[i] = 1. / nd[i];
	G.c.resize(np);
	for (int p = 0; p < np; p++) G.c[p] = r[p] > 0. ? 2. * std::abs(ker.dW(r[p])) / (r[p] * nd[G.pairs.pi[p]] * nd[G.pairs.pj[p]]) : 0.;
	return G;
}

// Particles on a curved surface in 3D (a bubble): the same densities and Brookshaw weights on the chord distances (the
// chord density reproduces the area element to second order, A1; the weights give the Laplace-Beltrami operator up to
// O((R / rho)^2), R the kernel radius, rho the radius of curvature). No periodicity; gravity along -z by default.
inline ThinFilmGraph Make_Surface_Graph(const std::vector<Vec3>& x, const Quintic2D& ker)
{
	ThinFilmGraph G;
	const int n = (int)x.size();
	G.N = n; G.x = x;
	Build_Pairs(x, ker.Radius(), G.pairs);
	G.pairs.Build_Incidence(n);
	const int np = G.pairs.Size();
	std::vector<real> r(np), nd(n, ker.W(0.));
	for (int p = 0; p < np; p++) r[p] = (x[G.pairs.pj[p]] - x[G.pairs.pi[p]]).norm();
	for (int p = 0; p < np; p++) { const real w = ker.W(r[p]); nd[G.pairs.pi[p]] += w; nd[G.pairs.pj[p]] += w; }
	G.a.resize(n);
	for (int i = 0; i < n; i++) G.a[i] = 1. / nd[i];
	G.c.resize(np);
	for (int p = 0; p < np; p++) G.c[p] = r[p] > 0. ? 2. * std::abs(ker.dW(r[p])) / (r[p] * nd[G.pairs.pi[p]] * nd[G.pairs.pj[p]]) : 0.;
	return G;
}

// Axisymmetric reference on a sphere of radius R centred at the origin: n zones of equal polar-angle width in
// [th0, th1] (0, pi: the whole sphere; th0 = 0 starts with the polar cap). Zone areas 2 pi R^2 (cos th_k - cos th_k+1);
// weights 2 pi sin(th_k+1/2) / dth between neighbouring zones (flux through the circle of circumference 2 pi R sin th
// per gradient / (R dth); the same weights give int |grad f|^2); height z = R (cos th_k + cos th_k+1) / 2, the zone's
// area-mean of z.
inline ThinFilmGraph Make_Sphere_Axisym(int n, real R, real th0 = 0., real th1 = kPi)
{
	ThinFilmGraph G;
	G.N = n;
	const real dth = (th1 - th0) / n;
	G.a.resize(n); G.height.resize(n);
	for (int k = 0; k < n; k++) {
		const real a = th0 + k * dth, b = a + dth, m = a + 0.5 * dth;
		G.a[k] = 2. * kPi * R * R * (std::cos(a) - std::cos(b));
		G.height[k] = 0.5 * R * (std::cos(a) + std::cos(b));
		G.x.push_back(Vec3(R * std::sin(m), 0., R * std::cos(m)));
	}
	for (int k = 0; k + 1 < n; k++) { G.pairs.pi.push_back(k); G.pairs.pj.push_back(k + 1); G.c.push_back(2. * kPi * std::sin(th0 + (k + 1) * dth) / dth); }
	G.pairs.Build_Incidence(n);
	return G;
}

struct ThinFilmParams {
	real sigma = 1.;        // film tension (both interfaces): squeeze energy sigma / 8 |grad h|^2
	real mu = 1.;           // liquid viscosity: mobility h^3 / (12 mu)
	Disjoining disj;
	real rho_g = 0.;        // rho g (gravity along -z)
	real Mobility(real h) const { return h * h * h / (12. * mu); }
};

// Energy, liquid pressure p = dE/dV and its Hessian on a graph
class ThinFilmEnergy {
public:
	ThinFilmEnergy(const ThinFilmGraph& G, const ThinFilmParams& p) : G_(G), p_(p) {}
	// Reservoir nodes (r[i] not NaN): a bath exchanging liquid at the fixed chemical potential r[i] (e.g. a meniscus:
	// its capillary suction plus the hydrostatic offset). Their energy is r[i] V_i, they take no part in the squeeze
	// energy (the film next to them has a free, zero-gradient edge for h), so p_i = r[i] whatever their volume.
	void Set_Reservoir(const std::vector<real>* r) { res_ = r; }
	bool Reservoir(int i) const { return res_ && !std::isnan((*res_)[i]); }
	// +inf if a thickness is not positive
	real Value(const std::vector<real>& V) const
	{
		const PairList& P = G_.pairs;
		real E = 0.;
		for (int i = 0; i < G_.N; i++) {
			if (!(V[i] > 0.)) return std::numeric_limits<real>::infinity();
			if (Reservoir(i)) { E += (*res_)[i] * V[i]; continue; }
			E += G_.a[i] * p_.disj.P(V[i] / G_.a[i]);
			if (p_.rho_g != 0.) E += p_.rho_g * V[i] * G_.Height(i);
		}
		real Es = 0.;
		for (int q = 0; q < P.Size(); q++) {
			if (Reservoir(P.pi[q]) || Reservoir(P.pj[q])) continue;
			const real d = V[P.pi[q]] / G_.a[P.pi[q]] - V[P.pj[q]] / G_.a[P.pj[q]]; Es += G_.c[q] * d * d;
		}
		return E + p_.sigma / 8. * Es;
	}
	// p_i = dE / dV_i = P'(h_i) + (sigma / 4) sum_q c_q (h_i - h_j) / a_i + rho_g z_i
	void Pressure(const std::vector<real>& V, std::vector<real>& pr) const
	{
		const PairList& P = G_.pairs;
		pr.resize(G_.N);
		Parallel_For(G_.N, [&](int i) {
			if (Reservoir(i)) { pr[i] = (*res_)[i]; return; }
			const real hi = V[i] / G_.a[i];
			real s = 0.;
			for (int e = P.inc_off[i]; e < P.inc_off[i + 1]; e++) {
				const int q = P.inc_pair[e], o = P.Other(e, i);
				if (!Reservoir(o)) s += G_.c[q] * (hi - V[o] / G_.a[o]);
			}
			pr[i] = p_.disj.dP(hi) + p_.sigma / 4. * s / G_.a[i] + (p_.rho_g != 0. ? p_.rho_g * G_.Height(i) : 0.);
		});
	}
	// (d2E / dV2) w at volumes V
	void Hess_Apply(const std::vector<real>& V, const std::vector<real>& w, std::vector<real>& out) const
	{
		const PairList& P = G_.pairs;
		out.resize(G_.N);
		Parallel_For(G_.N, [&](int i) {
			if (Reservoir(i)) { out[i] = 0.; return; }
			real s = 0.;
			for (int e = P.inc_off[i]; e < P.inc_off[i + 1]; e++) {
				const int q = P.inc_pair[e], o = P.Other(e, i);
				if (!Reservoir(o)) s += G_.c[q] * (w[i] / G_.a[i] - w[o] / G_.a[o]);
			}
			out[i] = p_.disj.d2P(V[i] / G_.a[i]) * w[i] / G_.a[i] + p_.sigma / 4. * s / G_.a[i];
		});
	}
	// diagonal of d2E / dV2
	real Hess_Diag(const std::vector<real>& V, int i) const
	{
		if (Reservoir(i)) return 0.;
		const PairList& P = G_.pairs;
		real s = 0.;
		for (int e = P.inc_off[i]; e < P.inc_off[i + 1]; e++) if (!Reservoir(P.Other(e, i))) s += G_.c[P.inc_pair[e]];
		return p_.disj.d2P(V[i] / G_.a[i]) / G_.a[i] + p_.sigma / 4. * s / (G_.a[i] * G_.a[i]);
	}
	// d2E / dV2 assembled: the disjoining diagonal plus the squeeze Laplacian (pair weights times pair_scale if given;
	// pairs with scale 0 left out: the sparsified preconditioner)
	void Hess_Triplets(const std::vector<real>& V, std::vector<Eigen::Triplet<real>>& t, const std::vector<real>* pair_scale = nullptr) const
	{
		const PairList& P = G_.pairs;
		for (int q = 0; q < P.Size(); q++) {
			const int i = P.pi[q], j = P.pj[q];
			const real ps = pair_scale ? (*pair_scale)[q] : 1.;
			if (ps == 0. || Reservoir(i) || Reservoir(j)) continue;
			const real cs = p_.sigma / 4. * G_.c[q] * ps;
			t.emplace_back(i, i, cs / (G_.a[i] * G_.a[i])); t.emplace_back(j, j, cs / (G_.a[j] * G_.a[j]));
			t.emplace_back(i, j, -cs / (G_.a[i] * G_.a[j])); t.emplace_back(j, i, -cs / (G_.a[i] * G_.a[j]));
		}
		for (int i = 0; i < G_.N; i++) if (!Reservoir(i)) t.emplace_back(i, i, p_.disj.d2P(V[i] / G_.a[i]) / G_.a[i]);
	}
	const ThinFilmGraph& Graph() const { return G_; }

private:
	const ThinFilmGraph& G_;
	const ThinFilmParams& p_;
	const std::vector<real>* res_ = nullptr;
};

// One backward-Euler step of a conserved quantity V on a graph (gradient flow in the H^-1 metric of the weighted graph
// Laplacian L, weights k_p): J(phi) = E(V0 - L phi) + phi^T L phi / (2 dt). Energy concept: Value(V) (+inf if
// inadmissible), Pressure(V, p) (p = dE/dV), Hess_Apply(V, w, out), Hess_Diag(V, i), Hess_Triplets(V, t) (assembled
// d2E/dV2, for the optional sparse preconditioner), Graph(). Drainage: V the liquid volumes; surface diffusion: V the
// surfactant amounts (ThinFilmMechanics.h).
template <class Energy>
class H1_Flux_Problem : public Problem {
public:
	H1_Flux_Problem(const Energy& E, const std::vector<real>& V0, const std::vector<real>& k, real dt)
		: E_(E), G_(E.Graph()), V0_(V0), k_(k), dt_(dt) {}
	int Dim() const override { return G_.N; }
	// out = L v
	void Lap(const VecX& v, std::vector<real>& out) const
	{
		const PairList& P = G_.pairs;
		out.resize(G_.N);
		Parallel_For(G_.N, [&](int i) {
			real s = 0.;
			for (int e = P.inc_off[i]; e < P.inc_off[i + 1]; e++) s += k_[P.inc_pair[e]] * (v[i] - v[P.Other(e, i)]);
			out[i] = s;
		});
	}
	void Volumes(const VecX& phi, std::vector<real>& V) const
	{
		Lap(phi, lp_);
		V.resize(G_.N);
		for (int i = 0; i < G_.N; i++) V[i] = V0_[i] - lp_[i];
	}
	real Value(const VecX& phi) override
	{
		Volumes(phi, V_);
		const real e = E_.Value(V_);
		if (!std::isfinite(e)) return e;
		real d = 0.;
		for (int i = 0; i < G_.N; i++) d += phi[i] * lp_[i];
		return e + 0.5 * d / dt_;
	}
	void Gradient(const VecX& phi, VecX& g) override
	{
		Volumes(phi, V_);
		E_.Pressure(V_, pr_);
		VecX t(G_.N);
		for (int i = 0; i < G_.N; i++) t[i] = phi[i] / dt_ - pr_[i];
		Lap(t, lt_);
		g.resize(G_.N);
		for (int i = 0; i < G_.N; i++) g[i] = lt_[i];
	}
	void Prepare_Hessian(const VecX& phi) override { Volumes(phi, Vh_); }
	// L d2E L v + L v / dt
	void Hess_Vec(const VecX& v, VecX& Hv) override
	{
		Lap(v, w_);
		E_.Hess_Apply(Vh_, w_, aw_);
		VecX t(G_.N);
		for (int i = 0; i < G_.N; i++) t[i] = aw_[i] + v[i] / dt_;
		Lap(t, lt_);
		Hv.resize(G_.N);
		for (int i = 0; i < G_.N; i++) Hv[i] = lt_[i];
	}
	// Optional sparse preconditioner: the exact Hessian L A L + L / dt assembled (A = d2E / dV2), with a tiny shift for
	// the constant mode of L. Exact and cheap on a 1D grid (7 diagonals); on particles the product couples neighbours
	// of neighbours of neighbours (dense rows), so it is off by default.
	bool sparse_precond = false;
	bool Preconditioner_Matrix(SpMat& M) override
	{
		if (!sparse_precond) return false;
		const PairList& P = G_.pairs;
		const int n = G_.N;
		std::vector<Eigen::Triplet<real>> tl, ta;
		for (int q = 0; q < P.Size(); q++) {
			const int i = P.pi[q], j = P.pj[q];
			tl.emplace_back(i, i, k_[q]); tl.emplace_back(j, j, k_[q]); tl.emplace_back(i, j, -k_[q]); tl.emplace_back(j, i, -k_[q]);
		}
		E_.Hess_Triplets(Vh_, ta);
		SpMat L(n, n), A(n, n);
		L.setFromTriplets(tl.begin(), tl.end());
		A.setFromTriplets(ta.begin(), ta.end());
		SpMat H = SpMat(L * A * L) + L / dt_;
		real dmean = 0.;
		for (int i = 0; i < n; i++) dmean += std::abs(H.coeff(i, i)) / n;
		SpMat S(n, n);
		S.setIdentity();
		M = H + (1e-10 * dmean) * S;
		return true;
	}
	// Optional factored preconditioner: the exact inverse of H = L A L + L / dt = L (A L + I / dt) from two sparse
	// factorizations with one-ring stencils, for 2D particle sets where H itself has dense rows:
	//   L u = r                    (u up to a constant: r has zero sum; L + tiny shift, Cholesky, once per step),
	//   (A L + I / dt) z = u       (as the block system [L, -I; I / dt, A] [z; w] = [0; u], sparse LU per Newton
	//                               iteration since A = d2E/dV2 changes with V).
	// Any constant in u only adds a constant to z ((A L + I / dt) 1 = 1 / dt), which H ignores.
	// With precond_radius > 0 both L and the pair part of A are built from the pairs shorter than that radius only,
	// rescaled so that each particle keeps its second moment sum_j c_ij r_ij^2 (the same continuum Laplacian): an
	// approximate inverse with a short stencil, whose factorizations stay cheap where the kernel's ~58 neighbours fill
	// in (2D particle sets); CG makes up the difference.
	bool factored_precond = false;
	real precond_radius = 0.;
	// Optional two-level preconditioner (MELP-like fine/coarse particle sets): prolong (n x nc, each row a partition of
	// unity over nearby coarse particles) maps coarse to fine. One symmetric V-cycle: Chebyshev smoothing of degree
	// cheb_degree (on D^-1 H, D the Jacobi diagonal), the Galerkin coarse correction P (P^T H P)^-1 P^T, smoothing again.
	// P^T H P = P^T L (A (L P)) + P^T (L P) / dt is assembled right to left, so the fine 3-ring matrix is never formed.
	const SpMat* prolong = nullptr;
	int cheb_degree = 3;
	bool Prepare_Preconditioner() override
	{
		if (prolong) return Prepare_Two_Level();
		if (!factored_precond) return false;
		const PairList& P = G_.pairs;
		const int n = G_.N, np = P.Size();
		if (ps_.empty()) {
			ps_.assign(np, 1.);
			if (precond_radius > 0.) {
				std::vector<real> m(n, 0.), mk(n, 0.), r2(np);
				for (int q = 0; q < np; q++) {
					r2[q] = (G_.x[P.pj[q]] - G_.x[P.pi[q]]).squaredNorm();
					const real w = G_.c[q] * r2[q];
					m[P.pi[q]] += w; m[P.pj[q]] += w;
					if (r2[q] < precond_radius * precond_radius) { mk[P.pi[q]] += w; mk[P.pj[q]] += w; }
				}
				for (int q = 0; q < np; q++) {
					const int i = P.pi[q], j = P.pj[q];
					ps_[q] = (r2[q] < precond_radius * precond_radius && mk[i] > 0. && mk[j] > 0.) ? 0.5 * (m[i] / mk[i] + m[j] / mk[j]) : 0.;
				}
			}
		}
		std::vector<Eigen::Triplet<real>> tl, tb;
		for (int q = 0; q < np; q++) {
			if (ps_[q] == 0.) continue;
			const int i = P.pi[q], j = P.pj[q];
			const real kq = k_[q] * ps_[q];
			tl.emplace_back(i, i, kq); tl.emplace_back(j, j, kq); tl.emplace_back(i, j, -kq); tl.emplace_back(j, i, -kq);
		}
		if (!lap_ready_) {
			SpMat L(n, n);
			L.setFromTriplets(tl.begin(), tl.end());
			real dmean = 0.;
			for (int i = 0; i < n; i++) dmean += L.coeff(i, i) / n;
			std::vector<Eigen::Triplet<real>> ts = tl;
			for (int i = 0; i < n; i++) ts.emplace_back(i, i, 1e-10 * dmean);
			SpMat Ls(n, n);
			Ls.setFromTriplets(ts.begin(), ts.end());
			lap_.compute(Ls);
			if (lap_.info() != Eigen::Success) return false;
			lap_ready_ = true;
		}
		tb = tl;
		std::vector<Eigen::Triplet<real>> ta;
		E_.Hess_Triplets(Vh_, ta, &ps_);
		for (const auto& t : ta) tb.emplace_back(n + t.row(), n + t.col(), t.value());
		for (int i = 0; i < n; i++) { tb.emplace_back(i, n + i, -1.); tb.emplace_back(n + i, i, 1. / dt_); }
		SpMat B(2 * n, 2 * n);
		B.setFromTriplets(tb.begin(), tb.end());
		if (!block_analyzed_) { block_.analyzePattern(B); block_analyzed_ = true; }
		block_.factorize(B);
		return block_.info() == Eigen::Success;
	}
	void Apply_Preconditioner(const VecX& r, VecX& z) override
	{
		if (prolong) { Apply_Two_Level(r, z); return; }
		const int n = G_.N;
		VecX r0 = r.array() - r.mean();
		const VecX u = lap_.solve(r0);
		VecX rhs = VecX::Zero(2 * n);
		rhs.tail(n) = u;
		const VecX s = block_.solve(rhs);
		z = s.head(n);
	}
	// diag(L A L) with A replaced by its diagonal, plus diag(L) / dt
	void Hess_Diag(VecX& d) override
	{
		const PairList& P = G_.pairs;
		std::vector<real> A(G_.N);
		for (int i = 0; i < G_.N; i++) A[i] = E_.Hess_Diag(Vh_, i);
		d.resize(G_.N);
		Parallel_For(G_.N, [&](int i) {
			real Lii = 0., s = 0.;
			for (int e = P.inc_off[i]; e < P.inc_off[i + 1]; e++) { const real kq = k_[P.inc_pair[e]]; Lii += kq; s += kq * kq * A[P.Other(e, i)]; }
			d[i] = std::max(Lii / dt_, Lii * Lii * A[i] + s + Lii / dt_);
		});
	}

	int two_level_cycles = 0;                              // V-cycles applied (diagnostics)
	double t_galerkin = 0., t_factor = 0., t_power = 0.; // setup times of the two-level preconditioner (s)

private:
	static double Since(std::chrono::steady_clock::time_point a) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count(); }
	bool Prepare_Two_Level()
	{
		const PairList& P = G_.pairs;
		const int n = G_.N, np = P.Size();
		auto tc = std::chrono::steady_clock::now();
		std::vector<Eigen::Triplet<real>> tl, ta;
		tl.reserve(4 * np);
		for (int q = 0; q < np; q++) {
			const int i = P.pi[q], j = P.pj[q];
			tl.emplace_back(i, i, k_[q]); tl.emplace_back(j, j, k_[q]); tl.emplace_back(i, j, -k_[q]); tl.emplace_back(j, i, -k_[q]);
		}
		E_.Hess_Triplets(Vh_, ta);
		SpMat L(n, n), A(n, n);
		L.setFromTriplets(tl.begin(), tl.end());
		A.setFromTriplets(ta.begin(), ta.end());
		const SpMat& Pr = *prolong;
		const SpMat LP = L * Pr;
		const SpMat ALP = A * LP;
		const SpMat LALP = L * ALP;
		SpMat Hc = SpMat(Pr.transpose() * LALP) + SpMat(Pr.transpose() * LP) / dt_;
		const int nc = (int)Hc.rows();
		real dmean = 0.;
		for (int c = 0; c < nc; c++) dmean += std::abs(Hc.coeff(c, c)) / nc;
		SpMat S(nc, nc);
		S.setIdentity();
		Hc += (1e-10 * dmean) * S;
		t_galerkin += Since(tc); tc = std::chrono::steady_clock::now();
		coarse_.compute(Hc);
		t_factor += Since(tc); tc = std::chrono::steady_clock::now();
		if (coarse_.info() != Eigen::Success) return false;
		// Jacobi diagonal and the largest eigenvalue of D^-1 H (power iteration) for the Chebyshev smoother
		Hess_Diag(dinv_);
		for (int i = 0; i < n; i++) dinv_[i] = dinv_[i] > 0. ? 1. / dinv_[i] : 0.;
		VecX v(n), Hv;
		for (int i = 0; i < n; i++) v[i] = std::sin(1.3 * i + 0.7); // generic start
		v.array() -= v.mean();
		real lam = 1.;
		for (int it = 0; it < 12; it++) {
			Hess_Vec(v, Hv);
			VecX w = dinv_.cwiseProduct(Hv);
			lam = w.norm() / std::max<real>(v.norm(), 1e-300);
			v = w / std::max<real>(w.norm(), 1e-300);
		}
		lmax_ = 1.2 * lam;
		lmin_ = lmax_ / 30.;
		t_power += Since(tc);
		return true;
	}
	// Chebyshev iteration for H x = b (preconditioned by D), cheb_degree steps from x
	void Chebyshev(const VecX& b, VecX& x)
	{
		const real theta = 0.5 * (lmax_ + lmin_), delta = 0.5 * (lmax_ - lmin_), sigma = theta / delta;
		VecX r, Hx, Hd;
		Hess_Vec(x, Hx);
		r = b - Hx;
		VecX d = dinv_.cwiseProduct(r) / theta;
		real rho = 1. / sigma;
		for (int k = 0; k < cheb_degree; k++) {
			x += d;
			Hess_Vec(d, Hd);
			r -= Hd;
			const real rho_new = 1. / (2. * sigma - rho);
			d = (rho_new * rho) * d + (2. * rho_new / delta) * dinv_.cwiseProduct(r);
			rho = rho_new;
		}
	}
	void Apply_Two_Level(const VecX& r, VecX& z)
	{
		two_level_cycles++;
		const SpMat& Pr = *prolong;
		const int n = G_.N;
		z = VecX::Zero(n);
		Chebyshev(r, z);
		VecX Hz;
		Hess_Vec(z, Hz);
		VecX rc = Pr.transpose() * (r - Hz);
		rc.array() -= rc.mean();
		const VecX ec = coarse_.solve(rc);
		z += Pr * ec;
		Chebyshev(r, z);
	}

	const Energy& E_;
	const ThinFilmGraph& G_;
	const std::vector<real>& V0_;
	const std::vector<real>& k_;
	real dt_;
	mutable std::vector<real> lp_, V_, pr_, lt_, w_, aw_, Vh_;
	Eigen::SimplicialLDLT<SpMat> lap_;
	Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> block_;
	bool lap_ready_ = false, block_analyzed_ = false;
	std::vector<real> ps_; // pair scales of the preconditioner (0: left out)
	Eigen::SimplicialLDLT<SpMat> coarse_;
	VecX dinv_;
	real lmax_ = 1., lmin_ = 0.03;
};

// Prolongation from coarse particles xc to fine particles xf (n_f x n_c): quintic weights within `radius` on the 3D
// (chord) distance, each row normalized to a partition of unity (constants are prolongated exactly); a fine particle
// with no coarse particle in range takes the nearest one.
inline SpMat Make_Prolongation(const std::vector<Vec3>& xf, const std::vector<Vec3>& xc, real radius)
{
	Quintic2D ker(radius);
	const real inv = 1. / radius;
	auto cell_of = [&](const Vec3& p) { return Eigen::Vector3i((int)std::floor(p[0] * inv), (int)std::floor(p[1] * inv), (int)std::floor(p[2] * inv)); };
	auto key_of = [](const Eigen::Vector3i& c) { const int64_t K = 1 << 20; return ((int64_t)(c[0] + K / 2) * K + (int64_t)(c[1] + K / 2)) * K + (int64_t)(c[2] + K / 2); };
	std::unordered_map<int64_t, std::vector<int>> grid;
	for (int c = 0; c < (int)xc.size(); c++) grid[key_of(cell_of(xc[c]))].push_back(c);
	const int nf = (int)xf.size();
	std::vector<std::vector<std::pair<int, real>>> rows(nf);
	Parallel_For(nf, [&](int i) {
		const Eigen::Vector3i c0 = cell_of(xf[i]);
		real sw = 0.;
		for (int dx = -1; dx <= 1; dx++) for (int dy = -1; dy <= 1; dy++) for (int dz = -1; dz <= 1; dz++) {
			auto it = grid.find(key_of(c0 + Eigen::Vector3i(dx, dy, dz)));
			if (it == grid.end()) continue;
			for (int c : it->second) {
				const real r = (xc[c] - xf[i]).norm();
				if (r < radius) { const real w = ker.W(r); rows[i].emplace_back(c, w); sw += w; }
			}
		}
		if (sw > 0.) { for (auto& e : rows[i]) e.second /= sw; }
		else {
			int best = 0; real bd = 1e300;
			for (int c = 0; c < (int)xc.size(); c++) { const real d = (xc[c] - xf[i]).squaredNorm(); if (d < bd) { bd = d; best = c; } }
			rows[i].emplace_back(best, 1.);
		}
	}, 64);
	std::vector<Eigen::Triplet<real>> t;
	for (int i = 0; i < nf; i++) for (const auto& e : rows[i]) t.emplace_back(i, e.first, e.second);
	SpMat P(nf, (int)xc.size());
	P.setFromTriplets(t.begin(), t.end());
	return P;
}

using Drainage_Problem = H1_Flux_Problem<ThinFilmEnergy>;

// Newton-CG solve of one H1_Flux_Problem step; V is overwritten with the new values. Stops when max|grad J| falls below
// rel_tol times its value at the start of the step, or when the Newton steps reach round-off (phi ~ dt p; the
// gradient's round-off floor rises with a stiff energy).
template <class Energy>
NewtonStats H1_Flux_Step(const Energy& E, std::vector<real>& V, const std::vector<real>& k, real dt, NewtonOptions o, real rel_tol, bool sparse_precond,
	bool factored_precond = false, real precond_radius = 0., const SpMat* prolong = nullptr, int cheb_degree = 3, std::vector<real>* flux = nullptr)
{
	const int n = E.Graph().N;
	const std::vector<real> V0 = V;
	H1_Flux_Problem<Energy> P(E, V0, k, dt);
	P.sparse_precond = sparse_precond;
	P.factored_precond = factored_precond;
	P.precond_radius = precond_radius;
	P.prolong = prolong;
	P.cheb_degree = cheb_degree;
	VecX phi = VecX::Zero(n), g0;
	P.Gradient(phi, g0);
	o.sparse_preconditioner = sparse_precond || factored_precond || prolong != nullptr;
	// the preconditioners follow the Hessian, which changes with V within the step (the two-level one too: while black
	// film forms, a stale one costs thousands of CG iterations, a rebuild about a second at 33k particles)
	o.precond_rebuild_every = o.sparse_preconditioner ? 1 : 0;
	{
		std::vector<real> pr;
		E.Pressure(V0, pr);
		const auto mm = std::minmax_element(pr.begin(), pr.end());
		o.step_tol = 1e-12 * dt * (*mm.second - *mm.first);
	}
	o.grad_tol = std::max<real>(std::numeric_limits<real>::min(), rel_tol * g0.lpNorm<Eigen::Infinity>());
	const NewtonStats st = Newton_CG_Minimize(P, phi, o);
	P.Volumes(phi, V);
	// the volume each pair passed in the step, from pi to pj (V = V0 - L phi): k_q (phi_i - phi_j)
	if (flux) {
		const PairList& Pl = E.Graph().pairs;
		flux->resize(Pl.Size());
		for (int q = 0; q < Pl.Size(); q++) (*flux)[q] = k[q] * (phi[Pl.pi[q]] - phi[Pl.pj[q]]);
	}
	return st;
}

struct Drainage_Stepper {
	ThinFilmParams prm;
	NewtonOptions newton;
	real rel_tol = 1e-9; // Newton stops when max|grad J| falls below rel_tol times its value at the start of the step (1e-10 is round-off)
	bool sparse_precond = false; // exact sparse Hessian as the CG preconditioner (1D grids)
	bool factored_precond = false; // factored inverse (2D particle sets)
	real precond_radius = 0.;      // > 0: the factored preconditioner on the pairs shorter than this (sparsified)
	const SpMat* prolong = nullptr; // two-level preconditioner: prolongation from a coarse particle set (MELP-like)
	int cheb_degree = 3;            // its Chebyshev smoothing degree
	const std::vector<real>* reservoir = nullptr; // per-node bath chemical potential, NaN off the bath (ThinFilmEnergy)

	// pair conductances k_p = M(mean h) c_p at the current volumes (lagged within the step)
	void Conductances(const ThinFilmGraph& G, const std::vector<real>& V, std::vector<real>& k) const
	{
		k.resize(G.pairs.Size());
		for (int q = 0; q < G.pairs.Size(); q++) {
			const int i = G.pairs.pi[q], j = G.pairs.pj[q];
			k[q] = prm.Mobility(0.5 * (V[i] / G.a[i] + V[j] / G.a[j])) * G.c[q];
		}
	}
	// flux (optional): the volume each pair passed, from pi to pj
	NewtonStats Step(const ThinFilmGraph& G, std::vector<real>& V, real dt, std::vector<real>* flux = nullptr) const
	{
		ThinFilmEnergy E(G, prm);
		E.Set_Reservoir(reservoir);
		std::vector<real> k;
		Conductances(G, V, k);
		return H1_Flux_Step(E, V, k, dt, newton, rel_tol, sparse_precond, factored_precond, precond_radius, prolong, cheb_degree, flux);
	}
};

} // namespace vf::thinfilm
