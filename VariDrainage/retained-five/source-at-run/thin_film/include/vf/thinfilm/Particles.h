//////////////////////////////////////////////////////////////////////////
// Particle helpers of the thin-film module: the quintic kernel, neighbor pairs, a Fibonacci sphere.
// Copies of film_particle's Kernel.h / SurfaceGeometry.h pieces, so that this module depends only on core.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include "vf/core/Parallel.h"
#include "vf/core/Types.h"

namespace vf::thinfilm {

// 2D smoothing kernels, support radius h (W(r) = 0 for r >= h, normalized: integral over the plane = 1).
// - Quintic (default): the quintic B-spline (MELP's UnitQUINTIC), W = a sum_k c_k (k - 3 r/h)_+^5, k = 3,2,1,
//   c = 1,-6,15, a = 63/(478 pi h^2);
// - WendlandC2 / WendlandC4: W = a psi(r/h), psi = (1-q)^4 (1+4q), a = 7/(pi h^2) / psi = (1-q)^6 (35q^2+18q+3),
//   a = 3/(pi h^2).
// Pairing: a particle set whose energy depends on kernel densities and that carries a positive pressure is unstable to
// particle pairing when the kernel's Fourier transform has negative parts and it has many neighbours (Dehnen & Aly,
// MNRAS 425, 1068 (2012); Wendland functions have a non-negative transform). Measured here (thin_film/README.md,
// tf_bench pcheck) with a physical surfactant pressure 2 Pi_s ~ 0.04 N/m on a perfect hexagonal lattice: the quintic
// pairs at h = 4 dx, not at 3 dx; Wendland C4 is marginal at 4 dx. Either way a lattice under that pressure is a
// crystal whose shear stiffness exceeds the film's (viscous) one by 10-70x at ms steps; the cure is the reference
// pressure (FilmMechParams::bath_pressure), not the kernel.
class Kernel2D {
public:
	enum Type { Quintic = 0, WendlandC2 = 1, WendlandC4 = 2 };
	explicit Kernel2D(real h = 1., Type t = Quintic) : type_(t) { Set_Radius(h); }
	void Set_Type(Type t) { type_ = t; Set_Radius(h_); }
	Type Kind() const { return type_; }
	void Set_Radius(real h)
	{
		h_ = h;
		a_ = (type_ == Quintic ? 63. / (478. * kPi) : type_ == WendlandC2 ? 7. / kPi : 3. / kPi) / (h * h);
	}
	real Radius() const { return h_; }

	real W(real r) const
	{
		if (type_ != Quintic) {
			const real q = r / h_;
			if (!(q < 1.)) return 0.;
			const real u = 1. - q, u2 = u * u;
			return type_ == WendlandC2 ? a_ * u2 * u2 * (1. + 4. * q) : a_ * u2 * u2 * u2 * (35. * q * q + 18. * q + 3.);
		}
		const real q3 = 3. * r / h_;
		real s = 0.;
		for (int m = 0; m < 3; m++) {
			const real u = kK[m] - q3;
			if (u > 0) { const real u2 = u * u; s += kC[m] * u2 * u2 * u; }
		}
		return a_ * s;
	}
	// dW/dr
	real dW(real r) const
	{
		if (type_ != Quintic) {
			const real q = r / h_;
			if (!(q < 1.)) return 0.;
			const real u = 1. - q, u2 = u * u;
			return (type_ == WendlandC2 ? -20. * q * u2 * u : -56. * q * u2 * u2 * u * (1. + 5. * q)) * a_ / h_;
		}
		const real q3 = 3. * r / h_;
		real s = 0.;
		for (int m = 0; m < 3; m++) {
			const real u = kK[m] - q3;
			if (u > 0) { const real u2 = u * u; s += kC[m] * u2 * u2; }
		}
		return a_ * s * (-15. / h_);
	}
	// d^2W/dr^2
	real d2W(real r) const
	{
		if (type_ != Quintic) {
			const real q = r / h_;
			if (!(q < 1.)) return 0.;
			const real u = 1. - q, u2 = u * u;
			return (type_ == WendlandC2 ? -20. * u2 * (1. - 4. * q) : -56. * u2 * u2 * (1. + 4. * q - 35. * q * q)) * a_ / (h_ * h_);
		}
		const real q3 = 3. * r / h_;
		real s = 0.;
		for (int m = 0; m < 3; m++) {
			const real u = kK[m] - q3;
			if (u > 0) s += kC[m] * u * u * u;
		}
		return a_ * s * (180. / (h_ * h_));
	}

private:
	static constexpr real kK[3] = { 3., 2., 1. };
	static constexpr real kC[3] = { 1., -6., 15. };
	Type type_ = Quintic;
	real h_ = 1., a_ = 1.;
};
using Quintic2D = Kernel2D; // the module's kernel type (named after its default)

// Unique pairs (i < j) within a radius, with the incidence list: for particle k, entries e in [off[k], off[k+1]) give
// the pair inc_pair[e] and whether k is its second particle (inc_side[e] = 1: k == pj).
struct PairList {
	std::vector<int> pi, pj;
	std::vector<int> inc_off, inc_pair;
	std::vector<char> inc_side;
	int Size() const { return (int)pi.size(); }
	int Other(int e, int k) const { const int p = inc_pair[e]; return inc_side[e] ? pi[p] : pj[p]; (void)k; }
	void Build_Incidence(int N)
	{
		inc_off.assign(N + 1, 0);
		for (int p = 0; p < Size(); p++) { inc_off[pi[p] + 1]++; inc_off[pj[p] + 1]++; }
		for (int k = 0; k < N; k++) inc_off[k + 1] += inc_off[k];
		inc_pair.resize(2 * Size()); inc_side.resize(2 * Size());
		std::vector<int> fill(inc_off.begin(), inc_off.end() - 1);
		for (int p = 0; p < Size(); p++) {
			inc_pair[fill[pi[p]]] = p; inc_side[fill[pi[p]]++] = 0;
			inc_pair[fill[pj[p]]] = p; inc_side[fill[pj[p]]++] = 1;
		}
	}
};

// All pairs closer than `radius` (cell grid; per-particle neighbor lists in parallel, then a deterministic order).
inline void Build_Pairs(const std::vector<Vec3>& x, real radius, PairList& pairs)
{
	pairs.pi.clear(); pairs.pj.clear();
	const int n = (int)x.size();
	const real inv = 1. / radius;
	auto cell_of = [&](const Vec3& p) { return Eigen::Vector3i((int)std::floor(p[0] * inv), (int)std::floor(p[1] * inv), (int)std::floor(p[2] * inv)); };
	auto key_of = [](const Eigen::Vector3i& c) {
		const int64_t K = 1 << 20;
		return ((int64_t)(c[0] + K / 2) * K + (int64_t)(c[1] + K / 2)) * K + (int64_t)(c[2] + K / 2);
	};
	std::unordered_map<int64_t, std::vector<int>> grid;
	grid.reserve(n);
	for (int i = 0; i < n; i++) grid[key_of(cell_of(x[i]))].push_back(i);
	const real r2 = radius * radius;
	std::vector<std::vector<int>> nb(n);
	Parallel_For(n, [&](int i) {
		const Eigen::Vector3i c = cell_of(x[i]);
		for (int dx = -1; dx <= 1; dx++) for (int dy = -1; dy <= 1; dy++) for (int dz = -1; dz <= 1; dz++) {
			auto it = grid.find(key_of(c + Eigen::Vector3i(dx, dy, dz)));
			if (it == grid.end()) continue;
			for (int j : it->second) if (j > i && (x[j] - x[i]).squaredNorm() < r2) nb[i].push_back(j);
		}
	}, 64);
	size_t total = 0;
	for (const auto& v : nb) total += v.size();
	pairs.pi.reserve(total); pairs.pj.reserve(total);
	for (int i = 0; i < n; i++) for (int j : nb[i]) { pairs.pi.push_back(i); pairs.pj.push_back(j); }
	pairs.Build_Incidence(n);
}

// N nearly uniform points on a sphere (golden-angle spiral).
inline std::vector<Vec3> Fibonacci_Sphere(int N, real R, const Vec3& center = Vec3::Zero())
{
	std::vector<Vec3> x(N);
	const real golden = kPi * (3. - std::sqrt(5.));
	for (int i = 0; i < N; i++) {
		const real z = 1. - 2. * (i + 0.5) / N;
		const real r = std::sqrt(std::max<real>(0., 1. - z * z));
		const real phi = golden * i;
		x[i] = center + R * Vec3(r * std::cos(phi), r * std::sin(phi), z);
	}
	return x;
}

} // namespace vf::thinfilm
