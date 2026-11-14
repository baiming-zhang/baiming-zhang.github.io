//////////////////////////////////////////////////////////////////////////
// Stage 2: mobile film surfaces (plug flow) on a flat film, thin_film/README.md.
//
// Material particles (positions x in the film plane, 2 DOFs each) carry a liquid volume V_i and a surfactant amount S_i
// per surface; the particle area a_i = 1 / N_i(x) (kernel density), so the thickness h_i = V_i N_i and the surface
// concentration Gamma_i = S_i N_i follow the particles' motion (plug flow: the liquid moves with the surfaces; the
// drainage relative to them is ThinFilm.h's step, applied in a split). One implicit (BDF2) step minimizes
//
//   E(x) = sum_i m_i / (2 hh^2) |x_i - xt_i|^2                         (inertia, m_i = rho V_i)
//        + sum_i [ 2 a_i g(Gamma_i) + a_i P(h_i) ]                     (surfactant on both surfaces, disjoining)
//        + (sigma_f / 8) sum_pairs c_p (h_i - h_j)^2                    (squeeze energy, weights frozen at the start)
//        + rho_g sum_i V_i (x_i . up)                                   (gravity in the plane of a vertical film)
//        + hh R((x - x_ref) / hh)                                       (Trouton viscosity, geometry frozen at the start)
//        + sum_i d_i . C_i d_i / (2 hh), d_i = x_i - x_ref,i            (friction against the frame, C_i = c_i I + e_i n_i n_i^T:
//                                                                         FilmMatState::drag, drag_n, drag_dir)
//
// - g(Gamma): the surfactant's free energy per surface without the clean tension (a uniform tension on material
//   particles is the kernel-scale tensile instability, docs/A2.md): Langmuir g = RT Gamma_inf [theta ln theta + (1 - theta)
//   ln(1 - theta)], theta = Gamma / Gamma_inf, surface pressure Gamma g' - g = -RT Gamma_inf ln(1 - theta) (ideal
//   RT Gamma for Gamma_inf -> infinity), film Gibbs elasticity 2 Gamma^2 g'' = 2 RT Gamma / (1 - theta);
// - R(v) = 2 sum_p mu hbar_p c_p ((v_j - v_i) . e_p)^2 + 1/2 sum_i mu a_i h_i (div_i v)^2 approximates
//   int mu h [ e(v):e(v) + (div v)^2 ] (Trouton's 4 mu h in extension): sum_p c_p ((v_j - v_i) . e_p)^2 ~ int [e:e / 2 +
//   (div v)^2 / 4] for the Brookshaw weights, the divergence term supplies the rest; div_i v = sum_j a_j (v_j - v_i) .
//   grad_i W_ij.
// Periodic in x and (optionally) in y; with walls in y, particles flagged `fixed` (a frame) are not DOFs.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <array>
#include <stdexcept>
#include <vector>
#include "vf/core/NewtonCG.h"
#include <Eigen/SparseLU>
#ifdef VF_HAVE_CHOLMOD
#include <Eigen/CholmodSupport>
#endif
#include "vf/core/Parallel.h"
#include "vf/thinfilm/Particles.h"
#include "vf/thinfilm/Solution.h"
#include "vf/thinfilm/ThinFilm.h"

namespace vf::thinfilm {

// Surfactant free energy per surface (without the clean tension): Langmuir, or ideal when Gamma_inf <= 0
struct Surfactant {
	real RT = 1., Gamma_inf = 0.;
	bool Valid(real G) const { return G > 0. && (Gamma_inf <= 0. || G < Gamma_inf); }
	real g(real G) const
	{
		if (Gamma_inf <= 0.) return RT * G * (std::log(G) - 1.);
		const real t = G / Gamma_inf;
		return RT * Gamma_inf * (t * std::log(t) + (1. - t) * std::log(1. - t));
	}
	real dg(real G) const { return Gamma_inf <= 0. ? RT * std::log(G) : RT * (std::log(G / Gamma_inf) - std::log(1. - G / Gamma_inf)); }
	real d2g(real G) const { return Gamma_inf <= 0. ? RT / G : RT / (G * (1. - G / Gamma_inf)); }
	real Pressure(real G) const { return G * dg(G) - g(G); } // = RT G (ideal), -RT Gamma_inf ln(1 - theta) (Langmuir)
	// the concentration with surface pressure Pi >= 0
	real Concentration(real Pi) const
	{
		Pi = std::max<real>(Pi, 0.);
		return Gamma_inf <= 0. ? Pi / RT : Gamma_inf * -std::expm1(-Pi / (RT * Gamma_inf));
	}
};

struct FilmMechParams {
	ThinFilmParams tf;          // sigma (squeeze), mu (viscosity), disjoining, rho_g (gravity along -up)
	Surfactant surf;
	real rho = 1.;              // liquid density (mass rho V)
	Vec3 up = Vec3(0., 1., 0.); // the film's up direction (in plane) for gravity
	bool bulk_viscosity = true; // the divergence term of the Trouton dissipation
	// Viscosity in the in-plane (Trouton) dissipation of the surfaces' motion, in code units; < 0: tf.mu. The same
	// liquid viscosity enters the drainage (h^3 / 12 mu) and the in-plane flow (4 mu h), but with thickness and
	// lateral units H and L the two code coefficients differ by (H / L)^2 (Units.h); equal only if H = L.
	real mu_s = -1.;
	real Mu_S() const { return mu_s >= 0. ? mu_s : tf.mu; }
	// Surface pressure Pi_m of a surfactant bath the film exchanges surface with (an open border): the energy is then
	// the grand-canonical E + 2 Pi_m sum_i a_i (the clean tension, left out of E, cancels against the bath's), so a
	// film at the bath's concentration is in equilibrium with it and a film stretched below it is under tension
	// 2 (Pi_m - Pi_s) (which carries its weight). 0: closed frames, no bath.
	real bath_pressure = 0.;
	// Soluble surfactant (Solution.h; off unless sol.On() and FilmMatState::B is set): each particle's amount
	// Q = 2 S + B is in local equilibrium between its two surfaces and its liquid, so its surfactant energy is
	// F(a) = 2 a g(theta) + V f(c) at the partition for its current area a (surf.RT must equal sol.RT, Gamma_inf 1).
	// Particles flagged FilmMatState::sbath have surfaces in exchange with a surfactant bath at sol.c_bath (a
	// meniscus's surface over its deep liquid): energy -2 Pi_bath a, their S and B left out.
	Solubility sol;
};

// Material state of a flat film
struct FilmMatState {
	std::vector<Vec3> x, v, x_prev, v_prev;
	std::vector<real> V, S;     // liquid volume, surfactant amount per surface
	std::vector<real> B;        // surfactant dissolved in the liquid (FilmMechParams::sol; empty: insoluble)
	std::vector<real> Vs;       // the liquid volume B is in equilibrium with where it is not V (> 0; empty or 0: V), e.g. a
	                            // meniscus surface over the meniscus's liquid (Meniscus.h)
	real Vsol(int i) const { return (int)Vs.size() == Size() && Vs[i] > 0. ? Vs[i] : V[i]; }
	std::vector<char> sbath;    // surfaces in exchange with the surfactant bath FilmMechParams::sol.c_bath (empty: none)
	std::vector<char> fixed;    // frame particles (not DOFs)
	std::vector<real> gscale;   // per-particle weight factor (empty: 1; 0: carried by a bath, e.g. a meniscus)
	std::vector<Vec3> fext;     // per-particle constant external force (potential -fext_i . x_i; empty: none)
	std::vector<real> bath_w;   // per-particle factor of the bath term 2 Pi_m a_i (empty: 1), e.g. a bath pressure that
	                            // acts only near the open border (its body force then goes into fext)
	std::vector<real> drag;     // per-particle friction coefficient c_i against the (fixed) frame: dissipation 1/2 c_i |v_i|^2
	                            // (empty: none), e.g. a meniscus's surface dragging the meniscus's liquid (Meniscus.h)
	std::vector<real> drag_n;   // extra friction e_i along the unit direction drag_dir_i: dissipation 1/2 e_i (drag_dir_i . v_i)^2
	std::vector<Vec3> drag_dir; // (empty: none), e.g. a meniscus surface moving across its junction (Meniscus.h)
	real dt_prev = 0.;
	bool has_prev = false;
	int Size() const { return (int)x.size(); }
};

// Periodic box with optional walls in y; pairs within the kernel radius times (1 + margin) with frozen periodic shifts
struct FilmBox {
	real Lx = 1., Ly = 1.;
	bool periodic_y = true;
	Vec3 Wrap(Vec3 d) const { d[0] -= Lx * std::round(d[0] / Lx); if (periodic_y) d[1] -= Ly * std::round(d[1] / Ly); return d; }
};

struct ShiftedPairs {
	PairList pairs;
	std::vector<Vec3> shift; // minimum-image shift: d_p = x_j - x_i + shift_p
};

inline void Build_Shifted_Pairs(const std::vector<Vec3>& x, const FilmBox& box, real radius, ShiftedPairs& sp)
{
	const int n = (int)x.size();
	const int ncx = std::max(1, (int)std::floor(box.Lx / radius)), ncy = std::max(1, (int)std::floor(box.Ly / radius));
	auto cell = [&](const Vec3& p) {
		int cx = (int)std::floor(p[0] / box.Lx * ncx), cy = (int)std::floor(p[1] / box.Ly * ncy);
		cx = ((cx % ncx) + ncx) % ncx; cy = box.periodic_y ? ((cy % ncy) + ncy) % ncy : std::clamp(cy, 0, ncy - 1);
		return cy * ncx + cx;
	};
	std::vector<std::vector<int>> cells(ncx * ncy);
	for (int i = 0; i < n; i++) cells[cell(x[i])].push_back(i);
	std::vector<std::vector<int>> nb(n);
	Parallel_For(n, [&](int i) {
		const int c0 = cell(x[i]), cx = c0 % ncx, cy = c0 / ncx;
		std::vector<int> seen;
		for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
			if (!box.periodic_y && (cy + dy < 0 || cy + dy >= ncy)) continue;
			const int c = (((cy + dy) % ncy + ncy) % ncy) * ncx + (((cx + dx) % ncx + ncx) % ncx);
			if (std::find(seen.begin(), seen.end(), c) != seen.end()) continue;
			seen.push_back(c);
			for (int j : cells[c]) if (j > i && box.Wrap(x[j] - x[i]).norm() < radius) nb[i].push_back(j);
		}
		std::sort(nb[i].begin(), nb[i].end());
	}, 64);
	sp.pairs.pi.clear(); sp.pairs.pj.clear(); sp.shift.clear();
	for (int i = 0; i < n; i++) for (int j : nb[i]) {
		sp.pairs.pi.push_back(i); sp.pairs.pj.push_back(j);
		sp.shift.push_back(box.Wrap(x[j] - x[i]) - (x[j] - x[i]));
	}
	sp.pairs.Build_Incidence(n);
}

// A fixed spherical film shape (a bubble whose shape the thin-film dynamics does not change). Positions on it are
// parametrized by the tangent plane at a start point and the radial projection x(u) = c + R u / |u|, u = x0 - c + T q.
struct SphereShape {
	Vec3 c = Vec3::Zero();
	real R = 1.;
	Vec3 Normal(const Vec3& x) const { return (x - c).normalized(); }
	Vec3 Project(const Vec3& x) const { return c + R * (x - c).normalized(); }
	void Basis(const Vec3& x, Vec3& t1, Vec3& t2) const
	{
		const Vec3 n = Normal(x);
		const Vec3 e = std::abs(n[2]) < 0.9 ? Vec3(0., 0., 1.) : Vec3(1., 0., 0.);
		t1 = e.cross(n).normalized(); t2 = n.cross(t1);
	}
	// dx/du w = (R / |u|) (w - u^ (u^ . w))
	Vec3 Jac(const Vec3& u, const Vec3& w) const { const real s = u.norm(); const Vec3 uh = u / s; return (R / s) * (w - uh * uh.dot(w)); }
	// G . d2x/du2 [a, b] as a vector in a: -(R / |u|^2) [u^ (G . P b) + (u^ . b) P G + (G . u^) P b], P = I - u^ u^T
	Vec3 Curv(const Vec3& u, const Vec3& G, const Vec3& b) const
	{
		const real s = u.norm();
		const Vec3 uh = u / s, Pb = b - uh * uh.dot(b), PG = G - uh * uh.dot(G);
		return -(R / (s * s)) * (uh * G.dot(Pb) + uh.dot(b) * PG + G.dot(uh) * Pb);
	}
};

class FilmMech_Problem : public Problem {
public:
	// x0: positions at the start of the step (the DOFs are offsets from them); xt: inertia targets; xref: dissipation
	// reference (velocity = (x - xref) / hh); everything frozen (pair weights, Trouton geometry) is evaluated at x0.
	// b1, b2 (optional): per-particle tangent basis of the DOFs (a curved film: the tangent plane at x0); default the
	// x and y axes of a flat film. sph (optional): the DOFs move the particles on this sphere exactly (x = projection of
	// x0 + T q), so the constraint reaction is radial at the actual position (no spurious torque).
	FilmMech_Problem(const Quintic2D& ker, const FilmMechParams& prm, const ShiftedPairs& sp, const FilmMatState& st,
		const std::vector<Vec3>& xt, const std::vector<Vec3>& xref, real hh,
		const std::vector<Vec3>* b1 = nullptr, const std::vector<Vec3>* b2 = nullptr, const SphereShape* sph = nullptr)
		: ker_(ker), p_(prm), sp_(sp), st_(st), xt_(xt), xref_(xref), hh_(hh), N_(st.Size()), b1_(b1), b2_(b2), sph_(sph)
	{
		dof_.assign(N_, -1);
		for (int i = 0; i < N_; i++) if (st.fixed.empty() || !st.fixed[i]) { dof_[i] = (int)free_.size(); free_.push_back(i); }
		const PairList& P = sp_.pairs;
		const int np = P.Size();
		const real h = ker_.Radius();
		// frozen quantities at x0
		std::vector<real> nd(N_, ker_.W(0.));
		std::vector<real> r0(np);
		for (int q = 0; q < np; q++) {
			r0[q] = D0(q).norm();
			if (r0[q] < h) { const real w = ker_.W(r0[q]); nd[P.pi[q]] += w; nd[P.pj[q]] += w; }
		}
		a0_.resize(N_); h0_.resize(N_);
		for (int i = 0; i < N_; i++) { a0_[i] = 1. / nd[i]; h0_[i] = st_.V[i] * nd[i]; }
		csq_.assign(np, 0.); cT_.assign(np, 0.); e_.assign(np, Vec3::Zero()); gb_.assign(np, Vec3::Zero());
		for (int q = 0; q < np; q++) {
			if (!(r0[q] < h) || !(r0[q] > 1e-12 * h)) continue;
			const int i = P.pi[q], j = P.pj[q];
			const real c = 2. * std::abs(ker_.dW(r0[q])) / (r0[q] * nd[i] * nd[j]);
			csq_[q] = c;
			cT_[q] = 2. * p_.Mu_S() * 0.5 * (h0_[i] + h0_[j]) * c;
			e_[q] = D0(q) / r0[q];
			gb_[q] = ker_.dW(r0[q]) * e_[q]; // dW(|x_j - x_i|)/dx_j
		}
		kb_.resize(N_);
		for (int i = 0; i < N_; i++) kb_[i] = p_.bulk_viscosity ? p_.Mu_S() * a0_[i] * h0_[i] : 0.;
		if (p_.sol.On()) {
			pib_ = p_.sol.Pi_Bath();
			lnc_.assign(N_, std::numeric_limits<real>::quiet_NaN());
		}
	}
	int Dim() const override { return 2 * (int)free_.size(); }
	const std::vector<int>& Free() const { return free_; }
	Vec3 B1(int i) const { return b1_ ? (*b1_)[i] : Vec3(1., 0., 0.); }
	Vec3 B2(int i) const { return b2_ ? (*b2_)[i] : Vec3(0., 1., 0.); }
	// positions at q; with a sphere also the unprojected points u = x0 - c + T q of the DOFs
	void To_X(const VecX& q, std::vector<Vec3>& x, std::vector<Vec3>* u = nullptr) const
	{
		x = st_.x;
		if (u) u->assign(N_, Vec3::Zero());
		for (size_t k = 0; k < free_.size(); k++) {
			const int i = free_[k];
			x[i] += q[2 * k] * B1(i) + q[2 * k + 1] * B2(i);
			if (sph_) { if (u) (*u)[i] = x[i] - sph_->c; x[i] = sph_->Project(x[i]); }
		}
	}

	real Value(const VecX& q) override
	{
		To_X(q, x_);
		Densities(x_, false);
		const real inv = 1. / (hh_ * hh_);
		real E = 0.;
		ei_.resize(N_);
		Parallel_For(N_, [&](int i) { ei_[i] = Ei(i, nd_[i]); }, 64);
		for (int i = 0; i < N_; i++) {
			const real e = ei_[i];
			if (!std::isfinite(e)) return e;
			E += e + p_.tf.rho_g * st_.V[i] * Weight(i) * p_.up.dot(x_[i]) - Fext(i).dot(x_[i]);
			if (dof_[i] >= 0) E += 0.5 * p_.rho * st_.V[i] * inv * (x_[i] - xt_[i]).squaredNorm() + 0.5 / hh_ * (x_[i] - xref_[i]).dot(Drag_F(i, x_[i] - xref_[i]));
		}
		E += Squeeze_Value();
		E += Visc_Value();
		return E;
	}
	void Gradient(const VecX& q, VecX& g) override
	{
		To_X(q, x_, sph_ ? &u_ : nullptr);
		Densities(x_, true);
		std::vector<Vec3> G;
		Coef_N(coef_);
		Full_Grad(x_, coef_, G);
		Gather(G, g, u_);
	}
	void Prepare_Hessian(const VecX& q) override
	{
		To_X(q, xh_, sph_ ? &uh_ : nullptr);
		Densities(xh_, true);
		Coef_N(coefh_);
		d2_.resize(N_);
		Parallel_For(N_, [&](int i) { d2_[i] = Ei_dd(i, nd_[i]); }, 64);
		ndh_ = nd_; gh_ = g_; Kd_ = d_; Kr_ = r_;
		if (sph_) Full_Grad(xh_, coefh_, Gh_); // the curvature term of the projection needs the x-gradient
	}
	void Hess_Vec(const VecX& vq, VecX& Hv) override
	{
		std::vector<Vec3> v(N_, Vec3::Zero()), H(N_, Vec3::Zero()), b(N_, Vec3::Zero());
		for (size_t k = 0; k < free_.size(); k++) {
			const int i = free_[k];
			b[i] = vq[2 * k] * B1(i) + vq[2 * k + 1] * B2(i);
			v[i] = sph_ ? sph_->Jac(uh_[i], b[i]) : b[i];
		}
		const PairList& P = sp_.pairs;
		const int np = P.Size();
		// s_i = grad N_i . v
		std::vector<real> s(N_, 0.), t(N_);
		for (int q = 0; q < np; q++) { const real dq = gh_[q].dot(v[P.pj[q]] - v[P.pi[q]]); s[P.pi[q]] += dq; s[P.pj[q]] += dq; }
		// t = (d2E/dN2) s: per-particle part and the squeeze Laplacian
		for (int i = 0; i < N_; i++) t[i] = d2_[i] * s[i];
		const real ks = p_.tf.sigma / 4.;
		for (int q = 0; q < np; q++) {
			if (csq_[q] == 0.) continue;
			const int i = P.pi[q], j = P.pj[q];
			const real w = ks * csq_[q] * (st_.V[i] * s[i] - st_.V[j] * s[j]);
			t[i] += st_.V[i] * w; t[j] -= st_.V[j] * w;
		}
		for (int q = 0; q < np; q++) {
			const int i = P.pi[q], j = P.pj[q];
			if (!(Kr_[q] < ker_.Radius())) continue;
			// J^T t and the grad^2 N term with the first-derivative coefficients
			const Vec3 dv = v[j] - v[i];
			Vec3 f = (t[i] + t[j]) * gh_[q];
			const real r = Kr_[q], w1 = ker_.dW(r), w2 = ker_.d2W(r);
			const Vec3 e = Kd_[q] / r;
			const Vec3 Kv = r > 1e-12 * ker_.Radius() ? Vec3(w2 * e * e.dot(dv) + (w1 / r) * (dv - e * e.dot(dv))) : Vec3(ker_.d2W(0.) * dv);
			f += (coefh_[i] + coefh_[j]) * Kv;
			H[j] += f; H[i] -= f;
		}
		const real inv = 1. / (hh_ * hh_);
		for (int i = 0; i < N_; i++) if (dof_[i] >= 0) H[i] += p_.rho * st_.V[i] * inv * v[i] + Drag_F(i, v[i]) / hh_;
		Add_Visc(v, H);
		if (sph_) {
			// chain rule through the projection: J^T H J b + G . d2x/du2 [., b]
			for (int i : free_) H[i] = sph_->Jac(uh_[i], H[i]) + sph_->Curv(uh_[i], Gh_[i], b[i]);
			Gather(H, Hv);
			return;
		}
		Gather(H, Hv);
	}
	void Hess_Diag(VecX& dg) override
	{
		const PairList& P = sp_.pairs;
		std::vector<Mat3> B(N_, Mat3::Zero());
		const real inv = 1. / (hh_ * hh_);
		for (int i = 0; i < N_; i++) B[i] += p_.rho * st_.V[i] * inv * Mat3::Identity() + Drag_M(i) / hh_;
		for (int q = 0; q < P.Size(); q++) {
			const int i = P.pi[q], j = P.pj[q];
			const Mat3 gg = gh_[q] * gh_[q].transpose();
			B[i] += (std::max<real>(0., d2_[i]) + std::max<real>(0., d2_[j])) * gg;
			B[j] += (std::max<real>(0., d2_[i]) + std::max<real>(0., d2_[j])) * gg;
			const Mat3 ee = e_[q] * e_[q].transpose();
			B[i] += 2. * cT_[q] / hh_ * ee; B[j] += 2. * cT_[q] / hh_ * ee;
		}
		dg.resize(Dim());
		for (size_t k = 0; k < free_.size(); k++) {
			const int i = free_[k];
			const Vec3 c1 = sph_ ? sph_->Jac(uh_[i], B1(i)) : B1(i), c2 = sph_ ? sph_->Jac(uh_[i], B2(i)) : B2(i);
			dg[2 * k] = std::max<real>(c1.dot(B[i] * c1), 1e-300); dg[2 * k + 1] = std::max<real>(c2.dot(B[i] * c2), 1e-300);
		}
	}

	// Factored preconditioner (factored_precond). The Hessian is
	//   H = A + J^T K J + D^T Kb D,
	// A: inertia, the pairwise Trouton term and the stress term sum_i dE/dN_i grad^2 N_i (1-ring); J = dN/dq, the
	// densities' response (1-ring); K = d2E/dN2, the surfactant's and the disjoining stiffness (diagonal) and the squeeze
	// term (a Laplacian on h); D the Trouton divergence (1-ring), Kb = kb / hh. Assembled as it stands it would couple
	// neighbours of neighbours (~230 particles with the 4 dx kernel). With physical surfactant (Gibbs elasticity
	// ~0.1 N/m) the compression stiffness J^T K J exceeds the shear stiffness (Trouton 4 mu h / dt, inertia) by 1e4-1e6
	// at ms steps: the surfaces are nearly incompressible, the Jacobi-preconditioned CG stalls, and no short-stencil
	// approximation of J serves (its divergence-free modes differ from J's, which the preconditioner then makes 1e5 too
	// stiff). So the preconditioner keeps J exactly, as the 1-ring saddle-point matrix
	//   [ A   J^T      D^T     ]
	//   [ J   -K^-1    0       ]
	//   [ D   0        -Kb^-1  ]
	// whose Schur complement on the first block is H: solving it with right-hand side (r, 0, 0) applies H^-1 r. K is
	// taken diagonal (the squeeze Laplacian by its diagonal: exact without squeeze energy, and the squeeze is
	// (H / L)^2-small in physical units) and floored at a small positive value (where the energy is not convex in N, as
	// in the disjoining spinodal, that particle's compression stiffness is left out). LDL^T (no pivoting) with AMD
	// ordering; the matrix has as many negative pivots as auxiliary rows exactly when the preconditioner is positive
	// definite (Sylvester): otherwise the stress term is projected onto its positive part per pair and the matrix
	// refactorized; failing that, Jacobi.
	// Built with CHOLMOD (VF_HAVE_CHOLMOD, precond_cholmod): H itself (the same parts, K and Kb diagonal) is assembled
	// - neighbours of neighbours, but symmetric positive definite - and factorized by CHOLMOD's supernodal Cholesky
	// (fill-reducing ordering, multithreaded BLAS); a failed factorization means not positive definite. The
	// LDL^T above is the fallback (single-threaded, simplicial: 2.7 s against 0.5 s for 2304 particles, 4 dx kernel).
	bool factored_precond = false;
	bool precond_cholmod = true;
	int precond_stress_projected = 0; // factorizations that needed the stress term projected (diagnostics)
	bool Prepare_Preconditioner() override
	{
		if (!factored_precond) return false;
		if (Assemble_Precond(false)) return true;
		precond_stress_projected++;
		return Assemble_Precond(true);
	}
	// diagnostics: factorize with or without the stress projection; negative pivots and auxiliary rows of the last one
	bool Factor_Precond(bool project_stress) { return Assemble_Precond(project_stress); }
	int precond_negatives = 0, precond_aux = 0;
	double precond_t_assemble = 0., precond_t_analyze = 0., precond_t_factor = 0.; // seconds
	long long precond_nnz_A = 0, precond_nnz_L = 0;     // lower triangle of the matrix, of its factor
	void Apply_Preconditioner(const VecX& r, VecX& z) override
	{
#ifdef VF_HAVE_CHOLMOD
		if (pc_use_chol_) { z = pc_chol_.solve(r); return; }
#endif
		VecX rhs = VecX::Zero(pc_n_);
		rhs.head(Dim()) = r;
		z = pc_ldlt_.solve(rhs).head(Dim());
	}

	// per-particle energy in terms of N (a = 1/N): surfactant on both surfaces (insoluble; soluble, in local
	// equilibrium with the particle's liquid; or a surfactant bath) and disjoining
	real Ei(int i, real N) const
	{
		const real h = st_.V[i] * N;
		if (!(h > 0.)) return std::numeric_limits<real>::infinity();
		return Surf_E(i, N, 0) + p_.tf.disj.P(h) / N + 2. * Bath(i) / N;
	}
	real Ei_d(int i, real N) const
	{
		const real h = st_.V[i] * N;
		return Surf_E(i, N, 1) + (h * p_.tf.disj.dP(h) - p_.tf.disj.P(h)) / (N * N) - 2. * Bath(i) / (N * N);
	}
	real Ei_dd(int i, real N) const
	{
		const real h = st_.V[i] * N;
		return Surf_E(i, N, 2) + (h * h * p_.tf.disj.d2P(h) - 2. * (h * p_.tf.disj.dP(h) - p_.tf.disj.P(h))) / (N * N * N) + 4. * Bath(i) / (N * N * N);
	}
	// the surfactant's part of Ei and its first and second derivatives in N
	real Surf_E(int i, real N, int order) const
	{
		if (SBath(i)) return order == 0 ? -2. * pib_ / N : order == 1 ? 2. * pib_ / (N * N) : -4. * pib_ / (N * N * N);
		if (Soluble(i)) {
			// F(a) at fixed V, Q: E(N) = F(1/N), E' = -F_a / N^2, E'' = F_aa / N^4 + 2 F_a / N^3
			real F, Fa, Faa;
			p_.sol.Energy(1. / N, st_.Vsol(i), 2. * st_.S[i] + st_.B[i], F, Fa, Faa, &lnc_[i]);
			return order == 0 ? F : order == 1 ? -Fa / (N * N) : Faa / (N * N * N * N) + 2. * Fa / (N * N * N);
		}
		if (!(st_.S[i] > 0.)) return 0.;
		const real G = st_.S[i] * N;
		if (order == 0) return p_.surf.Valid(G) ? 2. * p_.surf.g(G) / N : std::numeric_limits<real>::infinity();
		if (order == 1) return 2. * p_.surf.Pressure(G) / (N * N);
		return 2. * (G * G * p_.surf.d2g(G) - 2. * p_.surf.Pressure(G)) / (N * N * N);
	}
	bool SBath(int i) const { return p_.sol.On() && (int)st_.sbath.size() == N_ && st_.sbath[i]; }
	bool Soluble(int i) const { return p_.sol.On() && (int)st_.B.size() == N_ && 2. * st_.S[i] + st_.B[i] > 0. && st_.Vsol(i) > 0.; }
	real Weight(int i) const { return st_.gscale.empty() ? 1. : st_.gscale[i]; }
	Vec3 Fext(int i) const { return (int)st_.fext.size() == N_ ? st_.fext[i] : Vec3::Zero(); }
	real Bath(int i) const { return p_.bath_pressure * ((int)st_.bath_w.size() == N_ ? st_.bath_w[i] : 1.); }
	real Drag(int i) const { return (int)st_.drag.size() == N_ ? st_.drag[i] : 0.; }
	// the friction tensor c_i I + e_i n_i n_i^T (drag, drag_n, drag_dir) applied to v, and as a matrix
	Vec3 Drag_F(int i, const Vec3& v) const
	{
		Vec3 f = Drag(i) * v;
		if ((int)st_.drag_n.size() == N_ && st_.drag_n[i] != 0.) f += st_.drag_n[i] * st_.drag_dir[i].dot(v) * st_.drag_dir[i];
		return f;
	}
	Mat3 Drag_M(int i) const
	{
		Mat3 M = Drag(i) * Mat3::Identity();
		if ((int)st_.drag_n.size() == N_ && st_.drag_n[i] != 0.) M += st_.drag_n[i] * st_.drag_dir[i] * st_.drag_dir[i].transpose();
		return M;
	}
	const std::vector<real>& Density() const { return nd_; }

private:
	// The factored preconditioner at the prepared state; false if it is not positive definite. The parts: A (lower
	// triangle, DOF x DOF), J and D (particle x DOF), the diagonal K and Kb. With CHOLMOD H = A + J^T K J + D^T Kb D is
	// assembled and factorized (supernodal Cholesky); otherwise the saddle-point matrix (LDL^T).
	bool Assemble_Precond(bool project_stress)
	{
		const PairList& P = sp_.pairs;
		const int np = P.Size(), nf = (int)free_.size(), nq = 2 * nf;
		const bool bulk = p_.bulk_viscosity && p_.Mu_S() > 0.;
		using Trip = Eigen::Triplet<real>;
		std::vector<Trip> tA, tJ, tD;
		tA.reserve((size_t)np * 12 + 4 * nq); tJ.reserve((size_t)np * 8); if (bulk) tD.reserve((size_t)np * 8);
		auto add = [&](int r, int c, real v) { if (r >= c && v != 0.) tA.emplace_back(r, c, v); }; // lower triangle of A
		auto C = [&](int i, int a) -> Vec3 { const Vec3 b = a == 0 ? B1(i) : B2(i); return sph_ ? sph_->Jac(uh_[i], b) : b; };
		// symmetric 3x3 M acting on v_j - v_i (energy 1/2 dv^T M dv)
		auto pair_block = [&](int i, int j, const Mat3& M) {
			const int ki = dof_[i], kj = dof_[j];
			for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) {
				if (ki >= 0) add(2 * ki + a, 2 * ki + b, C(i, a).dot(M * C(i, b)));
				if (kj >= 0) add(2 * kj + a, 2 * kj + b, C(j, a).dot(M * C(j, b)));
				if (ki >= 0 && kj >= 0) { const real v = -C(i, a).dot(M * C(j, b)); add(2 * ki + a, 2 * kj + b, v); add(2 * kj + b, 2 * ki + a, v); }
			}
		};
		// row `row` of a 1-ring operator gets c . (v_j - v_i) (no explicit zeros: pairs beyond the kernel radius)
		auto diff_row = [&](std::vector<Trip>& t, int row, int i, int j, const Vec3& c) {
			const int ki = dof_[i], kj = dof_[j];
			for (int a = 0; a < 2; a++) {
				if (kj >= 0 && c.dot(C(j, a)) != 0.) t.emplace_back(row, 2 * kj + a, c.dot(C(j, a)));
				if (ki >= 0 && c.dot(C(i, a)) != 0.) t.emplace_back(row, 2 * ki + a, -c.dot(C(i, a)));
			}
		};
		const real inv = 1. / (hh_ * hh_), R = ker_.Radius();
		for (size_t k = 0; k < free_.size(); k++) {
			const int i = free_[k];
			for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) {
				real v = p_.rho * st_.V[i] * inv * C(i, a).dot(C(i, b)) + C(i, a).dot(Drag_M(i) * C(i, b)) / hh_;
				if (sph_) v += 0.5 * (B1B2(i, a).dot(sph_->Curv(uh_[i], Gh_[i], B1B2(i, b))) + B1B2(i, b).dot(sph_->Curv(uh_[i], Gh_[i], B1B2(i, a))));
				add(2 * (int)k + a, 2 * (int)k + b, v);
			}
		}
		std::vector<real> sqd(N_, 0.);
		const real ks = p_.tf.sigma / 4.;
		for (int q = 0; q < np; q++) {
			const int i = P.pi[q], j = P.pj[q];
			if (csq_[q] != 0.) { sqd[i] += ks * csq_[q] * st_.V[i] * st_.V[i]; sqd[j] += ks * csq_[q] * st_.V[j] * st_.V[j]; }
			Mat3 M = Mat3::Zero();
			if (p_.Mu_S() > 0. && cT_[q] != 0.) M += 2. * cT_[q] / hh_ * e_[q] * e_[q].transpose();
			if (Kr_[q] < R) {
				const real c = coefh_[i] + coefh_[j], r = Kr_[q];
				if (r > 1e-12 * R) {
					const Vec3 e = Kd_[q] / r;
					real kpar = c * ker_.d2W(r), kperp = c * ker_.dW(r) / r;
					if (project_stress) { kpar = std::max<real>(kpar, 0.); kperp = std::max<real>(kperp, 0.); }
					M += kpar * e * e.transpose() + kperp * (Mat3::Identity() - e * e.transpose());
				}
				else M += (project_stress ? std::max<real>(0., c * ker_.d2W(0.)) : c * ker_.d2W(0.)) * Mat3::Identity();
				diff_row(tJ, i, i, j, gh_[q]); // s_i, s_j += gh . (v_j - v_i)
				diff_row(tJ, j, i, j, gh_[q]);
			}
			pair_block(i, j, M);
			if (bulk) { // div_i -= a0_j gb . (v_j - v_i), div_j -= a0_i gb . (v_j - v_i)
				diff_row(tD, i, i, j, -a0_[j] * gb_[q]);
				diff_row(tD, j, i, j, -a0_[i] * gb_[q]);
			}
		}
		// K: the per-particle stiffness and the squeeze term's diagonal, without negative values (non-convex in N);
		// Kb = kb / hh
		VecX Kc(N_), Kb(N_);
		for (int i = 0; i < N_; i++) { Kc[i] = std::max<real>(d2_[i] + sqd[i], 0.); Kb[i] = bulk ? kb_[i] / hh_ : 0.; }
		precond_aux = 0; precond_negatives = 0;
		const auto t0 = std::chrono::steady_clock::now();
		auto sec = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };
#ifdef VF_HAVE_CHOLMOD
		if (precond_cholmod) {
			SpMat A(nq, nq), J(N_, nq), H;
			A.setFromTriplets(tA.begin(), tA.end());
			J.setFromTriplets(tJ.begin(), tJ.end());
			const SpMat KJ = Kc.asDiagonal() * J;
			const SpMat JtKJ = SpMat(J.transpose() * KJ).triangularView<Eigen::Lower>();
			H = A + JtKJ;
			if (bulk) {
				SpMat D(N_, nq);
				D.setFromTriplets(tD.begin(), tD.end());
				const SpMat KD = Kb.asDiagonal() * D;
				const SpMat DtKD = SpMat(D.transpose() * KD).triangularView<Eigen::Lower>();
				H += DtKD;
			}
			const auto t1 = std::chrono::steady_clock::now();
			pc_chol_.cholmod().print = 0;
			pc_chol_.analyzePattern(H); // the ordering (patterns can differ between the exact and projected stress terms)
			const auto t2 = std::chrono::steady_clock::now();
			pc_chol_.factorize(H);
			const auto t3 = std::chrono::steady_clock::now();
			precond_t_assemble += sec(t0, t1); precond_t_analyze += sec(t1, t2); precond_t_factor += sec(t2, t3);
			precond_nnz_A = (long long)H.nonZeros(); precond_nnz_L = (long long)pc_chol_.cholmod().lnz;
			pc_use_chol_ = true;
			return pc_chol_.info() == Eigen::Success;
		}
#endif
		pc_use_chol_ = false;
		// the saddle-point matrix [A J^T D^T; J -K^-1 0; D 0 -Kb^-1] (lower triangle); K floored at a small positive value
		pc_n_ = nq + N_ + (bulk ? N_ : 0);
		const int oK = nq, oB = nq + N_;
		std::vector<Trip>& t = tA;
		for (const Trip& x : tJ) t.emplace_back(oK + x.row(), x.col(), x.value());
		for (const Trip& x : tD) t.emplace_back(oB + x.row(), x.col(), x.value());
		const real kfloor = 1e-12 * std::max<real>(Kc.maxCoeff(), std::numeric_limits<real>::min());
		const real kbfloor = 1e-12 * std::max<real>(Kb.maxCoeff(), std::numeric_limits<real>::min());
		for (int i = 0; i < N_; i++) {
			t.emplace_back(oK + i, oK + i, -1. / std::max(Kc[i], kfloor));
			if (bulk) t.emplace_back(oB + i, oB + i, -1. / std::max(Kb[i], kbfloor));
		}
		SpMat M(pc_n_, pc_n_);
		M.setFromTriplets(t.begin(), t.end());
		const auto t1 = std::chrono::steady_clock::now();
		if (!pc_analyzed_) { pc_ldlt_.analyzePattern(M); pc_analyzed_ = true; }
		const auto t2 = std::chrono::steady_clock::now();
		pc_ldlt_.factorize(M);
		const auto t3 = std::chrono::steady_clock::now();
		precond_t_assemble += sec(t0, t1); precond_t_analyze += sec(t1, t2); precond_t_factor += sec(t2, t3);
		precond_nnz_A = (long long)M.nonZeros(); precond_nnz_L = (long long)pc_ldlt_.matrixL().nestedExpression().nonZeros();
		if (pc_ldlt_.info() != Eigen::Success) return false;
		const VecX& Dg = pc_ldlt_.vectorD();
		int neg = 0;
		for (int k = 0; k < pc_n_; k++) {
			if (!std::isfinite(Dg[k])) return false;
			neg += Dg[k] < 0.;
		}
		precond_negatives = neg; precond_aux = pc_n_ - nq;
		return neg == pc_n_ - nq;
	}
	Vec3 B1B2(int i, int a) const { return a == 0 ? B1(i) : B2(i); }
	Eigen::SimplicialLDLT<SpMat, Eigen::Lower, Eigen::AMDOrdering<int>> pc_ldlt_;
	bool pc_analyzed_ = false, pc_use_chol_ = false;
	int pc_n_ = 0;
#ifdef VF_HAVE_CHOLMOD
	Eigen::CholmodSupernodalLLT<SpMat, Eigen::Lower> pc_chol_;
#endif

	Vec3 D0(int q) const { return st_.x[sp_.pairs.pj[q]] - st_.x[sp_.pairs.pi[q]] + sp_.shift[q]; }
	void Densities(const std::vector<Vec3>& x, bool grad)
	{
		const PairList& P = sp_.pairs;
		const int np = P.Size();
		const real h = ker_.Radius();
		nd_.assign(N_, ker_.W(0.));
		d_.resize(np); r_.resize(np);
		if (grad) g_.resize(np);
		for (int q = 0; q < np; q++) {
			d_[q] = x[P.pj[q]] - x[P.pi[q]] + sp_.shift[q];
			r_[q] = d_[q].norm();
			if (r_[q] < h) {
				const real w = ker_.W(r_[q]);
				nd_[P.pi[q]] += w; nd_[P.pj[q]] += w;
				if (grad) g_[q] = r_[q] > 1e-12 * h ? Vec3(ker_.dW(r_[q]) / r_[q] * d_[q]) : Vec3::Zero();
			}
			else if (grad) g_[q] = Vec3::Zero();
		}
	}
	// dE/dN_i: per-particle term plus the squeeze term V_i (sigma / 4) sum_q c_q (h_i - h_j)
	void Coef_N(std::vector<real>& coef) const
	{
		const PairList& P = sp_.pairs;
		coef.resize(N_);
		Parallel_For(N_, [&](int i) { coef[i] = Ei_d(i, nd_[i]); }, 64);
		const real ks = p_.tf.sigma / 4.;
		for (int q = 0; q < P.Size(); q++) {
			if (csq_[q] == 0.) continue;
			const int i = P.pi[q], j = P.pj[q];
			const real w = ks * csq_[q] * (st_.V[i] * nd_[i] - st_.V[j] * nd_[j]);
			coef[i] += st_.V[i] * w; coef[j] -= st_.V[j] * w;
		}
	}
	void Add_Grad_N(const std::vector<real>& coef, std::vector<Vec3>& G) const
	{
		const PairList& P = sp_.pairs;
		for (int q = 0; q < P.Size(); q++) {
			const int i = P.pi[q], j = P.pj[q];
			const Vec3 f = (coef[i] + coef[j]) * g_[q];
			G[j] += f; G[i] -= f;
		}
	}
	real Squeeze_Value() const
	{
		const PairList& P = sp_.pairs;
		real E = 0.;
		for (int q = 0; q < P.Size(); q++) {
			if (csq_[q] == 0.) continue;
			const real dh = st_.V[P.pi[q]] * nd_[P.pi[q]] - st_.V[P.pj[q]] * nd_[P.pj[q]];
			E += csq_[q] * dh * dh;
		}
		return p_.tf.sigma / 8. * E;
	}
	// Trouton: (1/hh) [ sum_p cT_p ((d_j - d_i) . e_p)^2 + 1/2 sum_i kb_i div_i(d)^2 ], d = x - xref
	void Divergence(const std::vector<Vec3>& d, std::vector<real>& dv) const
	{
		const PairList& P = sp_.pairs;
		dv.assign(N_, 0.);
		for (int q = 0; q < P.Size(); q++) {
			const int i = P.pi[q], j = P.pj[q];
			const real u = (d[j] - d[i]).dot(gb_[q]);
			dv[i] -= a0_[j] * u; dv[j] -= a0_[i] * u;
		}
	}
	real Visc_Value() const
	{
		if (p_.Mu_S() <= 0.) return 0.;
		const PairList& P = sp_.pairs;
		std::vector<Vec3> d(N_);
		for (int i = 0; i < N_; i++) d[i] = x_[i] - xref_[i];
		real E = 0.;
		for (int q = 0; q < P.Size(); q++) { const real u = (d[P.pj[q]] - d[P.pi[q]]).dot(e_[q]); E += cT_[q] * u * u; }
		if (p_.bulk_viscosity) { std::vector<real> dv; Divergence(d, dv); for (int i = 0; i < N_; i++) E += 0.5 * kb_[i] * dv[i] * dv[i]; }
		return E / hh_;
	}
	// G += (1/hh) d/dd [ ... ] at d (also the Hessian-vector product, the form being quadratic)
	void Add_Visc(const std::vector<Vec3>& d, std::vector<Vec3>& G) const
	{
		if (p_.Mu_S() <= 0.) return;
		const PairList& P = sp_.pairs;
		for (int q = 0; q < P.Size(); q++) {
			const int i = P.pi[q], j = P.pj[q];
			const Vec3 f = 2. * cT_[q] / hh_ * (d[j] - d[i]).dot(e_[q]) * e_[q];
			G[j] += f; G[i] -= f;
		}
		if (!p_.bulk_viscosity) return;
		std::vector<real> dv;
		Divergence(d, dv);
		for (int q = 0; q < P.Size(); q++) {
			const int i = P.pi[q], j = P.pj[q];
			const real c = -(kb_[i] * dv[i] * a0_[j] + kb_[j] * dv[j] * a0_[i]) / hh_;
			G[j] += c * gb_[q]; G[i] -= c * gb_[q];
		}
	}
	// q-components of an x-space vector: T^T G (flat), T^T J^T G (sphere, J at u; J symmetric)
	void Gather(const std::vector<Vec3>& G, VecX& g, const std::vector<Vec3>& u = {}) const
	{
		g.resize(Dim());
		for (size_t k = 0; k < free_.size(); k++) {
			const int i = free_[k];
			const Vec3 Gi = (sph_ && !u.empty()) ? sph_->Jac(u[i], G[i]) : G[i];
			g[2 * k] = Gi.dot(B1(i)); g[2 * k + 1] = Gi.dot(B2(i));
		}
	}
	// x-gradient of the whole incremental potential at x (Densities(x, true) and Coef_N(coef) just evaluated)
	void Full_Grad(const std::vector<Vec3>& x, const std::vector<real>& coef, std::vector<Vec3>& G) const
	{
		G.assign(N_, Vec3::Zero());
		Add_Grad_N(coef, G);
		const real inv = 1. / (hh_ * hh_);
		for (int i = 0; i < N_; i++) {
			G[i] += p_.tf.rho_g * st_.V[i] * Weight(i) * p_.up - Fext(i);
			if (dof_[i] >= 0) G[i] += p_.rho * st_.V[i] * inv * (x[i] - xt_[i]) + Drag_F(i, x[i] - xref_[i]) / hh_;
		}
		std::vector<Vec3> d(N_);
		for (int i = 0; i < N_; i++) d[i] = x[i] - xref_[i];
		Add_Visc(d, G);
	}

	const Quintic2D& ker_;
	const FilmMechParams& p_;
	const ShiftedPairs& sp_;
	const FilmMatState& st_;
	const std::vector<Vec3>& xt_;
	const std::vector<Vec3>& xref_;
	real hh_;
	int N_;
	const std::vector<Vec3>* b1_;
	const std::vector<Vec3>* b2_;
	const SphereShape* sph_;
	std::vector<int> dof_, free_;
	std::vector<real> a0_, h0_, csq_, cT_, kb_;
	std::vector<Vec3> e_, gb_;
	std::vector<Vec3> x_, xh_, u_, uh_, Gh_;
	std::vector<real> nd_, ndh_, r_, Kr_, coef_, coefh_, d2_, ei_;
	real pib_ = 0.;                   // surface pressure of the surfactant bath (sol.c_bath)
	mutable std::vector<real> lnc_;  // per particle: ln c of its last partition (warm start; soluble surfactant)
	std::vector<Vec3> d_, g_, gh_, Kd_;
};

// One BDF2 step of the mobile surfaces (V and S fixed: the plug flow carries them). Flat: periodic box. On a sphere
// (on_sphere): 2 DOFs per particle in the tangent plane at its start position, mapped onto the sphere by the radial
// projection inside the energy (the constraint reaction is radial at the actual position); velocities kept tangential.
struct FilmMech_Stepper {
	Quintic2D kernel;
	FilmMechParams prm;
	FilmBox box;
	bool on_sphere = false;
	SphereShape sphere;
	NewtonOptions newton;
	real list_margin = 0.5;
	bool bdf2 = true;
	real rel_tol = 1e-9;
	// ghost_frame: for the mechanics the frame (fixed) particles take the thickness and concentration of the film next
	// to them (kernel-weighted means over their free neighbours), so the frame is mechanically neutral: it carries the
	// film's tension as a wire would. Without it the frame particles keep their own state, which drifts (their areas
	// follow the film's motion, they drain), and a frame under more tension than the film pulls film particles into it.
	// Their own V and S are restored after the step (the drainage of the film at the wire is unchanged).
	// 1: thickness and concentration. 2: thickness, and the concentration that gives each frame particle the film's
	// stress coefficient dE/dN = a^2 p (p = 2 Pi_s + h P' - P - 2 Pi_bath; kernel-weighted mean over its free
	// neighbours): the pair forces (c_i + c_j) |W'| are what holds the film, and a frame lattice denser than a
	// stretched film (smaller a) at the film's concentration pulls less than the film does, so under tension the top
	// rows come off the frame.
	int ghost_frame = 0;
	// max_halvings > 0: a step whose Newton solve fails, or that moves a particle by more than a quarter of the pair
	// list's margin (pairs that came into range during the step would be missing: particles could overlap), is redone
	// as two half steps, recursively up to max_halvings times (the half steps are BDF1, the history's dt differing)
	int max_halvings = 0;
	// factored_precond: the exact preconditioner of FilmMech_Problem (physical surfactant elasticity:
	// the surfaces are nearly incompressible and Jacobi stalls): CHOLMOD on the assembled Hessian (precond_cholmod, when built
	// with it) or the LDL^T of the 1-ring saddle-point form; refactorized every newton.precond_rebuild_every
	// iterations (0: once per solve). Diagnostics of the last step: precond_stress_projected.
	bool factored_precond = false;
	bool precond_cholmod = true; // CHOLMOD on the assembled Hessian when built with it (FilmMech_Problem)
	int precond_stress_projected = 0;

	NewtonStats Step(FilmMatState& st, real dt) { return Step_Adaptive(st, dt, max_halvings); }

private:
	NewtonStats Step_Adaptive(FilmMatState& st, real dt, int halvings)
	{
		if (halvings <= 0) return Step_Once(st, dt);
		const FilmMatState s0 = st;
		NewtonStats ns = Step_Once(st, dt);
		real dmax = 0.;
		for (int i = 0; i < st.Size(); i++) {
			if (!st.fixed.empty() && st.fixed[i]) continue;
			const Vec3 d = on_sphere ? Vec3(st.x[i] - s0.x[i]) : box.Wrap(st.x[i] - s0.x[i]);
			dmax = std::max(dmax, d.norm());
		}
		if (ns.converged && dmax <= 0.25 * list_margin * kernel.Radius()) return ns;
		st = s0;
		ns = Step_Adaptive(st, 0.5 * dt, halvings - 1);
		const NewtonStats n2 = Step_Adaptive(st, 0.5 * dt, halvings - 1);
		ns.converged = ns.converged && n2.converged;
		ns.iters += n2.iters; ns.total_cg_iters += n2.total_cg_iters;
		ns.grad_inf_history.insert(ns.grad_inf_history.end(), n2.grad_inf_history.begin(), n2.grad_inf_history.end());
		return ns;
	}

	NewtonStats Step_Once(FilmMatState& st, real dt)
	{
		const int N = st.Size();
		ShiftedPairs sp;
		std::vector<Vec3> b1, b2;
		if (on_sphere) {
			Build_Pairs(st.x, kernel.Radius() * (1. + list_margin), sp.pairs);
			sp.shift.assign(sp.pairs.Size(), Vec3::Zero());
			b1.resize(N); b2.resize(N);
			for (int i = 0; i < N; i++) sphere.Basis(st.x[i], b1[i], b2[i]);
		}
		else Build_Shifted_Pairs(st.x, box, kernel.Radius() * (1. + list_margin), sp);
		std::vector<real> V_own, S_own;
		if (ghost_frame && !st.fixed.empty()) { V_own = st.V; S_own = st.S; Set_Ghost_Frame(st, sp); }
		const bool use_bdf2 = bdf2 && st.has_prev && (int)st.x_prev.size() == N && std::abs(st.dt_prev - dt) <= 1e-12 * dt;
		const real hh = use_bdf2 ? 2. * dt / 3. : dt;
		std::vector<Vec3> xt(N), xref(N);
		for (int i = 0; i < N; i++) {
			xt[i] = use_bdf2 ? Vec3(4. / 3. * st.x[i] - 1. / 3. * st.x_prev[i] + 8. * dt / 9. * st.v[i] - 2. * dt / 9. * st.v_prev[i]) : Vec3(st.x[i] + dt * st.v[i]);
			xref[i] = use_bdf2 ? Vec3(4. / 3. * st.x[i] - 1. / 3. * st.x_prev[i]) : st.x[i];
		}
		FilmMech_Problem P(kernel, prm, sp, st, xt, xref, hh, on_sphere ? &b1 : nullptr, on_sphere ? &b2 : nullptr, on_sphere ? &sphere : nullptr);
		P.factored_precond = factored_precond;
		P.precond_cholmod = precond_cholmod;
		VecX q = VecX::Zero(P.Dim()), g0;
		for (size_t k = 0; k < P.Free().size(); k++) {
			const int i = P.Free()[k];
			const Vec3 d = xt[i] - st.x[i];
			q[2 * k] = d.dot(P.B1(i)); q[2 * k + 1] = d.dot(P.B2(i));
		}
		// start from the inertial prediction unless the current positions are better: after a fast transient the
		// prediction can carry the motion into the Langmuir or disjoining wall (energy +inf: no line search possible)
		{
			const real Eq = P.Value(q);
			const VecX q0 = VecX::Zero(P.Dim());
			const real E0 = P.Value(q0);
			if (!std::isfinite(Eq) || (std::isfinite(E0) && E0 < Eq)) q = q0;
		}
		P.Gradient(q, g0);
		NewtonOptions o = newton;
		o.sparse_preconditioner = factored_precond;
		o.grad_tol = std::max<real>(std::numeric_limits<real>::min(), rel_tol * g0.lpNorm<Eigen::Infinity>());
		o.step_tol = 1e-13 * kernel.Radius();
		o.max_step = std::min(o.max_step, 0.25 * list_margin * kernel.Radius());
		const NewtonStats ns = Newton_CG_Minimize(P, q, o);
		precond_stress_projected = P.precond_stress_projected;
		std::vector<Vec3> xn;
		P.To_X(q, xn);
		if (on_sphere) for (int i = 0; i < N; i++) if (st.fixed.empty() || !st.fixed[i]) xn[i] = sphere.Project(xn[i]);
		const std::vector<Vec3> x_old = st.x, v_old = st.v;
		for (int i = 0; i < N; i++) {
			if (!st.fixed.empty() && st.fixed[i]) { st.v[i].setZero(); continue; }
			st.v[i] = use_bdf2 ? Vec3(1.5 / dt * (xn[i] - 4. / 3. * st.x[i] + 1. / 3. * st.x_prev[i])) : Vec3((xn[i] - st.x[i]) / dt);
			st.x[i] = xn[i];
			if (on_sphere) { const Vec3 n = sphere.Normal(st.x[i]); st.v[i] -= st.v[i].dot(n) * n; continue; }
			// keep positions inside the periodic box
			st.x[i][0] -= box.Lx * std::floor(st.x[i][0] / box.Lx);
			if (box.periodic_y) st.x[i][1] -= box.Ly * std::floor(st.x[i][1] / box.Ly);
		}
		// history must follow the wrap: shift x_prev by the same lattice vectors
		st.x_prev = x_old;
		for (int i = 0; i < N; i++) {
			const Vec3 dmove = st.x[i] - xn[i]; // 0 or a lattice vector
			st.x_prev[i] += dmove;
		}
		st.v_prev = v_old; st.dt_prev = dt; st.has_prev = true;
		if (!V_own.empty()) for (int i = 0; i < N; i++) if (st.fixed[i]) { st.V[i] = V_own[i]; st.S[i] = S_own[i]; }
		return ns;
	}

private:
	// frame particles with free neighbours within the kernel radius: V, S such that h and Gamma (mode 1) or h and the
	// stress coefficient (mode 2) are the kernel-weighted means of the free neighbours' (densities at the current positions)
	void Set_Ghost_Frame(FilmMatState& st, const ShiftedPairs& sp) const
	{
		const int N = st.Size();
		const PairList& P = sp.pairs;
		const real R = kernel.Radius();
		std::vector<real> nd(N, kernel.W(0.)), w(P.Size(), 0.);
		for (int q = 0; q < P.Size(); q++) {
			const real r = (st.x[P.pj[q]] - st.x[P.pi[q]] + sp.shift[q]).norm();
			if (r < R) { w[q] = kernel.W(r); nd[P.pi[q]] += w[q]; nd[P.pj[q]] += w[q]; }
		}
		const Disjoining& D = prm.tf.disj;
		auto bath = [&](int j) { return prm.bath_pressure * ((int)st.bath_w.size() == N ? st.bath_w[j] : 1.); };
		auto coef = [&](int j) { // dE_j/dN_j = a^2 p
			const real h = st.V[j] * nd[j], G = st.S[j] * nd[j];
			const real p = (st.S[j] > 0. ? 2. * prm.surf.Pressure(G) : 0.) + h * D.dP(h) - D.P(h) - 2. * bath(j);
			return p / (nd[j] * nd[j]);
		};
		std::vector<real> ws(N, 0.), hs(N, 0.), gs(N, 0.);
		for (int q = 0; q < P.Size(); q++) {
			if (w[q] == 0.) continue;
			const int i = P.pi[q], j = P.pj[q];
			if (st.fixed[i] == st.fixed[j]) continue;
			const int f = st.fixed[i] ? i : j, m = st.fixed[i] ? j : i; // frame, film
			ws[f] += w[q]; hs[f] += w[q] * st.V[m] * nd[m];
			gs[f] += w[q] * (ghost_frame == 2 ? coef(m) : st.S[m] * nd[m]);
		}
		for (int i = 0; i < N; i++) {
			if (!st.fixed[i] || !(ws[i] > 0.)) continue;
			const real h = hs[i] / ws[i];
			st.V[i] = h / nd[i];
			if (ghost_frame != 2) { st.S[i] = gs[i] / ws[i] / nd[i]; continue; }
			// 2 Pi_s(Gamma) = c N^2 - (h P' - P) + 2 Pi_bath
			const real Pi = 0.5 * (gs[i] / ws[i] * nd[i] * nd[i] - (h * D.dP(h) - D.P(h)) + 2. * bath(i));
			st.S[i] = prm.surf.Concentration(std::max<real>(Pi, 1e-6 * prm.surf.RT)) / nd[i];
		}
	}
};

// Surfactant energy of the amounts S at fixed positions (graph areas a_i): E(S) = sum_i 2 a_i g(S_i / a_i);
// the chemical potential dE/dS_i = 2 g'(Gamma_i). Energy concept of H1_Flux_Problem.
class SurfactantEnergy {
public:
	SurfactantEnergy(const ThinFilmGraph& G, const Surfactant& s) : G_(G), s_(s) {}
	real Value(const std::vector<real>& S) const
	{
		real E = 0.;
		for (int i = 0; i < G_.N; i++) {
			const real Gm = S[i] / G_.a[i];
			if (!s_.Valid(Gm)) return std::numeric_limits<real>::infinity();
			E += 2. * G_.a[i] * s_.g(Gm);
		}
		return E;
	}
	void Pressure(const std::vector<real>& S, std::vector<real>& mu) const
	{
		mu.resize(G_.N);
		for (int i = 0; i < G_.N; i++) mu[i] = 2. * s_.dg(S[i] / G_.a[i]);
	}
	void Hess_Apply(const std::vector<real>& S, const std::vector<real>& w, std::vector<real>& out) const
	{
		out.resize(G_.N);
		for (int i = 0; i < G_.N; i++) out[i] = Hess_Diag(S, i) * w[i];
	}
	real Hess_Diag(const std::vector<real>& S, int i) const { return 2. * s_.d2g(S[i] / G_.a[i]) / G_.a[i]; }
	void Hess_Triplets(const std::vector<real>& S, std::vector<Eigen::Triplet<real>>& t, const std::vector<real>* = nullptr) const
	{
		for (int i = 0; i < G_.N; i++) t.emplace_back(i, i, Hess_Diag(S, i));
	}
	const ThinFilmGraph& Graph() const { return G_; }

private:
	const ThinFilmGraph& G_;
	const Surfactant& s_;
};

// Surface diffusion of the surfactant (positions fixed: a split step): the gradient flow of E(S) with the pair
// conductances k_p = M(Gamma_p) c_p, M = Ds / (2 g''(Gamma)) (= Ds Gamma (1 - theta) / (2 RT) for Langmuir), so that the
// flux is -M grad(2 g') = -Ds grad Gamma: Fick's law with diffusivity Ds for any g. Backward Euler; S conserved exactly.
struct Surface_Diffusion_Stepper {
	Surfactant surf;
	real Ds = 0.;
	NewtonOptions newton;
	real rel_tol = 1e-9;
	bool sparse_precond = false, factored_precond = false;
	real precond_radius = 0.;
	void Conductances(const ThinFilmGraph& G, const std::vector<real>& S, std::vector<real>& k) const
	{
		k.resize(G.pairs.Size());
		for (int q = 0; q < G.pairs.Size(); q++) {
			const int i = G.pairs.pi[q], j = G.pairs.pj[q];
			k[q] = Ds / (2. * surf.d2g(0.5 * (S[i] / G.a[i] + S[j] / G.a[j]))) * G.c[q];
		}
	}
	NewtonStats Step(const ThinFilmGraph& G, std::vector<real>& S, real dt) const
	{
		if (Ds <= 0.) { NewtonStats s; s.converged = true; return s; }
		SurfactantEnergy E(G, surf);
		std::vector<real> k;
		Conductances(G, S, k);
		return H1_Flux_Step(E, S, k, dt, newton, rel_tol, sparse_precond, factored_precond, precond_radius);
	}
};

// Soluble surfactant: each particle's amount Q = 2 S + B split again between its surfaces and its liquid at its
// current area (FilmMechParams::sol; Q conserved exactly: B = Vsol c, S = (Q - B) / 2). Bath particles (sbath) are left
// to their owner.
inline void Repartition(FilmMatState& st, const std::vector<real>& area, const Solubility& sol)
{
	const int n = st.Size();
	if (!sol.On() || (int)st.B.size() != n) return;
	Parallel_For(n, [&](int i) {
		if ((int)st.sbath.size() == n && st.sbath[i]) return;
		const real Q = 2. * st.S[i] + st.B[i];
		const real Vl = st.Vsol(i);
		if (!(Q > 0.) || !(Vl > 0.)) return;
		const Solubility::Part p = sol.Partition(area[i], Vl, Q);
		st.B[i] = Vl * p.c;
		st.S[i] = 0.5 * (Q - st.B[i]);
	}, 64);
}

// The dissolved surfactant carried by the drainage: F_q the volume pair q passed from pi to pj in the step
// (Drainage_Stepper::Step), V the volumes after it. Upwind (the liquid leaves at its donor's concentration), implicit
// in the concentrations: c_i V_i + sum_out F c_i - sum_in F c_j = B_i (before), an M-matrix whose columns sum to V > 0,
// so the concentrations stay positive and the amount is conserved exactly. Reservoir nodes (c_res not NaN: a meniscus)
// supply liquid at c_res and take what reaches them; their B is left alone. Returns the amount the reservoirs gave
// (+) or took (-). The solute's own driving force in the drainage (mixing: second order in the concentration
// differences) is left out, and so is diffusion along the film (D ~ 1e-10 - 1e-9 m^2/s: ~10 s per 100 um).
inline real Solute_Transport(const ThinFilmGraph& G, const std::vector<real>& flux, const std::vector<real>& V, std::vector<real>& B,
	const std::vector<real>* c_res = nullptr)
{
	const int n = G.N;
	const PairList& P = G.pairs;
	auto res = [&](int i) { return c_res && !std::isnan((*c_res)[i]); };
	std::vector<int> id(n, -1);
	int m = 0;
	for (int i = 0; i < n; i++) if (!res(i)) id[i] = m++;
	std::vector<Eigen::Triplet<real>> t;
	t.reserve(m + 2 * P.Size());
	VecX rhs(m);
	for (int i = 0; i < n; i++) if (id[i] >= 0) { t.emplace_back(id[i], id[i], V[i]); rhs[id[i]] = B[i]; }
	real from_res = 0.;
	for (int q = 0; q < P.Size(); q++) {
		const real F = flux[q];
		if (F == 0.) continue;
		const int d = F > 0. ? P.pi[q] : P.pj[q], r = F > 0. ? P.pj[q] : P.pi[q];
		const real f = std::abs(F);
		if (id[d] >= 0) t.emplace_back(id[d], id[d], f);
		if (id[r] < 0) continue;
		if (id[d] >= 0) t.emplace_back(id[r], id[d], -f);
		else { rhs[id[r]] += f * (*c_res)[d]; from_res += f * (*c_res)[d]; }
	}
	SpMat A(m, m);
	A.setFromTriplets(t.begin(), t.end());
	Eigen::SparseLU<SpMat> lu;
	lu.compute(A);
	const VecX c = lu.solve(rhs);
	for (int q = 0; q < P.Size(); q++) { // what reached the reservoirs
		const real F = flux[q];
		const int d = F > 0. ? P.pi[q] : P.pj[q], r = F > 0. ? P.pj[q] : P.pi[q];
		if (F != 0. && id[d] >= 0 && id[r] < 0) from_res -= std::abs(F) * c[id[d]];
	}
	for (int i = 0; i < n; i++) if (id[i] >= 0) B[i] = c[id[i]] * V[i];
	return from_res;
}

// Mobile film with drainage and surface diffusion, Lie splitting per step:
//   1. mechanics (FilmMech_Stepper): the surfaces move, V and S ride with them (plug flow);
//   2. drainage (Drainage_Stepper) at the new positions: Poiseuille flow of the liquid relative to the surfaces;
//   3. surface diffusion (Surface_Diffusion_Stepper) of the surfactant along the surfaces.
// Each substep is a minimizing movement of the same energy in its own variables (positions, V, S), so every substep
// is stable for any dt; the splitting error is first order in dt.
struct FilmCoupled_Stepper {
	FilmMech_Stepper mech;      // kernel, prm (shared by all substeps), box, Newton options
	Drainage_Stepper drain;     // prm is overwritten by mech.prm.tf
	Surface_Diffusion_Stepper diff; // surf is overwritten by mech.prm.surf; Ds as set
	struct Stats { NewtonStats mech, drain, diff; };
	ThinFilmGraph Graph(const FilmMatState& st) const
	{
		ThinFilmGraph G = mech.on_sphere ? Make_Surface_Graph(st.x, mech.kernel) : Make_Periodic_Patch(st.x, mech.box.Lx, mech.box.Ly, mech.kernel, mech.box.periodic_y);
		G.height.resize(G.N);
		for (int i = 0; i < G.N; i++) G.height[i] = mech.prm.up.dot(st.x[i]);
		return G;
	}
	bool Soluble(const FilmMatState& st) const { return mech.prm.sol.On() && (int)st.B.size() == st.Size(); }
	// soluble surfactant: the mechanics moves each particle's amount Q with it (the energy takes the partition at
	// every trial area); after it the split is updated, the drainage carries the dissolved part (Solute_Transport),
	// and the split is updated again at the new volumes
	Stats Step(FilmMatState& st, real dt)
	{
		Stats s;
		s.mech = mech.Step(st, dt);
		const ThinFilmGraph G = Graph(st);
		const bool sol = Soluble(st);
		if (sol) Repartition(st, G.a, mech.prm.sol);
		drain.prm = mech.prm.tf;
		std::vector<real> flux;
		s.drain = drain.Step(G, st.V, dt, sol ? &flux : nullptr);
		if (sol) { Solute_Transport(G, flux, st.V, st.B); Repartition(st, G.a, mech.prm.sol); }
		diff.surf = mech.prm.surf;
		if (!sol) s.diff = diff.Step(G, st.S, dt);
		else s.diff.converged = true; // surface diffusion of a soluble surfactant: not yet (Solute_Transport)
		return s;
	}
};

// Particle insertion and merging for a mobile film whose surfaces stretch or compress (the particles are material:
// a sliding film's stretched top loses resolution, its compressed bottom wastes it).
// The criteria are geometric (spacing), because kernel areas average over ~60 neighbours and barely respond to one
// split (an area criterion cascades: the children still read as sparse and split again, step after step):
// - insertion: a particle with a gap in some in-plane direction (of 8; the nearest neighbour inside a 45-degree cone,
//   or the kernel radius if none) wider than split_gap * dx0 gets one new particle at the middle of the gap, which
//   takes half of its V and S, its velocity and its BDF2 history shifted alike (the gap is halved at once);
// - merging: two particles closer than merge_dist * dx0 (closest pairs first, neither fixed nor already taken) become
//   one: V and S summed, position, velocity and history mass-weighted.
// Kernel areas are smoothed (over ~60 neighbours), so moving V and S alone would make the fields jump: a merged
// particle carries both shares in about one particle's area (h and Gamma nearly double, Gamma past saturation), a
// split one half a share in about the same area. So the fields are remapped conservatively: with the areas of the new
// particle set (make_graph), every affected particle (created, merged, absorbing, or within the kernel radius of a
// change) gets back the thickness and concentration its material had (its own, its parent's, or the area-weighted
// mean of the merged pair) times its new area, and the affected region is rescaled to its exact totals.
// Volume, surfactant and momentum are conserved exactly. Frame (fixed) particles never move. Flat films use the
// box's minimum image; on a sphere the new positions are projected back onto it.
struct RemeshParams {
	real dx0 = 1.;            // target particle spacing
	real split_gap = 1.9;     // insert into gaps wider than this (in dx0)
	real merge_dist = 0.45;   // merge pairs closer than this (in dx0)
	real kernel_radius = 4.;  // the remap region around each change
	bool absorb_fixed = true; // a free particle crowding a fixed (frame) one is absorbed by it; false: left alone (a
	                          // wire does not swallow film: absorption makes the frame a sink the film can be drawn into)
	bool remap_fixed = true;  // fixed particles within the kernel radius of a change are remapped with the free ones; false:
	                          // their V and S are left alone, so the free particles keep their totals exactly (a frame
	                          // that is a wall for an insoluble monolayer, Meniscus.h)
	std::function<bool(const Vec3&)> insert_ok; // where insertions are allowed (empty: everywhere), e.g. not below a bath
	// Optional numerical interface ownership. Pure remeshing then preserves V,S,B
	// in EACH region, not just their global totals. Physical motion/fluxes may still
	// cross the interface in other substeps. Requires no fixed-node absorption/remap.
	// This constrained redistribution is not a resolved meniscus or a physical law.
	std::function<int(const Vec3&)> remap_region;
};
struct RemeshStats { int splits = 0, merges = 0, remapped = 0; };

inline void Remap_Region_Amounts(const FilmMatState& before, FilmMatState& out,
	const ThinFilmGraph& G, const std::vector<char>& affected,
	const std::vector<real>& target_h, const std::vector<real>& target_s,
	const std::vector<real>& target_c, const std::function<int(const Vec3&)>& region)
{
	using Triple=std::array<real,3>;
	std::map<int,Triple> quota,target;
	const bool sol=(int)out.B.size()==out.Size();
	for(int i=0;i<before.Size();++i)if(before.fixed.empty()||!before.fixed[i]){
		auto& q=quota[region(before.x[i])];q[0]+=before.V[i];q[1]+=before.S[i];if(sol)q[2]+=before.B[i];
	}
	const auto original=quota;
	for(int i=0;i<out.Size();++i)if(out.fixed.empty()||!out.fixed[i]){
		const int label=region(out.x[i]);auto& q=quota[label];auto& t=target[label];
		if(!affected[i]){q[0]-=out.V[i];q[1]-=out.S[i];if(sol)q[2]-=out.B[i];}
		else{t[0]+=target_h[i]*G.a[i];t[1]+=target_s[i]*G.a[i];}
	}
	std::map<int,Triple> factors;
	for(const auto& entry:quota)if(!target.count(entry.first)){
		const auto old=original.find(entry.first);const Triple scale=old==original.end()?Triple{}:old->second;
		for(int j=0;j<3;++j)if(std::abs(entry.second[j])>1e-12*std::max(std::abs(scale[j]),1e-30))
			throw std::runtime_error("region inventory lost all remap support");
	}
	for(const auto& entry:target){
		const int label=entry.first;const auto& q=quota.at(label);auto& f=factors[label];
		const auto old=original.find(label);const Triple scale=old==original.end()?Triple{}:old->second;
		if(!(entry.second[0]>0.)){
			for(int j=0;j<3;++j)if(std::abs(q[j])>1e-12*std::max(std::abs(scale[j]),1e-30))throw std::runtime_error("empty region remap support");
			f={1.,1.,1.};continue;
		}
		for(int j=0;j<2;++j){
			if(q[j]<-1e-10*std::max(std::abs(entry.second[j]),1e-30))throw std::runtime_error("negative region remap quota");
			if(entry.second[j]<=0.&&q[j]>0.)throw std::runtime_error("no support for region remap quota");
			f[j]=entry.second[j]>0.?std::max<real>(0.,q[j])/entry.second[j]:1.;
		}
	}
	if(sol){
		for(int i=0;i<out.Size();++i)if(affected[i]&&(out.fixed.empty()||!out.fixed[i])){
			const int label=region(out.x[i]);const real V=factors.at(label)[0]*target_h[i]*G.a[i];
			const real Vs=(int)out.Vs.size()==out.Size()&&out.Vs[i]>0.?out.Vs[i]:V;
			target[label][2]+=target_c[i]*Vs;
		}
		for(auto& entry:factors){
			const auto& q=quota.at(entry.first);const real t=target.at(entry.first)[2];
			if(!(target.at(entry.first)[0]>0.))continue;
			if(q[2]<-1e-10*std::max(std::abs(t),1e-30)||(!(t>0.)&&q[2]>0.))throw std::runtime_error("infeasible dissolved region quota");
			entry.second[2]=t>0.?std::max<real>(0.,q[2])/t:1.;
		}
	}
	for(int i=0;i<out.Size();++i)if(affected[i]&&(out.fixed.empty()||!out.fixed[i])){
		const auto& f=factors.at(region(out.x[i]));const real V=f[0]*target_h[i]*G.a[i];
		if(!(V>0.))throw std::runtime_error("nonpositive regional remap volume");
		const real keep=out.V[i]/V;out.v[i]*=keep;
		if(out.has_prev&&(int)out.v_prev.size()==out.Size())out.v_prev[i]*=keep;
		out.V[i]=V;out.S[i]=f[1]*target_s[i]*G.a[i];
		if(sol)out.B[i]=f[2]*target_c[i]*out.Vsol(i);
	}
}

inline RemeshStats Remesh_Film(FilmMatState& st, const ThinFilmGraph& G, const RemeshParams& rp, const FilmBox* box, const SphereShape* sphere,
	const std::function<ThinFilmGraph(const std::vector<Vec3>&)>& make_graph)
{
	RemeshStats rs;
	if(rp.remap_region&&(rp.absorb_fixed||rp.remap_fixed))
		throw std::invalid_argument("regional remap requires fixed absorption/remapping disabled");
	const int n = st.Size();
	const PairList& P = G.pairs;
	auto disp = [&](int i, int j) { const Vec3 d = st.x[j] - st.x[i]; return (box && !sphere) ? box->Wrap(d) : d; };
	auto is_fixed = [&](int i) { return !st.fixed.empty() && st.fixed[i]; };
	std::vector<int> partner(n, -1), absorbed_by(n, -1);
	std::vector<char> taken(n, 0);
	const bool sol = (int)st.B.size() == n; // soluble surfactant: the dissolved amounts B go with the liquid
	const bool hasVs = sol && (int)st.Vs.size() == n;
	std::vector<real> addV(n, 0.), addS(n, 0.), addB(n, 0.);
	// merges: the closest pairs first; a free particle crowding a fixed (frame) one is absorbed by it (V and S go to
	// the frame particle, which stays put: the frame takes the momentum, as a wall)
	std::vector<std::pair<real, int>> close;
	for (int q = 0; q < P.Size(); q++) {
		const int i = P.pi[q], j = P.pj[q];
		if (is_fixed(i) && is_fixed(j)) continue;
		const real d = disp(i, j).norm();
		if (d < rp.merge_dist * rp.dx0) close.emplace_back(d, q);
	}
	std::sort(close.begin(), close.end());
	for (const auto& c : close) {
		int i = P.pi[c.second], j = P.pj[c.second];
		if (is_fixed(i)) std::swap(i, j); // i free
		if (taken[i]) continue;
		if (is_fixed(j)) {
			if (!rp.absorb_fixed) continue;
			absorbed_by[i] = j; taken[i] = 1; addV[j] += st.V[i]; addS[j] += st.S[i]; if (sol) addB[j] += st.B[i]; rs.merges++; continue;
		}
		if (taken[j]) continue;
		partner[i] = j; partner[j] = i; taken[i] = taken[j] = 1;
	}
	// insertions: the widest gap of each particle (one per particle, none closer than dx0 / 2 to another new one)
	std::vector<Vec3> gap_dir(n, Vec3::Zero());
	std::vector<real> gap_len(n, 0.);
	for (int i = 0; i < n; i++) {
		if (is_fixed(i) || taken[i]) continue;
		Vec3 t1(1., 0., 0.), t2(0., 1., 0.);
		if (sphere) sphere->Basis(st.x[i], t1, t2);
		for (int k = 0; k < 16; k++) {
			const real ang = 2. * kPi * k / 16.;
			const Vec3 e = std::cos(ang) * t1 + std::sin(ang) * t2;
			real dmin = G.pairs.Size() ? 1e300 : 0.;
			bool any = false;
			for (int q = P.inc_off[i]; q < P.inc_off[i + 1]; q++) {
				const Vec3 d = disp(i, P.Other(q, i));
				const real along = d.dot(e), dn = d.norm();
				if (dn > 0. && along > 0.7071 * dn) { dmin = std::min(dmin, dn); any = true; }
			}
			if (!any) continue; // nothing in that cone within the kernel: an edge (frame side) or a true hole; skip edges
			if (dmin > gap_len[i]) { gap_len[i] = dmin; gap_dir[i] = e; }
		}
	}
	std::vector<Vec3> inserted;
	std::vector<int> insert_from(n, 0);
	for (int i = 0; i < n; i++) {
		if (gap_len[i] <= rp.split_gap * rp.dx0) continue;
		const Vec3 off = 0.5 * gap_len[i] * gap_dir[i], p = st.x[i] + off;
		// the new particle must keep 0.7 dx0 from every neighbour of the parent (fixed ones included), not only from
		// those inside the cone that measured the gap: else it can land on a frame particle, get a tiny kernel area and
		// an impossible concentration
		bool ok = !rp.insert_ok || rp.insert_ok(p);
		for (int q = P.inc_off[i]; q < P.inc_off[i + 1] && ok; q++) if ((disp(i, P.Other(q, i)) - off).norm() < 0.7 * rp.dx0) ok = false;
		for (const Vec3& q : inserted) if (((box && !sphere) ? box->Wrap(q - p) : Vec3(q - p)).norm() < 0.5 * rp.dx0) { ok = false; break; }
		if (!ok) continue;
		inserted.push_back(p);
		insert_from[i] = 1;
	}
	FilmMatState out;
	out.dt_prev = st.dt_prev; out.has_prev = st.has_prev;
	const bool hist = st.has_prev && (int)st.x_prev.size() == n;
	std::vector<real> th, tg, tb; // target thickness, surface concentration and dissolved concentration of each new particle
	std::vector<char> changed;  // created, merged or absorbing
	std::vector<Vec3> sites;    // where something changed
	real cur_gs = 1.; // the gravity factor of the particle being emitted (st.gscale)
	auto push = [&](const Vec3& x, const Vec3& v, const Vec3& xp, const Vec3& vp, real V, real S, real B, char fx, real h, real g, real cb, char ch) {
		out.x.push_back(x); out.v.push_back(v); out.V.push_back(V); out.S.push_back(S); out.fixed.push_back(fx);
		if (sol) { out.B.push_back(B); tb.push_back(cb); }
		if (!st.gscale.empty()) out.gscale.push_back(cur_gs);
		if (hist) { out.x_prev.push_back(xp); out.v_prev.push_back(vp); }
		th.push_back(h); tg.push_back(g); changed.push_back(ch);
		if (ch) sites.push_back(x);
	};
	auto onto = [&](const Vec3& p) { return sphere ? sphere->Project(p) : p; };
	for (int i = 0; i < n; i++) {
		const Vec3 xp = hist ? st.x_prev[i] : Vec3::Zero(), vp = hist ? st.v_prev[i] : Vec3::Zero();
		const real hi = st.V[i] / G.a[i], gi = st.S[i] / G.a[i], Bi = sol ? st.B[i] : 0., ci = sol ? st.B[i] / st.Vsol(i) : 0.;
		const real Vsi = hasVs ? st.Vs[i] : 0.; // the partition volume, where it is not V
		cur_gs = st.gscale.empty() ? 1. : st.gscale[i];
		if (absorbed_by[i] >= 0) { sites.push_back(st.x[i]); continue; } // its V and S went to the frame particle
		if (partner[i] >= 0) {
			const int j = partner[i];
			if (j < i) continue; // emitted with the first of the pair
			const real m = st.V[i] + st.V[j], wi = st.V[i] / m, wj = st.V[j] / m;
			const Vec3 dj = disp(i, j), shift = st.x[i] + dj - st.x[j]; // j's image next to i
			const Vec3 x = onto(st.x[i] + wj * dj);
			const Vec3 v = wi * st.v[i] + wj * st.v[j];
			const Vec3 xpm = hist ? onto(wi * st.x_prev[i] + wj * (st.x_prev[j] + shift)) : Vec3::Zero();
			const Vec3 vpm = hist ? Vec3(wi * st.v_prev[i] + wj * st.v_prev[j]) : Vec3::Zero();
			const real ab = G.a[i] + G.a[j];
			const real Vlm = sol ? st.Vsol(i) + st.Vsol(j) : m;
			const real Bm = sol ? st.B[i] + st.B[j] : 0.;
			push(x, v, xpm, vpm, m, st.S[i] + st.S[j], Bm, 0, m / ab, (st.S[i] + st.S[j]) / ab, Bm / Vlm, 1);
			if (hasVs) out.Vs.push_back(Vsi > 0. || st.Vs[j] > 0. ? Vlm : 0.);
			rs.merges++;
			continue;
		}
		if (insert_from[i]) {
			// the parent stays, the new particle sits in the middle of the gap; half of V and S each (then remapped)
			const Vec3 off = 0.5 * gap_len[i] * gap_dir[i];
			push(st.x[i], st.v[i], xp, vp, 0.5 * st.V[i], 0.5 * st.S[i], 0.5 * Bi, 0, hi, gi, ci, 1);
			push(onto(st.x[i] + off), st.v[i], hist ? onto(xp + off) : xp, vp, 0.5 * st.V[i], 0.5 * st.S[i], 0.5 * Bi, 0, hi, gi, ci, 1);
			if (hasVs) { out.Vs.push_back(0.5 * Vsi); out.Vs.push_back(0.5 * Vsi); }
			rs.splits++;
			continue;
		}
		push(st.x[i], st.v[i], xp, vp, st.V[i] + addV[i], st.S[i] + addS[i], Bi + addB[i], st.fixed.empty() ? 0 : st.fixed[i], hi, gi, ci, addV[i] > 0. ? 1 : 0);
		if (hasVs) out.Vs.push_back(Vsi);
	}
	// keep flat positions inside the periodic box; the history moves by the same lattice vector
	if (box && !sphere) for (size_t k = 0; k < out.x.size(); k++) {
		Vec3 w = Vec3::Zero();
		w[0] = -box->Lx * std::floor(out.x[k][0] / box->Lx);
		if (box->periodic_y) w[1] = -box->Ly * std::floor(out.x[k][1] / box->Ly);
		out.x[k] += w;
		if (hist) out.x_prev[k] += w;
	}
	// conservative remap of the affected region
	if (!sites.empty()) {
		const ThinFilmGraph G1 = make_graph(out.x);
		const int m = (int)out.x.size();
		auto sep = [&](const Vec3& a, const Vec3& b) { const Vec3 d = b - a; return ((box && !sphere) ? box->Wrap(d) : d).norm(); };
		std::vector<char> aff(m, 0);
		Parallel_For(m, [&](int k) {
			if (changed[k]) { aff[k] = 1; return; }
			if (!rp.remap_fixed && out.fixed[k]) return;
			for (const Vec3& s : sites) if (sep(out.x[k], s) < rp.kernel_radius) { aff[k] = 1; return; }
		});
		if(rp.remap_region){
			Remap_Region_Amounts(st,out,G1,aff,th,tg,tb,rp.remap_region);
			for(int k=0;k<m;++k)if(aff[k]&&!out.fixed[k])++rs.remapped;
		}else{
		real sV = 0., sS = 0., nV = 0., nS = 0., sB = 0., nB = 0.;
		for (int k = 0; k < m; k++) if (aff[k]) { sV += out.V[k]; sS += out.S[k]; nV += th[k] * G1.a[k]; nS += tg[k] * G1.a[k]; }
		const real fV = nV > 0. ? sV / nV : 1., fS = nS > 0. ? sS / nS : 1.;
		// dissolved amounts: each particle gets its material's concentration times its new volume, rescaled to the total
		if (sol) for (int k = 0; k < m; k++) if (aff[k]) { sB += out.B[k]; nB += tb[k] * (hasVs && out.Vs[k] > 0. ? out.Vs[k] : fV * th[k] * G1.a[k]); }
		const real fB = nB > 0. ? sB / nB : 1.;
		for (int k = 0; k < m; k++) if (aff[k]) {
			const real Vnew = fV * th[k] * G1.a[k], keep = out.V[k] / Vnew; // each particle keeps its momentum
			out.v[k] *= keep;
			if (hist) out.v_prev[k] *= keep;
			out.V[k] = Vnew; out.S[k] = fS * tg[k] * G1.a[k]; rs.remapped++;
			if (sol) out.B[k] = fB * tb[k] * (hasVs && out.Vs[k] > 0. ? out.Vs[k] : Vnew);
		}
		}
	}
	st = std::move(out);
	return rs;
}

} // namespace vf::thinfilm
