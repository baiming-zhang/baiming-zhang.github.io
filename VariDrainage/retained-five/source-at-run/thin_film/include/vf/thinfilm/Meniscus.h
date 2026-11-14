//////////////////////////////////////////////////////////////////////////
// A film in a rectangular frame of wetted fibres, with a Plateau-border meniscus along each border
// (thin_film/README.md, "Framed film with menisci").
//
// Cross-section of a border: a fibre of radius a (the frame), the film in the plane through its axis, and on each
// face of the film a meniscus whose surface is a circular arc of radius r, tangent to the film plane (the junction)
// and to the fully wetted fibre. The arc turns through acos(r / (a + r)), so the meniscus surface between the
// junction and the fibre is w = r acos(r / (a + r)) wide per face (pi r / 2 for r << a).
//
// What the model takes from the meniscus (no thickness profile, no solution of the film's flow):
// - liquid: the menisci along the frame are taken in hydrostatic equilibrium (their liquid relaxes along the frame
//   much faster than the film drains into it): the liquid's chemical potential mu_m = p + rho g y is the same all
//   around, with p = -gamma / r(y) the Laplace pressure of the arc. One radius r_ref at one height y_ref (an input:
//   the frame's liquid content) fixes mu_m = -gamma / r_ref + rho g y_ref, and then r(y) = gamma / (rho g y - mu_m)
//   everywhere. The film exchanges liquid with the meniscus through the drainage: the meniscus's particles are
//   reservoir nodes at mu_m (ThinFilmEnergy::Set_Reservoir); the suction is the difference between mu_m and the
//   film's own potential. Not modelled: the meniscus's liquid content changing (mu_m is fixed).
// - surface: the meniscus surfaces carry the same insoluble monolayer as the film (the concentration is continuous
//   across the junction), and surface passes between film and meniscus with the flow. Unrolled into the film plane,
//   each meniscus surface is a strip of width w(y) between the film and the fibre, carried by particles like the
//   film's. They differ from film particles in three ways:
//   * weightless: the liquid under the strip is the meniscus's, held by the frame;
//   * friction: a meniscus surface moving at speed u shears the meniscus's liquid, which does not slip on the fibre:
//     a stress ~ mu u / depth. Along the border (the surface sliding parallel to the fibre: the liquid under it goes
//     along, sheared only towards the fibre) a closure, a modelling choice: zeta = c_drag mu / r per face, the
//     meniscus's depth taken as its radius; FilmMatState::drag = 2 zeta a_i. Across the junction (towards the fibre or
//     away from it) the liquid under the surface must flow back, in the wedge between the surface and the film's
//     mid-plane: lubrication, 3 mu / d(s) per face with d(s) the local depth (Meniscus_Section::Wedge_Integral),
//     largest at the junction (d = h / 2) - anisotropic friction FilmMatState::drag_n, drag_dir (MeniscusParams::wedge);
//   * thickness: the drainage treats them as the reservoir; for the mechanics (inertia, disjoining stress, squeeze
//     energy) and for the film they become when they cross the junction, they carry the thickness the film has at
//     the junction (kernel-weighted mean over the film particles within the kernel radius, extended into the strip
//     sweep by sweep): the meniscus's own thickness profile is not represented, its curvature is in mu_m.
//   With a soluble surfactant (FilmMechParams::sol, Solution.h) the surface is not a conserved monolayer: a meniscus
//   surface sits over the meniscus's deep liquid and exchanges with it within us-ms (it needs a sublayer only
//   Delta Gamma / c ~ 50 nm deep), so the strip is a surfactant bath at the menisci's concentration (sbath): constant
//   tension, its surfactant (and the liquid it carries into the film) at c_bath.
// - the fibre: fixed particles (the frame) outside the inner rectangle, at least 2 kernel radii thick. Its surface
//   is a wall for the meniscus surface (the fibre's wetting film is thin: its surface does not move). For the
//   mechanics they take the state of the strip next to them (FilmMech_Stepper::ghost_frame = 2), for the drainage
//   they are reservoir nodes at mu_m like the strip.
// A border without a meniscus (MeniscusParams::left etc. false) has the film attached directly to the frame.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "vf/thinfilm/ThinFilm.h"
#include "vf/thinfilm/ThinFilmMechanics.h"
#include "vf/thinfilm/Units.h"
#include "vf/thinfilm/InterfaceLedger.h"

namespace vf::thinfilm {

struct MeniscusParams {
	real fibre_radius = 125e-6; // m: the frame fibre's radius (measured: D / 2)
	real r_ref = 0.5e-3;        // m: the meniscus radius at the height y_ref (the frame's liquid content; an input)
	real y_ref = 0.;            // m above the frame's inner bottom edge
	real c_drag = 1.;           // friction closure zeta = c_drag mu / r per face, for motion along the border
	bool wedge = true;          // motion across the junction: lubrication friction 3 mu / d(s) in the meniscus wedge
	bool column = true;         // soluble surfactant: a meniscus surface in equilibrium with the liquid under it (false: a bath)
	bool left = true, right = true, bottom = true, top = true; // borders with a meniscus
};

// meniscus cross-section geometry (any consistent length unit): a fibre of radius a, an arc of radius r
struct Meniscus_Section {
	static real Arc_Angle(real a, real r) { return std::acos(r / (a + r)); }
	static real Width(real a, real r) { return r * Arc_Angle(a, r); }               // surface from junction to fibre
	static real Junction(real a, real r) { return std::sqrt(a * a + 2. * a * r); } // fibre axis to junction
	// liquid area between the film plane, the arc and the fibre, per face: the right triangle (fibre axis, junction,
	// arc centre) minus the arc's sector and the fibre's sector (their angles add to pi / 2)
	// Lubrication friction of a meniscus surface moving across its junction (towards the fibre or away from it): the
	// liquid under it, down to the film's mid-plane (a plane of symmetry), d(s) = h_j / 2 + r (1 - cos(s / r)) deep at
	// arc length s from the junction (h_j the film's thickness there), is dragged along by the surface and has to flow
	// back (the fibre closes the meniscus): u(z) = U (3 z^2 / d^2 - 1) / 2 (zero flux, zero shear at the mid-plane), a
	// stress 3 mu U / d per face. The integral of 1 / d over [s1, s2] (alpha = h_j / 2; lengths in one unit):
	static real Wedge_Integral(real r, real alpha, real s1, real s2)
	{
		const real k = std::sqrt((alpha + 2. * r) / alpha), c = 2. * r / std::sqrt(alpha * alpha + 2. * alpha * r);
		auto I = [&](real s) { return c * std::atan(k * std::tan(0.5 * s / r)); };
		return I(s2) - I(s1);
	}
	static real Liquid_Area(real a, real r)
	{
		const real xc = Junction(a, r), th = Arc_Angle(a, r);
		return 0.5 * xc * r - 0.5 * r * r * th - 0.5 * a * a * (0.5 * kPi - th);
	}
};

// Hex lattice filling the box, with the particles outside the inner rectangle [x0, x1] x [y0, y1] fixed (the frame,
// `frame` thick on every side). Periodic in x (the two side frames meet across it), walls in y.
struct Framed_Lattice {
	real x0 = 0., x1 = 0., y0 = 0., y1 = 0.;
	FilmBox box;
	std::vector<Vec3> x;
	std::vector<char> fixed;
	Framed_Lattice(real W, real H, real dx, real frame)
	{
		const real dy = dx * std::sqrt(3.) / 2.;
		const int nfx = (int)std::ceil(frame / dx), nfy = (int)std::ceil(frame / dy);
		const int nix = std::max(1, (int)std::lround(W / dx)), niy = std::max(1, (int)std::lround(H / dy));
		const int nx = nix + 2 * nfx, ny = niy + 2 * nfy;
		box.Lx = nx * dx; box.Ly = ny * dy; box.periodic_y = false;
		x0 = nfx * dx; x1 = x0 + nix * dx; y0 = nfy * dy; y1 = y0 + niy * dy;
		for (int r = 0; r < ny; r++) for (int c = 0; c < nx; c++) {
			const Vec3 p((c + 0.25 + 0.5 * (r & 1)) * dx, (r + 0.5) * dy, 0.);
			x.push_back(p);
			fixed.push_back(!(p[0] > x0 && p[0] < x1 && p[1] > y0 && p[1] < y1));
		}
	}
};

class Framed_Menisci {
public:
	// frame's inner rectangle [x0, x1] x [y0, y1] in code lateral units (the fibres' surfaces: inside is film or
	// meniscus); gravity along -y; heights in the code's absolute coordinate (as ThinFilmGraph::height)
	Framed_Menisci(const FilmMaterial& M, const CodeParams& cp, const MeniscusParams& mp, real x0, real x1, real y0, real y1)
		: mp_(mp), x0_(x0), x1_(x1), y0_(y0), y1_(y1)
	{
		const FilmScales& s = cp.sc;
		eps_ = s.H / s.L;
		gamma_ = M.gamma_ref() / s.S;
		rho_g_ = cp.tf.rho_g;
		a_ = mp.fibre_radius / s.L;
		mu_m_ = -gamma_ * eps_ / (mp.r_ref / s.L) + rho_g_ * (y0 + mp.y_ref / s.L);
		zcoef_ = 2. * mp.c_drag * M.mu * s.L / (s.S * s.T);
		mu_ = M.mu; Lm_ = s.L; Hm_ = s.H; zconv_ = s.L * s.L / (s.S * s.T);
	}
	// the menisci exist only below the atmosphere's pressure everywhere (p = mu_m - rho_g y < 0 up to the bottom)
	bool Valid() const { return mu_m_ - rho_g_ * y0_ < 0.; }
	real Mu() const { return mu_m_; }                       // chemical potential (code, as ThinFilmEnergy's pressure)
	real Pressure(real y) const { return mu_m_ - rho_g_ * Clamp_Y(y); } // Laplace pressure -gamma / r (code)
	real Radius(real y) const { return gamma_ * eps_ / -Pressure(y); }  // code lateral units
	real Width(real y) const { return Meniscus_Section::Width(a_, Radius(y)); }
	real Zeta(real y) const { return zcoef_ / Radius(y); }  // friction per area, both faces (code)
	real Fibre_Radius() const { return a_; }
	const MeniscusParams& Params() const { return mp_; }
	// 0: film, 1: meniscus surface (also a mobile particle pushed into the frame), 2: frame (fixed)
	int Zone(const Vec3& p, bool fixed) const
	{
		if (fixed) return 2;
		const real x = p[0], y = p[1];
		if (!(x > x0_ && x < x1_ && y > y0_ && y < y1_)) return 1;
		const real wy = Width(y);
		if ((mp_.left && x - x0_ < wy) || (mp_.right && x1_ - x < wy)) return 1;
		if ((mp_.bottom && y - y0_ < Width(y0_)) || (mp_.top && y1_ - y < Width(y1_))) return 1;
		return 0;
	}
	void Zones(const FilmMatState& st, std::vector<int>& z) const
	{
		z.resize(st.Size());
		for (int i = 0; i < st.Size(); i++) z[i] = Zone(st.x[i], !st.fixed.empty() && st.fixed[i]);
	}
	// Before a mechanics step and before the drainage: weights (gscale), friction (drag) and the meniscus's and
	// frame's thickness (V) from the film at the junction. G: the graph at the current positions. Returns the liquid
	// the meniscus gave to (+) or took from (-) the strip and frame particles by the thickness reset.
	// With a soluble surfactant (sol, and st.B set), Set_Soluble; dQ (if given) gets the surfactant the menisci gave (+)
	// to or took (-) from the free particles.
	real Prepare(FilmMatState& st, const ThinFilmGraph& G, const Kernel2D& ker, const FilmBox& box, const Solubility* sol = nullptr,
		real* dQ = nullptr) const
	{
		const int n = st.Size();
		std::vector<int> z;
		Zones(st, z);
		st.gscale.assign(n, 1.); st.drag.assign(n, 0.);
		const bool soluble = sol && sol->On() && (int)st.B.size() == n;
		std::vector<real> Vl0; // the liquid each particle's dissolved surfactant was in equilibrium with
		if (soluble) { Vl0.resize(n); for (int i = 0; i < n; i++) Vl0[i] = st.Vsol(i); }
		std::vector<real> h(n);
		std::vector<char> known(n);
		for (int i = 0; i < n; i++) {
			h[i] = st.V[i] / G.a[i];
			known[i] = z[i] == 0;
			if (z[i] == 0) continue;
			st.gscale[i] = 0.;
			if (z[i] == 1) st.drag[i] = Zeta(st.x[i][1]) * G.a[i];
		}
		// zero-gradient extension of the film's thickness: each sweep gives the particles next to known ones the
		// kernel-weighted mean of those
		const PairList& P = G.pairs;
		std::vector<real> w(P.Size());
		for (int q = 0; q < P.Size(); q++) w[q] = ker.W(box.Wrap(st.x[P.pj[q]] - st.x[P.pi[q]]).norm());
		std::vector<real> ws(n), hs(n);
		for (int sweep = 0; sweep < 64; sweep++) {
			std::fill(ws.begin(), ws.end(), 0.); std::fill(hs.begin(), hs.end(), 0.);
			for (int q = 0; q < P.Size(); q++) {
				const int i = P.pi[q], j = P.pj[q];
				if (known[i] == known[j] || w[q] <= 0.) continue;
				const int u = known[i] ? j : i, k = known[i] ? i : j;
				ws[u] += w[q]; hs[u] += w[q] * h[k];
			}
			bool any = false;
			for (int i = 0; i < n; i++) if (!known[i] && ws[i] > 0.) { h[i] = hs[i] / ws[i]; known[i] = 1; any = true; }
			if (!any) break;
		}
		real dV = 0.;
		for (int i = 0; i < n; i++) {
			if (z[i] == 0) continue;
			const real Vn = h[i] * G.a[i];
			dV += Vn - st.V[i]; st.V[i] = Vn;
		}
		const real dq = soluble ? Set_Soluble(st, G, z, h, Vl0, *sol) : 0.;
		if (dQ) *dQ = dq;
		if (mp_.wedge) Set_Wedge_Friction(st, G, z, h);
		return dV;
	}
	// The nearest junction whose strip the point is in: its distance s (towards the fibre; >= 0), the meniscus radius r
	// and strip width w there, and the junction's in-plane normal (code lateral units). False: in no strip.
	bool Strip_Geometry(const Vec3& p, real& s, real& r, real& w, Vec3& nrm) const
	{
		const real x = p[0], y = p[1], yc = Clamp_Y(y);
		s = 1e300; r = w = 0.; nrm = Vec3::Zero();
		auto cand = [&](bool on, real sb, real rb, real wb, const Vec3& nb) { if (on && sb >= 0. && sb < s) { s = sb; r = rb; w = wb; nrm = nb; } };
		cand(mp_.left, x0_ + Width(yc) - x, Radius(yc), Width(yc), Vec3(1., 0., 0.));
		cand(mp_.right, x - (x1_ - Width(yc)), Radius(yc), Width(yc), Vec3(1., 0., 0.));
		cand(mp_.bottom, y0_ + Width(y0_) - y, Radius(y0_), Width(y0_), Vec3(0., 1., 0.));
		cand(mp_.top, y - (y1_ - Width(y1_)), Radius(y1_), Width(y1_), Vec3(0., 1., 0.));
		return s < 1e300;
	}
	// Soluble surfactant. The frame (fibre) particles are a surfactant bath at sol.c_bath (sbath; a wall, their content
	// no one's). A meniscus surface is in local equilibrium with the meniscus's liquid under it, down to the film's
	// mid-plane: 2 d(s) per area, d(s) = h_j / 2 + r (1 - cos(s / r)) (MeniscusParams::column; Vs): film-like at the
	// junction (d = h_j / 2: a reservoir no deeper than the film's, so an elasticity like the film's), bath-like
	// towards the fibre (d ~ r: 1e2 times deeper). When a particle's column deepens (moving towards the fibre, or
	// entering the strip from the film) meniscus liquid at c_bath joins it; when it gets shallower its liquid leaves at
	// its own concentration. With column = false every meniscus particle is a bath (sbath; S = theta_b a, B = c_bath V).
	// Returns the surfactant the menisci gave to (+) or took from (-) the free particles.
	real Set_Soluble(FilmMatState& st, const ThinFilmGraph& G, const std::vector<int>& z, const std::vector<real>& h, const std::vector<real>& Vl0,
		const Solubility& sol) const
	{
		const int n = st.Size();
		st.sbath.assign(n, 0);
		if ((int)st.Vs.size() != n) st.Vs.assign(n, 0.);
		const real thb = sol.Theta_Bath();
		real dq = 0.;
		for (int i = 0; i < n; i++) {
			const bool fixed = !st.fixed.empty() && st.fixed[i];
			if (z[i] == 2 || (z[i] == 1 && !mp_.column)) {
				const real Q0 = 2. * st.S[i] + st.B[i];
				st.sbath[i] = 1; st.Vs[i] = 0.;
				st.S[i] = thb * G.a[i]; st.B[i] = sol.c_bath * st.V[i];
				if (!fixed) dq += 2. * st.S[i] + st.B[i] - Q0;
				continue;
			}
			real Vl = st.V[i];
			st.Vs[i] = 0.;
			real s, r, w;
			Vec3 nrm;
			if (z[i] == 1 && Strip_Geometry(st.x[i], s, r, w, nrm)) {
				const real d = 0.5 * h[i] * Hm_ + r * Lm_ * (1. - std::cos(std::min(s, w) / r)); // m
				st.Vs[i] = Vl = 2. * d / Hm_ * G.a[i];
			}
			const real B0 = st.B[i];
			if (!(Vl0[i] > 0.)) continue;
			const real Bn = Vl > Vl0[i] ? B0 + (Vl - Vl0[i]) * sol.c_bath : B0 / Vl0[i] * Vl;
			st.B[i] = Bn;
			if (!fixed) dq += Bn - B0;
		}
		return dq;
	}
	// The friction of the meniscus surfaces across their junction (Meniscus_Section::Wedge_Integral) on top of the
	// isotropic closure: each strip particle takes the wedge stress averaged over its span (sqrt(a) wide, centred at its
	// distance s from the nearest junction, clipped to the strip), along that junction's normal (drag_n, drag_dir).
	void Set_Wedge_Friction(FilmMatState& st, const ThinFilmGraph& G, const std::vector<int>& z, const std::vector<real>& h) const
	{
		const int n = st.Size();
		st.drag_n.assign(n, 0.); st.drag_dir.assign(n, Vec3::Zero());
		for (int i = 0; i < n; i++) {
			if (z[i] != 1) continue;
			real s, r, w;
			Vec3 nrm;
			if (!Strip_Geometry(st.x[i], s, r, w, nrm)) continue;
			const real span = std::sqrt(G.a[i]), s1 = std::max<real>(0., s - 0.5 * span), s2 = std::min(w, s + 0.5 * span);
			if (!(s2 > s1)) continue;
			const real alpha = 0.5 * h[i] * Hm_; // half the film's thickness at the junction (the strip carries it), m
			const real Iw = Meniscus_Section::Wedge_Integral(r * Lm_, alpha, s1 * Lm_, s2 * Lm_);
			const real zeta = 2. * 3. * mu_ * Iw / (span * Lm_) * zconv_; // both faces, averaged over the span, code units
			st.drag_n[i] = std::max<real>(0., zeta * G.a[i] - st.drag[i]);
			st.drag_dir[i] = nrm;
		}
	}
	// reservoir chemical potentials for the drainage: mu_m on the meniscus and the frame, NaN on the film
	void Reservoir(const FilmMatState& st, std::vector<real>& res) const
	{
		res.assign(st.Size(), std::numeric_limits<real>::quiet_NaN());
		for (int i = 0; i < st.Size(); i++) if (Zone(st.x[i], !st.fixed.empty() && st.fixed[i]) != 0) res[i] = mu_m_;
	}

private:
	real Clamp_Y(real y) const { return std::clamp(y, y0_, y1_); }
	MeniscusParams mp_;
	real x0_, x1_, y0_, y1_;
	real eps_ = 1., gamma_ = 0., rho_g_ = 0., a_ = 0., mu_m_ = 0., zcoef_ = 0.;
	real mu_ = 0., Lm_ = 1., Hm_ = 1., zconv_ = 1.; // SI viscosity; metres per lateral and per thickness unit; friction to code
};

// One step of a framed film with menisci (Lie splitting as FilmCoupled_Stepper): the menisci's setup, mechanics, the
// setup again at the new positions, drainage with the menisci and the frame as the reservoir (a step that does not
// converge is redone in 4, then 16 substeps), surface diffusion, remeshing (the frame neither absorbs nor is
// remapped: the free particles keep the monolayer's total exactly). With a soluble surfactant (cs.mech.prm.sol, st.B
// set) the meniscus surfaces are a surfactant bath (Framed_Menisci::Prepare), the drainage carries the dissolved
// surfactant into and out of the menisci (Solute_Transport), and each particle's amount is split again between its
// surfaces and its liquid after every substep (Repartition).
struct Framed_Film_Stepper {
	FilmCoupled_Stepper cs;
	RemeshParams rp;
	bool remesh = true;
	// Only for opt-in regional remapping: restore the same instantaneous local
	// equilibrium already assumed by the soluble mechanics after areas change.
	// This is not a finite-rate exchange model; each particle keeps 2S+B.
	bool equilibrate_region_remap = true;
	real drained = 0.;           // liquid the film zone gave to the menisci through the drainage (code units, cumulated)
	real exchanged = 0.;         // soluble surfactant the menisci gave (+) to the free particles: resets and transport
	int inserted = 0, merged = 0;
	Framed_Film_Stepper() { cs.mech.ghost_frame = 2; rp.absorb_fixed = false; rp.remap_fixed = false; }
	struct Stats { NewtonStats mech, drain, diff; };
	bool Soluble(const FilmMatState& st) const { return cs.Soluble(st); }
	Stats Step(FilmMatState& st, const Framed_Menisci& men, real dt, InterfaceLedger* audit = nullptr)
	{
		Stats s;
		const Solubility& so = cs.mech.prm.sol;
		const bool sol = Soluble(st);
		auto snapshot = [&](const char* stage) {
			if (audit) { std::vector<int> zones;men.Zones(st,zones);audit->Capture(stage,st,zones); }
		};
		if(audit){audit->Clear();snapshot("begin");}
		real dq = 0.;
		{
			const ThinFilmGraph G = cs.Graph(st);
			const real resetV=men.Prepare(st, G, cs.mech.kernel, cs.mech.box, &so, &dq);
			if(audit)audit->Scalar("prepare1_return",resetV,dq);
			snapshot("prepare1");
			exchanged += dq;
			if (sol) Repartition(st, G.a, so);
			snapshot("partition1");
		}
		if(audit){std::vector<int> zones;men.Zones(st,zones);audit->Begin_Material(st,zones);}
		s.mech = cs.mech.Step(st, dt);
		const ThinFilmGraph G = cs.Graph(st);
		std::vector<int> z;
		men.Zones(st, z);
		if(audit)audit->End_Material(z);
		snapshot("mechanics");
		const real resetV=men.Prepare(st, G, cs.mech.kernel, cs.mech.box, &so, &dq);
		if(audit)audit->Scalar("prepare2_return",resetV,dq);
		snapshot("prepare2");
		exchanged += dq;
		if (sol) Repartition(st, G.a, so);
		snapshot("partition2");
		auto film_volume = [&]() { real v = 0.; for (int i = 0; i < st.Size(); i++) if (z[i] == 0) v += st.V[i]; return v; };
		const real Vf0 = film_volume();
		std::vector<real> res;
		men.Reservoir(st, res);
		cs.drain.prm = cs.mech.prm.tf;
		cs.drain.reservoir = &res;
		const std::vector<real> V_start = st.V;
		ThinFilmEnergy drain_energy(G,cs.drain.prm);drain_energy.Set_Reservoir(&res);
		if(audit)audit->Scalar("drain_energy_before",0.,0.,drain_energy.Value(st.V));
		std::vector<real> flux, fsub;
		const bool want_flux=sol || audit;
		s.drain = cs.drain.Step(G, st.V, dt, want_flux ? &flux : nullptr);
		for (int sub = 4; !s.drain.converged && sub <= 16; sub *= 4) {
			st.V = V_start;
			if (want_flux) flux.assign(G.pairs.Size(), 0.);
			bool ok = true;
			for (int k = 0; k < sub && ok; k++) {
				s.drain = cs.drain.Step(G, st.V, dt / sub, want_flux ? &fsub : nullptr);
				ok = s.drain.converged;
				if (want_flux) for (size_t q = 0; q < flux.size(); q++) flux[q] += fsub[q];
			}
		}
		cs.drain.reservoir = nullptr;
		if(audit){
			audit->Drainage(G,flux,z);real work=0.;
			for(int i=0;i<G.N;++i)if(!std::isnan(res[i]))work+=res[i]*(st.V[i]-V_start[i]);
			audit->Scalar("drain_energy_after",0.,0.,drain_energy.Value(st.V),work);
		}
		snapshot("drainage");
		drained += Vf0 - film_volume();
		if (sol) {
			std::vector<real> cres(st.Size(), std::numeric_limits<real>::quiet_NaN());
			for (int i = 0; i < st.Size(); i++) if (!std::isnan(res[i])) cres[i] = so.c_bath;
			const real transportQ=Solute_Transport(G, flux, st.V, st.B, &cres);
			exchanged += transportQ;
			if(audit){audit->Transport(G,flux,st,z,cres);audit->Scalar("transport_return",0.,transportQ);}
			snapshot("transport");
			Repartition(st, G.a, so);
		}
		else snapshot("transport");
		snapshot("partition3");
		cs.diff.surf = cs.mech.prm.surf;
		if (!sol) s.diff = cs.diff.Step(G, st.S, dt);
		else s.diff.converged = true;
		snapshot("diffusion");
		if (remesh) {
			rp.kernel_radius = cs.mech.kernel.Radius();
			auto graph_of = [&](const std::vector<Vec3>& xs) { FilmMatState tmp; tmp.x = xs; return cs.Graph(tmp); };
			const ThinFilmGraph Gr = cs.Graph(st);
			const RemeshStats rs = Remesh_Film(st, Gr, rp, &cs.mech.box, nullptr, graph_of);
			inserted += rs.splits; merged += rs.merges;
			if(audit)audit->Remesh_Event(rs.splits,rs.merges);
		}
		snapshot("remesh_raw");
		if(remesh && rp.remap_region && sol && equilibrate_region_remap){
			const ThinFilmGraph post=cs.Graph(st);std::vector<int> zones;men.Zones(st,zones);
			st.sbath.resize(st.Size());
			for(int i=0;i<st.Size();++i)st.sbath[i]=zones[i]==2||(zones[i]==1&&!men.Params().column);
			Repartition(st,post.a,so);
			snapshot("remesh_equilibrium");
		}
		snapshot("remesh");
		return s;
	}
};

} // namespace vf::thinfilm
