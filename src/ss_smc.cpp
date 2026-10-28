// =============================================================================
// ss_smc.cpp -- marginal likelihood p(x | K) of the beta-binomial mixture by
// adaptive-tempering SMC (stepping-stone with resampling) along the prior path,
// and the same fixed-beta MCMC kernel run as ordinary chains (bbmm_mcmc_cpp).
// =============================================================================
//
// MODEL (K fixed; components, not non-empty clusters -- empty ones allowed)
//   pi        ~ Dirichlet(gamma, ..., gamma)          integrated out
//   z_i | pi  ~ Categorical(pi)
//   p_kq      ~ Beta(p_prior_a, p_prior_b)
//   nu_k      ~ Gamma(shape = nu_shape, rate = nu_rate), truncated to
//               nu_k >= nu_min (hard floor; nu_min = 0 is no floor)  one per component
//   or, with nu_pool (pooled dispersion, added 2026-10-16):
//   mu        ~ Normal(pool_m0, pool_s0^2)                          one per particle
//   log nu_k | mu ~ Normal(mu, pool_tau^2)                          every component, empties too
//               (no floor; the shared mu pulls every nu_k toward a common centre, so one
//               component cannot drift to a low-nu catch-all without paying for it)
//   X_iq | z_i = k ~ BetaBinomial(M_iq, nu_k p_kq, nu_k (1 - p_kq))
// Missing cells contribute nothing. The binomial coefficient is omitted, so
// log Z here differs from one with it included by a data-only constant that is
// the same for every K.
//
// LIKELIHOOD POWER zeta (coarsened posterior, Miller & Dunson 2019; the
// same power as the EM and MCMC fits in this package). The target is
//   prior(z, p, nu) * L(x | z, p, nu)^zeta,
// complete-data tempering with the priors untempered, so the evidence is
//   Z(zeta) = int prior * L^zeta.
// The SMC path below already passes through every power L^beta, so this is
// the same path stopped at beta = zeta. zeta = 1 (default) is the ordinary
// model.
//
// TEMPERING (option (a)): the complete-data likelihood is raised to beta,
//   pi_beta(z, p, nu) ∝ L(x | z, p, nu)^beta  p(z)  p(p)  p(nu),
//   L(x | z, p, nu) = prod_i f(x_i | p_{z_i}, nu_{z_i}).
// beta = 0 is the prior, which is sampled exactly.
//
// SMC (Del Moral, Doucet & Jasra 2006; adaptive schedule as in Zhou, Johansen
// & Aston 2016). With normalised weights W_j and log-likelihoods l_j:
//   next beta: largest step d <= 1 - beta with conditional ESS
//              (sum_j W_j w_j)^2 / sum_j W_j w_j^2 >= cess_target, w_j = exp(d l_j);
//   log Z    += log sum_j W_j exp(d l_j)          (the stepping-stone ratio);
//   W_j      ∝ W_j exp(d l_j); resample (systematic) if ESS < threshold * P;
//   move every particle with an MCMC kernel invariant for pi_beta.
// The sum of the log increments is the stepping-stone / SMC estimate of log Z.
// The weighted mean and variance of l after each move are returned so R can
// form TI, corrected TI and generalised TI from the same run.
//
// MCMC KERNEL at fixed beta (n_sweeps times):
//   1. empty components: (p_k, nu_k) redrawn from the prior (their exact
//      conditional, since they carry no data);
//   2. n_sm split-merge proposals (fixed-K port of the MFM sampler's
//      Jain-Neal-style move, profiles drawn from the prior);
//   3. one Gibbs sweep of z with pi collapsed:
//      P(z_i = k | ...) ∝ (n_{-i,k} + gamma) f(x_i | k)^beta;
//   4. random-walk Metropolis on logit(p_kq), every (k, q) with data;
//   5. random-walk Metropolis on log(nu_k), every non-empty k;
//   6. (nu_pool only) blocked Gibbs for (mu, lambda of the empty components):
//      mu | lambda_occupied ~ Normal (conjugate; the empties are integrated
//      out, as they carry no data), then each empty lambda_k | mu ~ N(mu, tau^2).
// With nu_pool every prior draw of nu (beta = 0, empties, split-merge) is
// conditional on the particle's current mu, and every prior density of nu is
// N(log nu; mu, tau^2); mu is never changed by steps 1-5, so they stay
// invariant for the joint target conditional on mu.
// Proposal sds are global, adapted between SMC steps from the previous
// step's acceptance rates.
//
// SPLIT-MERGE AT FIXED K. Pick an ordered pair i != j uniformly.
//   z_i == z_j = c (split): needs an empty label; if none, the move is
//     rejected. Choose e uniformly among the E empty labels; draw theta_A,
//     theta_B from the prior; i -> A, j -> B, every other member of c is
//     allocated in index order with probability ∝ (n_side + gamma) f^beta.
//     A keeps label c, B takes label e (e's old parameters are discarded).
//     log alpha = beta [l(A) + l(B) - l(c)] + log p(z') - log p(z) + log E - log q_alloc.
//   z_i = a != z_j = b (merge): every member of b moves to a; theta_a and a
//     fresh theta_b are drawn from the prior. The reverse split chooses b with
//     probability 1/E' (E' = empties after the merge) and allocates with the
//     OLD theta_a, theta_b.
//     log alpha = beta [l(M) - l(a) - l(b)] + log p(z') - log p(z) + log q_rev - log E'.
//   Prior densities of drawn parameters cancel against the proposal in both.
//   log p(z) = sum_k lgamma(n_k + gamma) + const.
//
// RANDOM NUMBERS. Each particle's moves at SMC step t use a private
// mt19937_64 seeded from (seed, t, particle index), so results do not depend
// on the number of threads. Resampling uses its own stream.
// =============================================================================

#include <Rcpp.h>
#include <random>
#include <chrono>
#include "bbcore.h"

using namespace Rcpp;

namespace {

// ---------------------------------------------------------------------------
// Random numbers.
// ---------------------------------------------------------------------------
inline uint64_t splitmix64(uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

struct Rng {
  std::mt19937_64 g;
  std::normal_distribution<double> nd{0.0, 1.0};
  explicit Rng(uint64_t s) : g(s) {}
  double unif() { return (double)(g() >> 11) * 0x1.0p-53; }          // [0, 1)
  double unif_pos() { double u; do u = unif(); while (u <= 0.0); return u; }
  double norm() { return nd(g); }
  double gamma(double shape) { std::gamma_distribution<double> d(shape, 1.0); return d(g); }
  double beta(double a, double b) { const double x = gamma(a), y = gamma(b); return x / (x + y); }
  int below(int n) { return (int)std::floor(unif() * n); }
};

uint64_t stream_seed(uint64_t seed, uint64_t step, uint64_t index) {
  return splitmix64(splitmix64(splitmix64(seed) ^ step) ^ (index * 0x632BE59BD9B4E019ULL));
}

inline double log_sum_exp(const double* x, int n) {
  double m = -INFINITY;
  for (int k = 0; k < n; ++k) m = std::max(m, x[k]);
  if (!std::isfinite(m)) return m;
  double s = 0.0;
  for (int k = 0; k < n; ++k) s += std::exp(x[k] - m);
  return m + std::log(s);
}

inline double inv_logit(double eta) {
  return eta >= 0 ? 1.0 / (1.0 + std::exp(-eta)) : std::exp(eta) / (1.0 + std::exp(eta));
}
inline double logit(double p) { return std::log(p) - std::log1p(-p); }

// ---------------------------------------------------------------------------
// Settings.
// ---------------------------------------------------------------------------
struct Settings {
  int K, P, n_sweeps, n_sm, n_threads, max_steps;
  double gamma, pa, pb, nu_shape, nu_rate;
  double nu_min = 0.0;        // hard floor on nu_k: prior truncated to [nu_min, inf)
  double log_nu_tail = 0.0;   // log P(nu >= nu_min) under the untruncated Gamma
  double cess_target, resample_threshold, p_sd, nu_sd;
  bool adapt;
  bool sm_informed;          // data-informed split-merge profiles (else prior draws)
  double sm_conc0, sm_nu_ref;
  int sm_method;             // 0 prior, 1 anchor (informed), 2 Jain-Neal launch state
  int n_launch, n_p_launch;
  double sm_p_rw_sd, sm_nu_rw_sd, delta_sd, delta_nu_sd;
  double var_eta_prior = 0, var_lam_prior = 0;
  bool debug_sm = false;
  int launch_method = 1;                          // 0 "mh" (MFM port), 1 "pooled"
  int max_rounds = 1;
  bool nu_pool = false;                           // hierarchical log-normal prior on nu with shared mu
  double pool_tau = 0.25, pool_m0 = 0.0, pool_s0 = 1.0;
  double zeta = 1.0;                              // likelihood power of the target (coarsened posterior)
  double rho_target = 0.0;                        // SMC: stop the kernel rounds once corr(l before, l after) < rho_target
  uint64_t seed;
};

// nu from its prior, Gamma(nu_shape, nu_rate) truncated to [nu_min, inf), by
// rejection (at the defaults, nu_min = 20, about 1% of draws are rejected).
// With nu_pool: log nu ~ N(mu, pool_tau^2), mu the particle's current centre.
inline double draw_nu(Rng& r, const Settings& s, double mu) {
  if (s.nu_pool) return std::exp(mu + s.pool_tau * r.norm());
  const double lo = std::max(s.nu_min, 1e-8);
  for (;;) {
    const double nu = r.gamma(s.nu_shape) / s.nu_rate;
    if (nu >= lo && std::isfinite(nu)) return nu;
  }
}

// ---------------------------------------------------------------------------
// Data plus the per-component nu term. bbcore's UnitConstants assume one nu
// shared by all components; here nu varies by component, so each component
// carries its own table over the distinct trial counts of each question:
//   mt[off_m[q] + r] = lgamma(mv[q][r] + nu_k) - lgamma(nu_k).
// ---------------------------------------------------------------------------
struct Ext {
  bbcore::Data d;
  std::vector<int> off_m;     // offsets of question blocks in an m table
  int Tm = 0;
  std::vector<int> cell_gm;   // absolute m-table position of each cell
};

Ext make_ext(const int* X, const int* M, int N, int Q) {
  Ext e;
  e.d = bbcore::make_data(X, M, N, Q);
  e.off_m.resize(Q);
  for (int q = 0; q < Q; ++q) { e.off_m[q] = e.Tm; e.Tm += (int)e.d.mv[q].size(); }
  e.cell_gm.resize(e.d.n_cells());
  for (int j = 0; j < e.d.n_cells(); ++j) e.cell_gm[j] = e.off_m[e.d.cell_q[j]] + e.d.cell_mi[j];
  return e;
}

// One component's cached tables at (p_k, nu_k).
struct Comp {
  bbcore::Tables tab;
  std::vector<double> mt;
  void build(const Ext& e, const double* p, double nu) {
    const int Q = e.d.Q;
    std::vector<double> nuv(Q, nu);
    tab.build(e.d, p, nuv.data());
    mt.assign(e.Tm, 0.0);
    for (int q = 0; q < Q; ++q)
      bbcore::shifted_lgamma(nu, e.d.mv[q].data(), (int)e.d.mv[q].size(), mt.data() + e.off_m[q]);
  }
};

inline double unit_ll(const Ext& e, const Comp& c, int i) {
  const bbcore::Data& d = e.d;
  const double* t = c.tab.t.data();
  const double* m = c.mt.data();
  double s = 0.0;
  for (int j = d.cell_start[i]; j < d.cell_start[i + 1]; ++j)
    s += t[d.cell_gx[j]] + t[d.cell_gy[j]] - m[e.cell_gm[j]];
  return s;
}

// A freshly drawn parameter set, evaluated directly (no tables).
struct Proposal;
struct Fresh {
  std::vector<double> p, a, b, lga, lgb;
  double nu = 1.0, lgnu = 0.0;
  // nu from its prior; p_q ~ Beta(c_q r_q, c_q (1 - r_q)).
  void draw(Rng& r, const Settings& s, double mu, const std::vector<double>& rr, const std::vector<double>& cc) {
    const int Q = (int)rr.size();
    p.resize(Q); a.resize(Q); b.resize(Q); lga.resize(Q); lgb.resize(Q);
    nu = draw_nu(r, s, mu);
    lgnu = std::lgamma(nu);
    for (int q = 0; q < Q; ++q) {
      double v = r.beta(cc[q] * rr[q], cc[q] * (1.0 - rr[q]));
      v = std::min(std::max(v, 1e-12), 1.0 - 1e-12);
      p[q] = v; a[q] = nu * v; b[q] = nu * (1.0 - v);
      lga[q] = std::lgamma(a[q]); lgb[q] = std::lgamma(b[q]);
    }
  }
};

// ---------------------------------------------------------------------------
// Split-merge profile proposals. Each question's profile is drawn from
//   Beta(c_q r_q, c_q (1 - r_q)).
// Prior mode: r_q = a / (a + b), c_q = a + b, i.e. exactly the Beta prior, so
// the prior and proposal densities cancel in the acceptance ratio.
// Informed mode:
//   anchor (one unit i):  r_q = (x_iq + 0.5) / (m_iq + 1)
//   pool   (a set S):     r_q = (sum_S x + 0.5) / (sum_S m + 1)
//   c_q = sm_conc0 + beta * sum over the unit(s) of m nu_ref / (m + nu_ref),
// the tempered effective sample size of a beta-binomial observation with
// dispersion nu_ref. Questions with no data keep c_q = sm_conc0, r_q = 0.5.
// nu is always proposed from its prior, so its prior and proposal cancel.
// ---------------------------------------------------------------------------
struct Proposal { std::vector<double> r, c; };

Proposal prior_proposal(const Settings& s, int Q) {
  return Proposal{std::vector<double>(Q, s.pa / (s.pa + s.pb)), std::vector<double>(Q, s.pa + s.pb)};
}

// Accumulate sufficient statistics of a set of units, then form the proposal.
struct PoolAcc {
  std::vector<double> sx, sm, ess;
  explicit PoolAcc(int Q) : sx(Q, 0.0), sm(Q, 0.0), ess(Q, 0.0) {}
  void add(const Ext& e, int i, double nu_ref) {
    const bbcore::Data& d = e.d;
    for (int j = d.cell_start[i]; j < d.cell_start[i + 1]; ++j) {
      const int q = d.cell_q[j];
      const double m = d.cell_x[j] + d.cell_y[j];
      sx[q] += d.cell_x[j]; sm[q] += m; ess[q] += m * nu_ref / (m + nu_ref);
    }
  }
  Proposal make(const Settings& s, double beta) const {
    const int Q = (int)sx.size();
    if (!s.sm_informed) return prior_proposal(s, Q);
    Proposal pr{std::vector<double>(Q), std::vector<double>(Q)};
    for (int q = 0; q < Q; ++q) {
      if (sm[q] > 0) { pr.r[q] = (sx[q] + 0.5) / (sm[q] + 1.0); pr.c[q] = s.sm_conc0 + beta * ess[q]; }
      else { pr.r[q] = 0.5; pr.c[q] = s.sm_conc0; }
    }
    return pr;
  }
};

double log_dbeta(double p, double a, double b) {
  return (a - 1.0) * std::log(p) + (b - 1.0) * std::log1p(-p) - (std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
}
double log_q(const Proposal& pr, const double* p) {
  double s = 0.0;
  for (size_t q = 0; q < pr.r.size(); ++q) s += log_dbeta(p[q], pr.c[q] * pr.r[q], pr.c[q] * (1.0 - pr.r[q]));
  return s;
}
double log_prior_p(const Settings& s, const double* p, int Q) {
  double t = 0.0;
  for (int q = 0; q < Q; ++q) t += log_dbeta(p[q], s.pa, s.pb);
  return t;
}

inline double unit_ll_fresh(const Ext& e, const Fresh& f, int i) {
  const bbcore::Data& d = e.d;
  double s = 0.0;
  for (int j = d.cell_start[i]; j < d.cell_start[i + 1]; ++j) {
    const int q = d.cell_q[j], x = d.cell_x[j], y = d.cell_y[j];
    s += std::lgamma(x + f.a[q]) - f.lga[q] + std::lgamma(y + f.b[q]) - f.lgb[q]
         - (std::lgamma(x + y + f.nu) - f.lgnu);
  }
  return s;
}

// ---------------------------------------------------------------------------
// Particle and per-thread workspace.
// ---------------------------------------------------------------------------
struct Particle {
  std::vector<int> z;        // 0-based labels
  std::vector<double> p;     // K x Q, row-major: p[k * Q + q]
  std::vector<double> nu;    // K
  double mu = 0.0;           // nu_pool: centre of log nu_k (unused otherwise)
  double ll = 0.0;           // log L(x | z, p, nu) at the current state
};

struct MoveStats {
  double p_acc = 0, p_try = 0, nu_acc = 0, nu_try = 0;
  double split_acc = 0, split_try = 0, merge_acc = 0, merge_try = 0, split_noempty = 0;
};

struct Work {
  std::vector<Comp> comps;
  std::vector<double> L;                 // N x K, L[i + N k]
  std::vector<int> nk;
  std::vector<bbcore::Histogram> H;      // success/failure histograms
  std::vector<std::vector<int>> hm;      // per-component trial-count histograms (length Tm)
  std::vector<std::vector<int>> cells;   // member cells per (k, q): cells[k * Q + q]
  std::vector<double> scratch;
  std::vector<double> lp;
  Work(const Ext& e, int K) : comps(K), L((size_t)e.d.N * K), nk(K), H(K), hm(K),
    cells((size_t)K * e.d.Q), lp(K) {
    int mx = 0;
    for (int q = 0; q < e.d.Q; ++q)
      mx = std::max(mx, std::max(e.d.Dx(q) + e.d.Dy(q), (int)e.d.mv[q].size()));
    scratch.resize(mx + 1);
  }
};


void draw_component_prior(Particle& pt, int k, const Settings& s, Rng& r, int Q) {
  pt.nu[k] = draw_nu(r, s, pt.mu);
  for (int q = 0; q < Q; ++q)
    pt.p[(size_t)k * Q + q] = std::min(std::max(r.beta(s.pa, s.pb), 1e-12), 1.0 - 1e-12);
}

void count_nk(const Particle& pt, Work& w, int K) {
  std::fill(w.nk.begin(), w.nk.end(), 0);
  for (int zi : pt.z) w.nk[zi]++;
  (void)K;
}

void build_comp(const Ext& e, const Particle& pt, Work& w, int k) {
  const int Q = e.d.Q, N = e.d.N;
  w.comps[k].build(e, pt.p.data() + (size_t)k * Q, pt.nu[k]);
  double* col = w.L.data() + (size_t)N * k;
  for (int i = 0; i < N; ++i) col[i] = unit_ll(e, w.comps[k], i);
}

void build_all(const Ext& e, const Particle& pt, Work& w, int K) {
  for (int k = 0; k < K; ++k) build_comp(e, pt, w, k);
}

double complete_ll(const Ext& e, const Particle& pt, const Work& w) {
  const int N = e.d.N;
  double s = 0.0;
  for (int i = 0; i < N; ++i) s += w.L[i + (size_t)N * pt.z[i]];
  return s;
}

void build_members(const Ext& e, const Particle& pt, Work& w, int K) {
  const bbcore::Data& d = e.d;
  const int Q = d.Q;
  for (int k = 0; k < K; ++k) {
    w.H[k].reset(d);
    w.hm[k].assign(e.Tm, 0);
    for (int q = 0; q < Q; ++q) w.cells[(size_t)k * Q + q].clear();
  }
  for (int i = 0; i < d.N; ++i) {
    const int k = pt.z[i];
    w.H[k].update(d, i, +1);
    for (int j = d.cell_start[i]; j < d.cell_start[i + 1]; ++j) {
      w.hm[k][e.cell_gm[j]]++;
      w.cells[(size_t)k * Q + d.cell_q[j]].push_back(j);
    }
  }
}

// Profile-dependent part of component k's log-likelihood in question q at (p, nu).
double col_sum(const Ext& e, Work& w, int k, int q, double p, double nu) {
  const int Q = e.d.Q;
  const auto& cq = w.cells[(size_t)k * Q + q];
  if (cq.empty()) return 0.0;
  if (bbcore::prefer_histogram(e.d, q, (int)cq.size()))
    return bbcore::column_sum_hist(e.d, w.H[k], q, p, nu, w.scratch.data());
  return bbcore::column_sum_direct(e.d, cq, p, nu);
}

// nu-only part of component k's log-likelihood: -sum over member cells of
// [lgamma(m + nu) - lgamma(nu)].
double m_part(const Ext& e, Work& w, int k, double nu) {
  const bbcore::Data& d = e.d;
  const int Q = d.Q;
  double s = 0.0;
  const double lgnu = std::lgamma(nu);
  for (int q = 0; q < Q; ++q) {
    const auto& cq = w.cells[(size_t)k * Q + q];
    if (cq.empty()) continue;
    const int Dm = (int)d.mv[q].size();
    if ((double)Dm < 2.5 * cq.size()) {
      bbcore::shifted_lgamma(nu, d.mv[q].data(), Dm, w.scratch.data());
      const int* h = w.hm[k].data() + e.off_m[q];
      for (int r = 0; r < Dm; ++r) if (h[r]) s -= h[r] * w.scratch[r];
    } else {
      for (int j : cq) s -= std::lgamma(d.cell_x[j] + d.cell_y[j] + nu) - lgnu;
    }
  }
  return s;
}

// ---------------------------------------------------------------------------
// Split-merge (see header). Requires w.L, w.nk current; keeps them current.
// Profiles come from the proposals above; with sm_informed = false they are
// prior draws and the prior/proposal terms below cancel exactly.
// ---------------------------------------------------------------------------
void split_merge(const Ext& e, Particle& pt, Work& w, const Settings& s, double beta,
                 Rng& r, MoveStats& st) {
  const int N = e.d.N, Q = e.d.Q, K = s.K;
  const double g = s.gamma;
  if (N < 2) return;
  const int i = r.below(N);
  int j = r.below(N - 1); if (j >= i) ++j;
  const int ci = pt.z[i], cj = pt.z[j];

  // Anchor proposals for the two sides of a split (same in both directions).
  PoolAcc ai(Q), aj(Q);
  ai.add(e, i, s.sm_nu_ref); aj.add(e, j, s.sm_nu_ref);
  const Proposal qA = ai.make(s, beta), qB = aj.make(s, beta);

  if (ci == cj) {                                       // ---------------- split
    st.split_try += 1;
    std::vector<int> empties;
    for (int k = 0; k < K; ++k) if (w.nk[k] == 0) empties.push_back(k);
    if (empties.empty()) { st.split_noempty += 1; return; }
    const int E = (int)empties.size();
    const int el = empties[r.below(E)];
    Fresh A, B; A.draw(r, s, pt.mu, qA.r, qA.c); B.draw(r, s, pt.mu, qB.r, qB.c);
    PoolAcc pool(Q);
    double sumA = unit_ll_fresh(e, A, i), sumB = unit_ll_fresh(e, B, j);
    double ll_old = w.L[i + (size_t)N * ci] + w.L[j + (size_t)N * ci];
    pool.add(e, i, s.sm_nu_ref); pool.add(e, j, s.sm_nu_ref);
    int nA = 1, nB = 1;
    double logq = 0.0;
    std::vector<int> toB; toB.push_back(j);
    for (int l = 0; l < N; ++l) {
      if (l == i || l == j || pt.z[l] != ci) continue;
      ll_old += w.L[l + (size_t)N * ci];
      pool.add(e, l, s.sm_nu_ref);
      const double la = unit_ll_fresh(e, A, l), lb = unit_ll_fresh(e, B, l);
      const double wa = std::log(nA + g) + beta * la, wb = std::log(nB + g) + beta * lb;
      const double mx = std::max(wa, wb);
      const double norm = mx + std::log(std::exp(wa - mx) + std::exp(wb - mx));
      if (r.unif() < std::exp(wa - norm)) { ++nA; sumA += la; logq += wa - norm; }
      else { ++nB; sumB += lb; logq += wb - norm; toB.push_back(l); }
    }
    const int nc = nA + nB;
    const double* pc = pt.p.data() + (size_t)ci * Q;
    const Proposal qM = pool.make(s, beta);             // reverse merge's proposal for theta_c
    const double lprop = log_prior_p(s, A.p.data(), Q) + log_prior_p(s, B.p.data(), Q)
      - log_prior_p(s, pc, Q) - log_q(qA, A.p.data()) - log_q(qB, B.p.data()) + log_q(qM, pc);
    const double lpz = std::lgamma(nA + g) + std::lgamma(nB + g) - std::lgamma(nc + g) - std::lgamma(g);
    const double la = beta * (sumA + sumB - ll_old) + lpz + std::log((double)E) - logq + lprop;
    if (std::log(r.unif_pos()) < la) {
      for (int l : toB) pt.z[l] = el;
      for (int q = 0; q < Q; ++q) { pt.p[(size_t)ci * Q + q] = A.p[q]; pt.p[(size_t)el * Q + q] = B.p[q]; }
      pt.nu[ci] = A.nu; pt.nu[el] = B.nu;
      w.nk[ci] = nA; w.nk[el] = nB;
      build_comp(e, pt, w, ci); build_comp(e, pt, w, el);
      st.split_acc += 1;
    }
    return;
  }

  // ---------------------------------------------------------------- merge
  st.merge_try += 1;
  const int a = ci, b = cj;
  PoolAcc pool(Q);
  for (int l = 0; l < N; ++l) if (pt.z[l] == a || pt.z[l] == b) pool.add(e, l, s.sm_nu_ref);
  const Proposal qM = pool.make(s, beta);
  const Proposal q0 = prior_proposal(s, Q);
  Fresh Mg, Bn; Mg.draw(r, s, pt.mu, qM.r, qM.c); Bn.draw(r, s, pt.mu, q0.r, q0.c);
  double ll_new = 0.0, ll_old = 0.0, logq_rev = 0.0;
  int nA = 1, nB = 1;
  for (int l = 0; l < N; ++l) {
    const int zl = pt.z[l];
    if (zl != a && zl != b) continue;
    ll_new += unit_ll_fresh(e, Mg, l);
    const double La = w.L[l + (size_t)N * a], Lb = w.L[l + (size_t)N * b];
    ll_old += (zl == a) ? La : Lb;
    if (l == i || l == j) continue;
    const double wa = std::log(nA + g) + beta * La, wb = std::log(nB + g) + beta * Lb;
    const double mx = std::max(wa, wb);
    const double norm = mx + std::log(std::exp(wa - mx) + std::exp(wb - mx));
    if (zl == a) { logq_rev += wa - norm; ++nA; } else { logq_rev += wb - norm; ++nB; }
  }
  int E_after = 1;
  for (int k = 0; k < K; ++k) if (w.nk[k] == 0) ++E_after;
  const int na = w.nk[a], nb = w.nk[b];
  const double* pa_ = pt.p.data() + (size_t)a * Q;
  const double* pb_ = pt.p.data() + (size_t)b * Q;
  const double lprop = log_prior_p(s, Mg.p.data(), Q) - log_prior_p(s, pa_, Q) - log_prior_p(s, pb_, Q)
    + log_q(qA, pa_) + log_q(qB, pb_) - log_q(qM, Mg.p.data());
  const double lpz = std::lgamma(na + nb + g) + std::lgamma(g) - std::lgamma(na + g) - std::lgamma(nb + g);
  const double la = beta * (ll_new - ll_old) + lpz + logq_rev - std::log((double)E_after) + lprop;
  if (std::log(r.unif_pos()) < la) {
    for (int l = 0; l < N; ++l) if (pt.z[l] == b) pt.z[l] = a;
    for (int q = 0; q < Q; ++q) { pt.p[(size_t)a * Q + q] = Mg.p[q]; pt.p[(size_t)b * Q + q] = Bn.p[q]; }
    pt.nu[a] = Mg.nu; pt.nu[b] = Bn.nu;
    w.nk[a] = na + nb; w.nk[b] = 0;
    build_comp(e, pt, w, a); build_comp(e, pt, w, b);
    st.merge_acc += 1;
  }
}

// ---------------------------------------------------------------------------
// Jain-Neal launch-state split-merge at fixed K (sm_method = "jain_neal").
// A fixed-K, tempered port of the authors' MFM sampler's split-merge move; the
// derivation there carries over with these changes:
//   * target(S, theta) = beta * loglik(S; p, nu) + log prior density of
//     theta = (logit p_1..logit p_Q, log nu) -- normalised, since a split has
//     two prior terms against one;
//   * the launch state G(L | U, i, j) uses tempered weights (n + gamma) f^beta
//     and tempered MH steps; it still depends only on (U, i, j, beta), so it is
//     the same in both directions and cancels;
//   * K is fixed: a split moves side B into an empty label e chosen uniformly
//     among the E empties (rejected if none); the emptied label's old
//     parameters are discarded, and the reverse merge redraws them from the
//     prior, so that density cancels. Hence
//       split:  log a = tA + tB - tC + log p(z')/p(z) + log E - log q(A,B|L) - log g(delta)
//       merge:  log a = tU - tA - tB + log p(z')/p(z) - log E' + log q(A,B|L) + log g(delta)
//     with E' the number of empty labels after the merge;
//   * launch random-walk sds are divided by sqrt(1 + beta n_side); the delta
//     sds come from jn_delta_sds (information at the launch state combined
//     with the prior). Both depend only on (A, B, L, beta).
// ---------------------------------------------------------------------------
inline double log_prior_eta(const Settings& s, double p) {
  return log_dbeta(p, s.pa, s.pb) + std::log(p) + std::log1p(-p);
}
inline double log_dnorm(double x, double mu, double sd);
inline double log_prior_lambda(const Settings& s, double nu, double mu) {
  // Normalised truncated density in lambda = log nu: the truncation constant
  // matters, since a split carries two prior terms against one.
  // With nu_pool: N(lambda; mu, tau^2), normalised.
  if (s.nu_pool) return log_dnorm(std::log(nu), mu, s.pool_tau);
  if (nu < s.nu_min) return -INFINITY;
  return s.nu_shape * std::log(s.nu_rate) - std::lgamma(s.nu_shape)
         + (s.nu_shape - 1.0) * std::log(nu) - s.nu_rate * nu + std::log(nu) - s.log_nu_tail;
}
inline double log_dnorm(double x, double mu, double sd) {
  const double z = (x - mu) / sd;
  return -0.5 * z * z - std::log(sd) - 0.91893853320467274178;   // log sqrt(2 pi)
}
inline bool valid_p(double p) { return p > 0.0 && p < 1.0 && std::isfinite(p); }
inline bool valid_nu(double nu) { return nu > 1e-8 && std::isfinite(nu); }

// log-likelihood of unit i at (p, nu) with precomputed a, b, lgamma terms
inline double unit_ll_params(const Ext& e, int i, const double* a, const double* b,
                             const double* lga, const double* lgb, double nu, double lgnu) {
  const bbcore::Data& d = e.d;
  double s = 0.0;
  for (int j = d.cell_start[i]; j < d.cell_start[i + 1]; ++j) {
    const int q = d.cell_q[j], x = d.cell_x[j], y = d.cell_y[j];
    s += std::lgamma(x + a[q]) - lga[q] + std::lgamma(y + b[q]) - lgb[q] - (std::lgamma(x + y + nu) - lgnu);
  }
  return s;
}

struct ParamCache {
  std::vector<double> a, b, lga, lgb; double nu = 1, lgnu = 0;
  void set(const std::vector<double>& p, double nu_) {
    const int Q = (int)p.size();
    a.resize(Q); b.resize(Q); lga.resize(Q); lgb.resize(Q);
    nu = nu_; lgnu = std::lgamma(nu);
    for (int q = 0; q < Q; ++q) { a[q] = nu * p[q]; b[q] = nu * (1 - p[q]); lga[q] = std::lgamma(a[q]); lgb[q] = std::lgamma(b[q]); }
  }
};

// Sum of log-likelihoods of a set of units at (p, nu); -Inf if invalid.
double set_ll(const Ext& e, const std::vector<int>& S, const std::vector<double>& p, double nu) {
  if (!valid_nu(nu)) return -INFINITY;
  for (double v : p) if (!valid_p(v)) return -INFINITY;
  ParamCache c; c.set(p, nu);
  double s = 0.0;
  for (int i : S) s += unit_ll_params(e, i, c.a.data(), c.b.data(), c.lga.data(), c.lgb.data(), c.nu, c.lgnu);
  return s;
}

double theta_log_prior(const Settings& s, const std::vector<double>& p, double nu, double mu) {
  if (!valid_nu(nu)) return -INFINITY;
  double t = log_prior_lambda(s, nu, mu);
  for (double v : p) { if (!valid_p(v)) return -INFINITY; t += log_prior_eta(s, v); }
  return t;
}

struct JNSide {
  std::vector<double> p; double nu = 1;
  std::vector<double> ll;                 // per U position, untempered
};

struct JN {
  std::vector<int> U;                     // ascending
  std::vector<char> inA;
  int posA = -1, posB = -1;
  JNSide A, B;
};

void side_eval(const Ext& e, const JN& w, JNSide& sd) {
  ParamCache c; c.set(sd.p, sd.nu);
  sd.ll.resize(w.U.size());
  for (size_t pos = 0; pos < w.U.size(); ++pos)
    sd.ll[pos] = unit_ll_params(e, w.U[pos], c.a.data(), c.b.data(), c.lga.data(), c.lgb.data(), c.nu, c.lgnu);
}

enum class Sweep { LAUNCH, FINAL, FORCED };

// Restricted scan over the non-anchor units of U; FINAL/FORCED return log q.
double jn_sweep(const JN& w0, JN& w, const Settings& s, double beta, Sweep mode,
                const std::vector<char>* target, Rng& r) {
  (void)w0;
  std::vector<int> order;
  for (int pos = 0; pos < (int)w.U.size(); ++pos) if (pos != w.posA && pos != w.posB) order.push_back(pos);
  if (mode == Sweep::LAUNCH)
    for (int h = (int)order.size() - 1; h > 0; --h) std::swap(order[h], order[r.below(h + 1)]);
  std::vector<char> local;
  std::vector<char>& inA = (mode == Sweep::FORCED) ? (local = w.inA) : w.inA;
  int nA = 0; for (char v : inA) nA += v;
  int nB = (int)inA.size() - nA;
  double log_q = 0.0;
  for (int pos : order) {
    if (inA[pos]) nA--; else nB--;
    const double la = std::log(nA + s.gamma) + beta * w.A.ll[pos];
    const double lb = std::log(nB + s.gamma) + beta * w.B.ll[pos];
    const double mx = std::max(la, lb);
    const double ln = mx + std::log(std::exp(la - mx) + std::exp(lb - mx));
    bool toA;
    if (mode == Sweep::FORCED) toA = (*target)[pos];
    else toA = r.unif() < std::exp(la - ln);
    if (mode != Sweep::LAUNCH) log_q += (toA ? la : lb) - ln;
    if (toA) nA++; else nB++;
    inA[pos] = toA;
  }
  return log_q;
}

// Tempered MH refinement of one launch side, using its current members.
void jn_side_update(const Ext& e, JN& w, JNSide& sd, bool sideA, const Settings& s, double beta, double mu, Rng& r) {
  const bbcore::Data& d = e.d;
  const int Q = d.Q;
  std::vector<int> mem;
  for (size_t pos = 0; pos < w.U.size(); ++pos) if ((bool)w.inA[pos] == sideA) mem.push_back(w.U[pos]);
  std::vector<std::vector<int>> cells(Q);
  for (int i : mem) for (int j = d.cell_start[i]; j < d.cell_start[i + 1]; ++j) cells[d.cell_q[j]].push_back(j);
  const double n = (double)mem.size();
  const double sd_p = s.sm_p_rw_sd / std::sqrt(1.0 + beta * n);
  const double sd_nu = s.sm_nu_rw_sd / std::sqrt(1.0 + beta * n);
  auto msum = [&](double nu) {
    const double lg = std::lgamma(nu); double t = 0;
    for (int q = 0; q < Q; ++q) for (int j : cells[q]) t -= std::lgamma(d.cell_x[j] + d.cell_y[j] + nu) - lg;
    return t;
  };
  std::vector<double> S(Q);
  for (int q = 0; q < Q; ++q) S[q] = cells[q].empty() ? 0.0 : bbcore::column_sum_direct(d, cells[q], sd.p[q], sd.nu);
  for (int pass = 0; pass < s.n_p_launch; ++pass) {
    for (int q = 0; q < Q; ++q) {
      const double p1 = inv_logit(logit(sd.p[q]) + sd_p * r.norm());
      const double u = r.unif_pos();
      if (!valid_p(p1)) continue;
      const double S1 = cells[q].empty() ? 0.0 : bbcore::column_sum_direct(d, cells[q], p1, sd.nu);
      if (std::log(u) < beta * (S1 - S[q]) + log_prior_eta(s, p1) - log_prior_eta(s, sd.p[q])) { sd.p[q] = p1; S[q] = S1; }
    }
    const double nu1 = sd.nu * std::exp(sd_nu * r.norm());
    const double u = r.unif_pos();
    if (!valid_nu(nu1)) continue;
    double l0 = msum(sd.nu), l1 = msum(nu1);
    std::vector<double> S1(Q);
    for (int q = 0; q < Q; ++q) {
      S1[q] = cells[q].empty() ? 0.0 : bbcore::column_sum_direct(d, cells[q], sd.p[q], nu1);
      l0 += S[q]; l1 += S1[q];
    }
    if (std::log(u) < beta * (l1 - l0) + log_prior_lambda(s, nu1, mu) - log_prior_lambda(s, sd.nu, mu)) { sd.nu = nu1; S = S1; }
  }
  side_eval(e, w, sd);
}

// Proposal sds for delta = theta_A - theta_B (logit p per question, log nu).
// Per question, the tempered posterior variance of eta_A - eta_B is about
//   (1 / I_A + 1 / I_B) / beta,  I_S = sum over S's cells in q of
//   [m nu_L / (m + nu_L)] p_L (1 - p_L)
// (the Fisher information of logit p for a beta-binomial cell, evaluated at
// the launch profiles), and the prior variance of a difference of two prior
// draws is 2 Var(logit p) = 2 [trigamma(a) + trigamma(b)]. The proposal
// variance combines the two as precisions, times delta_sd^2. For log nu the
// data part is (1/n_A + 1/n_B) / beta and the prior part 2 trigamma(shape).
// Everything depends only on (A, B, L, beta), so g is the same density in the
// split and the reverse merge.
//
// The proposal mean is the launch difference shrunk toward 0 by the data's
// share of the precision, w = prec_data / (prec_data + prec_prior): at small
// beta the launch profiles are near-independent prior draws and their
// difference is noise, while the prior difference has mean 0.
void jn_delta_sds(const Ext& e, const JN& w, const std::vector<int>& A, const std::vector<int>& B,
                  const Settings& s, double beta, std::vector<double>& sdp, double& sdn,
                  std::vector<double>& ctrp, double& ctrn) {
  const bbcore::Data& d = e.d;
  const int Q = d.Q;
  std::vector<double> IA(Q, 0.0), IB(Q, 0.0);
  auto info = [&](const std::vector<int>& S, const JNSide& sd, std::vector<double>& I) {
    for (int i : S)
      for (int j = d.cell_start[i]; j < d.cell_start[i + 1]; ++j) {
        const int q = d.cell_q[j];
        const double m = d.cell_x[j] + d.cell_y[j];
        I[q] += m * sd.nu / (m + sd.nu) * sd.p[q] * (1.0 - sd.p[q]);
      }
  };
  info(A, w.A, IA); info(B, w.B, IB);
  const double prior_eta = 2.0 * s.var_eta_prior;
  sdp.assign(Q, 0.0); ctrp.assign(Q, 0.0);
  for (int q = 0; q < Q; ++q) {
    double pd = 0.0;
    if (IA[q] > 0 && IB[q] > 0 && beta > 0) pd = beta / (1.0 / IA[q] + 1.0 / IB[q]);
    const double prec = 1.0 / prior_eta + pd;
    sdp[q] = s.delta_sd / std::sqrt(prec);
    ctrp[q] = (pd / prec) * (logit(w.A.p[q]) - logit(w.B.p[q]));
  }
  const double prior_lam = 2.0 * s.var_lam_prior;
  const double pdn = beta / (1.0 / A.size() + 1.0 / B.size());
  const double prec_n = 1.0 / prior_lam + pdn;
  sdn = s.delta_nu_sd / std::sqrt(prec_n);
  ctrn = (pdn / prec_n) * (std::log(w.A.nu) - std::log(w.B.nu));
}

// Pooled launch update: each side's profile becomes its members' pooled
// proportions, shrunk by the Beta prior, (sum x + a) / (sum m + a + b); nu
// stays at sm_nu_ref. Deterministic given the members.
void jn_side_pooled(const Ext& e, JN& w, JNSide& sd, bool sideA, const Settings& s) {
  const bbcore::Data& d = e.d;
  const int Q = d.Q;
  std::vector<double> sx(Q, 0.0), sm(Q, 0.0);
  for (size_t pos = 0; pos < w.U.size(); ++pos) {
    if ((bool)w.inA[pos] != sideA) continue;
    const int u = w.U[pos];
    for (int c = d.cell_start[u]; c < d.cell_start[u + 1]; ++c) {
      sx[d.cell_q[c]] += d.cell_x[c]; sm[d.cell_q[c]] += d.cell_x[c] + d.cell_y[c];
    }
  }
  for (int q = 0; q < Q; ++q) sd.p[q] = (sx[q] + s.pa) / (sm[q] + s.pa + s.pb);
  sd.nu = s.sm_nu_ref;
  side_eval(e, w, sd);
}

// L ~ G(L | U, i, j, beta). G may be any procedure that uses only (U, i, j,
// beta) and fresh randomness -- never the current partition or parameters.
//   launch_method 0 ("mh", the MFM port): sides start from prior draws and a
//     random half-split, then n_launch rounds of [restricted sweep; tempered
//     MH refinement of each side].
//   launch_method 1 ("pooled"): sides start from the anchors' own smoothed
//     proportions, (x + a) / (m + a + b), nu = sm_nu_ref; then n_launch rounds
//     of [restricted sweep; pooled-proportion update of each side]. A
//     k-means-like launch: far closer to the data than prior draws.
void jn_launch(const Ext& e, JN& w, int i, int j, const Settings& s, double beta, double mu, Rng& r) {
  const int Q = e.d.Q;
  const bbcore::Data& d = e.d;
  if (s.launch_method == 1) {
    for (int side = 0; side < 2; ++side) {
      JNSide& sd = side == 0 ? w.A : w.B;
      const int u = side == 0 ? i : j;
      sd.p.assign(Q, s.pa / (s.pa + s.pb));
      for (int c = d.cell_start[u]; c < d.cell_start[u + 1]; ++c)
        sd.p[d.cell_q[c]] = (d.cell_x[c] + s.pa) / (d.cell_x[c] + d.cell_y[c] + s.pa + s.pb);
      sd.nu = s.sm_nu_ref;
    }
  } else {
    const Proposal q0 = prior_proposal(s, Q);
    for (JNSide* sd : {&w.A, &w.B}) {
      Fresh f; f.draw(r, s, mu, q0.r, q0.c);
      sd->p = f.p; sd->nu = f.nu;
    }
  }
  w.inA.assign(w.U.size(), 0);
  for (size_t pos = 0; pos < w.U.size(); ++pos) {
    const int u = w.U[pos];
    if (u == i) { w.inA[pos] = 1; w.posA = (int)pos; }
    else if (u == j) { w.inA[pos] = 0; w.posB = (int)pos; }
    else w.inA[pos] = r.unif() < 0.5;
  }
  side_eval(e, w, w.A); side_eval(e, w, w.B);
  for (int l = 0; l < s.n_launch; ++l) {
    jn_sweep(w, w, s, beta, Sweep::LAUNCH, nullptr, r);
    if (s.launch_method == 1) {
      jn_side_pooled(e, w, w.A, true, s);
      jn_side_pooled(e, w, w.B, false, s);
    } else {
      jn_side_update(e, w, w.A, true, s, beta, mu, r);
      jn_side_update(e, w, w.B, false, s, beta, mu, r);
    }
  }
}

void split_merge_jn(const Ext& e, Particle& pt, Work& w, const Settings& s, double beta,
                    Rng& r, MoveStats& st) {
  const int N = e.d.N, Q = e.d.Q, K = s.K;
  const double g = s.gamma;
  if (N < 2) return;
  const int i = r.below(N);
  int j = r.below(N - 1); if (j >= i) ++j;
  const int ki = pt.z[i], kj = pt.z[j];
  JN jn;

  if (ki == kj) {                                         // ------------ split
    st.split_try += 1;
    std::vector<int> empties;
    for (int k = 0; k < K; ++k) if (w.nk[k] == 0) empties.push_back(k);
    if (empties.empty()) { st.split_noempty += 1; return; }
    const int E = (int)empties.size();
    const int el = empties[r.below(E)];
    const int c = ki;
    for (int l = 0; l < N; ++l) if (pt.z[l] == c) jn.U.push_back(l);
    const int n = (int)jn.U.size();
    jn_launch(e, jn, i, j, s, beta, pt.mu, r);
    const double log_q = jn_sweep(jn, jn, s, beta, Sweep::FINAL, nullptr, r);
    std::vector<int> A, B;
    for (int pos = 0; pos < n; ++pos) (jn.inA[pos] ? A : B).push_back(jn.U[pos]);
    const int nA = (int)A.size(), nB = (int)B.size();
    std::vector<double> sdp, ctrp; double sdn, ctrn;
    jn_delta_sds(e, jn, A, B, s, beta, sdp, sdn, ctrp, ctrn);
    std::vector<double> pA(Q), pB(Q);
    double log_g = 0.0;
    const double* pc = pt.p.data() + (size_t)c * Q;
    for (int q = 0; q < Q; ++q) {
      const double ctr = ctrp[q];
      const double dl = ctr + sdp[q] * r.norm();
      log_g += log_dnorm(dl, ctr, sdp[q]);
      const double etaC = logit(pc[q]);
      pA[q] = inv_logit(etaC + ((double)nB / n) * dl);
      pB[q] = inv_logit(etaC - ((double)nA / n) * dl);
    }
    const double dn = ctrn + sdn * r.norm();
    log_g += log_dnorm(dn, ctrn, sdn);
    const double lamC = std::log(pt.nu[c]);
    const double nuA = std::exp(lamC + ((double)nB / n) * dn), nuB = std::exp(lamC - ((double)nA / n) * dn);
    std::vector<double> pcv(pc, pc + Q);
    double llC = 0.0;
    for (int l : jn.U) llC += w.L[l + (size_t)N * c];
    const double tA = beta * set_ll(e, A, pA, nuA) + theta_log_prior(s, pA, nuA, pt.mu);
    const double tB = beta * set_ll(e, B, pB, nuB) + theta_log_prior(s, pB, nuB, pt.mu);
    const double tC = beta * llC + theta_log_prior(s, pcv, pt.nu[c], pt.mu);
    const double lpz = std::lgamma(nA + g) + std::lgamma(nB + g) - std::lgamma(n + g) - std::lgamma(g);
    double la = tA + tB - tC + lpz + std::log((double)E) - log_q - log_g;
    if (s.debug_sm) Rprintf("split beta=%.2e n=%d nA=%d nB=%d tA=%.2f tB=%.2f tC=%.2f lpz=%.2f logE=%.2f logq=%.2f logg=%.2f -> %.2f\n",
                            beta, n, nA, nB, tA, tB, tC, lpz, std::log((double)E), log_q, log_g, la);
    if (std::isnan(la)) la = -INFINITY;
    if (std::log(r.unif_pos()) < la) {
      for (int l : B) pt.z[l] = el;
      for (int q = 0; q < Q; ++q) { pt.p[(size_t)c * Q + q] = pA[q]; pt.p[(size_t)el * Q + q] = pB[q]; }
      pt.nu[c] = nuA; pt.nu[el] = nuB;
      w.nk[c] = nA; w.nk[el] = nB;
      build_comp(e, pt, w, c); build_comp(e, pt, w, el);
      st.split_acc += 1;
    }
    return;
  }

  // -------------------------------------------------------------- merge
  st.merge_try += 1;
  std::vector<int> A, B;
  for (int l = 0; l < N; ++l) {
    if (pt.z[l] == ki) { jn.U.push_back(l); A.push_back(l); }
    else if (pt.z[l] == kj) { jn.U.push_back(l); B.push_back(l); }
  }
  const int nA = (int)A.size(), nB = (int)B.size(), n = nA + nB;
  jn_launch(e, jn, i, j, s, beta, pt.mu, r);
  std::vector<char> target(n);
  for (int pos = 0; pos < n; ++pos) target[pos] = (pt.z[jn.U[pos]] == ki);
  const double log_q = jn_sweep(jn, jn, s, beta, Sweep::FORCED, &target, r);
  std::vector<double> sdp, ctrp; double sdn, ctrn;
  jn_delta_sds(e, jn, A, B, s, beta, sdp, sdn, ctrp, ctrn);
  const double* pa = pt.p.data() + (size_t)ki * Q;
  const double* pb = pt.p.data() + (size_t)kj * Q;
  std::vector<double> pC(Q), pav(pa, pa + Q), pbv(pb, pb + Q);
  double log_g = 0.0;
  for (int q = 0; q < Q; ++q) {
    const double ea = logit(pa[q]), eb = logit(pb[q]);
    pC[q] = inv_logit((nA * ea + nB * eb) / n);
    log_g += log_dnorm(ea - eb, ctrp[q], sdp[q]);
  }
  const double la_ = std::log(pt.nu[ki]), lb_ = std::log(pt.nu[kj]);
  const double nuC = std::exp((nA * la_ + nB * lb_) / n);
  log_g += log_dnorm(la_ - lb_, ctrn, sdn);
  double llA = 0.0, llB = 0.0;
  for (int l : A) llA += w.L[l + (size_t)N * ki];
  for (int l : B) llB += w.L[l + (size_t)N * kj];
  const double tU = beta * set_ll(e, jn.U, pC, nuC) + theta_log_prior(s, pC, nuC, pt.mu);
  const double tA = beta * llA + theta_log_prior(s, pav, pt.nu[ki], pt.mu);
  const double tB = beta * llB + theta_log_prior(s, pbv, pt.nu[kj], pt.mu);
  int E_after = 1;
  for (int k = 0; k < K; ++k) if (w.nk[k] == 0) ++E_after;
  const double lpz = std::lgamma(n + g) + std::lgamma(g) - std::lgamma(nA + g) - std::lgamma(nB + g);
  double la = tU - tA - tB + lpz - std::log((double)E_after) + log_q + log_g;
  if (s.debug_sm) Rprintf("merge beta=%.2e nA=%d nB=%d tU=%.2f tA=%.2f tB=%.2f lpz=%.2f logE=%.2f logq=%.2f logg=%.2f -> %.2f\n",
                          beta, nA, nB, tU, tA, tB, lpz, std::log((double)E_after), log_q, log_g, la);
  if (std::isnan(la)) la = -INFINITY;
  if (std::log(r.unif_pos()) < la) {
    const Proposal q0 = prior_proposal(s, Q);
    Fresh fb; fb.draw(r, s, pt.mu, q0.r, q0.c);     // emptied label: fresh prior draw
    for (int l : B) pt.z[l] = ki;
    for (int q = 0; q < Q; ++q) { pt.p[(size_t)ki * Q + q] = pC[q]; pt.p[(size_t)kj * Q + q] = fb.p[q]; }
    pt.nu[ki] = nuC; pt.nu[kj] = fb.nu;
    w.nk[ki] = n; w.nk[kj] = 0;
    build_comp(e, pt, w, ki); build_comp(e, pt, w, kj);
    st.merge_acc += 1;
  }
}

// ---------------------------------------------------------------------------
// The MCMC kernel at fixed beta. Leaves pt.ll current.
// ---------------------------------------------------------------------------
void kernel(const Ext& e, Particle& pt, Work& w, const Settings& s, double beta,
            double p_sd, double nu_sd, Rng& r, MoveStats& st) {
  const int N = e.d.N, Q = e.d.Q, K = s.K;
  const double g = s.gamma;
  for (int sweep = 0; sweep < s.n_sweeps; ++sweep) {
    // 1. empty components from the prior
    count_nk(pt, w, K);
    for (int k = 0; k < K; ++k) if (w.nk[k] == 0) draw_component_prior(pt, k, s, r, Q);
    build_all(e, pt, w, K);

    // 2. split-merge
    for (int m = 0; m < s.n_sm; ++m) {
      if (s.sm_method == 2) split_merge_jn(e, pt, w, s, beta, r, st);
      else split_merge(e, pt, w, s, beta, r, st);
    }

    // 3. allocations, pi collapsed
    for (int i = 0; i < N; ++i) {
      w.nk[pt.z[i]]--;
      for (int k = 0; k < K; ++k) w.lp[k] = std::log(w.nk[k] + g) + beta * w.L[i + (size_t)N * k];
      const double lse = log_sum_exp(w.lp.data(), K);
      double u = r.unif(), cum = 0.0;
      int pick = K - 1;
      for (int k = 0; k < K; ++k) { cum += std::exp(w.lp[k] - lse); if (u < cum) { pick = k; break; } }
      pt.z[i] = pick;
      w.nk[pick]++;
    }
    build_members(e, pt, w, K);

    // 4. profiles: random walk on logit(p_kq)
    for (int k = 0; k < K; ++k) {
      if (w.nk[k] == 0) continue;
      const double nu = pt.nu[k];
      for (int q = 0; q < Q; ++q) {
        if (w.cells[(size_t)k * Q + q].empty()) continue;  // only the prior: exact draw
        const double p0 = pt.p[(size_t)k * Q + q];
        const double eta1 = logit(p0) + p_sd * r.norm();
        const double p1 = std::min(std::max(inv_logit(eta1), 1e-12), 1.0 - 1e-12);
        const double s0 = col_sum(e, w, k, q, p0, nu), s1 = col_sum(e, w, k, q, p1, nu);
        // log Beta prior + log Jacobian p(1-p) of the logit transform
        const double lp0 = s.pa * std::log(p0) + s.pb * std::log1p(-p0);
        const double lp1 = s.pa * std::log(p1) + s.pb * std::log1p(-p1);
        st.p_try += 1;
        if (std::log(r.unif_pos()) < beta * (s1 - s0) + lp1 - lp0) {
          pt.p[(size_t)k * Q + q] = p1;
          st.p_acc += 1;
        }
      }
      // a component's questions with no member data: exact prior draw
      for (int q = 0; q < Q; ++q)
        if (w.cells[(size_t)k * Q + q].empty())
          pt.p[(size_t)k * Q + q] = std::min(std::max(r.beta(s.pa, s.pb), 1e-12), 1.0 - 1e-12);
    }

    // 5. dispersion: random walk on log(nu_k)
    for (int k = 0; k < K; ++k) {
      if (w.nk[k] == 0) continue;
      const double nu0 = pt.nu[k];
      const double nu1 = nu0 * std::exp(nu_sd * r.norm());
      if (!(nu1 > 1e-8) || !std::isfinite(nu1)) continue;
      if (nu1 < s.nu_min) { st.nu_try += 1; continue; }   // outside the truncated prior: reject
      double l0 = m_part(e, w, k, nu0), l1 = m_part(e, w, k, nu1);
      for (int q = 0; q < Q; ++q) {
        const double pq = pt.p[(size_t)k * Q + q];
        l0 += col_sum(e, w, k, q, pq, nu0);
        l1 += col_sum(e, w, k, q, pq, nu1);
      }
      // Prior density of lambda = log nu: Gamma with the log-scale Jacobian,
      // shape * log(nu) - rate * nu, or N(lambda; mu, tau^2) with nu_pool.
      const double lp0 = s.nu_pool ? log_dnorm(std::log(nu0), pt.mu, s.pool_tau) : s.nu_shape * std::log(nu0) - s.nu_rate * nu0;
      const double lp1 = s.nu_pool ? log_dnorm(std::log(nu1), pt.mu, s.pool_tau) : s.nu_shape * std::log(nu1) - s.nu_rate * nu1;
      st.nu_try += 1;
      if (std::log(r.unif_pos()) < beta * (l1 - l0) + lp1 - lp0) { pt.nu[k] = nu1; st.nu_acc += 1; }
    }

    // 6. pooled dispersion: mu | occupied lambdas (empties integrated out),
    //    then the empties' lambdas | mu. Exact blocked Gibbs; the likelihood
    //    does not involve mu.
    if (s.nu_pool) {
      double sl = 0.0; int T = 0;
      for (int k = 0; k < K; ++k) if (w.nk[k] > 0) { sl += std::log(pt.nu[k]); ++T; }
      const double t2 = s.pool_tau * s.pool_tau, s02 = s.pool_s0 * s.pool_s0;
      const double prec = 1.0 / s02 + T / t2;
      const double mean = (s.pool_m0 / s02 + sl / t2) / prec;
      pt.mu = mean + r.norm() / std::sqrt(prec);
      for (int k = 0; k < K; ++k) if (w.nk[k] == 0) pt.nu[k] = draw_nu(r, s, pt.mu);
    }
  }
  build_all(e, pt, w, K);
  pt.ll = complete_ll(e, pt, w);
}

// Exact draw from the prior (beta = 0).
void draw_prior(const Ext& e, Particle& pt, Work& w, const Settings& s, Rng& r) {
  const int N = e.d.N, Q = e.d.Q, K = s.K;
  pt.z.assign(N, 0); pt.p.assign((size_t)K * Q, 0.5); pt.nu.assign(K, 1.0);
  std::vector<double> pi(K);
  double tot = 0.0;
  for (int k = 0; k < K; ++k) { pi[k] = r.gamma(s.gamma); tot += pi[k]; }
  for (int k = 0; k < K; ++k) pi[k] /= tot;
  for (int i = 0; i < N; ++i) {
    double u = r.unif(), cum = 0.0;
    int pick = K - 1;
    for (int k = 0; k < K; ++k) { cum += pi[k]; if (u < cum) { pick = k; break; } }
    pt.z[i] = pick;
  }
  if (s.nu_pool) pt.mu = s.pool_m0 + s.pool_s0 * r.norm();
  for (int k = 0; k < K; ++k) draw_component_prior(pt, k, s, r, Q);
  build_all(e, pt, w, K);
  pt.ll = complete_ll(e, pt, w);
}

// Pearson correlation of two equal-length vectors (NaN if either is constant).
double corr(const std::vector<double>& a, const std::vector<double>& b) {
  const int n = (int)a.size();
  double ma = 0, mb = 0;
  for (int i = 0; i < n; ++i) { ma += a[i]; mb += b[i]; }
  ma /= n; mb /= n;
  double sab = 0, saa = 0, sbb = 0;
  for (int i = 0; i < n; ++i) {
    sab += (a[i] - ma) * (b[i] - mb); saa += (a[i] - ma) * (a[i] - ma); sbb += (b[i] - mb) * (b[i] - mb);
  }
  return sab / std::sqrt(saa * sbb);
}

// Conditional ESS (as a fraction of P) of incremental weights exp(d * l_j).
double cess(const std::vector<double>& logW, const std::vector<double>& ll, double d) {
  const int P = (int)ll.size();
  double mx = -INFINITY;
  for (int j = 0; j < P; ++j) mx = std::max(mx, d * ll[j]);
  double s1 = 0.0, s2 = 0.0;
  for (int j = 0; j < P; ++j) {
    const double W = std::exp(logW[j]);
    const double v = std::exp(d * ll[j] - mx);
    s1 += W * v; s2 += W * v * v;
  }
  return s1 * s1 / s2;
}

}  // namespace

Settings read_settings(const List& control) {
  Settings s;
  s.K = as<int>(control["K"]);
  s.P = as<int>(control["n_particles"]);
  s.n_sweeps = as<int>(control["n_sweeps"]);
  s.n_sm = as<int>(control["n_sm"]);
  s.n_threads = as<int>(control["n_threads"]);
  s.max_steps = as<int>(control["max_steps"]);
  s.gamma = as<double>(control["gamma"]);
  s.pa = as<double>(control["p_prior_a"]);
  s.pb = as<double>(control["p_prior_b"]);
  s.nu_shape = as<double>(control["nu_shape"]);
  s.nu_rate = as<double>(control["nu_rate"]);
  s.nu_min = as<double>(control["nu_min"]);
  if (!(s.nu_min >= 0.0)) stop("nu_min must be >= 0.");
  s.log_nu_tail = R::pgamma(s.nu_min, s.nu_shape, 1.0 / s.nu_rate, 0, 1);   // log P(nu >= nu_min)
  s.cess_target = as<double>(control["cess_target"]);
  s.resample_threshold = as<double>(control["resample_threshold"]);
  s.p_sd = as<double>(control["p_sd"]);
  s.nu_sd = as<double>(control["nu_sd"]);
  s.adapt = as<bool>(control["adapt"]);
  s.sm_method = as<int>(control["sm_method"]);
  s.sm_informed = (s.sm_method == 1);
  s.n_launch = as<int>(control["n_launch"]);
  s.n_p_launch = as<int>(control["n_p_launch"]);
  s.sm_p_rw_sd = as<double>(control["sm_p_rw_sd"]);
  s.sm_nu_rw_sd = as<double>(control["sm_nu_rw_sd"]);
  s.delta_sd = as<double>(control["delta_sd"]);
  s.delta_nu_sd = as<double>(control["delta_nu_sd"]);
  s.var_eta_prior = R::trigamma(s.pa) + R::trigamma(s.pb);
  if (control.containsElementNamed("debug_sm")) s.debug_sm = as<bool>(control["debug_sm"]);
  s.max_rounds = as<int>(control["max_rounds"]);
  s.zeta = as<double>(control["zeta"]);
  if (!(s.zeta > 0.0 && s.zeta <= 1.0)) stop("zeta must lie in (0, 1].");
  s.launch_method = as<int>(control["launch_method"]);
  s.rho_target = as<double>(control["rho_target"]);
  s.var_lam_prior = R::trigamma(s.nu_shape);
  if (control.containsElementNamed("nu_pool")) s.nu_pool = as<bool>(control["nu_pool"]);
  if (s.nu_pool) {
    s.pool_tau = as<double>(control["pool_tau"]);
    s.pool_m0 = as<double>(control["pool_m0"]);
    s.pool_s0 = as<double>(control["pool_s0"]);
    if (!(s.pool_tau > 0 && s.pool_s0 > 0)) stop("pool_tau and pool_s0 must be > 0.");
    if (s.nu_min > 0) stop("nu_pool replaces the nu floor: set nu_min = 0.");
    s.var_lam_prior = s.pool_tau * s.pool_tau;   // Jain-Neal delta sds: Var(lambda | mu)
  }
  s.sm_conc0 = as<double>(control["sm_conc0"]);
  s.sm_nu_ref = as<double>(control["sm_nu_ref"]);
  s.seed = (uint64_t)as<double>(control["seed"]);
  return s;
}

// [[Rcpp::export]]
List ss_smc_cpp(IntegerMatrix X, IntegerMatrix M, List control) {
  Settings s = read_settings(control);
  const bool progress = as<bool>(control["progress"]);
  const auto t_start = std::chrono::steady_clock::now();   // for the progress lines only
  // Fixed schedule: increasing betas ending at 1 (empty = adaptive).
  const NumericVector beta_grid = as<NumericVector>(control["beta_grid"]);
  const bool fixed_grid = beta_grid.size() > 0;

  const int N = X.nrow(), Q = X.ncol(), K = s.K, P = s.P;
  const Ext e = make_ext(X.begin(), M.begin(), N, Q);

  std::vector<Particle> parts(P);
  std::vector<MoveStats> stats(P);
  std::vector<double> ll(P), logW(P, -std::log((double)P));
  bool threw = false;

  // ---- beta = 0: exact prior draws
  bbcore::parallel_for(P, s.n_threads, [&](int b, int en) {
    Work w(e, K);
    for (int j = b; j < en; ++j) {
      Rng r(stream_seed(s.seed, 0, j));
      draw_prior(e, parts[j], w, s, r);
    }
  });
  for (int j = 0; j < P; ++j) ll[j] = parts[j].ll;

  auto wmean_var = [&](double& mean, double& var) {
    double m = 0.0;
    for (int j = 0; j < P; ++j) m += std::exp(logW[j]) * ll[j];
    double v = 0.0;
    for (int j = 0; j < P; ++j) v += std::exp(logW[j]) * (ll[j] - m) * (ll[j] - m);
    mean = m; var = v;
  };

  std::vector<double> tr_beta, tr_logz_inc, tr_ess, tr_mean, tr_var, tr_pacc, tr_nuacc,
    tr_split, tr_merge, tr_noempty, tr_psd, tr_nusd, tr_resampled, tr_rounds, tr_rho;
  double mean0, var0;
  wmean_var(mean0, var0);
  tr_beta.push_back(0.0); tr_logz_inc.push_back(0.0); tr_ess.push_back(P);
  tr_mean.push_back(mean0); tr_var.push_back(var0);
  tr_pacc.push_back(NA_REAL); tr_nuacc.push_back(NA_REAL); tr_split.push_back(NA_REAL);
  tr_merge.push_back(NA_REAL); tr_noempty.push_back(NA_REAL);
  tr_psd.push_back(s.p_sd); tr_nusd.push_back(s.nu_sd); tr_resampled.push_back(0);
  tr_rounds.push_back(0); tr_rho.push_back(NA_REAL);

  double beta = 0.0, logZ = 0.0, p_sd = s.p_sd, nu_sd = s.nu_sd;
  Rng rs(stream_seed(s.seed, 0xFFFFFFFFULL, 0));   // resampling stream
  int step = 0;

  while (beta < s.zeta) {
    if (++step > s.max_steps) { threw = true; break; }
    // ---- next beta: from the fixed grid, or by bisection on the conditional ESS
    double d_hi = s.zeta - beta, d;
    if (fixed_grid) d = (step - 1 < beta_grid.size() - 1) ? beta_grid[step - 1] - beta : d_hi;
    else if (cess(logW, ll, d_hi) >= s.cess_target) d = d_hi;
    else {
      double lo = 0.0, hi = d_hi;
      for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (cess(logW, ll, mid) >= s.cess_target) lo = mid; else hi = mid;
      }
      d = lo > 0.0 ? lo : hi;
    }
    const double beta_new = (d == d_hi) ? s.zeta : beta + d;

    // ---- stepping-stone increment and reweighting
    std::vector<double> lw(P);
    for (int j = 0; j < P; ++j) lw[j] = logW[j] + d * ll[j];
    const double lse = log_sum_exp(lw.data(), P);
    logZ += lse;                                  // logW was normalised
    for (int j = 0; j < P; ++j) logW[j] = lw[j] - lse;
    double ess_inv = 0.0;
    for (int j = 0; j < P; ++j) ess_inv += std::exp(2.0 * logW[j]);
    const double ess = 1.0 / ess_inv;

    // ---- systematic resampling
    bool resampled = false;
    if (ess < s.resample_threshold * P) {
      std::vector<Particle> np(P);
      const double u0 = rs.unif() / P;
      double cum = std::exp(logW[0]);
      int src = 0;
      for (int j = 0; j < P; ++j) {
        const double u = u0 + (double)j / P;
        while (u > cum && src < P - 1) { ++src; cum += std::exp(logW[src]); }
        np[j] = parts[src];
      }
      parts.swap(np);
      std::fill(logW.begin(), logW.end(), -std::log((double)P));
      resampled = true;
    }

    // ---- move every particle at beta_new
    beta = beta_new;
    std::fill(stats.begin(), stats.end(), MoveStats());
    // Rounds of n_sweeps each, until the particles' log-likelihoods have
    // decorrelated from their values before the moves (corr < rho_target),
    // or max_rounds is reached (Buchholz, Chopin & Jacob 2021). With
    // max_rounds = 1 this is a fixed n_sweeps per step.
    const std::vector<double> ll_before = ll;
    int rounds = 0;
    double rho = NA_REAL;
    while (rounds < s.max_rounds) {
      const int rnd = rounds;
      bbcore::parallel_for(P, s.n_threads, [&](int b, int en) {
        Work w(e, K);
        for (int j = b; j < en; ++j) {
          Rng r(stream_seed(s.seed, (uint64_t)step * 4096ULL + (uint64_t)rnd, j));
          kernel(e, parts[j], w, s, beta, p_sd, nu_sd, r, stats[j]);
        }
      });
      for (int j = 0; j < P; ++j) ll[j] = parts[j].ll;
      ++rounds;
      rho = corr(ll_before, ll);
      if (!(rho >= s.rho_target)) break;          // also stops when rho is NaN (no spread)
    }
    tr_rounds.push_back(rounds); tr_rho.push_back(rho);

    MoveStats tot;
    for (const auto& x : stats) {
      tot.p_acc += x.p_acc; tot.p_try += x.p_try; tot.nu_acc += x.nu_acc; tot.nu_try += x.nu_try;
      tot.split_acc += x.split_acc; tot.split_try += x.split_try;
      tot.merge_acc += x.merge_acc; tot.merge_try += x.merge_try; tot.split_noempty += x.split_noempty;
    }
    const double pacc = tot.p_try > 0 ? tot.p_acc / tot.p_try : NA_REAL;
    const double nuacc = tot.nu_try > 0 ? tot.nu_acc / tot.nu_try : NA_REAL;

    double mean, var;
    wmean_var(mean, var);
    tr_beta.push_back(beta); tr_logz_inc.push_back(lse); tr_ess.push_back(ess);
    tr_mean.push_back(mean); tr_var.push_back(var);
    tr_pacc.push_back(pacc); tr_nuacc.push_back(nuacc);
    tr_split.push_back(tot.split_try > 0 ? tot.split_acc / tot.split_try : NA_REAL);
    tr_merge.push_back(tot.merge_try > 0 ? tot.merge_acc / tot.merge_try : NA_REAL);
    tr_noempty.push_back(tot.split_try > 0 ? tot.split_noempty / tot.split_try : NA_REAL);
    tr_psd.push_back(p_sd); tr_nusd.push_back(nu_sd); tr_resampled.push_back(resampled);

    // ---- adapt proposal sds for the next step (aim ~0.3 acceptance)
    if (s.adapt) {
      auto scale = [](double acc) {
        if (!std::isfinite(acc)) return 1.0;
        if (acc < 0.10) return 0.5;
        if (acc < 0.20) return 0.8;
        if (acc > 0.50) return 1.6;
        if (acc > 0.40) return 1.2;
        return 1.0;
      };
      p_sd = std::min(std::max(p_sd * scale(pacc), 1e-4), 3.0);
      nu_sd = std::min(std::max(nu_sd * scale(nuacc), 1e-4), 3.0);
    }

    if (progress && (step % 25 == 0 || beta >= s.zeta)) {
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
      Rcpp::Rcout << "step " << step << "  sec " << (long) secs << "  beta " << beta << "  logZ " << logZ
                  << "  ESS " << ess << "  acc p/nu " << pacc << "/" << nuacc << "\n";
      Rcpp::checkUserInterrupt();
    }
  }
  if (threw) stop("SMC exceeded max_steps before reaching beta = zeta.");

  // ---- final particles
  IntegerMatrix Z(P, N);
  NumericVector nu_out(P * K), p_out(P * K * Q), W(P), LL(P), MU(P);
  for (int j = 0; j < P; ++j) {
    MU[j] = s.nu_pool ? parts[j].mu : NA_REAL;
    for (int i = 0; i < N; ++i) Z(j, i) = parts[j].z[i] + 1;
    for (int k = 0; k < K; ++k) nu_out[j + P * k] = parts[j].nu[k];
    for (int k = 0; k < K; ++k)
      for (int q = 0; q < Q; ++q) p_out[j + P * (k + K * q)] = parts[j].p[(size_t)k * Q + q];
    W[j] = std::exp(logW[j]);
    LL[j] = parts[j].ll;
  }
  nu_out.attr("dim") = IntegerVector::create(P, K);
  p_out.attr("dim") = IntegerVector::create(P, K, Q);

  DataFrame steps = DataFrame::create(
    Named("beta") = tr_beta, Named("log_increment") = tr_logz_inc, Named("ess") = tr_ess,
    Named("mean_ll") = tr_mean, Named("var_ll") = tr_var, Named("resampled") = tr_resampled,
    Named("acc_p") = tr_pacc, Named("acc_nu") = tr_nuacc, Named("acc_split") = tr_split,
    Named("acc_merge") = tr_merge, Named("split_no_empty") = tr_noempty,
    Named("p_sd") = tr_psd, Named("nu_sd") = tr_nusd, Named("rounds") = tr_rounds,
    Named("rho") = tr_rho);

  return List::create(Named("logZ") = logZ, Named("steps") = steps, Named("z") = Z,
                      Named("p") = p_out, Named("nu") = nu_out, Named("weights") = W,
                      Named("ll") = LL, Named("mu") = MU);
}

// Complete-data log-likelihood of one state, for tests: z 1-based, p K x Q, nu length K.
// [[Rcpp::export]]
double ss_complete_ll_cpp(IntegerMatrix X, IntegerMatrix M, IntegerVector z, NumericMatrix p,
                          NumericVector nu) {
  const int N = X.nrow(), Q = X.ncol(), K = p.nrow();
  const Ext e = make_ext(X.begin(), M.begin(), N, Q);
  Particle pt;
  pt.z.resize(N); for (int i = 0; i < N; ++i) pt.z[i] = z[i] - 1;
  pt.p.resize((size_t)K * Q);
  for (int k = 0; k < K; ++k) for (int q = 0; q < Q; ++q) pt.p[(size_t)k * Q + q] = p(k, q);
  pt.nu.assign(nu.begin(), nu.end());
  Work w(e, K);
  build_all(e, pt, w, K);
  return complete_ll(e, pt, w);
}

// Run the fixed-beta kernel as an MCMC chain from z_init (1-based), for tests
// of invariance against the exact posterior. Returns the allocations after
// each iteration (n_iter x N, 1-based).
// [[Rcpp::export]]
IntegerMatrix ss_kernel_chain_cpp(IntegerMatrix X, IntegerMatrix M, List control, double beta,
                                  int n_iter, IntegerVector z_init) {
  Settings s = read_settings(control);
  const int N = X.nrow(), Q = X.ncol(), K = s.K;
  const Ext e = make_ext(X.begin(), M.begin(), N, Q);
  Work w(e, K);
  Rng r(stream_seed(s.seed, 0, 0));
  Particle pt;
  draw_prior(e, pt, w, s, r);
  for (int i = 0; i < N; ++i) pt.z[i] = z_init[i] - 1;
  MoveStats st;
  IntegerMatrix out(n_iter, N);
  for (int it = 0; it < n_iter; ++it) {
    kernel(e, pt, w, s, beta, s.p_sd, s.nu_sd, r, st);
    for (int i = 0; i < N; ++i) out(it, i) = pt.z[i] + 1;
  }
  return out;
}

// Final rejuvenation: draw n_final particles from a finished run's weighted
// population (systematic resampling), then move each one n_rounds times with
// the fixed-beta kernel on its own random stream, so that duplicates
// separate. The kernel leaves pi_beta invariant, so the result is still an
// approximate sample from pi_beta; the evidence is not touched.
//   P3: the run's p array (P x K x Q, R layout); NU: P x K; W: weights.
// [[Rcpp::export]]
List ss_rejuvenate_cpp(IntegerMatrix X, IntegerMatrix M, List control, double beta,
                       IntegerMatrix Z, NumericVector P3, NumericMatrix NU, NumericVector W,
                       int n_final, int n_rounds, double p_sd, double nu_sd, NumericVector MU) {
  Settings s = read_settings(control);
  const int N = X.nrow(), Q = X.ncol(), K = s.K, P = Z.nrow();
  const Ext e = make_ext(X.begin(), M.begin(), N, Q);
  // systematic resampling of n_final indices
  Rng rs(stream_seed(s.seed, 0xFEEDULL, 0));
  double tot = 0.0; for (int j = 0; j < P; ++j) tot += W[j];
  std::vector<int> src(n_final);
  {
    const double u0 = rs.unif() / n_final;
    double cum = W[0] / tot; int k = 0;
    for (int j = 0; j < n_final; ++j) {
      const double u = u0 + (double)j / n_final;
      while (u > cum && k < P - 1) { ++k; cum += W[k] / tot; }
      src[j] = k;
    }
  }
  std::vector<Particle> parts(n_final);
  for (int j = 0; j < n_final; ++j) {
    const int a = src[j];
    Particle& pt = parts[j];
    pt.z.resize(N); for (int i = 0; i < N; ++i) pt.z[i] = Z(a, i) - 1;
    pt.p.resize((size_t)K * Q);
    for (int k = 0; k < K; ++k) for (int q = 0; q < Q; ++q) pt.p[(size_t)k * Q + q] = P3[a + P * (k + K * q)];
    pt.nu.resize(K); for (int k = 0; k < K; ++k) pt.nu[k] = NU(a, k);
    if (s.nu_pool) { if (MU.size() != P) stop("nu_pool rejuvenation needs the run's mu."); pt.mu = MU[a]; }
  }
  std::vector<MoveStats> stats(n_final);
  bbcore::parallel_for(n_final, s.n_threads, [&](int b, int en) {
    Work w(e, K);
    for (int j = b; j < en; ++j) {
      Rng r(stream_seed(s.seed, 0xBEEF0000ULL, j));
      for (int rd = 0; rd < n_rounds; ++rd) kernel(e, parts[j], w, s, beta, p_sd, nu_sd, r, stats[j]);
    }
  });
  IntegerMatrix Zo(n_final, N);
  NumericVector nu_out(n_final * K), p_out(n_final * K * Q), LL(n_final), MUo(n_final);
  double pa = 0, pt_ = 0, na = 0, nt = 0, sa = 0, st_ = 0, ma = 0, mt = 0;
  for (int j = 0; j < n_final; ++j) {
    for (int i = 0; i < N; ++i) Zo(j, i) = parts[j].z[i] + 1;
    for (int k = 0; k < K; ++k) nu_out[j + n_final * k] = parts[j].nu[k];
    for (int k = 0; k < K; ++k)
      for (int q = 0; q < Q; ++q) p_out[j + n_final * (k + K * q)] = parts[j].p[(size_t)k * Q + q];
    LL[j] = parts[j].ll;
    MUo[j] = s.nu_pool ? parts[j].mu : NA_REAL;
    pa += stats[j].p_acc; pt_ += stats[j].p_try; na += stats[j].nu_acc; nt += stats[j].nu_try;
    sa += stats[j].split_acc; st_ += stats[j].split_try; ma += stats[j].merge_acc; mt += stats[j].merge_try;
  }
  nu_out.attr("dim") = IntegerVector::create(n_final, K);
  p_out.attr("dim") = IntegerVector::create(n_final, K, Q);
  return List::create(Named("z") = Zo, Named("p") = p_out, Named("nu") = nu_out, Named("ll") = LL, Named("mu") = MUo,
                      Named("weights") = NumericVector(n_final, 1.0 / n_final), Named("source") = IntegerVector(src.begin(), src.end()) + 1,
                      Named("acceptance") = NumericVector::create(Named("p") = pa / std::max(pt_, 1.0),
                        Named("nu") = na / std::max(nt, 1.0), Named("split") = sa / std::max(st_, 1.0),
                        Named("merge") = ma / std::max(mt, 1.0)));
}

// Fixed-K MCMC: the SMC's fixed-beta kernel (split-merge, collapsed z Gibbs,
// logit-p and log-nu random walks) run as C independent chains at beta (the
// likelihood power zeta), each from a given starting state.
//   Z0: C x N starting allocations (1-based); P0: C x K x Q profiles (R layout);
//   NU0: C x K dispersions; MU0: C centres of log nu (nu_pool only, else unused).
// During the first n_burn iterations each chain's random-walk sds are tuned
// every `block` iterations toward ~0.3 acceptance; they are frozen afterwards,
// so the saved draws come from a fixed, pi_beta-invariant kernel. Every chain
// has its own random stream, seeded from (seed, chain), so results do not
// depend on n_threads. Draws are saved every `thin` iterations after burn-in;
// row r = c * n_save + s of each output holds chain c's s-th saved draw.
// [[Rcpp::export]]
List bbmm_mcmc_cpp(IntegerMatrix X, IntegerMatrix M, List control, double beta,
                   int n_iter, int n_burn, int thin, IntegerMatrix Z0, NumericVector P0,
                   NumericMatrix NU0, NumericVector MU0, bool progress) {
  Settings s = read_settings(control);
  const int N = X.nrow(), Q = X.ncol(), K = s.K, C = Z0.nrow();
  if (Z0.ncol() != N || NU0.nrow() != C || NU0.ncol() != K || P0.size() != (R_xlen_t)C * K * Q)
    stop("Starting state has the wrong dimensions.");
  if (thin < 1 || n_iter < thin || n_burn < 0) stop("Need thin >= 1, n_iter >= thin and n_burn >= 0.");
  const Ext e = make_ext(X.begin(), M.begin(), N, Q);
  const int n_save = n_iter / thin, n_total = n_burn + n_iter, block = 50;

  std::vector<Particle> parts(C);
  std::vector<Rng> rngs;
  rngs.reserve(C);
  for (int c = 0; c < C; ++c) {
    Particle& pt = parts[c];
    pt.z.resize(N);
    for (int i = 0; i < N; ++i) {
      if (Z0(c, i) < 1 || Z0(c, i) > K) stop("Starting allocations must lie in 1..K.");
      pt.z[i] = Z0(c, i) - 1;
    }
    pt.p.resize((size_t)K * Q);
    for (int k = 0; k < K; ++k)
      for (int q = 0; q < Q; ++q)
        pt.p[(size_t)k * Q + q] = std::min(std::max(P0[c + C * (k + K * q)], 1e-12), 1.0 - 1e-12);
    pt.nu.resize(K);
    for (int k = 0; k < K; ++k) pt.nu[k] = std::max(NU0(c, k), std::max(s.nu_min, 1e-8));
    if (s.nu_pool) pt.mu = MU0.size() == C ? MU0[c] : s.pool_m0;
    rngs.emplace_back(stream_seed(s.seed, 0xC4A1ULL, c));
  }

  std::vector<double> p_sd(C, s.p_sd), nu_sd(C, s.nu_sd);
  std::vector<MoveStats> blk(C), kept(C);
  std::vector<int> zb((size_t)C * n_save * N);
  std::vector<double> pb((size_t)C * n_save * K * Q), nub((size_t)C * n_save * K), llb((size_t)C * n_save),
    mub((size_t)C * n_save);
  auto scale = [](double acc) {
    if (!std::isfinite(acc)) return 1.0;
    if (acc < 0.10) return 0.5;
    if (acc < 0.20) return 0.8;
    if (acc > 0.50) return 1.6;
    if (acc > 0.40) return 1.2;
    return 1.0;
  };
  const auto t_start = std::chrono::steady_clock::now();
  int next_report = n_total / 10;

  for (int it0 = 0; it0 < n_total; it0 += block) {
    const int it1 = std::min(n_total, it0 + block);
    std::fill(blk.begin(), blk.end(), MoveStats());
    bbcore::parallel_for(C, s.n_threads, [&](int b, int en) {
      Work w(e, K);
      for (int c = b; c < en; ++c) {
        Particle& pt = parts[c];
        for (int it = it0; it < it1; ++it) {
          MoveStats& st = it < n_burn ? blk[c] : kept[c];
          kernel(e, pt, w, s, beta, p_sd[c], nu_sd[c], rngs[c], st);
          if (it < n_burn) continue;
          const int a = it - n_burn + 1;
          if (a % thin != 0) continue;
          const size_t r = (size_t)c * n_save + (a / thin - 1);
          for (int i = 0; i < N; ++i) zb[r * N + i] = pt.z[i] + 1;
          for (int k = 0; k < K; ++k) {
            nub[r * K + k] = pt.nu[k];
            for (int q = 0; q < Q; ++q) pb[(r * K + k) * Q + q] = pt.p[(size_t)k * Q + q];
          }
          llb[r] = pt.ll;
          mub[r] = s.nu_pool ? pt.mu : NA_REAL;
        }
      }
    });
    if (s.adapt && it1 <= n_burn)
      for (int c = 0; c < C; ++c) {
        const double pa = blk[c].p_try > 0 ? blk[c].p_acc / blk[c].p_try : NA_REAL;
        const double na = blk[c].nu_try > 0 ? blk[c].nu_acc / blk[c].nu_try : NA_REAL;
        p_sd[c] = std::min(std::max(p_sd[c] * scale(pa), 1e-4), 3.0);
        nu_sd[c] = std::min(std::max(nu_sd[c] * scale(na), 1e-4), 3.0);
      }
    if (progress && it1 >= next_report) {
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
      Rcpp::Rcout << "iteration " << it1 << " / " << n_total << (it1 <= n_burn ? " (burn-in)" : "")
                  << "  sec " << (long)secs << "\n";
      next_report += std::max(1, n_total / 10);
    }
    Rcpp::checkUserInterrupt();
  }

  const int R = C * n_save;
  IntegerMatrix Z(R, N);
  NumericVector p_out((R_xlen_t)R * K * Q), LL(R), MU(R);
  NumericMatrix NU(R, K);
  IntegerVector chain(R), iter(R);
  for (int r = 0; r < R; ++r) {
    for (int i = 0; i < N; ++i) Z(r, i) = zb[(size_t)r * N + i];
    for (int k = 0; k < K; ++k) {
      NU(r, k) = nub[(size_t)r * K + k];
      for (int q = 0; q < Q; ++q) p_out[r + (R_xlen_t)R * (k + K * q)] = pb[((size_t)r * K + k) * Q + q];
    }
    LL[r] = llb[r]; MU[r] = mub[r];
    chain[r] = r / n_save + 1; iter[r] = n_burn + (r % n_save + 1) * thin;
  }
  p_out.attr("dim") = IntegerVector::create(R, K, Q);
  NumericMatrix acc(C, 4);
  for (int c = 0; c < C; ++c) {
    const MoveStats& st = kept[c];
    acc(c, 0) = st.p_try > 0 ? st.p_acc / st.p_try : NA_REAL;
    acc(c, 1) = st.nu_try > 0 ? st.nu_acc / st.nu_try : NA_REAL;
    acc(c, 2) = st.split_try > 0 ? st.split_acc / st.split_try : NA_REAL;
    acc(c, 3) = st.merge_try > 0 ? st.merge_acc / st.merge_try : NA_REAL;
  }
  colnames(acc) = CharacterVector::create("p", "nu", "split", "merge");
  return List::create(Named("z") = Z, Named("p") = p_out, Named("nu") = NU, Named("ll") = LL, Named("mu") = MU,
                      Named("chain") = chain, Named("iteration") = iter, Named("acceptance") = acc,
                      Named("p_sd") = NumericVector(p_sd.begin(), p_sd.end()),
                      Named("nu_sd") = NumericVector(nu_sd.begin(), nu_sd.end()));
}
