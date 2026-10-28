// bbcore.h -- shared, cached beta-binomial likelihood for the BBMM samplers.
//
// Model (all three samplers):
//   X_iq | z_i = k  ~  BetaBinomial(M_iq, alpha = nu_q * p_kq, beta = nu_q * (1 - p_kq))
// with x = X_iq successes, y = M_iq - X_iq failures, m = M_iq trials. Up to the
// binomial coefficient (which depends only on the data and is always omitted),
//
//   log f(x | p, nu) =  [lgamma(x + a) - lgamma(a)]           a = nu p
//                     + [lgamma(y + b) - lgamma(b)]           b = nu (1 - p)
//                     - [lgamma(m + nu) - lgamma(nu)].
//
// WHY THIS FILE EXISTS
// Every sampler spends almost all of its time evaluating the first two
// brackets for many (unit, question, component) triples. Two facts make that
// cacheable without any approximation:
//
//  1. Within question q only a limited set of distinct counts occurs (at most
//     the number of units, often far fewer on big datasets). For a component
//     with profile p_kq we evaluate lgamma(v + a) - lgamma(a) once per DISTINCT
//     value v, walking the sorted values with the exact recurrence
//     lgamma(z + 1) = lgamma(z) + log(z). A unit's log-likelihood is then a sum
//     of table lookups. ("Tables" below.)
//
//  2. A component's total log-likelihood in column q depends on its members
//     only through how many members show each distinct value. Keeping that
//     histogram up to date ("Histogram" below) makes the cost of a profile
//     update proportional to the number of distinct values, not the number of
//     members -- the right trade on large datasets, where components are big.
//
// The term in m depends only on nu, so it is cached per unit ("UnitConstants").
//
// Everything here is exact up to floating-point rounding (the recurrence is
// restarted from std::lgamma every RESTART steps; tests in test_core.R check
// agreement with direct evaluation to 1e-9). No sampler logic lives here.
//
// Layout conventions
//   Matrices passed in from R are column-major (R's layout); internally the
//   data are stored per unit (row-major "cells") because allocation updates
//   visit one unit at a time.
//   Missing cells (NA, or M = 0) are simply absent: they contribute zero.

#ifndef BBCORE_H
#define BBCORE_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <thread>
#include <vector>

namespace bbcore {

// ---------------------------------------------------------------------------
// Tuning constants (measured on an Apple M4 Pro: lgamma ~5.3 ns, log ~1.8 ns).
// ---------------------------------------------------------------------------

// Largest gap between consecutive distinct values bridged by the product
// recurrence (one multiply per step and a single log) instead of a fresh
// lgamma call.
constexpr int GAP_MAX = 16;

// Recompute exactly with std::lgamma at least this often along a table, so
// rounding error from the recurrence cannot accumulate.
constexpr int RESTART = 256;

// Matrices with fewer cells x components than this are filled on one thread;
// below it, thread start-up costs more than it saves.
constexpr double PARALLEL_MIN_WORK = 2.0e5;

// ---------------------------------------------------------------------------
// shifted_lgamma: out[r] = lgamma(a + v[r]) - lgamma(a) for sorted distinct
// values v[0] < v[1] < ... (v[0] may be 0). Exact recurrence, see header.
// ---------------------------------------------------------------------------
inline void shifted_lgamma(double a, const int* v, int D, double* out) {
  const double lga = std::lgamma(a);
  double cur = 0.0;      // lgamma(a + prev) - lgamma(a)
  int prev = 0;          // value at which `cur` is exact
  int since_restart = 0;
  for (int r = 0; r < D; ++r) {
    const int gap = v[r] - prev;
    if (gap > 0) {
      if (gap <= GAP_MAX && since_restart < RESTART) {
        // lgamma(a+v) = lgamma(a+prev) + log( prod_{j=prev}^{v-1} (a+j) )
        double prod = 1.0;
        for (int j = prev; j < v[r]; ++j) prod *= (a + j);
        cur += std::log(prod);
        since_restart += gap;
      } else {
        cur = std::lgamma(a + v[r]) - lga;
        since_restart = 0;
      }
      prev = v[r];
    }
    out[r] = cur;
  }
}

// ---------------------------------------------------------------------------
// Data: the count matrices reorganised for fast access.
// ---------------------------------------------------------------------------
struct Data {
  int N = 0, Q = 0;

  // Per question q: sorted distinct observed successes, failures and trials.
  std::vector<std::vector<int>> xv, yv, mv;

  // Offsets of question q's success / failure blocks inside a component's
  // likelihood table (see Tables). T = total table length.
  std::vector<int> off_x, off_y;
  int T = 0;

  // Observed cells of unit i are cells [cell_start[i], cell_start[i+1]).
  // For each cell: its question, and the table positions of its success and
  // failure values (absolute positions, so a lookup is table[gx] + table[gy]).
  std::vector<int> cell_start, cell_q, cell_gx, cell_gy, cell_mi;

  // Raw counts per cell (needed for direct evaluation with a fresh profile).
  std::vector<int> cell_x, cell_y;

  // For the nu terms: how many observed cells in question q have trials mv[q][r].
  std::vector<std::vector<int>> m_count;

  // Number of observed cells in question q.
  std::vector<int> n_obs;

  int n_cells() const { return (int)cell_q.size(); }
  int Dx(int q) const { return (int)xv[q].size(); }
  int Dy(int q) const { return (int)yv[q].size(); }
};

// Build Data from column-major N x Q integer matrices. A cell is missing when
// either entry is negative (R's NA_INTEGER is INT_MIN) or M == 0.
inline Data make_data(const int* X, const int* M, int N, int Q) {
  Data d;
  d.N = N; d.Q = Q;
  d.xv.resize(Q); d.yv.resize(Q); d.mv.resize(Q); d.m_count.resize(Q);
  d.n_obs.assign(Q, 0);
  auto missing = [&](int i, int q) {
    const int x = X[i + (size_t)q * N], m = M[i + (size_t)q * N];
    return x < 0 || m <= 0 || x > m;
  };
  // Distinct values per question.
  for (int q = 0; q < Q; ++q) {
    std::vector<int> xs, ys, ms;
    for (int i = 0; i < N; ++i) {
      if (missing(i, q)) continue;
      const int x = X[i + (size_t)q * N], m = M[i + (size_t)q * N];
      xs.push_back(x); ys.push_back(m - x); ms.push_back(m);
    }
    auto uniq = [](std::vector<int>& v) {
      std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end());
    };
    d.n_obs[q] = (int)ms.size();
    std::vector<int> ms_all = ms;
    uniq(xs); uniq(ys); uniq(ms);
    d.xv[q] = xs; d.yv[q] = ys; d.mv[q] = ms;
    d.m_count[q].assign(ms.size(), 0);
    for (int m : ms_all)
      d.m_count[q][std::lower_bound(ms.begin(), ms.end(), m) - ms.begin()]++;
  }
  // Table offsets: [q=0 successes][q=0 failures][q=1 successes]...
  d.off_x.resize(Q); d.off_y.resize(Q);
  int pos = 0;
  for (int q = 0; q < Q; ++q) {
    d.off_x[q] = pos; pos += d.Dx(q);
    d.off_y[q] = pos; pos += d.Dy(q);
  }
  d.T = pos;
  // Cells, unit by unit.
  d.cell_start.assign(N + 1, 0);
  for (int i = 0; i < N; ++i) {
    d.cell_start[i] = d.n_cells();
    for (int q = 0; q < Q; ++q) {
      if (missing(i, q)) continue;
      const int x = X[i + (size_t)q * N], m = M[i + (size_t)q * N], y = m - x;
      const int rx = std::lower_bound(d.xv[q].begin(), d.xv[q].end(), x) - d.xv[q].begin();
      const int ry = std::lower_bound(d.yv[q].begin(), d.yv[q].end(), y) - d.yv[q].begin();
      const int rm = std::lower_bound(d.mv[q].begin(), d.mv[q].end(), m) - d.mv[q].begin();
      d.cell_q.push_back(q);
      d.cell_gx.push_back(d.off_x[q] + rx);
      d.cell_gy.push_back(d.off_y[q] + ry);
      d.cell_mi.push_back(rm);
      d.cell_x.push_back(x);
      d.cell_y.push_back(y);
    }
  }
  d.cell_start[N] = d.n_cells();
  return d;
}

// ---------------------------------------------------------------------------
// Tables: one component's cached likelihood terms.
//   t[off_x[q] + r] = lgamma(xv[q][r] + a_q) - lgamma(a_q)
//   t[off_y[q] + r] = lgamma(yv[q][r] + b_q) - lgamma(b_q)
// ---------------------------------------------------------------------------
struct Tables {
  std::vector<double> t;

  // (Re)build the block for question q from profile value p and dispersion nu.
  void build_column(const Data& d, int q, double p, double nu) {
    if ((int)t.size() != d.T) t.assign(d.T, 0.0);
    shifted_lgamma(nu * p, d.xv[q].data(), d.Dx(q), t.data() + d.off_x[q]);
    shifted_lgamma(nu * (1.0 - p), d.yv[q].data(), d.Dy(q), t.data() + d.off_y[q]);
  }

  // Build every question from a profile vector p[0..Q-1] and nu[0..Q-1].
  void build(const Data& d, const double* p, const double* nu) {
    t.assign(d.T, 0.0);
    for (int q = 0; q < d.Q; ++q) build_column(d, q, p[q], nu[q]);
  }
};

// ---------------------------------------------------------------------------
// UnitConstants: the nu-only term of each unit's log-likelihood,
//   c[i] = - sum over observed q of [lgamma(m_iq + nu_q) - lgamma(nu_q)],
// computed once per value of nu from per-question tables over distinct m.
// ---------------------------------------------------------------------------
struct UnitConstants {
  std::vector<double> c;

  void build(const Data& d, const double* nu) {
    std::vector<std::vector<double>> lm(d.Q);
    for (int q = 0; q < d.Q; ++q) {
      lm[q].resize(d.mv[q].size());
      shifted_lgamma(nu[q], d.mv[q].data(), (int)d.mv[q].size(), lm[q].data());
    }
    c.assign(d.N, 0.0);
    for (int i = 0; i < d.N; ++i) {
      double s = 0.0;
      for (int j = d.cell_start[i]; j < d.cell_start[i + 1]; ++j)
        s -= lm[d.cell_q[j]][d.cell_mi[j]];
      c[i] = s;
    }
  }
};

// Sum over all units of the nu-only term, from the distinct-m histogram
// (cost: number of distinct trial counts, independent of N).
inline double total_m_term(const Data& d, const double* nu) {
  double s = 0.0;
  std::vector<double> lm;
  for (int q = 0; q < d.Q; ++q) {
    lm.resize(d.mv[q].size());
    shifted_lgamma(nu[q], d.mv[q].data(), (int)d.mv[q].size(), lm.data());
    for (size_t r = 0; r < lm.size(); ++r) s -= d.m_count[q][r] * lm[r];
  }
  return s;
}

// ---------------------------------------------------------------------------
// Unit log-likelihoods.
// ---------------------------------------------------------------------------

// From a component's tables: exact log f(x_i | component), binomial
// coefficient omitted. `c_i` is UnitConstants::c[i].
inline double unit_loglik(const Data& d, const Tables& tab, int i, double c_i) {
  const double* t = tab.t.data();
  double s = c_i;
  for (int j = d.cell_start[i]; j < d.cell_start[i + 1]; ++j)
    s += t[d.cell_gx[j]] + t[d.cell_gy[j]];
  return s;
}

// Direct evaluation for a profile that has no tables (e.g. a freshly drawn
// auxiliary or birth component). lga[q], lgb[q] must hold lgamma(a_q),
// lgamma(b_q) for this profile; pass them in so a caller evaluating many units
// against the same fresh profile computes them once.
inline double unit_loglik_direct(const Data& d, int i, const double* a, const double* b,
                                 const double* lga, const double* lgb, double c_i) {
  double s = c_i;
  for (int j = d.cell_start[i]; j < d.cell_start[i + 1]; ++j) {
    const int q = d.cell_q[j];
    s += std::lgamma(d.cell_x[j] + a[q]) - lga[q] + std::lgamma(d.cell_y[j] + b[q]) - lgb[q];
  }
  return s;
}

// Run body(begin, end) over [0, n) in contiguous blocks on up to n_threads
// threads. Each block writes only its own outputs, so no locking is needed.
inline void parallel_for(int n, int n_threads, const std::function<void(int, int)>& body) {
  if (n_threads <= 1 || n < 2) { body(0, n); return; }
  n_threads = std::min(n_threads, n);
  std::vector<std::thread> pool;
  const int chunk = (n + n_threads - 1) / n_threads;
  for (int t = 0; t < n_threads; ++t) {
    const int b = t * chunk, e = std::min(n, b + chunk);
    if (b < e) pool.emplace_back(body, b, e);
  }
  for (auto& th : pool) th.join();
}

// Fill out[i + N*k] = log f(x_i | component k) for all units and the K
// components whose tables are given. Uses threads only when the work is large.
inline void loglik_matrix(const Data& d, const std::vector<const Tables*>& tabs,
                          const UnitConstants& uc, double* out, int n_threads) {
  const int K = (int)tabs.size(), N = d.N;
  const double work = (double)d.n_cells() * K;
  const int nt = work < PARALLEL_MIN_WORK ? 1 : n_threads;
  parallel_for(N, nt, [&](int b, int e) {
    for (int k = 0; k < K; ++k)
      for (int i = b; i < e; ++i)
        out[i + (size_t)N * k] = unit_loglik(d, *tabs[k], i, uc.c[i]);
  });
}

// ---------------------------------------------------------------------------
// Histogram: one component's membership summarised per distinct value.
//   h[g]   = number of member cells whose success (or failure) value sits at
//            table position g
//   n_q[q] = number of member cells observed in question q
// ---------------------------------------------------------------------------
struct Histogram {
  std::vector<int> h;
  std::vector<int> n_q;
  int size = 0;          // number of member units

  void reset(const Data& d) { h.assign(d.T, 0); n_q.assign(d.Q, 0); size = 0; }

  // Add (sign = +1) or remove (sign = -1) unit i.
  void update(const Data& d, int i, int sign) {
    if (h.empty()) reset(d);
    for (int j = d.cell_start[i]; j < d.cell_start[i + 1]; ++j) {
      h[d.cell_gx[j]] += sign;
      h[d.cell_gy[j]] += sign;
      n_q[d.cell_q[j]] += sign;
    }
    size += sign;
  }
};

// ---------------------------------------------------------------------------
// Component column sums: the profile-dependent part of a component's
// log-likelihood in question q,
//   S_kq(p) = sum over members i observed in q of
//             [lgamma(x_iq + a) - lgamma(a)] + [lgamma(y_iq + b) - lgamma(b)].
// Two exact routes; column_sum() picks the cheaper one.
// ---------------------------------------------------------------------------

// Route 1 (large components): rebuild the column's tables at p into `scratch`
// (length >= Dx(q) + Dy(q)) and weight by the histogram. Cost ~ Dx + Dy.
// On return `scratch` holds the column's tables at p, ready to copy into the
// component's Tables if the caller accepts p.
inline double column_sum_hist(const Data& d, const Histogram& H, int q, double p,
                              double nu, double* scratch) {
  const int Dx = d.Dx(q), Dy = d.Dy(q);
  shifted_lgamma(nu * p, d.xv[q].data(), Dx, scratch);
  shifted_lgamma(nu * (1.0 - p), d.yv[q].data(), Dy, scratch + Dx);
  const int* hx = H.h.data() + d.off_x[q];
  const int* hy = H.h.data() + d.off_y[q];
  double s = 0.0;
  for (int r = 0; r < Dx; ++r) if (hx[r]) s += hx[r] * scratch[r];
  for (int r = 0; r < Dy; ++r) if (hy[r]) s += hy[r] * scratch[Dx + r];
  return s;
}

// Route 2 (small components): loop over the member cells of question q.
// `cells_q` lists, for each member observed in q, its cell index.
inline double column_sum_direct(const Data& d, const std::vector<int>& cells_q, double p,
                                double nu) {
  const double a = nu * p, b = nu * (1.0 - p);
  const double lga = std::lgamma(a), lgb = std::lgamma(b);
  double s = 0.0;
  for (int j : cells_q)
    s += std::lgamma(d.cell_x[j] + a) - lga + std::lgamma(d.cell_y[j] + b) - lgb;
  return s;
}

// True when the histogram route is expected to be cheaper for a component
// with n member cells in question q (direct: 2 lgamma per member; histogram:
// about one log per distinct value plus the weighted sum).
inline bool prefer_histogram(const Data& d, int q, int n_member_cells) {
  return (double)(d.Dx(q) + d.Dy(q)) < 2.5 * n_member_cells;
}

}  // namespace bbcore

#endif  // BBCORE_H
