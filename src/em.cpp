// em.cpp -- the two kernels of the tempered EM (bbmm_em): the beta-binomial
// log-likelihood of every unit under every component at a common, fixed nu,
// and the Newton M-step for the profiles p_kq.
//
// Both evaluate special functions directly, once per (unit, component,
// question), so the cost does not grow with the size of the vote counts.
// The log density is that of VGAM::dbetabinom (binomial coefficient
// included). Missing cells must arrive as x = n = 0; they contribute nothing.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// log f(x_iq | p_kq, nu) summed over questions: an N x K matrix.
//   log f = log C(n, x) + lgamma(x + a) - lgamma(a) + lgamma(y + b) - lgamma(b)
//           - lgamma(n + nu) + lgamma(nu),   a = nu p, b = nu (1 - p), y = n - x.
// [[Rcpp::export]]
NumericMatrix bb_loglik_em_cpp(NumericMatrix x, NumericMatrix n, NumericMatrix p, double nu) {
  const int N = x.nrow(), Q = x.ncol(), K = p.nrow();
  NumericMatrix out(N, K);
  const double lg_nu = std::lgamma(nu);
  std::vector<double> base(N, 0.0);              // terms that do not depend on the component
  for (int q = 0; q < Q; ++q)
    for (int i = 0; i < N; ++i) {
      const double nn = n(i, q);
      if (nn <= 0) continue;
      const double xx = x(i, q), yy = nn - xx;
      base[i] += std::lgamma(nn + 1) - std::lgamma(xx + 1) - std::lgamma(yy + 1) - std::lgamma(nn + nu) + lg_nu;
    }
  for (int k = 0; k < K; ++k) {
    double* o = &out(0, k);
    for (int i = 0; i < N; ++i) o[i] = base[i];
    for (int q = 0; q < Q; ++q) {
      const double a = nu * p(k, q), b = nu * (1 - p(k, q)), lga = std::lgamma(a), lgb = std::lgamma(b);
      for (int i = 0; i < N; ++i) {
        const double nn = n(i, q);
        if (nn <= 0) continue;
        const double xx = x(i, q);
        o[i] += std::lgamma(xx + a) - lga + std::lgamma(nn - xx + b) - lgb;
      }
    }
  }
  return out;
}

// M-step for p: for each (k, q), maximise
//   sum_i w_ik log f(x_iq | p, nu) + (prior_a - 1) log p + (prior_b - 1) log(1 - p)
// by safeguarded Newton (bisection fallback inside a shrinking bracket).
// The weights w are the E-step responsibilities times zeta, so the Beta prior
// is not tempered. Units with negligible weight are skipped.
// [[Rcpp::export]]
NumericMatrix bb_newton_em_cpp(NumericMatrix w, NumericMatrix x, NumericMatrix n, double nu,
                               NumericMatrix initial, double prior_a, double prior_b,
                               int maxit, double tol, double eps) {
  const int N = x.nrow(), Q = x.ncol(), K = w.ncol();
  NumericMatrix p = clone(initial);
  std::vector<int> idx;
  idx.reserve(N);
  for (int k = 0; k < K; ++k) {
    double wmax = 0;
    for (int i = 0; i < N; ++i) wmax = std::max(wmax, w(i, k));
    for (int q = 0; q < Q; ++q) {
      idx.clear();
      for (int i = 0; i < N; ++i)
        if (n(i, q) > 0 && w(i, k) > 1e-12 * wmax) idx.push_back(i);
      double z = std::min(std::max(p(k, q), eps), 1 - eps), lo = eps, hi = 1 - eps;
      bool done = false;
      for (int it = 0; it < maxit; ++it) {
        const double a = nu * z, b = nu * (1 - z);
        const double da = R::digamma(a), db = R::digamma(b), ta = R::trigamma(a), tb = R::trigamma(b);
        double g = 0, h = 0;
        for (int i : idx) {
          const double wi = w(i, k), xx = x(i, q), yy = n(i, q) - xx;
          g += wi * (R::digamma(xx + a) - da - R::digamma(yy + b) + db);
          h += wi * (R::trigamma(xx + a) - ta + R::trigamma(yy + b) - tb);
        }
        g = nu * g + (prior_a - 1) / z - (prior_b - 1) / (1 - z);
        h = nu * nu * h - (prior_a - 1) / (z * z) - (prior_b - 1) / ((1 - z) * (1 - z));
        if (g > 0) lo = z; else hi = z;
        const double step = g / h;
        double cand = z - step;
        if (std::isfinite(step) && std::abs(step) < tol * (1 + std::abs(z))) { done = true; break; }
        if (!std::isfinite(cand) || cand <= lo || cand >= hi) cand = (lo + hi) / 2;
        z = cand;
      }
      if (!done && hi - lo > 1e-6) stop("M-step for p did not converge.");
      p(k, q) = z;
    }
  }
  return p;
}
