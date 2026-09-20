// Determinants of inefficiency, shared by the two samplers.
//
// A one-sided term's scale may be made a function of covariates,
//
//   sigma_u,j^2 = sigma_u^2 * exp(z_j' gamma)     (half-normal, truncated
//                                                  normal)
//   1 / lambda_j = (1 / lambda) * exp(z_j' gamma) (exponential)
//
// so that gamma = 0 returns the homoskedastic model exactly and the baseline
// parameter keeps the meaning, and the prior, it has there. The determinants
// carry no intercept for that reason: the baseline parameter is the value the
// scale takes where every determinant is zero.
//
// The truncated normal family adds a pre-truncation mean,
//
//   u_j ~ N+(z_j' delta, sigma_u,j^2),
//
// which is the constant of Stevenson (1980) when the only determinant is an
// intercept and the specification of Battese and Coelli (1995) otherwise.
//
// Neither gamma nor delta has a conjugate full conditional. For gamma the
// obstruction is the exponential itself; for delta and, in the truncated
// normal, for the baseline scale, it is the normalising constant
// Phi(z' delta / sigma), which depends on the very parameters being drawn.
// Each is therefore updated by a random walk Metropolis step inside the Gibbs
// sweep. The sweep remains a valid MCMC kernel: every other block is drawn
// from its exact full conditional, and a Metropolis step leaves the target
// invariant on its own.

#ifndef BSFA_DETERMINANTS_H
#define BSFA_DETERMINANTS_H

#include <RcppArmadillo.h>

// Log of the standard normal distribution function, computed on the log scale
// throughout so that the far left tail, where Phi underflows to zero and the
// truncated normal's normaliser would become infinite, is still representable.
inline double bsfa_log_pnorm(const double x) {
  return ::Rf_pnorm5(x, 0.0, 1.0, 1, 1);
}

// A random walk Metropolis block with an adapted proposal.
//
// The proposal shape starts at the prior covariance, which is the only scale
// information available before the chain has moved, and is replaced during
// burn-in by the empirical covariance of the draws so far, scaled by the
// 2.38^2 / q of Roberts, Gelman and Gilks (1997). A scalar step size is
// adapted alongside it by Robbins-Monro towards an acceptance rate of 0.234.
//
// All adaptation stops at the end of burn-in. Adapting on past draws makes
// the chain non-Markovian, and while the diminishing adaptation conditions of
// Roberts and Rosenthal (2007) would permit it to continue, stopping is the
// simpler guarantee that the retained draws come from a kernel with the right
// invariant distribution.
struct bsfa_mh_block {

  arma::vec value;      // current state
  arma::vec prior_mean;
  arma::mat prior_prec;

  arma::mat shape;      // Cholesky factor of the proposal covariance
  double log_step;
  int accepted;
  int proposed;

  // Running moments of the draws, for the adapted proposal covariance.
  arma::vec sum;
  arma::mat sum_outer;
  int count;
  bool took;            // whether the last update accepted
  bool shaped;          // whether the empirical covariance has been adopted

  bsfa_mh_block()
    : log_step(0.0), accepted(0), proposed(0), count(0), took(false),
      shaped(false) {}

  void setup(const arma::vec& init, const arma::vec& mean,
             const arma::mat& prec) {
    value = init;
    prior_mean = mean;
    prior_prec = prec;
    const arma::uword q = value.n_elem;
    arma::mat cov;
    if (!arma::inv_sympd(cov, prec)) {
      cov = arma::eye<arma::mat>(q, q) * 0.01;
    }
    if (!arma::chol(shape, cov, "lower")) {
      shape = arma::eye<arma::mat>(q, q) * 0.1;
    }
    // A starting step that moves a fraction of a prior standard deviation.
    log_step = std::log(2.38 / std::sqrt(static_cast<double>(q))) - 1.0;
    sum.zeros(q);
    sum_outer.zeros(q, q);
  }

  arma::uword size() const { return value.n_elem; }

  double log_prior(const arma::vec& x) const {
    const arma::vec d = x - prior_mean;
    return -0.5 * arma::dot(d, prior_prec * d);
  }

  arma::vec propose() const {
    arma::vec e(value.n_elem);
    for (arma::uword j = 0; j < e.n_elem; j++) {
      e(j) = ::norm_rand();
    }
    return value + std::exp(log_step) * (shape * e);
  }

  // Accept or reject a proposal whose log target difference is 'log_ratio'.
  bool decide(const arma::vec& candidate, const double log_ratio) {
    proposed++;
    took = std::isfinite(log_ratio) && (std::log(::unif_rand()) < log_ratio);
    if (took) {
      value = candidate;
      accepted++;
    }
    return took;
  }

  // Called once per sweep, after the block has been updated.
  void adapt(const int iter, const int burnin) {
    if (iter >= burnin) {
      return;
    }
    const double target = 0.234;
    // Robbins-Monro on the log step, with a decaying gain.
    const double gain = 1.0 / std::sqrt(static_cast<double>(iter + 1));
    log_step += gain * ((took ? 1.0 : 0.0) - target);
    if (log_step > 5.0) log_step = 5.0;
    if (log_step < -20.0) log_step = -20.0;

    // The moments are collected only over the second half of burn-in. The
    // first half is the chain travelling from its starting value to the bulk
    // of the posterior, and a covariance that includes that journey is far
    // wider than the posterior's, which is what drives the acceptance rate
    // down rather than towards its target.
    if (iter < burnin / 2) {
      return;
    }
    sum += value;
    sum_outer += value * value.t();
    count++;

    // Refresh the proposal shape occasionally, once there is enough of a
    // trace to estimate a covariance from. The step size is reset to the
    // scaling of Roberts, Gelman and Gilks only when the empirical
    // covariance is first adopted; resetting it at every refresh would
    // discard the Robbins-Monro correction each time and leave the step
    // permanently untuned.
    if (count >= 100 && (iter + 1) % 100 == 0) {
      const double c = static_cast<double>(count);
      const arma::vec mn = sum / c;
      arma::mat cov = sum_outer / c - mn * mn.t();
      const arma::uword q = value.n_elem;
      cov += arma::eye<arma::mat>(q, q) * 1e-10;
      cov *= 2.38 * 2.38 / static_cast<double>(q);
      arma::mat L;
      if (arma::chol(L, cov, "lower")) {
        shape = L;
        if (!shaped) {
          log_step = 0.0;
          shaped = true;
        }
      }
    }
  }

  double acceptance() const {
    return (proposed > 0) ?
      static_cast<double>(accepted) / static_cast<double>(proposed) : 0.0;
  }
};

// The scale multiplier exp(Z gamma) of every term, or a vector of ones when
// the model has no scale determinants.
inline arma::vec bsfa_scale_factor(const arma::mat& Z, const arma::vec& gamma,
                                   const arma::uword n) {
  if (Z.n_cols == 0) {
    return arma::vec(n, arma::fill::ones);
  }
  return arma::exp(Z * gamma);
}

// The pre-truncation mean Z delta, or a vector of zeros for a family that has
// none.
inline arma::vec bsfa_mean_shift(const arma::mat& Z, const arma::vec& delta,
                                 const arma::uword n) {
  if (Z.n_cols == 0) {
    return arma::vec(n, arma::fill::zeros);
  }
  return Z * delta;
}

// Log kernel of the one-sided terms given the scale and, for the truncated
// normal, the pre-truncation mean. This is the part of the sweep's target
// that the Metropolis blocks move, written once so that the three of them
// cannot drift apart.
//
// ineff: 0 half-normal, 1 exponential, 2 truncated normal.
inline double bsfa_log_oneside(const arma::vec& u, const arma::vec& fac,
                               const arma::vec& m, const double base,
                               const int ineff) {
  const arma::uword n = u.n_elem;
  double out = 0.0;
  if (ineff == 1) {
    // Exponential with rate lambda_j = base * exp(-z_j' gamma).
    for (arma::uword j = 0; j < n; j++) {
      const double rate = base / fac(j);
      out += std::log(rate) - rate * u(j);
    }
    return out;
  }
  // Half-normal (m == 0) and truncated normal share a density; for the
  // half-normal the normaliser is Phi(0) = 1/2 for every observation and
  // cancels in every ratio, but it is cheap and keeps one code path.
  for (arma::uword j = 0; j < n; j++) {
    const double s2 = base * fac(j);
    const double sd = std::sqrt(s2);
    const double d = u(j) - m(j);
    out += -0.5 * std::log(s2) - 0.5 * d * d / s2 - bsfa_log_pnorm(m(j) / sd);
  }
  return out;
}

// One one-sided term, with its family, its baseline scale, its determinants
// and the Metropolis blocks they need.
//
// The two-component sampler has one of these and the four-component sampler
// has two, the persistent term and the transient one, which differ only in
// how many of them there are and in what they are indexed by. Holding the
// whole apparatus in one place is what keeps the two samplers, and the two
// terms within the second of them, from drifting apart.
struct bsfa_oneside {

  int ineff;            // 0 half-normal, 1 exponential, 2 truncated normal
  arma::mat Zs;         // determinants of the scale
  arma::mat Zm;         // determinants of the pre-truncation mean
  bsfa_mh_block mh_gamma;
  bsfa_mh_block mh_delta;
  bsfa_mh_block mh_scale;
  arma::vec fac;        // exp(Zs gamma), one per term
  arma::vec m;          // Zm delta, one per term
  double base;          // sigma_u^2 for 0 and 2, the rate lambda for 1
  double a, b;          // shape and rate of the prior on the baseline
  arma::uword n;        // number of terms

  arma::uword q_s() const { return Zs.n_cols; }
  arma::uword q_m() const { return Zm.n_cols; }

  void setup(const int ineff_, const arma::mat& Zs_, const arma::mat& Zm_,
             const arma::vec& g0, const arma::mat& G0i,
             const arma::vec& d0, const arma::mat& D0i,
             const arma::vec& gamma_init, const arma::vec& delta_init,
             const double base_init, const double a_, const double b_,
             const arma::uword n_) {
    ineff = ineff_;
    Zs = Zs_;
    Zm = Zm_;
    base = base_init;
    a = a_;
    b = b_;
    n = n_;
    if (q_s() > 0) {
      mh_gamma.setup(gamma_init, g0, G0i);
    }
    if (q_m() > 0) {
      mh_delta.setup(delta_init, d0, D0i);
    }
    if (ineff == 2) {
      mh_scale.setup(arma::vec(1).fill(std::log(base)),
                     arma::vec(1, arma::fill::zeros),
                     arma::mat(1, 1, arma::fill::eye));
    }
    fac = bsfa_scale_factor(Zs, mh_gamma.value, n);
    m = bsfa_mean_shift(Zm, mh_delta.value, n);
  }

  // The prior contribution of this term to the truncated normal full
  // conditional of the term itself, given the likelihood's precision and
  // scaled mean so far.
  void contribute(const arma::uword j, double& prec, double& mean) const {
    if (ineff == 1) {
      mean -= base / fac(j);
    } else {
      const double s2 = base * fac(j);
      prec += 1.0 / s2;
      // Zero for the half-normal, whose Zm has no columns.
      mean += m(j) / s2;
    }
  }

  // The baseline scale. Conjugate for the half-normal and the exponential
  // once each term is divided by its own multiplier; a Metropolis step on
  // the log scale for the truncated normal, whose normaliser carries it.
  void draw_base(const arma::vec& u) {
    if (ineff == 0) {
      double ss = 0.0;
      for (arma::uword j = 0; j < n; j++) {
        ss += u(j) * u(j) / fac(j);
      }
      base = 1.0 / ::Rf_rgamma(a + 0.5 * static_cast<double>(n),
                               1.0 / (b + 0.5 * ss));
    } else if (ineff == 1) {
      double su = 0.0;
      for (arma::uword j = 0; j < n; j++) {
        su += u(j) / fac(j);
      }
      base = ::Rf_rgamma(a + static_cast<double>(n), 1.0 / (b + su));
    } else {
      const arma::vec cand = mh_scale.propose();
      const double cand_base = std::exp(cand(0));
      const double cur = std::log(base);
      // The inverse gamma prior, plus the log Jacobian of the move to the
      // log scale, which is the log of the parameter itself.
      const double lp_cur = bsfa_log_oneside(u, fac, m, base, 2) -
        a * cur - b / base;
      const double lp_cand = bsfa_log_oneside(u, fac, m, cand_base, 2) -
        a * cand(0) - b / cand_base;
      if (mh_scale.decide(cand, lp_cand - lp_cur)) {
        base = cand_base;
      }
    }
  }

  void draw_determinants(const arma::vec& u) {
    if (q_s() > 0) {
      const arma::vec cand = mh_gamma.propose();
      const arma::vec fac_cand = arma::exp(Zs * cand);
      const double lr =
        (bsfa_log_oneside(u, fac_cand, m, base, ineff) +
           mh_gamma.log_prior(cand)) -
        (bsfa_log_oneside(u, fac, m, base, ineff) +
           mh_gamma.log_prior(mh_gamma.value));
      if (mh_gamma.decide(cand, lr)) {
        fac = fac_cand;
      }
    }
    if (q_m() > 0) {
      const arma::vec cand = mh_delta.propose();
      const arma::vec m_cand = Zm * cand;
      const double lr =
        (bsfa_log_oneside(u, fac, m_cand, base, ineff) +
           mh_delta.log_prior(cand)) -
        (bsfa_log_oneside(u, fac, m, base, ineff) +
           mh_delta.log_prior(mh_delta.value));
      if (mh_delta.decide(cand, lr)) {
        m = m_cand;
      }
    }
  }

  void adapt(const int iter, const int burnin) {
    if (q_s() > 0) mh_gamma.adapt(iter, burnin);
    if (q_m() > 0) mh_delta.adapt(iter, burnin);
    if (ineff == 2) mh_scale.adapt(iter, burnin);
  }

  // What add_posterior_coefficients() stores as the term's scalar parameter:
  // the standard deviation for the two normal families, the rate for the
  // exponential.
  double reported() const {
    return (ineff == 1) ? base : std::sqrt(base);
  }
};

#endif  // BSFA_DETERMINANTS_H
