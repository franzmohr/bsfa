// Gibbs sampler with data augmentation for the four-component stochastic
// frontier model of Kumbhakar, Lien and Hjalmarsson (2014), Colombi et al.
// (2014) and Tsionas and Kumbhakar (2014).
//
// The sampled model is
//
//   y_it = x_it' beta + mu_i + s * eta_i + s * u_it + v_it,
//
// with v_it ~ N(0, sigma_v^2), mu_i ~ N(0, sigma_mu^2) an unrestricted unit
// effect, eta_i >= 0 persistent inefficiency and u_it >= 0 transient
// inefficiency. s = -1 for a production frontier and s = +1 for a cost
// frontier. The two one-sided terms follow the same family, half-normal or
// exponential, with separate parameters.
//
// Splitting the unit-specific part into mu_i and eta_i is the point of the
// model: the two-component panel sampler in gibbs_sf.cpp has only eta_i, so
// every persistent difference between units is booked as inefficiency.

// As in gibbs_sf.cpp, coefficients may be placed under stochastic search
// variable selection (George, Sun and Ni, 2008) or under the Bayesian variable
// selection of Korobilis (2013).

// Each of the two one-sided terms may carry determinants of its scale and,
// under the truncated normal family, of its pre-truncation mean. The
// persistent term is one per unit and the transient one is one per
// observation, so their determinants sit at different levels and are given
// separately. Both are handled by the shared machinery in determinants.h.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "determinants.h"

using namespace Rcpp;

// Draw the inclusion indicators of the coefficients under stochastic search
// variable selection and write the implied prior precisions into the diagonal
// of 'Bi', as in gibbs_sf.cpp. The two mixture components are compared on the
// log scale, where the tight one does not underflow.
static void draw_inclusion4(arma::vec& inc,
                            arma::mat& Bi,
                            const arma::vec& beta,
                            const arma::uvec& idx,
                            const arma::vec& tau0,
                            const arma::vec& tau1,
                            const arma::vec& prob) {
  for (arma::uword j = 0; j < idx.n_elem; j++) {
    const double b = beta(idx(j));
    const double l1 = -std::log(tau1(j)) - 0.5 * b * b / (tau1(j) * tau1(j)) +
      std::log(prob(j));
    const double l0 = -std::log(tau0(j)) - 0.5 * b * b / (tau0(j) * tau0(j)) +
      std::log(1.0 - prob(j));
    const double p = 1.0 / (1.0 + std::exp(l0 - l1));
    const bool in = ::unif_rand() < p;
    inc(j) = in ? 1.0 : 0.0;
    const double tau = in ? tau1(j) : tau0(j);
    Bi(idx(j), idx(j)) = 1.0 / (tau * tau);
  }
}


// Draw the inclusion indicators under the Bayesian variable selection of
// Korobilis (2013), where the regressor itself is switched on and off rather
// than its prior tightened. The frontier is x_i' Lambda beta with Lambda the
// diagonal matrix of indicators, so an excluded regressor leaves the
// likelihood altogether and its coefficient falls back on its prior.
//
// Each indicator is drawn from its full conditional, a Bernoulli whose odds
// are the prior odds times the ratio of the two likelihoods, and the ratio
// only involves the column being switched:
//
//   log L(in) - log L(out) = -( ||x_j b_j||^2 - 2 r' x_j b_j ) / (2 sigma^2),
//
// with r the residual of the frontier that leaves column j out. The fitted
// values are carried through the loop and updated one column at a time, so
// the whole sweep costs one pass over the data per selected coefficient.
static void draw_bvs_inclusion4(arma::vec& inc,
                                arma::vec& keep,
                                arma::vec& fit,
                                const arma::vec& beta,
                                const arma::vec& y_tilde,
                                const arma::mat& X,
                                const arma::uvec& idx,
                                const arma::vec& prob,
                                const double sigma_v2) {
  for (arma::uword j = 0; j < idx.n_elem; j++) {
    const arma::uword p = idx(j);
    const arma::vec xj = X.col(p) * beta(p);
    const arma::vec base = fit - keep(p) * xj;
    const arma::vec r = y_tilde - base;
    const double d = arma::dot(xj, xj) - 2.0 * arma::dot(r, xj);
    const double odds = std::log(prob(j)) - std::log(1.0 - prob(j)) -
      0.5 * d / sigma_v2;
    const double p_in = 1.0 / (1.0 + std::exp(-odds));
    keep(p) = (::unif_rand() < p_in) ? 1.0 : 0.0;
    inc(j) = keep(p);
    fit = base + keep(p) * xj;
  }
}

// Draw from the standard normal distribution truncated to [alpha, Inf), as in
// gibbs_sf.cpp: naive rejection for moderate alpha, and the exponential
// accept-reject scheme of Robert (1995) in the right tail.
static double rstdnorm_lb4(const double alpha) {
  if (alpha < 0.45) {
    double z;
    do {
      z = ::norm_rand();
    } while (z < alpha);
    return z;
  }
  const double a_star = 0.5 * (alpha + std::sqrt(alpha * alpha + 4.0));
  double z, rho;
  do {
    z = alpha + ::exp_rand() / a_star;
    rho = std::exp(-0.5 * (z - a_star) * (z - a_star));
  } while (::unif_rand() > rho);
  return z;
}

// Draw from N(mu, sd^2) truncated to [0, Inf).
static double rtnorm_pos4(const double mu, const double sd) {
  return mu + sd * rstdnorm_lb4(-mu / sd);
}

//' Gibbs sampler for the four-component stochastic frontier model
//'
//' Workhorse behind \code{\link{add_posterior_coefficients}} for objects of
//' class \code{"sfmodel4_exp"} and \code{"sfmodel4_hn"}. Not intended to be
//' called directly, since it performs no input checking.
//'
//' @param y vector of observations on the dependent variable.
//' @param X matrix of regressors.
//' @param g zero-based integer vector mapping each observation to its unit.
//' @param n_units number of distinct units in \code{g}.
//' @param b0 prior mean of the frontier coefficients.
//' @param B0i prior precision matrix of the frontier coefficients.
//' @param ssvs_idx zero-based positions of the coefficients placed under
//'   stochastic search variable selection; empty for no variable selection.
//' @param tau0 prior standard deviations of those coefficients when they are
//'   excluded from the frontier.
//' @param tau1 prior standard deviations of those coefficients when they are
//'   included.
//' @param prob_prior prior inclusion probabilities of those coefficients.
//' @param varsel 0 for no variable selection, 1 for stochastic search variable
//'   selection and 2 for the Bayesian variable selection of Korobilis (2013).
//' @param a_v shape of the gamma prior on the error precision.
//' @param b_v rate of the gamma prior on the error precision.
//' @param a_mu shape of the gamma prior on the precision of the unit effect.
//' @param b_mu rate of the gamma prior on the precision of the unit effect.
//' @param a_eta shape of the gamma prior on the persistent inefficiency
//'   parameter.
//' @param b_eta rate of the gamma prior on the persistent inefficiency
//'   parameter.
//' @param a_u shape of the gamma prior on the transient inefficiency
//'   parameter.
//' @param b_u rate of the gamma prior on the transient inefficiency parameter.
//' @param beta_init starting values of the frontier coefficients.
//' @param sigma_v2_init starting value of the error variance.
//' @param sigma_mu2_init starting value of the variance of the unit effect.
//' @param par_eta_init starting value of the persistent inefficiency
//'   parameter.
//' @param par_u_init starting value of the transient inefficiency parameter.
//' @param mu_init starting values of the unit effects.
//' @param eta_init starting values of the persistent inefficiency terms.
//' @param u_init starting values of the transient inefficiency terms.
//' @param Zs_eta,Zm_eta determinants of the scale and of the pre-truncation
//'   mean of the persistent term, one row per unit; zero columns for none.
//' @param Zs_u,Zm_u the same for the transient term, one row per observation.
//' @param g0_eta,G0i_eta,d0_eta,D0i_eta prior mean and precision of the
//'   persistent term's scale and mean determinant coefficients.
//' @param g0_u,G0i_u,d0_u,D0i_u the same for the transient term.
//' @param gamma_eta_init,delta_eta_init,gamma_u_init,delta_u_init starting
//'   values of those coefficients.
//' @param ineff 0 for half-normal, 1 for exponential, 2 for truncated normal
//'   inefficiency.
//' @param s -1 for a production frontier, 1 for a cost frontier.
//' @param draws number of retained iterations before thinning.
//' @param burnin number of discarded iterations.
//' @param thin thinning interval.
//' @param u_thin 0 not to store the augmented terms, otherwise the interval at
//'   which the retained draws are stored.
//' @param verbose how often to report progress; 0 for no reporting.
//'
//' @return A named list of draw matrices.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List gibbs_sf4(const arma::vec& y,
                     const arma::mat& X,
                     const arma::uvec& g,
                     const int n_units,
                     const arma::vec& b0,
                     const arma::mat& B0i,
                     const arma::uvec& ssvs_idx,
                     const arma::vec& tau0,
                     const arma::vec& tau1,
                     const arma::vec& prob_prior,
                     const int varsel,
                     const double a_v, const double b_v,
                     const double a_mu, const double b_mu,
                     const double a_eta, const double b_eta,
                     const double a_u, const double b_u,
                     const arma::vec& beta_init,
                     const double sigma_v2_init,
                     const double sigma_mu2_init,
                     const double par_eta_init,
                     const double par_u_init,
                     const arma::vec& mu_init,
                     const arma::vec& eta_init,
                     const arma::vec& u_init,
                     const arma::mat& Zs_eta,
                     const arma::mat& Zm_eta,
                     const arma::mat& Zs_u,
                     const arma::mat& Zm_u,
                     const arma::vec& g0_eta,
                     const arma::mat& G0i_eta,
                     const arma::vec& d0_eta,
                     const arma::mat& D0i_eta,
                     const arma::vec& g0_u,
                     const arma::mat& G0i_u,
                     const arma::vec& d0_u,
                     const arma::mat& D0i_u,
                     const arma::vec& gamma_eta_init,
                     const arma::vec& delta_eta_init,
                     const arma::vec& gamma_u_init,
                     const arma::vec& delta_u_init,
                     const int ineff,
                     const double s,
                     const int draws, const int burnin, const int thin,
                     const int u_thin,
                     const int verbose) {

  const arma::uword n = y.n_elem;
  const arma::uword k = X.n_cols;
  const int n_keep = draws / thin;
  // The transient terms alone are one column per observation and per retained
  // draw, so on a panel of any size they dominate everything else the sampler
  // returns. They carry a thinning interval of their own.
  const int n_keep_u = (u_thin > 0) ? n_keep / u_thin : 0;

  // The prior precision changes across sweeps under variable selection, but
  // only in the diagonal entries of the coefficients it selects on, whose
  // prior mean add_priors() has pinned at zero; B0i * b0 is therefore the same
  // in every sweep either way.
  const arma::mat XtX = X.t() * X;
  const arma::vec B0ib0 = B0i * b0;
  const arma::uword n_sel = ssvs_idx.n_elem;
  arma::mat Bi = B0i;
  arma::vec inc(n_sel, arma::fill::ones);
  // As in gibbs_sf.cpp: the Korobilis indicators multiply the coefficients, so
  // the frontier is built from beta_eff. Without that algorithm the multiplier
  // stays at one.
  arma::vec keep(k, arma::fill::ones);
  arma::vec beta_eff = beta_init;

  arma::vec unit_size(n_units, arma::fill::zeros);
  for (arma::uword i = 0; i < n; i++) {
    unit_size(g(i)) += 1.0;
  }

  // Starting values, supplied by add_initial_values().
  arma::vec beta = beta_init;
  double sigma_v2 = sigma_v2_init;
  double sigma_mu2 = sigma_mu2_init;
  arma::vec mu = mu_init;
  arma::vec eta = eta_init;
  arma::vec u = u_init;

  // The two one-sided terms. They differ in how many there are of them and
  // in the level their determinants sit at, and in nothing else.
  const arma::uword n_u = static_cast<arma::uword>(n_units);
  bsfa_oneside term_eta;
  bsfa_oneside term_u;
  term_eta.setup(ineff, Zs_eta, Zm_eta, g0_eta, G0i_eta, d0_eta, D0i_eta,
                 gamma_eta_init, delta_eta_init,
                 (ineff == 1) ? par_eta_init : par_eta_init * par_eta_init,
                 a_eta, b_eta, n_u);
  term_u.setup(ineff, Zs_u, Zm_u, g0_u, G0i_u, d0_u, D0i_u,
               gamma_u_init, delta_u_init,
               (ineff == 1) ? par_u_init : par_u_init * par_u_init,
               a_u, b_u, n);

  arma::mat beta_store(k, n_keep, arma::fill::zeros);
  arma::vec sigma_v_store(n_keep, arma::fill::zeros);
  arma::vec sigma_mu_store(n_keep, arma::fill::zeros);
  arma::vec par_eta_store(n_keep, arma::fill::zeros);
  arma::vec par_u_store(n_keep, arma::fill::zeros);
  arma::mat mu_store(n_keep_u > 0 ? n_units : 0, n_keep_u, arma::fill::zeros);
  arma::mat eta_store(n_keep_u > 0 ? n_units : 0, n_keep_u, arma::fill::zeros);
  arma::mat u_store(n_keep_u > 0 ? n : 0, n_keep_u, arma::fill::zeros);
  arma::mat inc_store(n_sel, n_keep, arma::fill::zeros);
  arma::mat gamma_eta_store(term_eta.q_s(), n_keep, arma::fill::zeros);
  arma::mat delta_eta_store(term_eta.q_m(), n_keep, arma::fill::zeros);
  arma::mat gamma_u_store(term_u.q_s(), n_keep, arma::fill::zeros);
  arma::mat delta_u_store(term_u.q_m(), n_keep, arma::fill::zeros);

  const int n_iter = burnin + draws;
  int store = 0;
  int store_u = 0;

  arma::vec xb(n, arma::fill::zeros);
  arma::vec unit_sum(n_units, arma::fill::zeros);

  for (int iter = 0; iter < n_iter; iter++) {

    xb = X * beta_eff;

    // --- Transient inefficiency, one per observation -----------------------
    for (arma::uword i = 0; i < n; i++) {
      const double r = y(i) - xb(i) - mu(g(i)) - s * eta(g(i));
      double prec = 1.0 / sigma_v2;
      double mean = s * r / sigma_v2;
      term_u.contribute(i, prec, mean);
      u(i) = rtnorm_pos4(mean / prec, std::sqrt(1.0 / prec));
    }

    // --- Persistent inefficiency, one per unit -----------------------------
    unit_sum.zeros();
    for (arma::uword i = 0; i < n; i++) {
      unit_sum(g(i)) += y(i) - xb(i) - mu(g(i)) - s * u(i);
    }
    for (int j = 0; j < n_units; j++) {
      double prec = unit_size(j) / sigma_v2;
      double mean = s * unit_sum(j) / sigma_v2;
      term_eta.contribute(static_cast<arma::uword>(j), prec, mean);
      eta(j) = rtnorm_pos4(mean / prec, std::sqrt(1.0 / prec));
    }

    // --- Unit effects, one per unit, unrestricted --------------------------
    unit_sum.zeros();
    for (arma::uword i = 0; i < n; i++) {
      unit_sum(g(i)) += y(i) - xb(i) - s * eta(g(i)) - s * u(i);
    }
    for (int j = 0; j < n_units; j++) {
      const double prec = unit_size(j) / sigma_v2 + 1.0 / sigma_mu2;
      const double mean = (unit_sum(j) / sigma_v2) / prec;
      mu(j) = mean + std::sqrt(1.0 / prec) * ::norm_rand();
    }

    // --- Frontier coefficients ---------------------------------------------
    arma::vec y_tilde(n);
    for (arma::uword i = 0; i < n; i++) {
      y_tilde(i) = y(i) - mu(g(i)) - s * eta(g(i)) - s * u(i);
    }

    // A coefficient the Korobilis indicators have switched off contributes
    // nothing to the likelihood, so its row and column of the cross-product
    // and its entry of the cross-moment are zeroed and it is drawn from its
    // prior.
    arma::mat XtXg = XtX;
    arma::vec Xty = X.t() * y_tilde;
    if (varsel == 2) {
      XtXg.each_col() %= keep;
      XtXg.each_row() %= keep.t();
      Xty %= keep;
    }

    const arma::mat V = arma::inv_sympd(Bi + XtXg / sigma_v2);
    const arma::vec m = V * (B0ib0 + Xty / sigma_v2);
    arma::vec z(k);
    for (arma::uword j = 0; j < k; j++) {
      z(j) = ::norm_rand();
    }
    beta = m + arma::chol(V, "lower") * z;

    // --- Inclusion indicators ----------------------------------------------
    // Drawn after the coefficients, so that the pair stored in a retained
    // sweep is a draw from their joint posterior as that sweep left it.
    if (varsel == 1) {
      draw_inclusion4(inc, Bi, beta, ssvs_idx, tau0, tau1, prob_prior);
      beta_eff = beta;
    } else if (varsel == 2) {
      arma::vec fit = X * (keep % beta);
      draw_bvs_inclusion4(inc, keep, fit, beta, y_tilde, X, ssvs_idx,
                          prob_prior, sigma_v2);
      beta_eff = keep % beta;
    } else {
      beta_eff = beta;
    }

    // --- Variance parameters -----------------------------------------------
    const arma::vec eps = y_tilde - X * beta_eff;
    sigma_v2 = 1.0 / ::Rf_rgamma(a_v + 0.5 * static_cast<double>(n),
                                 1.0 / (b_v + 0.5 * arma::dot(eps, eps)));

    sigma_mu2 = 1.0 / ::Rf_rgamma(a_mu + 0.5 * static_cast<double>(n_units),
                                  1.0 / (b_mu + 0.5 * arma::dot(mu, mu)));

    term_eta.draw_base(eta);
    term_eta.draw_determinants(eta);
    term_eta.adapt(iter, burnin);

    term_u.draw_base(u);
    term_u.draw_determinants(u);
    term_u.adapt(iter, burnin);

    // --- Store ---------------------------------------------------------------
    // The draw kept is the last sweep of each thinning block, so that the
    // retained draws carry the iteration index that .mcmc_draws() attaches
    // to them in R.
    if (iter >= burnin && ((iter - burnin + 1) % thin == 0) && store < n_keep) {
      // The coefficient kept is the one the frontier was built from, so that
      // a regressor the Korobilis indicators switched off enters the stored
      // draw as the zero it was.
      beta_store.col(store) = beta_eff;
      if (n_sel > 0) {
        inc_store.col(store) = inc;
      }
      sigma_v_store(store) = std::sqrt(sigma_v2);
      sigma_mu_store(store) = std::sqrt(sigma_mu2);
      par_eta_store(store) = term_eta.reported();
      par_u_store(store) = term_u.reported();
      if (term_eta.q_s() > 0) {
        gamma_eta_store.col(store) = term_eta.mh_gamma.value;
      }
      if (term_eta.q_m() > 0) {
        delta_eta_store.col(store) = term_eta.mh_delta.value;
      }
      if (term_u.q_s() > 0) {
        gamma_u_store.col(store) = term_u.mh_gamma.value;
      }
      if (term_u.q_m() > 0) {
        delta_u_store.col(store) = term_u.mh_delta.value;
      }
      if (n_keep_u > 0 && ((store + 1) % u_thin == 0) && store_u < n_keep_u) {
        mu_store.col(store_u) = mu;
        eta_store.col(store_u) = eta;
        u_store.col(store_u) = u;
        store_u++;
      }
      store++;
    }

    if (verbose > 0 && ((iter + 1) % verbose == 0)) {
      Rcpp::Rcout << "Iteration " << (iter + 1) << " of " << n_iter << "\n";
    }
    if ((iter + 1) % 256 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }

  // Acceptance rates of every Metropolis block that was active, named so
  // that a block which barely moves can be found without counting.
  Rcpp::NumericVector acc;
  Rcpp::CharacterVector acc_names;
  if (term_eta.q_s() > 0) {
    acc.push_back(term_eta.mh_gamma.acceptance());
    acc_names.push_back("scale_eta");
  }
  if (term_eta.q_m() > 0) {
    acc.push_back(term_eta.mh_delta.acceptance());
    acc_names.push_back("mean_eta");
  }
  if (term_u.q_s() > 0) {
    acc.push_back(term_u.mh_gamma.acceptance());
    acc_names.push_back("scale_u");
  }
  if (term_u.q_m() > 0) {
    acc.push_back(term_u.mh_delta.acceptance());
    acc_names.push_back("mean_u");
  }
  if (ineff == 2) {
    acc.push_back(term_eta.mh_scale.acceptance());
    acc_names.push_back("par_eta");
    acc.push_back(term_u.mh_scale.acceptance());
    acc_names.push_back("par_u");
  }
  if (acc.size() > 0) {
    acc.names() = acc_names;
  }

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta_store.t(),
    Rcpp::Named("sigma_v") = sigma_v_store,
    Rcpp::Named("sigma_mu") = sigma_mu_store,
    Rcpp::Named("par_eta") = par_eta_store,
    Rcpp::Named("par_u") = par_u_store,
    Rcpp::Named("mu") = n_keep_u > 0 ? Rcpp::wrap(mu_store.t()) : R_NilValue,
    Rcpp::Named("eta") = n_keep_u > 0 ? Rcpp::wrap(eta_store.t()) : R_NilValue,
    Rcpp::Named("u") = n_keep_u > 0 ? Rcpp::wrap(u_store.t()) : R_NilValue,
    Rcpp::Named("inclusion") = n_sel > 0 ? Rcpp::wrap(inc_store.t()) :
                                           R_NilValue,
    Rcpp::Named("scale_eta") = term_eta.q_s() > 0 ?
      Rcpp::wrap(gamma_eta_store.t()) : R_NilValue,
    Rcpp::Named("mean_eta") = term_eta.q_m() > 0 ?
      Rcpp::wrap(delta_eta_store.t()) : R_NilValue,
    Rcpp::Named("scale_u") = term_u.q_s() > 0 ?
      Rcpp::wrap(gamma_u_store.t()) : R_NilValue,
    Rcpp::Named("mean_u") = term_u.q_m() > 0 ?
      Rcpp::wrap(delta_u_store.t()) : R_NilValue,
    Rcpp::Named("acceptance") = acc.size() > 0 ? Rcpp::wrap(acc) : R_NilValue);
}
