// Gibbs sampler with data augmentation for composed-error stochastic
// frontier models, following van den Broeck, Koop, Osiewalski and Steel
// (1994) and Koop (2003, ch. 7).
//
// The sampled model is
//
//   y_i = x_i' beta + s * u_{g(i)} + v_i,   v_i ~ N(0, sigma_v^2),
//
// where s = -1 for a production frontier and s = +1 for a cost frontier,
// and g(i) maps observation i to the unit that owns its inefficiency term.
// With one unit per observation this is the cross-sectional model; with one
// unit per panel individual it is the time-invariant model of Pitt and Lee
// (1981), which is the specification used in the reference MATLAB code for
// exercise 14.13 of Koop, Poirier and Tobias (2007).

// Coefficients may be placed under stochastic search variable selection
// (George, Sun and Ni, 2008): a coefficient under selection carries a mixture
// of two normal priors centred on zero, a tight one of standard deviation
// tau0 standing for its absence from the frontier and a loose one of tau1
// standing for its presence, and an inclusion indicator is drawn for it in
// every sweep.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// Draw the inclusion indicators of the coefficients under stochastic search
// variable selection and write the implied prior precisions into the diagonal
// of 'Bi'.
//
// The two mixture components are compared on the log scale. Formed directly,
// the tight component's density underflows to zero as soon as a coefficient is
// a few of its own standard deviations from zero, which for a tau0 of the size
// this prior calls for happens routinely, and the ratio of the two densities
// is then 0/0 rather than the inclusion probability of one.
static void draw_inclusion(arma::vec& inc,
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

// Draw from the standard normal distribution truncated to [alpha, Inf).
// Naive rejection is used for moderate alpha; in the right tail it degenerates,
// so the exponential accept-reject scheme of Robert (1995) takes over.
static double rstdnorm_lb(const double alpha) {
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
static double rtnorm_pos(const double mu, const double sd) {
  return mu + sd * rstdnorm_lb(-mu / sd);
}

//' Gibbs sampler for the stochastic frontier model
//'
//' Workhorse behind \code{\link{add_posterior_coefficients}}. Not intended to
//' be called directly, since it performs no input checking.
//'
//' @param y vector of observations on the dependent variable.
//' @param X matrix of regressors.
//' @param g zero-based integer vector mapping each observation to the unit that
//'   owns its inefficiency term.
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
//' @param a_v shape of the gamma prior on the error precision.
//' @param b_v rate of the gamma prior on the error precision.
//' @param a_u shape of the gamma prior on the inefficiency parameter.
//' @param b_u rate of the gamma prior on the inefficiency parameter. For
//'   \code{ineff = 0} these refer to the precision of the half-normal
//'   distribution, for \code{ineff = 1} to the exponential rate.
//' @param beta_init starting values of the frontier coefficients.
//' @param sigma_v2_init starting value of the error variance.
//' @param par_u_init starting value of the inefficiency parameter, the scale
//'   for \code{ineff = 0} and the rate for \code{ineff = 1}.
//' @param u_init starting values of the inefficiency terms.
//' @param ineff 0 for half-normal, 1 for exponential inefficiency.
//' @param s -1 for a production frontier, 1 for a cost frontier.
//' @param draws number of retained iterations before thinning.
//' @param burnin number of discarded iterations.
//' @param thin thinning interval.
//' @param u_thin 0 not to store the augmented inefficiency draws, otherwise
//'   the interval at which the retained draws are stored.
//' @param verbose how often to report progress; 0 for no reporting.
//'
//' @return A named list of draw matrices.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List gibbs_sf(const arma::vec& y,
                    const arma::mat& X,
                    const arma::uvec& g,
                    const int n_units,
                    const arma::vec& b0,
                    const arma::mat& B0i,
                    const arma::uvec& ssvs_idx,
                    const arma::vec& tau0,
                    const arma::vec& tau1,
                    const arma::vec& prob_prior,
                    const double a_v, const double b_v,
                    const double a_u, const double b_u,
                    const arma::vec& beta_init,
                    const double sigma_v2_init,
                    const double par_u_init,
                    const arma::vec& u_init,
                    const int ineff,
                    const double s,
                    const int draws, const int burnin, const int thin,
                    const int u_thin,
                    const int verbose) {

  const arma::uword n = y.n_elem;
  const arma::uword k = X.n_cols;
  const int n_keep = draws / thin;
  // The augmented terms are one column per unit and per retained draw, which
  // is the largest thing the sampler produces. They carry a thinning interval
  // of their own so that a long chain need not store every one of them.
  const int n_keep_u = (u_thin > 0) ? n_keep / u_thin : 0;

  // Quantities that do not change across iterations. The prior precision does
  // change under variable selection, but only in the diagonal entries of the
  // coefficients it selects on, whose prior mean add_priors() has pinned at
  // zero; B0i * b0 is therefore the same in every sweep either way.
  const arma::mat XtX = X.t() * X;
  const arma::vec B0ib0 = B0i * b0;
  const arma::uword n_sel = ssvs_idx.n_elem;
  arma::mat Bi = B0i;
  arma::vec inc(n_sel, arma::fill::ones);

  arma::vec unit_size(n_units, arma::fill::zeros);
  for (arma::uword i = 0; i < n; i++) {
    unit_size(g(i)) += 1.0;
  }

  // Starting values, supplied by add_initial_values().
  arma::vec beta = beta_init;
  double sigma_v2 = sigma_v2_init;
  double sigma_u2 = par_u_init * par_u_init;  // half-normal scale
  double lambda = par_u_init;                 // exponential rate
  arma::vec u = u_init;

  // Storage.
  arma::mat beta_store(k, n_keep, arma::fill::zeros);
  arma::vec sigma_v_store(n_keep, arma::fill::zeros);
  arma::vec par_u_store(n_keep, arma::fill::zeros);
  arma::mat u_store(n_keep_u > 0 ? n_units : 0, n_keep_u, arma::fill::zeros);
  arma::mat inc_store(n_sel, n_keep, arma::fill::zeros);

  const int n_iter = burnin + draws;
  int store = 0;
  int store_u = 0;

  for (int iter = 0; iter < n_iter; iter++) {

    // --- Inefficiency terms, one per unit, truncated normal ---------------
    const arma::vec e = y - X * beta;
    arma::vec unit_sum(n_units, arma::fill::zeros);
    for (arma::uword i = 0; i < n; i++) {
      unit_sum(g(i)) += e(i);
    }

    for (int j = 0; j < n_units; j++) {
      double prec = unit_size(j) / sigma_v2;
      double mean = s * unit_sum(j) / sigma_v2;
      if (ineff == 0) {           // half-normal prior on u
        prec += 1.0 / sigma_u2;
      } else {                    // exponential prior on u
        mean -= lambda;
      }
      u(j) = rtnorm_pos(mean / prec, std::sqrt(1.0 / prec));
    }

    // --- Frontier coefficients --------------------------------------------
    arma::vec u_long(n);
    for (arma::uword i = 0; i < n; i++) {
      u_long(i) = u(g(i));
    }
    const arma::vec y_tilde = y - s * u_long;

    const arma::mat V = arma::inv_sympd(Bi + XtX / sigma_v2);
    const arma::vec m = V * (B0ib0 + X.t() * y_tilde / sigma_v2);
    arma::vec z(k);
    for (arma::uword j = 0; j < k; j++) {
      z(j) = ::norm_rand();
    }
    beta = m + arma::chol(V, "lower") * z;

    // --- Inclusion indicators ----------------------------------------------
    // Drawn after the coefficients rather than before them, so that the pair
    // stored in a retained sweep is a draw from their joint posterior as that
    // sweep left it. The precision they write is what the next sweep draws
    // the coefficients with.
    if (n_sel > 0) {
      draw_inclusion(inc, Bi, beta, ssvs_idx, tau0, tau1, prob_prior);
    }

    // --- Error variance ----------------------------------------------------
    const arma::vec eps = y_tilde - X * beta;
    sigma_v2 = 1.0 / ::Rf_rgamma(a_v + 0.5 * static_cast<double>(n),
                                 1.0 / (b_v + 0.5 * arma::dot(eps, eps)));

    // --- Inefficiency hyperparameter ---------------------------------------
    if (ineff == 0) {
      sigma_u2 = 1.0 / ::Rf_rgamma(a_u + 0.5 * static_cast<double>(n_units),
                                   1.0 / (b_u + 0.5 * arma::dot(u, u)));
    } else {
      lambda = ::Rf_rgamma(a_u + static_cast<double>(n_units),
                           1.0 / (b_u + arma::accu(u)));
    }

    // --- Store --------------------------------------------------------------
    // The draw kept is the last sweep of each thinning block, so that the
    // retained draws carry the iteration index that .mcmc_draws() attaches
    // to them in R.
    if (iter >= burnin && ((iter - burnin + 1) % thin == 0) && store < n_keep) {
      beta_store.col(store) = beta;
      if (n_sel > 0) {
        inc_store.col(store) = inc;
      }
      sigma_v_store(store) = std::sqrt(sigma_v2);
      par_u_store(store) = (ineff == 0) ? std::sqrt(sigma_u2) : lambda;
      if (n_keep_u > 0 && ((store + 1) % u_thin == 0) && store_u < n_keep_u) {
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

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta_store.t(),
    Rcpp::Named("sigma_v") = sigma_v_store,
    Rcpp::Named("par_u") = par_u_store,
    Rcpp::Named("u") = n_keep_u > 0 ? Rcpp::wrap(u_store.t()) : R_NilValue,
    Rcpp::Named("inclusion") = n_sel > 0 ? Rcpp::wrap(inc_store.t()) :
                                           R_NilValue);
}
