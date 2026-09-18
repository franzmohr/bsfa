#' Posterior efficiency scores
#'
#' Summarises the posterior distribution of the efficiency scores
#' \eqn{r_j = \exp(-u_j)} implied by the augmented draws of an estimated model.
#'
#' Because the sampler treats \eqn{u} as a latent variable, every retained
#' iteration carries a complete set of efficiency scores. The summaries below
#' are therefore posterior summaries in the ordinary sense, and the intervals
#' have their usual interpretation. This is the practical advantage of the
#' Bayesian approach here: the maximum likelihood counterpart has to report the
#' Jondrow et al. (1982) conditional mean, which is a point predictor of an
#' unobserved quantity and carries no comparable measure of uncertainty.
#'
#' @param object an object of class \code{"sfmodel_exp"} or
#'   \code{"sfmodel_hn"}, estimated with \code{keep_u = TRUE}.
#' @param probs quantiles of the posterior to report.
#' @param ... unused, for compatibility with the generic.
#'
#' @return A data frame with one row per unit, giving the posterior mean,
#'   standard deviation and the requested quantiles of the efficiency score.
#'
#' @references
#' Jondrow, J., Lovell, C. A. K., Materov, I. S., & Schmidt, P. (1982). On the
#' estimation of technical inefficiency in the stochastic frontier production
#' function model. \emph{Journal of Econometrics}, 19(2--3), 233--238.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 100, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
#'
#' model <- create_sfmodel_exp(y ~ x1, data = d,
#'                             iterations = 500, burnin = 200)
#' model <- add_posterior_coefficients(add_priors(model))
#'
#' head(efficiency(model))
#'
#' @export
efficiency <- function(object, ...) {
  UseMethod("efficiency")
}

#' @rdname efficiency
#' @export
efficiency.sfmodel <- function(object, probs = c(0.05, 0.5, 0.95), ...) {

  if (is.null(object$posterior)) {
    stop("Object does not contain posterior draws. ",
         "See ?add_posterior_coefficients.")
  }
  if (is.null(object$posterior$u$coeffs)) {
    stop("No inefficiency draws stored. Re-run add_posterior_coefficients() ",
         "with keep_u = TRUE.")
  }

  r <- exp(-as.matrix(object$posterior$u$coeffs))

  efficiency_table(r, object$data$unit_labels, probs)
}
