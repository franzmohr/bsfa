#' Draw from the posterior of a stochastic frontier model
#'
#' Runs the Gibbs sampler on a model object prepared by
#' \code{\link{create_sfmodel_exp}} or \code{\link{create_sfmodel_hn}} and the
#' \code{add_*} functions.
#'
#' Each sweep draws the inefficiency terms from their truncated normal full
#' conditionals, the frontier coefficients from a normal conditioning on
#' \eqn{y + u} rather than \eqn{y}, and the variance parameters from gamma
#' distributions. This is the data augmentation scheme of van den Broeck, Koop,
#' Osiewalski and Steel (1994).
#'
#' Treating \eqn{u} as a latent variable rather than integrating it out means
#' the efficiency scores come out of the sampler as draws, so
#' \code{\link{efficiency}} reports genuine posterior distributions. Maximum
#' likelihood has to fall back on the Jondrow et al. (1982) conditional mean, a
#' point predictor of an unobserved quantity with no comparable measure of
#' uncertainty.
#'
#' Priors must have been added. Initial values are added with the default
#' method if they are missing, since least squares starting values are a safe
#' choice; the seed is left alone if none was set.
#'
#' @param object an object of class \code{"sfmodel_exp"} or
#'   \code{"sfmodel_hn"}.
#' @param keep_u whether to store the inefficiency draws. Required by
#'   \code{\link{efficiency}}, but note that the stored object grows with the
#'   number of units times the number of retained draws.
#' @param keep_ll whether to store pointwise log-likelihood contributions, for
#'   use with information criteria. Only available for a model without
#'   \code{id}, since the marginal likelihood of a panel unit does not
#'   factorise over its observations.
#' @param verbose either \code{FALSE}, \code{TRUE} for progress at every ten
#'   per cent of iterations, or an integer reporting interval.
#' @param ... unused, for compatibility with the generic.
#'
#' @return An object of class \code{"bsfa_exp"} or \code{"bsfa_hn"}, both
#'   inheriting from \code{"bsfa"}, with elements
#'   \item{draws}{a list with the posterior draws of \code{beta},
#'     \code{sigma_v} and the inefficiency parameter, named \code{lambda} or
#'     \code{sigma_u} according to the model.}
#'   \item{u}{matrix of inefficiency draws, or \code{NULL}.}
#'   \item{log_lik}{matrix of pointwise log-likelihood draws, or \code{NULL}.}
#'   and the specification the draws came from.
#'
#' @references
#' Jondrow, J., Lovell, C. A. K., Materov, I. S., & Schmidt, P. (1982). On the
#' estimation of technical inefficiency in the stochastic frontier production
#' function model. \emph{Journal of Econometrics}, 19(2--3), 233--238.
#'
#' van den Broeck, J., Koop, G., Osiewalski, J., & Steel, M. F. J. (1994).
#' Stochastic frontier models: A Bayesian perspective. \emph{Journal of
#' Econometrics}, 61(2), 273--303.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
#'
#' est <- create_sfmodel_exp(y ~ x1 + x2, data = d,
#'                           iterations = 500, burnin = 200)
#' est <- add_priors(est)
#' est <- add_initial_values(est)
#' est <- draw_posterior(est)
#' summary(est)
#'
#' @export
draw_posterior <- function(object, ...) {
  UseMethod("draw_posterior")
}

#' @rdname draw_posterior
#' @export
draw_posterior.sfmodel_exp <- function(object, keep_u = TRUE, keep_ll = FALSE,
                                       verbose = FALSE, ...) {
  draw_posterior_sf(object, keep_u = keep_u, keep_ll = keep_ll,
                    verbose = verbose, result_class = "bsfa_exp")
}

#' @rdname draw_posterior
#' @export
draw_posterior.sfmodel_hn <- function(object, keep_u = TRUE, keep_ll = FALSE,
                                      verbose = FALSE, ...) {
  draw_posterior_sf(object, keep_u = keep_u, keep_ll = keep_ll,
                    verbose = verbose, result_class = "bsfa_hn")
}

#' Run the sampler for a prepared model object
#'
#' @param object a model object with priors attached.
#' @param keep_u whether to store the inefficiency draws.
#' @param keep_ll whether to store pointwise log-likelihood contributions.
#' @param verbose progress reporting.
#' @param result_class the class to give the result, alongside \code{"bsfa"}.
#'
#' @return An object of class \code{"bsfa"}.
#'
#' @keywords internal
draw_posterior_sf <- function(object, keep_u, keep_ll, verbose,
                              result_class) {

  if (is.null(object$priors)) {
    stop("Add priors before drawing from the posterior. See ?add_priors.")
  }
  if (is.null(object$initial)) {
    object <- add_initial_values(object)
  }

  if (keep_ll && object$model$panel) {
    warning("'keep_ll' is not available for panel models; setting it to ",
            "FALSE.")
    keep_ll <- FALSE
  }

  if (!is.null(object$seed)) {
    set.seed(object$seed)
  }

  n_iter <- object$burnin + object$iterations
  verbose_int <- if (isTRUE(verbose)) max(1L, floor(n_iter / 10)) else
    if (isFALSE(verbose)) 0L else as.integer(verbose)

  out <- gibbs_sf(y = object$data$y,
                  X = object$data$X,
                  g = as.integer(object$data$g) - 1L,
                  n_units = object$data$n_units,
                  b0 = object$priors$b0,
                  B0i = object$priors$B0i,
                  a_v = object$priors$shape_v,
                  b_v = object$priors$rate_v,
                  a_u = object$priors$shape_u,
                  b_u = object$priors$rate_u,
                  beta_init = object$initial$beta,
                  sigma_v2_init = object$initial$sigma_v2,
                  par_u_init = object$initial$par_u,
                  u_init = object$initial$u,
                  ineff = if (object$model$ineff == "halfnormal") 0L else 1L,
                  s = if (object$model$type == "production") -1 else 1,
                  draws = as.integer(object$iterations),
                  burnin = as.integer(object$burnin),
                  thin = as.integer(object$thin),
                  keep_u = keep_u,
                  keep_ll = keep_ll,
                  verbose = verbose_int)

  colnames(out$beta) <- colnames(object$data$X)
  if (!is.null(out$u)) {
    colnames(out$u) <- object$data$unit_labels
  }

  draws <- list(beta = out$beta, sigma_v = as.numeric(out$sigma_v))
  draws[[object$model$par_u_name]] <- as.numeric(out$par_u)

  structure(
    list(draws = draws,
         u = out$u,
         log_lik = out$log_lik,
         units = object$data$unit_labels,
         model = object$model,
         formula = object$formula,
         terms = object$terms,
         priors = object$priors,
         initial = object$initial,
         seed = object$seed,
         n = object$n,
         k = object$k,
         n_units = object$data$n_units,
         mcmc = list(iterations = object$iterations,
                     burnin = object$burnin,
                     thin = object$thin),
         call = object$call),
    class = c(result_class, "bsfa"))
}
