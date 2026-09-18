#' Add log-likelihood
#'
#' Calculates and adds the pointwise posterior log-likelihood to a stochastic
#' frontier model that already carries posterior draws.
#'
#' The log-likelihood is evaluated at the composed error, with the inefficiency
#' term integrated out rather than conditioned on, so that it measures the fit
#' of the model rather than the fit of one particular set of augmented draws.
#' The two specifications give different densities, so there is one method per
#' class.
#'
#' For the half-normal model the composed error \eqn{v - u} is skew normal,
#' \deqn{f(\varepsilon) = \frac{2}{\sigma} \phi(\varepsilon / \sigma)
#'       \Phi(-\lambda \varepsilon / \sigma),}
#' with \eqn{\sigma^2 = \sigma_u^2 + \sigma_v^2} and
#' \eqn{\lambda = \sigma_u / \sigma_v}. For the exponential model it is
#' \deqn{f(\varepsilon) = \lambda \exp(\lambda \varepsilon +
#'       \lambda^2 \sigma_v^2 / 2)
#'       \Phi(-\varepsilon / \sigma_v - \lambda \sigma_v).}
#' In both cases the signs of the arguments are reversed for a cost frontier.
#'
#' Only models without an \code{id} are supported, because the marginal
#' likelihood of a panel unit does not factorise over its observations: the
#' units share one inefficiency term, so integrating it out couples the
#' observations that belong to the same unit.
#'
#' @param object an object of class \code{"sfmodel_exp"} or
#'   \code{"sfmodel_hn"}, usually the result of a call to
#'   \code{\link{add_posterior_coefficients}}.
#' @param ... additional arguments.
#'
#' @return The object in \code{object} with \code{posterior$loglik} added, an
#'   \code{\link[coda]{mcmc}} object with one row per draw and one column per
#'   observation.
#'
#' @seealso \code{\link{add_posterior_coefficients}},
#'   \code{\link{selection_criteria}}
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
#'
#' model <- create_sfmodel_exp(y ~ x1 + x2, data = d,
#'                             iterations = 500, burnin = 200)
#' model <- add_posterior_coefficients(add_priors(model))
#' model <- add_posterior_loglik(model)
#'
#' selection_criteria(model)
#'
#' @family posterior simulation
#' @export
add_posterior_loglik <- function(object, ...) {
  UseMethod("add_posterior_loglik")
}

#' @rdname add_posterior_loglik
#' @export
add_posterior_loglik.sfmodel_exp <- function(object, ...) {
  loglik_over_draws(object)
}

#' @rdname add_posterior_loglik
#' @export
add_posterior_loglik.sfmodel_hn <- function(object, ...) {
  loglik_over_draws(object)
}

#' Evaluate the pointwise log-likelihood at every draw
#'
#' @param object a model object carrying posterior draws.
#'
#' @return The model object with \code{posterior$loglik} attached.
#'
#' @keywords internal
loglik_over_draws <- function(object) {

  check_posterior_blocks(object, c("beta", sf_scalar_blocks(object)))

  if (object$model$panel) {
    stop("The pointwise log-likelihood is not available for panel models, ",
         "because the marginal likelihood of a unit does not factorise over ",
         "its observations.")
  }

  beta <- as.matrix(object$posterior$beta$coeffs)
  sigma_v <- as.numeric(object$posterior$sigma_v$coeffs)
  par_u <- as.numeric(object$posterior[[object$model$par_u_name]]$coeffs)

  ll <- matrix(NA_real_, nrow = nrow(beta), ncol = object$n)
  for (d in seq_len(nrow(beta))) {
    ll[d, ] <- sf_loglik_point(object, beta[d, ], sigma_v[d], par_u[d])
  }

  object$posterior$loglik <- .mcmc_draws(object, ll)
  object
}

#' Pointwise log-likelihood at one parameter value
#'
#' The mathematical kernel behind \code{\link{add_posterior_loglik}} and
#' \code{\link{selection_criteria}}. The two densities are kept side by side
#' here, where they can be compared, rather than split across methods; the
#' dispatch that distinguishes the models lives at the level of the exported
#' functions.
#'
#' @param object a model object.
#' @param beta the frontier coefficients.
#' @param sigma_v the standard deviation of the symmetric error.
#' @param par_u the inefficiency parameter, the rate for the exponential model
#'   and the scale for the half-normal one.
#'
#' @return A numeric vector with one entry per observation.
#'
#' @keywords internal
sf_loglik_point <- function(object, beta, sigma_v, par_u) {

  e <- as.numeric(object$data$y - object$data$X %*% beta)
  s <- if (object$model$type == "production") -1 else 1

  if (object$model$ineff == "exponential") {
    # Written out, this density adds lambda^2 sigma_v^2 / 2 to a log
    # distribution function of about the same size and opposite sign, and the
    # two cancel: at lambda = 10^4 the sum has already lost nine digits, and by
    # 10^9 it returns powers of two. Cancelling them by hand leaves
    #   log(lambda) - e^2 / (2 sigma_v^2) + M(s e / sigma_v - lambda sigma_v),
    # with M the function of log_phi_ratio(), which is the same expression
    # evaluated where nothing large is subtracted from anything large.
    z <- s * e / sigma_v
    log(par_u) - 0.5 * z^2 + log_phi_ratio(z - par_u * sigma_v)
  } else {
    sig <- sqrt(par_u^2 + sigma_v^2)
    lam <- par_u / sigma_v
    log(2) - log(sig) + stats::dnorm(e / sig, log = TRUE) +
      stats::pnorm(s * lam * e / sig, log.p = TRUE)
  }
}
