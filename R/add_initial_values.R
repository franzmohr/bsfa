#' Add initial values to a stochastic frontier model
#'
#' Attaches the starting values of the Gibbs sampler to a model object. Priors
#' must be added first, since both methods draw on them.
#'
#' With \code{method = "ols"} the coefficients start at their least squares
#' estimates and the error variance at the least squares residual variance,
#' which ignores the one-sided term but is close enough to be a sensible
#' starting point. The inefficiency parameter starts at its prior mean and the
#' inefficiency terms themselves at zero.
#'
#' With \code{method = "prior"} every parameter is drawn from its prior instead.
#' This is the more honest choice when checking convergence, because running
#' several chains from prior draws disperses the starting points rather than
#' concentrating them all on the same least squares fit.
#'
#' How wide that dispersion is depends on the prior, and under a vague one it
#' is very wide indeed: with the default shape and rate of 0.01 on the error
#' precision, the median starting value for \eqn{\sigma_v^2} is of the order of
#' \eqn{10^{28}}. Chains do come back from there within a few dozen sweeps, but
#' a draw that overflows to \code{Inf} is not a starting point at all, so such
#' a draw is repeated. An improper prior on the coefficients has no draw either
#' and is reported rather than left to fail in the linear algebra.
#'
#' The function also stores the seed of the posterior simulation as element
#' \code{seed} of \code{object$model}, unless the model has one already. It is
#' drawn from R's random number generator, so \code{set.seed()} before this
#' call makes it reproducible, and \code{\link{add_seed}} replaces it
#' afterwards. A \code{set.seed()} call placed after this one therefore does
#' not change the draws.
#'
#' @param object an object of class \code{"sfmodel_exp"} or
#'   \code{"sfmodel_hn"}, with priors already added.
#' @param method either \code{"ols"} or \code{"prior"}. See details.
#' @param ... unused, for compatibility with the generic.
#'
#' @return The model object with the element \code{initial} attached, and the
#'   element \code{seed} of \code{model} if it did not have one.
#'
#' @seealso \code{\link{add_priors}}
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
#'
#' model <- create_sfmodel_exp(y ~ x1 + x2, data = d,
#'                             iterations = 500, burnin = 200)
#' model <- add_priors(model)
#' model <- add_initial_values(model)
#' model$initial$beta
#'
#' @export
add_initial_values <- function(object, ...) {
  UseMethod("add_initial_values")
}

#' @rdname add_initial_values
#' @export
add_initial_values.sfmodel_exp <- function(object, method = "ols", ...) {

  init <- initial_common(object, method)

  # Prior mean of the gamma prior on the exponential rate.
  init$par_u <- switch(
    method,
    ols = object$priors$shape_u / object$priors$rate_u,
    prior = draw_finite(
      function() stats::rgamma(1, shape = object$priors$shape_u,
                               rate = object$priors$rate_u), "lambda"))

  attach_initial(object, init)
}

#' @rdname add_initial_values
#' @export
add_initial_values.sfmodel_hn <- function(object, method = "ols", ...) {

  init <- initial_common(object, method)

  # The gamma prior is on the precision, so the prior mean of sigma_u^2 is
  # rate / (shape - 1) and the starting value is its square root.
  init$par_u <- switch(
    method,
    ols = sqrt(object$priors$rate_u / (object$priors$shape_u - 1)),
    prior = draw_finite(
      function() 1 / sqrt(stats::rgamma(1, shape = object$priors$shape_u,
                                        rate = object$priors$rate_u)),
      "sigma_u"))

  attach_initial(object, init)
}

#' Starting values shared by both inefficiency distributions
#'
#' @param object a model object with priors attached.
#' @param method either \code{"ols"} or \code{"prior"}.
#'
#' @return A list with \code{beta}, \code{sigma_v2}, \code{u} and
#'   \code{method}.
#'
#' @keywords internal
initial_common <- function(object, method) {

  method <- match.arg(method, c("ols", "prior"))

  if (is.null(object$priors)) {
    stop("Add priors before initial values.")
  }

  y <- object$data$y
  X <- object$data$X

  if (method == "ols") {
    beta <- as.numeric(stats::.lm.fit(X, y)$coefficients)
    resid <- y - X %*% beta
    sigma_v2 <- sum(resid^2) / (object$n - object$k)
  } else {
    # A flat prior on the coefficients has no draw, so a singular precision is
    # reported as the improper prior it is rather than as a LAPACK failure.
    # Factorising the precision also avoids inverting it: with B0i = U'U, the
    # vector U^-1 z has covariance B0i^-1.
    U <- tryCatch(chol(object$priors$B0i), error = function(e) {
      stop("Cannot draw starting values from an improper prior on the ",
           "coefficients: 'coef$v_i' is singular, so the prior has no ",
           "covariance to draw from. Use method = \"ols\", or give the ",
           "coefficients a proper prior.")
    })
    beta <- as.numeric(object$priors$b0 +
                         backsolve(U, stats::rnorm(object$k)))
    sigma_v2 <- draw_finite(
      function() 1 / stats::rgamma(1, shape = object$priors$shape_v,
                                   rate = object$priors$rate_v),
      "sigma_v2")
  }

  list(beta = beta,
       sigma_v2 = sigma_v2,
       u = rep(0, object$data$n_units),
       method = method)
}

#' Attach the shared starting values and a seed
#'
#' @param object a model object with priors attached.
#' @param init the starting values built by \code{initial_common}.
#'
#' @return The model object with \code{initial} and, if it had none,
#'   \code{model$seed} attached.
#'
#' @keywords internal
attach_initial <- function(object, init) {

  object$initial <- init

  # The seed of the posterior simulation, unless the model has one already. It
  # is drawn from R's generator, so set.seed() before this call makes it
  # reproducible, and add_seed() replaces it afterwards.
  if (is.null(object$model$seed)) {
    object$model$seed <- .draw_model_seed()
  }

  drop_stale_posterior(object)
}
