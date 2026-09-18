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
#' @param object an object of class \code{"sfmodel_exp"} or
#'   \code{"sfmodel_hn"}, with priors already added.
#' @param method either \code{"ols"} or \code{"prior"}. See details.
#' @param ... unused, for compatibility with the generic.
#'
#' @return The model object with the element \code{initial} attached.
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
  init$par_u <- switch(method,
                       ols = object$priors$shape_u / object$priors$rate_u,
                       prior = stats::rgamma(1, shape = object$priors$shape_u,
                                             rate = object$priors$rate_u))

  object$initial <- init
  object
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
    prior = 1 / sqrt(stats::rgamma(1, shape = object$priors$shape_u,
                                   rate = object$priors$rate_u)))

  object$initial <- init
  object
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
    beta <- as.numeric(object$priors$b0 +
                         t(chol(solve(object$priors$B0i))) %*%
                         stats::rnorm(object$k))
    sigma_v2 <- 1 / stats::rgamma(1, shape = object$priors$shape_v,
                                  rate = object$priors$rate_v)
  }

  list(beta = beta,
       sigma_v2 = sigma_v2,
       u = rep(0, object$data$n_units),
       method = method)
}
