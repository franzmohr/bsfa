#' Add priors to a stochastic frontier model
#'
#' Attaches the prior specification to a model object produced by
#' \code{\link{create_sfmodel_exp}} or \code{\link{create_sfmodel_hn}}.
#'
#' The frontier coefficients receive a normal prior with mean \code{mu} and
#' precision \code{v_i}, and the error precision \eqn{\sigma_v^{-2}} a gamma
#' prior with the given shape and rate.
#'
#' The inefficiency term is hard to elicit directly, so its prior is anchored on
#' \code{r_star}, the efficiency score the researcher considers equally likely
#' to be exceeded as not, a priori. How that anchor maps into a prior differs
#' between the two models, which is why there is one method per class.
#'
#' For the exponential model the mapping is exact. Placing a gamma prior of
#' shape 1 and rate \eqn{c = -\log(r^*)} on the rate parameter \eqn{\lambda}
#' implies the marginal prior density
#' \deqn{p(u) = \int_0^\infty \lambda e^{-\lambda u} c e^{-c\lambda} d\lambda
#'            = c / (u + c)^2,}
#' whose median is exactly \eqn{c}, so that the prior median of \eqn{\exp(-u)}
#' is \eqn{r^*}. This is the elicitation of van den Broeck, Koop, Osiewalski and
#' Steel (1994), and the default shape of 1 corresponds to the two degrees of
#' freedom they use.
#'
#' For the half-normal model no such exact result is available, and the anchor
#' is matched only in expectation: the prior mean of \eqn{\sigma_u^2} is set to
#' the value that would place the median of \eqn{u} at \eqn{-\log(r^*)}, namely
#' \eqn{(-\log(r^*) / \Phi^{-1}(0.75))^2}. The shape must exceed one for that
#' mean to exist, hence the higher default.
#'
#' @param object an object of class \code{"sfmodel_exp"} or
#'   \code{"sfmodel_hn"}.
#' @param coef a named list of prior specifications for the frontier
#'   coefficients, with elements \code{mu}, the prior mean, and \code{v_i}, the
#'   prior precision. Scalars are expanded to a common mean and a diagonal
#'   precision matrix; a vector for \code{mu} and a matrix for \code{v_i} may be
#'   given instead.
#' @param sigma a named list of prior specifications for the error precision,
#'   with elements \code{shape} and \code{rate}.
#' @param lambda a named list of prior specifications for the rate of the
#'   exponential inefficiency distribution, with elements \code{r_star}, the
#'   prior median efficiency, and \code{shape}.
#' @param sigma_u a named list of prior specifications for the scale of the
#'   half-normal inefficiency distribution, with elements \code{r_star} and
#'   \code{shape}.
#' @param ... unused, for compatibility with the generic.
#'
#' @return The model object with the element \code{priors} attached.
#'
#' @references
#' van den Broeck, J., Koop, G., Osiewalski, J., & Steel, M. F. J. (1994).
#' Stochastic frontier models: A Bayesian perspective. \emph{Journal of
#' Econometrics}, 61(2), 273--303.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
#'
#' model <- create_sfmodel_exp(y ~ x1 + x2, data = d,
#'                             iterations = 500, burnin = 200)
#' model <- add_priors(model, lambda = list(r_star = 0.85))
#' model$priors$rate_u
#'
#' @export
add_priors <- function(object, ...) {
  UseMethod("add_priors")
}

#' @rdname add_priors
#' @export
add_priors.sfmodel_exp <- function(object,
                                   coef = list(mu = 0, v_i = 0.01),
                                   sigma = list(shape = 0.01, rate = 0.01),
                                   lambda = list(r_star = 0.75, shape = 1),
                                   ...) {

  lambda <- merge_prior_list(lambda, list(r_star = 0.75, shape = 1), "lambda")
  check_r_star(lambda$r_star)

  object$priors <- c(prior_coef_sigma(object, coef, sigma),
                     list(shape_u = lambda$shape,
                          rate_u = -log(lambda$r_star),
                          r_star = lambda$r_star))
  object
}

#' @rdname add_priors
#' @export
add_priors.sfmodel_hn <- function(object,
                                  coef = list(mu = 0, v_i = 0.01),
                                  sigma = list(shape = 0.01, rate = 0.01),
                                  sigma_u = list(r_star = 0.75, shape = 2.5),
                                  ...) {

  sigma_u <- merge_prior_list(sigma_u, list(r_star = 0.75, shape = 2.5),
                              "sigma_u")
  check_r_star(sigma_u$r_star)
  if (sigma_u$shape <= 1) {
    stop("The shape of the prior on sigma_u must exceed 1, so that the prior ",
         "mean of sigma_u^2 exists.")
  }

  # Match the prior mean of sigma_u^2 to the scale that places the median of a
  # half-normal variate at -log(r_star). See the details section.
  scale_u <- (-log(sigma_u$r_star) / stats::qnorm(0.75))^2

  object$priors <- c(prior_coef_sigma(object, coef, sigma),
                     list(shape_u = sigma_u$shape,
                          rate_u = scale_u * (sigma_u$shape - 1),
                          r_star = sigma_u$r_star))
  object
}

#' Build the prior blocks shared by both inefficiency distributions
#'
#' @param object a model object.
#' @param coef a named list with elements \code{mu} and \code{v_i}.
#' @param sigma a named list with elements \code{shape} and \code{rate}.
#'
#' @return A list with \code{b0}, \code{B0i}, \code{shape_v} and \code{rate_v}.
#'
#' @keywords internal
prior_coef_sigma <- function(object, coef, sigma) {

  k <- object$k

  coef <- merge_prior_list(coef, list(mu = 0, v_i = 0.01), "coef")
  sigma <- merge_prior_list(sigma, list(shape = 0.01, rate = 0.01), "sigma")

  b0 <- coef$mu
  if (length(b0) == 1L) {
    b0 <- rep(b0, k)
  }
  if (length(b0) != k) {
    stop("'coef$mu' must be of length 1 or ", k, ".")
  }

  B0i <- coef$v_i
  if (length(B0i) == 1L) {
    B0i <- diag(B0i, k)
  }
  B0i <- as.matrix(B0i)
  if (!identical(dim(B0i), c(as.integer(k), as.integer(k)))) {
    stop("'coef$v_i' must be a scalar or a ", k, " x ", k, " matrix.")
  }

  # Zeros are allowed: they give the improper limiting prior p(h) proportional
  # to 1/h, which is the non-informative choice used in much of the literature
  # and still leaves a proper posterior here.
  if (sigma$shape < 0 || sigma$rate < 0) {
    stop("The shape and rate of the prior on the error precision must be ",
         "non-negative.")
  }

  list(b0 = b0, B0i = B0i, shape_v = sigma$shape, rate_v = sigma$rate)
}

#' Merge a user-supplied prior list into its defaults
#'
#' Unlike \code{modifyList}, this reports unknown element names rather than
#' silently adding them, which otherwise makes typos in a list-based interface
#' very hard to spot.
#'
#' @param user the list supplied by the user, or \code{NULL}.
#' @param default the list of default values.
#' @param what the argument name, used in error messages.
#'
#' @return The merged list.
#'
#' @keywords internal
merge_prior_list <- function(user, default, what) {

  if (is.null(user)) {
    return(default)
  }
  if (!is.list(user)) {
    stop("'", what, "' must be a named list.")
  }

  unknown <- setdiff(names(user), names(default))
  if (length(unknown) > 0) {
    stop("Unknown element(s) in '", what, "': ",
         paste(unknown, collapse = ", "), ". Expected: ",
         paste(names(default), collapse = ", "), ".")
  }

  default[names(user)] <- user
  default
}

#' Check that a prior median efficiency is a valid probability
#'
#' @param r_star the value to check.
#'
#' @return Invisibly \code{TRUE}; called for the error it raises.
#'
#' @keywords internal
check_r_star <- function(r_star) {
  if (length(r_star) != 1L || is.na(r_star) || r_star <= 0 || r_star >= 1) {
    stop("'r_star' must be a single number strictly between 0 and 1.")
  }
  invisible(TRUE)
}
