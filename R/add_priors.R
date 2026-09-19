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
#' Because \code{r_star} is an efficiency, it fixes the scale of \eqn{u} in the
#' units of the response. This is only meaningful if the response is on a
#' logarithmic scale, where \eqn{u} is a proportional shortfall. A model fitted
#' to output in levels will return efficiency scores that depend on the units
#' the output happens to be measured in.
#'
#' In both models the anchor is matched on the marginal prior of \eqn{u},
#' that is after the parameter of the inefficiency distribution has been
#' integrated out, so that the prior median of \eqn{\exp(-u)} is \eqn{r^*}
#' exactly, whatever shape is used.
#'
#' For the exponential model, a gamma prior of shape \eqn{a} and rate \eqn{b}
#' on \eqn{\lambda} implies the marginal prior density
#' \deqn{p(u) = \int_0^\infty \lambda e^{-\lambda u}
#'              \frac{b^a}{\Gamma(a)} \lambda^{a - 1} e^{-b\lambda} d\lambda
#'            = \frac{a b^a}{(u + b)^{a + 1}},}
#' a Lomax density whose median is \eqn{b (2^{1/a} - 1)}, so the rate is set to
#' \eqn{-\log(r^*) / (2^{1/a} - 1)}. At the default shape of 1 this is simply
#' \eqn{-\log(r^*)} and the density is \eqn{c / (u + c)^2}, which is the
#' elicitation of van den Broeck, Koop, Osiewalski and Steel (1994); their two
#' degrees of freedom are that default.
#'
#' For the half-normal model the prior sits on \eqn{\sigma_u^2} and is inverse
#' gamma, so \eqn{u} is marginally half t on \eqn{2a} degrees of freedom with
#' scale \eqn{\sqrt{b / a}}, and the rate is set to
#' \eqn{a \left( -\log(r^*) / t_{2a}^{-1}(0.75) \right)^2}. The shape must
#' still exceed one, because the starting values that \code{method = "ols"}
#' derives use the prior mean of \eqn{\sigma_u^2}, hence the higher default.
#'
#' Earlier versions matched the half-normal anchor through the prior mean of
#' \eqn{\sigma_u^2} rather than through this median, which left it optimistic
#' by a factor that depended only on the shape: at the default of 2.5 the
#' implied prior median efficiency was \eqn{(r^*)^{0.834}}, so that asking for
#' 0.75 gave 0.787 and asking for 0.5 gave 0.561. The exponential model was
#' matched exactly at shape 1 but not at any other, where the error was larger
#' still: at shape 2, \eqn{r^* = 0.75} implied 0.888.
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
  el <- elicit_exp(lambda, "lambda")

  object$priors <- c(prior_coef_sigma(object, coef, sigma),
                     list(shape_u = el$shape, rate_u = el$rate,
                          r_star = lambda$r_star))
  drop_stale_posterior(object)
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
  el <- elicit_hn(sigma_u, "sigma_u")

  object$priors <- c(prior_coef_sigma(object, coef, sigma),
                     list(shape_u = el$shape, rate_u = el$rate,
                          r_star = sigma_u$r_star))
  drop_stale_posterior(object)
}

#' Elicit the gamma prior on an exponential inefficiency rate
#'
#' Sets the rate so that the median of the marginal prior on \eqn{u}, after
#' the rate parameter has been integrated out, falls at \code{-log(r_star)}.
#' That marginal is a Lomax distribution of shape \code{shape} and scale
#' \code{rate}, whose median is \code{rate * (2^(1/shape) - 1)}, so the
#' scale has to be divided by that factor. At the default shape of one the
#' factor is one and the rate is \code{-log(r_star)}, which is the
#' elicitation of van den Broeck, Koop, Osiewalski and Steel (1994).
#'
#' @param spec a list with elements \code{r_star} and \code{shape}.
#' @param what the name of the argument the specification came from, used in
#'   error messages.
#'
#' @return A list with \code{shape} and \code{rate}.
#'
#' @keywords internal
elicit_exp <- function(spec, what) {
  check_probability(spec$r_star, paste0(what, "$r_star"))
  check_shape(spec$shape, what, min = 0)
  list(shape = spec$shape,
       rate = -log(spec$r_star) / (2^(1 / spec$shape) - 1))
}

#' Elicit the gamma prior on a half-normal inefficiency scale
#'
#' Sets the rate so that the median of the marginal prior on \eqn{u}, after
#' the scale has been integrated out, falls at \code{-log(r_star)}. With an
#' inverse gamma prior on \eqn{\sigma_u^2} that marginal is a half t on
#' \code{2 * shape} degrees of freedom and scale \code{sqrt(rate / shape)},
#' so the rate follows from its median.
#'
#' @param spec a list with elements \code{r_star} and \code{shape}.
#' @param what the name of the argument the specification came from, used in
#'   error messages.
#'
#' @return A list with \code{shape} and \code{rate}.
#'
#' @keywords internal
elicit_hn <- function(spec, what) {
  check_probability(spec$r_star, paste0(what, "$r_star"))
  check_shape(spec$shape, what, min = 1)
  scale2 <- (-log(spec$r_star) / stats::qt(0.75, df = 2 * spec$shape))^2
  list(shape = spec$shape, rate = spec$shape * scale2)
}

#' Check the shape of a gamma prior on an inefficiency parameter
#'
#' The exponential model needs a positive shape for the prior to be proper.
#' The half-normal model needs one above 1 as well, though no longer for the
#' reason it once did: the anchor used to be matched through the prior mean
#' of the squared scale, which is \code{rate / (shape - 1)} and exists only
#' then, and is now matched through the median of the marginal prior, which
#' exists for any positive shape. The bound stays because the starting value
#' \code{method = "ols"} derives is still that mean.
#'
#' @param shape the value to check.
#' @param what the name of the argument it came from.
#' @param min the value the shape has to exceed, 0 or 1.
#'
#' @return Invisibly \code{TRUE}; called for the error it raises.
#'
#' @keywords internal
check_shape <- function(shape, what, min) {
  if (length(shape) != 1L || !is.numeric(shape) || !is.finite(shape) ||
      shape <= min) {
    stop("The shape of a ",
         if (min == 0) "gamma prior on an exponential rate must be positive"
         else paste("half-normal scale prior must exceed 1, so that the",
                    "prior mean of its square exists and can be used as a",
                    "starting value"),
         "; '", what, "$shape' is not.")
  }
  invisible(TRUE)
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
  if (!is.numeric(b0) || !all(is.finite(b0))) {
    stop("'coef$mu' must be finite numbers.")
  }
  if (length(b0) == 1L) {
    b0 <- rep(b0, k)
  }
  if (length(b0) != k) {
    stop("'coef$mu' must be of length 1 or ", k, ".")
  }

  B0i <- coef$v_i
  if (!is.numeric(B0i) || !all(is.finite(B0i))) {
    stop("'coef$v_i' must be finite numbers.")
  }
  if (length(B0i) == 1L) {
    B0i <- diag(B0i, k)
  }
  B0i <- as.matrix(B0i)
  if (!identical(dim(B0i), c(as.integer(k), as.integer(k)))) {
    stop("'coef$v_i' must be a scalar or a ", k, " x ", k, " matrix.")
  }
  check_precision(B0i)

  # Zeros are allowed: they give the improper limiting prior p(h) proportional
  # to 1/h, which is the non-informative choice used in much of the literature
  # and still leaves a proper posterior here.
  if (length(sigma$shape) != 1L || length(sigma$rate) != 1L ||
      !is.numeric(sigma$shape) || !is.numeric(sigma$rate) ||
      !is.finite(sigma$shape) || !is.finite(sigma$rate) ||
      sigma$shape < 0 || sigma$rate < 0) {
    stop("The shape and rate of the prior on the error precision must be ",
         "non-negative.")
  }

  list(b0 = b0, B0i = B0i, shape_v = sigma$shape, rate_v = sigma$rate)
}

#' Check that a prior precision matrix is one
#'
#' A precision matrix has to be symmetric and positive semi-definite. Neither
#' is checked by the sampler, which reaches \code{inv_sympd()} with whatever it
#' is given: an asymmetric matrix draws one Armadillo warning per sweep and
#' then samples from a prior that is not the one asked for, and an indefinite
#' one aborts part way through a run with a message about the inverse rather
#' than about the prior.
#'
#' A matrix of zeros passes, since it is the flat limiting prior on the
#' coefficients, the counterpart of the improper prior the error precision
#' already allows.
#'
#' @param B0i the matrix to check.
#'
#' @return Invisibly \code{TRUE}; called for the error it raises.
#'
#' @keywords internal
check_precision <- function(B0i) {

  if (!isSymmetric(unname(B0i))) {
    stop("'coef$v_i' must be symmetric, since it is a precision matrix.")
  }

  ev <- eigen(B0i, symmetric = TRUE, only.values = TRUE)$values
  tol <- sqrt(.Machine$double.eps) * max(1, max(abs(ev)))
  if (min(ev) < -tol) {
    stop("'coef$v_i' must be positive semi-definite, since it is a precision ",
         "matrix; its smallest eigenvalue is ", format(min(ev)), ".")
  }

  invisible(TRUE)
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
