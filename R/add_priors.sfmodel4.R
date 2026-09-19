#' Add priors to a four-component stochastic frontier model
#'
#' Attaches the prior specification to a model object produced by
#' \code{\link{create_sfmodel4_exp}} or \code{\link{create_sfmodel4_hn}}.
#'
#' The blocks for the frontier coefficients and the error precision are those of
#' \code{\link{add_priors}} for the two-component models. Two further blocks are
#' needed here.
#'
#' The unit effect \eqn{\mu_i \sim N(0, \sigma_\mu^2)} takes a gamma prior on
#' its precision. With few units the posterior for \eqn{\sigma_\mu} is sensitive
#' to that choice, and the vague default is not a neutral one in the sense that
#' Gelman (2006) discusses; a sensitivity check is worth running before reading
#' much into the split between \eqn{\mu_i} and \eqn{\eta_i}.
#'
#' The two one-sided terms are elicited separately, each from its own prior
#' median efficiency, exactly as in the two-component case. Note that the
#' anchors apply per component: with \code{r_star = 0.9} for each, the prior
#' median of the overall efficiency \eqn{\exp(-\eta_i - u_{it})} is near 0.81
#' rather than 0.9. The default of 0.9 is chosen with that in mind.
#'
#' @param object an object of class \code{"sfmodel4_exp"} or
#'   \code{"sfmodel4_hn"}.
#' @param coef a named list of prior specifications for the frontier
#'   coefficients, with elements \code{mu}, the prior mean, and \code{v_i}, the
#'   prior precision.
#' @param sigma a named list of prior specifications for the error precision,
#'   with elements \code{shape} and \code{rate}.
#' @param sigma_mu a named list of prior specifications for the precision of
#'   the unit effect, with elements \code{shape} and \code{rate}.
#' @param lambda_eta a named list for the rate of the exponential distribution
#'   of persistent inefficiency, with elements \code{r_star} and \code{shape}.
#' @param lambda_u a named list for the rate of the exponential distribution of
#'   transient inefficiency, with elements \code{r_star} and \code{shape}.
#' @param sigma_eta a named list for the scale of the half-normal distribution
#'   of persistent inefficiency, with elements \code{r_star} and \code{shape}.
#' @param sigma_u a named list for the scale of the half-normal distribution of
#'   transient inefficiency, with elements \code{r_star} and \code{shape}.
#' @param varsel a named list of prior specifications for the variable
#'   selection algorithm. Required if the model was created with a
#'   \code{varsel} algorithm, and not allowed otherwise. It takes the same
#'   elements as in \code{\link{add_priors}} for the two-component models, and
#'   applies to the frontier coefficients only.
#' @param ... unused, for compatibility with the generic.
#'
#' @return The model object with the element \code{priors} attached.
#'
#' @references
#' Gelman, A. (2006). Prior distributions for variance parameters in
#' hierarchical models. \emph{Bayesian Analysis}, 1(3), 515--534.
#'
#' @seealso \code{\link{add_priors}} for the two-component models.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf4(n = 60, n_time = 8, beta = c(1, 0.5, 0.3))
#'
#' model <- create_sfmodel4_exp(y ~ x1 + x2, data = d, id = "id",
#'                              iterations = 500, burnin = 200)
#' model <- add_priors(model,
#'                     lambda_eta = list(r_star = 0.95),
#'                     lambda_u = list(r_star = 0.9))
#' model$priors$r_star
#'
#' @export
add_priors.sfmodel4_exp <- function(object,
                                    coef = list(mu = 0, v_i = 0.01),
                                    sigma = list(shape = 0.01, rate = 0.01),
                                    sigma_mu = list(shape = 0.01,
                                                    rate = 0.01),
                                    lambda_eta = list(r_star = 0.9, shape = 1),
                                    lambda_u = list(r_star = 0.9, shape = 1),
                                    varsel = NULL,
                                    ...) {

  defaults <- list(r_star = 0.9, shape = 1)
  eta <- merge_prior_list(lambda_eta, defaults, "lambda_eta")
  u <- merge_prior_list(lambda_u, defaults, "lambda_u")

  prior_four(object, coef, sigma, sigma_mu,
             elicit_exp(eta, "lambda_eta"), elicit_exp(u, "lambda_u"),
             c(persistent = eta$r_star, transient = u$r_star), varsel)
}

#' @rdname add_priors.sfmodel4_exp
#' @export
add_priors.sfmodel4_hn <- function(object,
                                   coef = list(mu = 0, v_i = 0.01),
                                   sigma = list(shape = 0.01, rate = 0.01),
                                   sigma_mu = list(shape = 0.01, rate = 0.01),
                                   sigma_eta = list(r_star = 0.9,
                                                    shape = 2.5),
                                   sigma_u = list(r_star = 0.9, shape = 2.5),
                                   varsel = NULL,
                                   ...) {

  defaults <- list(r_star = 0.9, shape = 2.5)
  eta <- merge_prior_list(sigma_eta, defaults, "sigma_eta")
  u <- merge_prior_list(sigma_u, defaults, "sigma_u")

  prior_four(object, coef, sigma, sigma_mu,
             elicit_hn(eta, "sigma_eta"), elicit_hn(u, "sigma_u"),
             c(persistent = eta$r_star, transient = u$r_star), varsel)
}

#' Assemble the prior blocks of a four-component model
#'
#' @param object a model object.
#' @param coef a named list with elements \code{mu} and \code{v_i}.
#' @param sigma a named list with elements \code{shape} and \code{rate}.
#' @param sigma_mu a named list with elements \code{shape} and \code{rate}.
#' @param eta the elicited prior on the persistent inefficiency parameter.
#' @param u the elicited prior on the transient inefficiency parameter.
#' @param r_star the two prior median efficiencies, named.
#' @param varsel the variable selection specification, or \code{NULL}.
#'
#' @return The model object with \code{priors} attached.
#'
#' @keywords internal
prior_four <- function(object, coef, sigma, sigma_mu, eta, u, r_star,
                       varsel = NULL) {

  sigma_mu <- merge_prior_list(sigma_mu, list(shape = 0.01, rate = 0.01),
                               "sigma_mu")
  if (length(sigma_mu$shape) != 1L || length(sigma_mu$rate) != 1L ||
      !is.numeric(sigma_mu$shape) || !is.numeric(sigma_mu$rate) ||
      !is.finite(sigma_mu$shape) || !is.finite(sigma_mu$rate) ||
      sigma_mu$shape <= 0 || sigma_mu$rate <= 0) {
    stop("The shape and rate of the prior on the precision of the unit ",
         "effect must be positive. Unlike the error precision, this one is ",
         "not identified without a proper prior when the number of units is ",
         "small.")
  }

  object$priors <- c(prior_coef_sigma(object, coef, sigma, varsel),
                     list(shape_mu = sigma_mu$shape,
                          rate_mu = sigma_mu$rate,
                          shape_eta = eta$shape,
                          rate_eta = eta$rate,
                          shape_u = u$shape,
                          rate_u = u$rate,
                          r_star = r_star))
  drop_stale_posterior(object)
}
