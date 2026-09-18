#' Add initial values to a four-component stochastic frontier model
#'
#' Attaches the starting values of the Gibbs sampler. Priors must be added
#' first, since both methods draw on them.
#'
#' With \code{method = "ols"} the coefficients start at their least squares
#' estimates. The least squares residual variance has to cover the symmetric
#' error, the unit effect and the two one-sided terms at once, so it is split
#' rather than handed to any one of them: half goes to \eqn{\sigma_v^2} and a
#' quarter to \eqn{\sigma_\mu^2}. The inefficiency parameters start at their
#' prior means and all four latent terms at zero.
#'
#' With \code{method = "prior"} the coefficients and the variances are drawn
#' from their priors instead, which disperses the starting points across chains.
#' That matters more here than in the two-component models, because the split
#' between the unit effect and persistent inefficiency is the part of this model
#' most likely to mix slowly.
#'
#' The function also stores the seed of the posterior simulation as element
#' \code{seed} of \code{object$model}, unless the model has one already.
#'
#' @param object an object of class \code{"sfmodel4_exp"} or
#'   \code{"sfmodel4_hn"}, with priors already added.
#' @param method either \code{"ols"} or \code{"prior"}. See details.
#' @param ... unused, for compatibility with the generic.
#'
#' @return The model object with the element \code{initial} attached.
#'
#' @seealso \code{\link{add_initial_values}} for the two-component models.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf4(n = 60, n_time = 8, beta = c(1, 0.5, 0.3))
#'
#' model <- create_sfmodel4_exp(y ~ x1 + x2, data = d, id = "id",
#'                              iterations = 500, burnin = 200)
#' model <- add_initial_values(add_priors(model))
#' model$initial$sigma_mu2
#'
#' @export
add_initial_values.sfmodel4_exp <- function(object, method = "ols", ...) {
  attach_initial(object, initial_four(object, method, "exponential"))
}

#' @rdname add_initial_values.sfmodel4_exp
#' @export
add_initial_values.sfmodel4_hn <- function(object, method = "ols", ...) {
  attach_initial(object, initial_four(object, method, "halfnormal"))
}

#' Starting values for a four-component model
#'
#' @param object a model object with priors attached.
#' @param method either \code{"ols"} or \code{"prior"}.
#' @param ineff either \code{"exponential"} or \code{"halfnormal"}.
#'
#' @return A list of starting values.
#'
#' @keywords internal
initial_four <- function(object, method, ineff) {

  init <- initial_common(object, method)

  # The least squares residual variance has to cover four sources of variation
  # at once, so it is split rather than assigned to any one of them.
  if (method == "ols") {
    init$sigma_mu2 <- init$sigma_v2 / 4
    init$sigma_v2 <- init$sigma_v2 / 2
  } else {
    init$sigma_mu2 <- 1 / stats::rgamma(1, shape = object$priors$shape_mu,
                                        rate = object$priors$rate_mu)
  }

  from_prior <- function(shape, rate) {
    if (ineff == "exponential") {
      if (method == "ols") shape / rate else
        stats::rgamma(1, shape = shape, rate = rate)
    } else {
      # The gamma prior sits on the precision, so the scale is the square root
      # of the mean of its inverse.
      if (method == "ols") sqrt(rate / (shape - 1)) else
        1 / sqrt(stats::rgamma(1, shape = shape, rate = rate))
    }
  }

  init$par_eta <- from_prior(object$priors$shape_eta, object$priors$rate_eta)
  init$par_u <- from_prior(object$priors$shape_u, object$priors$rate_u)

  init$mu <- rep(0, object$data$n_units)
  init$eta <- rep(0, object$data$n_units)
  init$u <- rep(0, object$n)

  init
}
