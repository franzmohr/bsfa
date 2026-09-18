#' Model selection criteria
#'
#' Generic function used to calculate selection criteria.
#'
#' @param object an object with suitable input data passed forward to method.
#' @param ci a numeric between 0 and 1 specifying the probability of the
#'   credible band. Defaults to 0.95.
#' @param ... arguments passed forward to method.
#'
#' @details The criteria are obtained from the draws of the log-likelihood that
#' \code{\link{add_posterior_loglik}} adds to a model. With \eqn{ll_t^{(i)}} the
#' log-likelihood of observation \eqn{t} in draw \eqn{i} of \eqn{R} draws,
#' \eqn{LL^{(i)} = \sum_{t = 1}^{T} ll_t^{(i)}} and \eqn{\kappa} the number of
#' estimated parameters, these are
#' \itemize{
#'  \item \code{"LL"}: the log-likelihood \eqn{LL^{(i)}}, summarised by the mean,
#' the median and the bounds of the credible band of its draws.
#'  \item \code{"AIC"}: \eqn{D + 2 \kappa};
#'  \item \code{"BIC"}: \eqn{D + \ln(T) \kappa};
#'  \item \code{"HQ"}: \eqn{D + 2 \ln(\ln(T)) \kappa};
#'  \item \code{"WAIC"}: \eqn{-2 \sum_{t = 1}^{T} \left( \ln \left(
#' \frac{1}{R} \sum_{i = 1}^{R} \exp(ll_t^{(i)}) \right) -
#' \mathrm{Var}_i \left[ ll_t^{(i)} \right] \right)}.
#' }
#'
#' \eqn{D} is the deviance of the model at its point estimate, obtained by
#' evaluating the log-likelihood once at the posterior mean of the parameters.
#' Since AIC, BIC and HQ correct the deviance of a fitted model for the optimism
#' of having fitted it, it is the deviance at the point estimate that their
#' penalties belong to. The mean of the deviance over the posterior is the
#' larger quantity, by about the effective number of parameters, so adding a
#' penalty to it would charge the complexity of the model a second time and tilt
#' every comparison towards the smaller model.
#'
#' \eqn{\kappa} counts the frontier coefficients, the standard deviation of the
#' symmetric error and the parameter of the inefficiency distribution. The
#' inefficiency terms themselves are latent variables rather than parameters and
#' are not counted, which is the convention but is worth stating: a model with
#' one \eqn{u} per observation has as many latent variables as observations, and
#' a count that included them would be meaningless.
#'
#' AIC, BIC and HQ are point estimates rather than quantities with a posterior
#' distribution, so their credible bands are \code{NA}. WAIC is a point estimate
#' as well, and its band is a normal interval built from its standard error,
#' which is the usual way it is reported.
#'
#' LOOIC is not provided. The Pareto smoothed importance sampling it needs is a
#' piece of machinery in its own right and deserves its own testing, so it is
#' left out rather than added hastily.
#'
#' The criteria require a model without an \code{id}, since that is the
#' condition under which the pointwise log-likelihood exists; see
#' \code{\link{add_posterior_loglik}}.
#'
#' @return A list of class \code{"selcrit"}, which also inherits the class of
#'   the model, with the element \code{model} and one data frame per criterion.
#'   Each has the columns \code{mean}, \code{median}, \code{qlower} and
#'   \code{qupper}, where bands that do not apply are \code{NA}.
#'
#' @seealso \code{\link{add_posterior_loglik}}
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
#'
#' # Compare the two inefficiency distributions on the same data
#' fit <- function(model) {
#'   add_posterior_loglik(add_posterior_coefficients(add_priors(model)))
#' }
#' m_exp <- fit(create_sfmodel_exp(y ~ x1 + x2, data = d,
#'                                 iterations = 500, burnin = 200))
#' m_hn <- fit(create_sfmodel_hn(y ~ x1 + x2, data = d,
#'                               iterations = 500, burnin = 200))
#'
#' selection_criteria(m_exp)$WAIC
#' selection_criteria(m_hn)$WAIC
#'
#' @export
selection_criteria <- function(object, ...) {
  UseMethod("selection_criteria")
}

#' @rdname selection_criteria
#' @export
selection_criteria.sfmodel4 <- function(object, ...) {
  stop("Selection criteria are not available for the four-component model, ",
       "because the pointwise log-likelihood they are built from is not. ",
       "Integrating out the unit effect and the persistent inefficiency term ",
       "couples the observations of a unit, so the likelihood does not ",
       "factorise over them. See ?add_posterior_loglik.")
}

#' @rdname selection_criteria
#' @export
selection_criteria.sfmodel <- function(object, ci = 0.95, ...) {

  if (is.null(object$posterior$loglik)) {
    stop("Object does not contain draws of the log-likelihood. ",
         if (object$model$panel)
           paste("A panel model cannot have one: the units share an",
                 "inefficiency term, so the likelihood does not factorise",
                 "over a unit's observations.") else
           "See ?add_posterior_loglik.")
  }
  check_probability(ci, "ci")

  ll <- as.matrix(object$posterior$loglik)
  tt <- ncol(ll)
  kappa <- object$k + 2

  probs <- c((1 - ci) / 2, 1 - (1 - ci) / 2)

  # Deviance at the posterior mean of the parameters, which is the point the
  # penalties of AIC, BIC and HQ belong to. See the details section.
  beta_bar <- colMeans(as.matrix(object$posterior$beta$coeffs))
  sigma_v_bar <- mean(object$posterior$sigma_v$coeffs)
  par_u_bar <- mean(object$posterior[[object$model$par_u_name]]$coeffs)
  deviance <- -2 * sum(sf_loglik_point(object, beta_bar, sigma_v_bar,
                                       par_u_bar))

  ll_draws <- rowSums(ll)
  ll_q <- stats::quantile(ll_draws, probs = probs)

  # The log of the mean of the likelihood, computed on the log scale so that a
  # badly fitting draw cannot overflow or underflow the exponential.
  lppd <- apply(ll, 2, function(z) {
    m <- max(z)
    m + log(mean(exp(z - m)))
  })
  p_waic <- apply(ll, 2, stats::var)
  elpd <- lppd - p_waic
  waic <- -2 * sum(elpd)
  waic_se <- 2 * sqrt(tt * stats::var(elpd))
  waic_band <- waic + stats::qnorm(probs) * waic_se

  point <- function(value, lower = NA_real_, upper = NA_real_) {
    data.frame(mean = value, median = NA_real_,
               qlower = lower, qupper = upper)
  }

  structure(
    list(model = object$model,
         LL = data.frame(mean = mean(ll_draws),
                         median = stats::median(ll_draws),
                         qlower = unname(ll_q[1]),
                         qupper = unname(ll_q[2])),
         AIC = point(deviance + 2 * kappa),
         BIC = point(deviance + log(tt) * kappa),
         HQ = point(deviance + 2 * log(log(tt)) * kappa),
         WAIC = point(waic, waic_band[1], waic_band[2])),
    class = c("selcrit", class(object)))
}

#' @export
print.selcrit <- function(x, digits = 4, ...) {

  crit <- x[setdiff(names(x), "model")]
  tab <- do.call(rbind, crit)
  rownames(tab) <- names(crit)

  cat("Selection criteria for a Bayesian stochastic frontier model\n\n")
  cat("Frontier:     ", x$model$type, "\n", sep = "")
  cat("Inefficiency: ", x$model$ineff, "\n\n", sep = "")
  print(round(tab, digits))

  invisible(x)
}
