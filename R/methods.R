#' @export
print.bsfa <- function(x, ...) {
  cat("Bayesian stochastic frontier model\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nFrontier: ", x$model$type,
      "\nInefficiency: ", x$model$ineff,
      "\nObservations: ", x$n,
      "\nInefficiency terms: ", x$n_units, "\n", sep = "")
  cat("\nPosterior means of the coefficients:\n")
  print(round(colMeans(x$draws$beta), 4))
  invisible(x)
}

#' Summarise a Bayesian stochastic frontier model
#'
#' @param object an object of class \code{"bsfa"}.
#' @param probs quantiles of the posterior to report.
#' @param ... unused, for compatibility with the generic.
#'
#' @return An object of class \code{"summary.bsfa"}.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
#'
#' est <- create_sfmodel_exp(y ~ x1 + x2, data = d,
#'                           iterations = 500, burnin = 200)
#' est <- draw_posterior(add_priors(est))
#' summary(est)
#'
#' @export
summary.bsfa <- function(object, probs = c(0.05, 0.5, 0.95), ...) {

  par_u_name <- object$model$par_u_name
  par_u <- object$draws[[par_u_name]]

  pars <- cbind(object$draws$beta, sigma_v = object$draws$sigma_v, par_u)
  colnames(pars)[ncol(pars)] <- par_u_name

  # A signal-to-noise ratio is only defined for the half-normal model, where
  # both components are scale parameters of the same kind.
  if (object$model$ineff == "halfnormal") {
    pars <- cbind(pars, lambda = par_u / object$draws$sigma_v)
  }

  qs <- t(apply(pars, 2, stats::quantile, probs = probs))
  colnames(qs) <- paste0(format(100 * probs, trim = TRUE), "%")

  tab <- cbind(mean = colMeans(pars),
               sd = apply(pars, 2, stats::sd),
               qs)

  eff <- if (is.null(object$u)) NULL else {
    r <- exp(-object$u)
    c(mean = mean(r), min = min(colMeans(r)), max = max(colMeans(r)))
  }

  structure(list(coefficients = tab,
                 efficiency = eff,
                 call = object$call,
                 model = object$model,
                 n = object$n,
                 n_units = object$n_units,
                 mcmc = object$mcmc),
            class = "summary.bsfa")
}

#' @export
print.summary.bsfa <- function(x, digits = 4, ...) {
  cat("Bayesian stochastic frontier model\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nFrontier: ", x$model$type,
      "\nInefficiency: ", x$model$ineff,
      "\nObservations: ", x$n,
      "\nInefficiency terms: ", x$n_units,
      "\nDraws: ", x$mcmc$iterations, " retained after ", x$mcmc$burnin,
      " burn-in, thinning ", x$mcmc$thin, "\n", sep = "")
  cat("\nPosterior summary:\n")
  print(round(x$coefficients, digits))
  if (!is.null(x$efficiency)) {
    cat("\nPosterior mean efficiency: ", round(x$efficiency[["mean"]], digits),
        "\nRange of unit means: ", round(x$efficiency[["min"]], digits),
        " to ", round(x$efficiency[["max"]], digits), "\n", sep = "")
  }
  invisible(x)
}
