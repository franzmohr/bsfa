#' Summarising Bayesian stochastic frontier models
#'
#' Summary method for class \code{"sfmodel"}.
#'
#' @param object an object of class \code{"sfmodel_exp"} or
#'   \code{"sfmodel_hn"}, usually the result of a call to
#'   \code{\link{add_posterior_coefficients}}.
#' @param ci a numeric between 0 and 1 specifying the probability of the
#'   credible band. Defaults to 0.95.
#' @param x an object of class \code{"summary.sfmodel"}, usually the result of
#'   a call to \code{\link{summary.sfmodel}}.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{summary.sfmodel} returns a list of class
#'   \code{"summary.sfmodel"}, which contains the following components:
#'   \item{coefficients}{summary statistics of the posterior draws of the
#'     frontier coefficients and the two error components.}
#'   \item{efficiency}{summary statistics of the posterior mean efficiency
#'     across units, or \code{NULL} if the inefficiency draws were not kept.}
#'   \item{specifications}{a list containing information on the model
#'     specification.}
#'
#' @seealso \code{\link{efficiency}} for the unit-level scores.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
#'
#' model <- create_sfmodel_exp(y ~ x1 + x2, data = d,
#'                             iterations = 500, burnin = 200)
#' model <- add_posterior_coefficients(add_priors(model))
#'
#' summary(model)
#'
#' @export
summary.sfmodel <- function(object, ci = 0.95, ...) {

  check_probability(ci, "ci")

  pars <- sf_par_draws(object)

  probs <- c((1 - ci) / 2, 0.5, 1 - (1 - ci) / 2)
  qs <- t(apply(pars, 2, stats::quantile, probs = probs))
  colnames(qs) <- c(paste0(format(100 * probs[1], trim = TRUE), "%"),
                    "median",
                    paste0(format(100 * probs[3], trim = TRUE), "%"))

  tab <- cbind(mean = colMeans(pars),
               sd = apply(pars, 2, stats::sd),
               qs)

  # For a four-component model the efficiency method returns the overall score,
  # which is the one a two-component model would report; the split between the
  # persistent and transient parts is left to efficiency() itself.
  eff <- NULL
  if (!is.null(object$posterior$u$coeffs)) {
    r <- efficiency(object)$mean
    eff <- c(mean = mean(r), stats::quantile(r, probs = c(0, 0.5, 1)))
    names(eff) <- c("mean", "min", "median", "max")
  }

  structure(
    list(coefficients = tab,
         efficiency = eff,
         specifications = list(type = object$model$type,
                               ineff = object$model$ineff,
                               panel = object$model$panel,
                               components = object$model$components,
                               n = object$n,
                               k = object$k,
                               n_units = object$data$n_units,
                               iterations = object$iterations,
                               burnin = object$burnin,
                               thin = object$thin,
                               ci = ci,
                               call = object$call)),
    class = "summary.sfmodel")
}

#' @rdname summary.sfmodel
#' @export
print.summary.sfmodel <- function(x, digits = 4, ...) {

  spec <- x$specifications

  four <- identical(spec$components, 4L)

  cat(if (four) "Four-component Bayesian stochastic frontier model\n\n" else
        "Bayesian stochastic frontier model\n\n")
  cat("Call:\n")
  print(spec$call)
  cat("\nFrontier:           ", spec$type,
      "\nInefficiency:       ", spec$ineff,
      if (four) ", persistent and transient" else "",
      "\nObservations:       ", spec$n,
      "\nUnits:              ", spec$n_units,
      "\nDraws:              ", spec$iterations, " after ", spec$burnin,
      " burn-in, thinning ", spec$thin, "\n", sep = "")

  cat("\nPosterior summary, ", format(100 * spec$ci), "% credible bands:\n",
      sep = "")
  print(round(x$coefficients, digits))

  if (!is.null(x$efficiency)) {
    cat("\nPosterior mean ", if (four) "overall " else "",
        "efficiency, per ",
        if (four || !x$specifications$panel) "observation" else "unit", ":\n",
        sep = "")
    print(round(x$efficiency, digits))
  }

  invisible(x)
}
