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
#'     frontier coefficients and the two error components, ending in
#'     \code{ESS}, the effective sample size of each block, and, for a model
#'     with variable selection, in \code{PIP}, the posterior probability that
#'     the regressor belongs in the frontier. It is the number of
#'     independent draws the chain is worth, so the Monte Carlo error of a
#'     posterior mean is about its standard deviation divided by the square
#'     root of it. A value far below the number of draws is the sign that the
#'     chain needs to be longer, and in the four-component model it usually
#'     is.}
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

  # The effective sample size sits beside the band it qualifies. The draws are
  # autocorrelated, and in the four-component model severely so, where the unit
  # effect and persistent inefficiency trade off against each other from sweep
  # to sweep; without this the summary gives no sign of it.
  tab <- cbind(mean = colMeans(pars),
               sd = apply(pars, 2, stats::sd),
               qs,
               ESS = sf_ess(pars))

  # Under variable selection the posterior inclusion probability sits beside
  # the coefficient it belongs to. It is the share of draws in which the
  # regressor was in the frontier, and the mean and band on the same row are
  # averages over both states, so a coefficient with a low probability here has
  # a posterior concentrated near zero because most of its draws came from the
  # excluded half of the prior.
  tab <- add_pip_column(object, tab)

  varsel <- if (is.null(object$posterior$inclusion$coeffs)) NULL else
    object$model$varsel

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
                               varsel = varsel,
                               determinants = lapply(
                                 object$model$determinants,
                                 function(z) if (is.null(z)) NULL else
                                   z$names),
                               acceptance = object$acceptance,
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
      if (is.null(spec$varsel)) "" else
        paste0("\nVariable selection: ", toupper(spec$varsel)),
      "\nInefficiency:       ", spec$ineff,
      if (four) ", persistent and transient" else "",
      "\nObservations:       ", spec$n,
      "\nUnits:              ", spec$n_units,
      "\nDraws:              ", spec$iterations, " after ", spec$burnin,
      " burn-in, thinning ", spec$thin, "\n", sep = "")

  # The determinants are named beside the specification they belong to, so
  # that a row called 'z1' in the table below says which term it is
  # explaining rather than leaving that to be inferred from the formula.
  det <- spec$determinants
  det <- det[!vapply(det, is.null, logical(1))]
  if (length(det) > 0) {
    cat("Determinants:       ",
        paste(paste0(names(det), ": ",
                     vapply(det, paste, character(1), collapse = ", ")),
              collapse = "\n                     "), "\n", sep = "")
  }

  cat("\nPosterior summary, ", format(100 * spec$ci), "% credible bands:\n",
      sep = "")
  # Printed as a data frame so that the effective sample size, which is a count
  # of draws, is not formatted to the same decimals as a standard deviation.
  stat <- setdiff(colnames(x$coefficients), "ESS")
  tab <- as.data.frame(round(x$coefficients[, stat, drop = FALSE], digits))
  tab$ESS <- round(x$coefficients[, "ESS"])
  print(tab)

  # A Metropolis block that rarely moves has an effective sample size near
  # zero, and the table above would show that as a narrow band rather than as
  # a chain that never went anywhere. The rate is reported so that the two
  # cannot be confused.
  if (!is.null(spec$acceptance)) {
    cat("
Metropolis acceptance: ",
        paste(paste0(names(spec$acceptance), " ",
                     format(round(spec$acceptance, 3))), collapse = ", "),
        "
", sep = "")
    if (any(spec$acceptance < 0.05)) {
      cat("  A rate this low means the block barely moved; read its draws ",
          "with care.
", sep = "")
    }
  }

  if (!is.null(x$efficiency)) {
    cat("\nPosterior mean ", if (four) "overall " else "",
        "efficiency, per ",
        if (four || !x$specifications$panel) "observation" else "unit", ":\n",
        sep = "")
    print(round(x$efficiency, digits))
  }

  invisible(x)
}
