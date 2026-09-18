#' How far do the efficiency scores move with the prior?
#'
#' Refits a model across a grid of prior median efficiencies and reports how
#' much the scores change, which is a direct measurement of the thing
#' \code{\link{skewness_test}} tests a proxy for.
#'
#' The inefficiency term is the part of a frontier model the data speak to
#' least, and its prior is expressed as an anchor, \code{r_star}, in the units
#' of the response. Where the sample is informative the posterior barely
#' notices the anchor; where it is not, the scores follow it. Since that
#' difference decides whether the scores can be reported as findings, it is
#' worth measuring rather than assuming, and it costs one refit per grid point.
#'
#' This is the Bayesian form of the question. \code{\link{skewness_test}} asks
#' whether the residuals carry the asymmetry a one-sided term would produce,
#' which is a property of the data and needs no sampling, but is a test, with a
#' test's error rate and its asymptotics. This asks what the anchor does to the
#' answer, which is the quantity of interest itself and has no null
#' distribution to get wrong. On a sample simulated with no inefficiency at
#' all, the mean score moved by 0.036 as the anchor ran from 0.5 to 0.9,
#' against 0.005 for a sample that had plenty; the skewness test called that
#' same empty sample informative.
#'
#' Every fit uses one seed, the model's own where it has one, so that the table
#' shows the prior moving the answer rather than the sampler doing it. The rest
#' of the prior specification is carried over unchanged, so a model that
#' already has priors keeps them and varies only the anchor. For a
#' four-component model both anchors, persistent and transient, move together.
#'
#' The frontier coefficients are reported alongside, and are usually the
#' reassurance: they hold still while the scores move, which is the pattern to
#' expect and says the fitted frontier is not in question even where the
#' distance to it is.
#'
#' @param object an object of class \code{"sfmodel"}. Priors need not have been
#'   added; where they have, everything but the anchor is kept.
#' @param r_star the prior median efficiencies to try, at least two.
#' @param keep_u passed to \code{\link{add_posterior_coefficients}}. The scores
#'   need the draws, but an interval may be given to thin them.
#' @param verbose whether to report progress, one message per fit.
#'
#' @return A list of class \code{"sfsens"}, with \code{efficiency}, a data
#'   frame of the score summaries by anchor; \code{coefficients}, the posterior
#'   means of the frontier coefficients by anchor; \code{spread}, the range of
#'   the mean score; \code{coef_spread}; and the \code{seed} used throughout.
#'
#' @seealso \code{\link{skewness_test}} for the check that needs no sampling,
#'   \code{\link{add_priors}} for what the anchor means.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 200, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
#'
#' model <- create_sfmodel_exp(y ~ x1, data = d,
#'                             iterations = 500, burnin = 200)
#' prior_sensitivity(model, r_star = c(0.6, 0.9))
#'
#' @export
prior_sensitivity <- function(object, r_star = c(0.5, 0.75, 0.9),
                              keep_u = TRUE, verbose = FALSE) {

  if (!inherits(object, "sfmodel")) {
    stop("'object' must be a stochastic frontier model.")
  }
  if (length(r_star) < 2L) {
    stop("'r_star' needs at least two values: the point is to vary it.")
  }
  for (rs in r_star) {
    check_probability(rs, "r_star")
  }
  r_star <- sort(unique(as.numeric(r_star)))

  # One seed for every fit, so that what the table shows is the prior moving
  # the answer and not the sampler.
  seed <- object$model$seed
  if (is.null(seed)) {
    seed <- .draw_model_seed()
  }

  eff <- NULL
  coefs <- NULL
  for (rs in r_star) {
    if (isTRUE(verbose)) {
      message("Fitting at r_star = ", format(rs), " ...")
    }
    fit <- add_posterior_coefficients(
      add_seed(sensitivity_priors(object, rs), seed), keep_u = keep_u)
    s <- efficiency(fit)$mean
    eff <- rbind(eff, data.frame(r_star = rs, mean = mean(s),
                                 sd = stats::sd(s), min = min(s),
                                 max = max(s)))
    coefs <- rbind(coefs, coef.sfmodel(fit))
  }
  rownames(eff) <- NULL
  rownames(coefs) <- format(r_star)

  structure(list(efficiency = eff, coefficients = coefs,
                 spread = max(eff$mean) - min(eff$mean),
                 coef_spread = apply(coefs, 2, function(z) max(z) - min(z)),
                 seed = seed, model = object$model),
            class = "sfsens")
}

#' Rebuild a prior specification with a different anchor
#'
#' The anchor is one element of a specification that may have several, so it is
#' replaced rather than the rest discarded: the coefficient and error priors
#' are read back off the object and passed again, and the shape of the
#' inefficiency prior is kept, so that only the median efficiency moves.
#'
#' @param object a model object, with or without priors.
#' @param rs the prior median efficiency to use.
#'
#' @return The model object with priors attached.
#'
#' @keywords internal
sensitivity_priors <- function(object, rs) {

  p <- object$priors
  four <- identical(object$model$components, 4L)
  args <- list(object = object)

  if (!is.null(p)) {
    args$coef <- list(mu = p$b0, v_i = p$B0i)
    args$sigma <- list(shape = p$shape_v, rate = p$rate_v)
    if (four) {
      args$sigma_mu <- list(shape = p$shape_mu, rate = p$rate_mu)
    }
  }

  # The argument the anchor belongs to is named for the parameter it governs,
  # which is what the class decides.
  nm <- if (four) c(object$model$par_eta_name, object$model$par_u_name) else
    object$model$par_u_name
  shapes <- if (is.null(p)) rep(list(NULL), length(nm)) else
    if (four) list(p$shape_eta, p$shape_u) else list(p$shape_u)

  for (i in seq_along(nm)) {
    args[[nm[i]]] <- c(list(r_star = rs),
                       if (is.null(shapes[[i]])) NULL else
                         list(shape = shapes[[i]]))
  }

  suppressMessages(do.call(add_priors, args))
}

#' @export
print.sfsens <- function(x, digits = 4, ...) {

  cat("Sensitivity of the efficiency scores to the prior\n\n")
  cat("Frontier:     ", x$model$type, "\n", sep = "")
  cat("Inefficiency: ", x$model$ineff,
      if (identical(x$model$components, 4L))
        ", both anchors moved together" else "", "\n", sep = "")
  cat("Seed:         ", x$seed, " for every fit\n\n", sep = "")

  cat("Mean efficiency by prior median efficiency:\n")
  print(round(x$efficiency, digits))

  cat("\nFrontier coefficients:\n")
  print(round(x$coefficients, digits))

  cat("\n  scores move by  ", format(round(x$spread, digits)),
      " across the grid\n", sep = "")
  cat("  coefficients by ", format(round(max(x$coef_spread), digits)),
      " at most\n\n", sep = "")

  cat(strwrap(if (x$spread > 0.02) paste(
    "The scores follow the anchor. What they say about the level of",
    "efficiency is largely what the prior was told to say, and belongs in a",
    "report as a prior sensitivity rather than as a finding. The ranking of",
    "units may still be informative where the level is not.")
    else paste(
      "The scores barely move with the anchor, so they are estimated from the",
      "data rather than carried over from the prior.")), sep = "\n")
  cat("\n")

  invisible(x)
}
