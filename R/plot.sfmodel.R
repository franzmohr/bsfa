#' Plotting draws of a Bayesian stochastic frontier model
#'
#' A plot function for objects of class \code{"sfmodel"}.
#'
#' @param x an object of class \code{"sfmodel_exp"} or \code{"sfmodel_hn"},
#'   usually the result of a call to \code{\link{add_posterior_coefficients}}.
#' @param type either \code{"hist"} (default) for histograms, \code{"trace"}
#'   for trace plots, \code{"boxplot"} for boxplots, or \code{"efficiency"} for
#'   the distribution of the unit-level efficiency scores.
#' @param ci interval used to calculate the credible bands drawn with
#'   \code{type = "efficiency"}.
#' @param units for \code{type = "efficiency"}, the number of units shown in the
#'   caterpillar plot, taken evenly across the ranked scores. Defaults to 50.
#'   \code{NULL} or \code{Inf} shows all of them.
#' @param max_cols an integer of the maximum number of panels per row. Defaults
#'   to 3.
#' @param ... further graphical parameters.
#'
#' @details With \code{type = "efficiency"} the function draws two panels: a
#'   histogram of the posterior mean efficiency across units, which is the
#'   summary usually reported in applied work, and a caterpillar plot of the
#'   ranked scores with their credible bands. The second panel is worth looking
#'   at, because the first hides how little the data often say about any
#'   individual unit. A cross-sectional model in particular has one observation
#'   per inefficiency term, so the bands can be wide enough to cover most of the
#'   range the histogram spreads over.
#'
#' @return A plot of the posterior draws.
#'
#' @seealso \code{\link{efficiency}} for the underlying numbers.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
#'
#' model <- create_sfmodel_exp(y ~ x1 + x2, data = d,
#'                             iterations = 500, burnin = 200)
#' model <- add_posterior_coefficients(add_priors(model))
#'
#' plot(model)
#' plot(model, type = "trace")
#' plot(model, type = "efficiency")
#'
#' @export
plot.sfmodel <- function(x, type = c("hist", "trace", "boxplot", "efficiency"),
                         ci = 0.95, units = 50, max_cols = 3, ...) {

  type <- match.arg(type)

  if (type == "efficiency") {
    return(invisible(plot_sf_efficiency(x, ci = ci, units = units, ...)))
  }

  draws <- sf_par_draws(x)
  nms <- colnames(draws)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = grid_dim(ncol(draws), max_cols),
                mar = c(4, 4, 2.5, 1))

  for (j in seq_along(nms)) {
    switch(
      type,
      hist = graphics::hist(draws[, j], main = nms[j], xlab = NULL,
                            col = grDevices::grey(0.85), breaks = 25, ...),
      trace = graphics::plot(stats::time(x$posterior$beta$coeffs),
                             draws[, j], type = "l", main = nms[j],
                             xlab = "Iteration", ylab = NULL, ...),
      boxplot = graphics::boxplot(draws[, j], main = nms[j], ...))
  }

  invisible(x)
}

#' Collect the parameter draws of a stochastic frontier model
#'
#' Assembles the frontier coefficients and the scalar parameters into one
#' matrix: the error standard deviation and the inefficiency parameter for a
#' two-component model, and in addition the standard deviation of the unit
#' effect and the second inefficiency parameter for a four-component one.
#'
#' For the two-component half-normal model the signal-to-noise ratio is
#' appended. It is only defined there, because both components are then scale
#' parameters of the same kind, and because the four-component model has two
#' one-sided terms with no single such ratio between them.
#'
#' @param object a model object carrying posterior draws.
#'
#' @return A matrix with one row per draw and one column per parameter.
#'
#' @keywords internal
sf_par_draws <- function(object) {

  four <- identical(object$model$components, 4L)
  blocks <- sf_scalar_blocks(object)
  check_posterior_blocks(object, c("beta", blocks))

  out <- as.matrix(object$posterior$beta$coeffs)
  for (b in blocks) {
    out <- cbind(out, as.matrix(object$posterior[[b]]$coeffs))
  }

  if (object$model$ineff == "halfnormal" && !four) {
    out <- cbind(out,
                 lambda = out[, object$model$par_u_name] / out[, "sigma_v"])
  }

  out
}

#' Panel layout for a given number of plots
#'
#' @param n the number of panels.
#' @param max_cols the maximum number of panels per row.
#'
#' @return A vector of length two for \code{par(mfrow = )}.
#'
#' @keywords internal
grid_dim <- function(n, max_cols) {
  cols <- min(n, max_cols)
  c(ceiling(n / cols), cols)
}

#' Draw the efficiency panels
#'
#' @param object a model object carrying inefficiency draws.
#' @param ci the probability of the credible band.
#' @param units the number of units shown in the caterpillar plot.
#' @param ... further graphical parameters.
#'
#' @return Invisibly, the efficiency summary that was plotted.
#'
#' @keywords internal
plot_sf_efficiency <- function(object, ci, units, ...) {

  check_probability(ci, "ci")

  probs <- c((1 - ci) / 2, 0.5, 1 - (1 - ci) / 2)
  # For a four-component model this is the overall efficiency, the default of
  # its efficiency method and the quantity a two-component model would report.
  eff <- efficiency(object, probs = probs)

  band <- paste0(format(100 * probs[c(1, 3)], trim = TRUE), "%")
  low <- eff[[band[1]]]
  high <- eff[[band[2]]]

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))

  graphics::hist(eff$mean, breaks = 25, col = grDevices::grey(0.85),
                 main = "Posterior means", xlab = "Efficiency", ...)

  ord <- order(eff$mean)
  # Thinning the ranked units keeps the caterpillar readable for large panels
  # while still spanning the whole range of scores.
  if (!is.null(units) && is.finite(units) && units < length(ord)) {
    ord <- ord[round(seq(1, length(ord), length.out = units))]
  }
  pos <- seq_along(ord)

  graphics::plot(eff$mean[ord], pos, type = "n",
                 xlim = range(c(low[ord], high[ord])),
                 main = paste0(format(100 * ci), "% credible bands"),
                 xlab = "Efficiency", ylab = "Unit, ranked")
  graphics::segments(low[ord], pos, high[ord], pos,
                     col = grDevices::grey(0.6))
  graphics::points(eff$mean[ord], pos, pch = 16, cex = 0.5)

  invisible(eff)
}
