#' Is there any inefficiency in the data?
#'
#' Tests the least squares residuals for skewness in the direction the
#' frontier implies, which is the standard diagnostic for whether the data
#' support a one-sided term at all.
#'
#' A production frontier puts the composed error at \eqn{v_i - u_i}, which is
#' skewed to the left because \eqn{u} is non-negative; a cost frontier reverses
#' it. If the residuals of an ordinary least squares fit carry no such skew,
#' the sample contains no evidence of inefficiency. Schmidt and Lin (1984) call
#' this the wrong skew problem, and Waldman (1982) showed that it is exactly
#' the case in which the maximum likelihood estimate sits at the boundary,
#' \eqn{\sigma_u = 0}, so that the frontier collapses to the least squares fit.
#'
#' It matters more here than it does for a maximum likelihood estimator, not
#' less. A boundary estimate is conspicuous: the fitted model simply announces
#' that there is no inefficiency. A posterior cannot do that. The prior holds
#' the inefficiency term away from zero, so the sampler returns efficiency
#' scores of entirely reasonable appearance that are a reading of the prior
#' rather than of the data. On a sample simulated with no inefficiency
#' whatsoever, this package reports a mean efficiency near 0.89 and says
#' nothing. This test is the warning that is otherwise missing.
#'
#' The statistic is that of Coelli (1995). With \eqn{m_2} and \eqn{m_3} the
#' second and third moments of the residuals,
#' \deqn{M3T = m_3 / \sqrt{6 m_2^3 / n},}
#' which is standard normal when the residuals are symmetric. The p-value is
#' one sided, in the direction the frontier implies, so a small one is evidence
#' that there is inefficiency to estimate.
#'
#' The test needs no posterior draws and can be run on a model as soon as it is
#' created. A wrong or absent skew is not a reason to abandon the model: it is
#' a reason to treat the efficiency scores as prior driven, to check that the
#' response really is on a logarithmic scale, and to say as much when reporting
#' them.
#'
#' @param object an object of class \code{"sfmodel"}.
#'
#' @return A list of class \code{"sfskew"}, with the skewness of the residuals,
#'   the statistic, its one sided p-value, the direction the frontier implies
#'   and a verdict of \code{"informative"}, \code{"weak"} or \code{"wrong"}.
#'
#' @references
#' Coelli, T. (1995). Estimators and hypothesis tests for a stochastic frontier
#' function: A Monte Carlo analysis. \emph{Journal of Productivity Analysis},
#' 6(3), 247--268.
#'
#' Schmidt, P., & Lin, T.-F. (1984). Simple tests of alternative specifications
#' in stochastic frontier models. \emph{Journal of Econometrics}, 24(3),
#' 349--361.
#'
#' Waldman, D. M. (1982). A stationary point for the stochastic frontier
#' likelihood. \emph{Journal of Econometrics}, 18(2), 275--279.
#'
#' @seealso \code{\link{efficiency}} for the scores this qualifies,
#'   \code{\link{add_priors}} for the prior they fall back on.
#'
#' @examples
#' set.seed(1234)
#'
#' # A sample with inefficiency in it.
#' d <- sim_sf(n = 300, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
#' skewness_test(create_sfmodel_exp(y ~ x1, data = d))
#'
#' # And one without, where the efficiency scores would be the prior.
#' d0 <- data.frame(x1 = rnorm(300))
#' d0$y <- 1 + 0.5 * d0$x1 + rnorm(300, sd = 0.3)
#' skewness_test(create_sfmodel_exp(y ~ x1, data = d0))
#'
#' @export
skewness_test <- function(object) {

  if (!inherits(object, "sfmodel")) {
    stop("'object' must be a stochastic frontier model.")
  }

  r <- as.numeric(stats::.lm.fit(object$data$X, object$data$y)$residuals)
  n <- length(r)
  m2 <- mean(r^2)
  m3 <- mean(r^3)

  if (m2 <= 0) {
    stop("The least squares residuals have no variation, so there is no ",
         "skewness to measure.")
  }

  skew <- m3 / m2^1.5
  m3t <- m3 / sqrt(6 * m2^3 / n)

  # The frontier fixes which tail the evidence has to be in: a production
  # frontier subtracts the one-sided term and so skews the residuals left.
  s <- if (object$model$type == "production") -1 else 1
  p <- if (s < 0) stats::pnorm(m3t) else stats::pnorm(m3t, lower.tail = FALSE)

  verdict <- if (sign(m3) != s && m3 != 0) "wrong" else
    if (p > 0.1) "weak" else "informative"

  structure(list(skewness = skew, statistic = m3t, p.value = p,
                 n = n, expected = if (s < 0) "negative" else "positive",
                 type = object$model$type, verdict = verdict),
            class = "sfskew")
}

#' @export
print.sfskew <- function(x, digits = 4, ...) {

  cat("Skewness of the least squares residuals\n\n")
  cat("Frontier:     ", x$type, " (implies ", x$expected, " skewness)\n",
      sep = "")
  cat("Observations: ", x$n, "\n\n", sep = "")
  cat("  skewness       ", format(round(x$skewness, digits)), "\n", sep = "")
  cat("  M3T statistic  ", format(round(x$statistic, digits)), "\n", sep = "")
  cat("  p-value        ", format(round(x$p.value, digits)), "\n\n", sep = "")

  cat(strwrap(switch(
    x$verdict,
    informative = paste(
      "The residuals are skewed as the frontier implies, so the sample",
      "carries information about the one-sided term and the efficiency",
      "scores are estimated from it."),
    weak = paste(
      "The residuals are barely skewed. The sample says little about the",
      "one-sided term, so the efficiency scores will follow the prior more",
      "than the data. Check that the response is on a logarithmic scale,",
      "and vary r_star to see how far the scores move with it."),
    wrong = paste(
      "The residuals are skewed the wrong way, which is the case in which a",
      "maximum likelihood fit would collapse onto least squares. There is no",
      "evidence of inefficiency in this sample. The sampler will still return",
      "efficiency scores, but they will be a reading of the prior rather than",
      "of the data."))), sep = "\n")
  cat("\n")

  invisible(x)
}
