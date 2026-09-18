#' Posterior efficiency scores of a four-component model
#'
#' Summarises the posterior distribution of the efficiency scores implied by
#' the augmented draws. The four-component model distinguishes three of them,
#' which is the main reason to fit it.
#'
#' \describe{
#'   \item{\code{"persistent"}}{\eqn{\exp(-\eta_i)}, one score per unit. This is
#'     the part of inefficiency that does not move: it reflects how the unit is
#'     set up rather than how it is run, and is what a policy or an investment
#'     would have to change.}
#'   \item{\code{"transient"}}{\eqn{\exp(-u_{it})}, one score per observation.
#'     The part that varies within a unit, and the one that answers whether
#'     efficiency improved over the sample.}
#'   \item{\code{"overall"}}{\eqn{\exp(-\eta_i - u_{it})}, the product of the
#'     two, one score per observation. This is what a two-component model
#'     reports, and reporting it alone throws away the split.}
#' }
#'
#' The unit effect \eqn{\mu_i} enters none of them. That is the point of the
#' model: persistent differences between units that are not inefficiency stay
#' out of the efficiency scores instead of depressing them.
#'
#' @param object an object of class \code{"sfmodel4_exp"} or
#'   \code{"sfmodel4_hn"}, estimated with \code{keep_u = TRUE}.
#' @param type one of \code{"overall"}, \code{"persistent"} or
#'   \code{"transient"}. See details.
#' @param probs quantiles of the posterior to report.
#' @param ... unused, for compatibility with the generic.
#'
#' @return A data frame with the posterior mean, standard deviation and the
#'   requested quantiles of the efficiency score. Its first column gives the
#'   unit, which for the observation-level types repeats.
#'
#' @seealso \code{\link{efficiency}} for the two-component models.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf4(n = 60, n_time = 8, beta = c(1, 0.5, 0.3))
#'
#' model <- create_sfmodel4_exp(y ~ x1 + x2, data = d, id = "id",
#'                              iterations = 500, burnin = 200)
#' model <- add_posterior_coefficients(add_priors(model))
#'
#' head(efficiency(model, type = "persistent"))
#' head(efficiency(model, type = "transient"))
#'
#' @export
efficiency.sfmodel4 <- function(object,
                                type = c("overall", "persistent",
                                         "transient"),
                                probs = c(0.05, 0.5, 0.95), ...) {

  type <- match.arg(type)

  if (is.null(object$posterior$eta$coeffs)) {
    stop("No inefficiency draws stored. Re-run add_posterior_coefficients() ",
         "with keep_u = TRUE.")
  }

  eta <- as.matrix(object$posterior$eta$coeffs)
  g <- object$data$g

  if (type == "persistent") {
    r <- exp(-eta)
    unit <- object$data$unit_labels
  } else {
    u <- as.matrix(object$posterior$u$coeffs)
    # The persistent term is one per unit, so it is expanded to the
    # observations of that unit before being combined with the transient one.
    r <- if (type == "transient") exp(-u) else exp(-(eta[, g, drop = FALSE] + u))
    unit <- object$data$unit_labels[g]
  }

  qs <- t(apply(r, 2, stats::quantile, probs = probs))
  colnames(qs) <- paste0(format(100 * probs, trim = TRUE), "%")

  data.frame(unit = unit,
             mean = colMeans(r),
             sd = apply(r, 2, stats::sd),
             qs,
             row.names = NULL,
             check.names = FALSE,
             stringsAsFactors = FALSE)
}
