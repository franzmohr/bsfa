#' Simulate data from a stochastic frontier model
#'
#' Generates artificial data from the model estimated by \code{\link{bsfa}},
#' for use in examples, tests and prior predictive checks.
#'
#' @param n number of units. With \code{n_time > 1} the returned data frame has
#'   \code{n * n_time} rows and the inefficiency term is held fixed within a
#'   unit, as in the time-invariant panel model.
#' @param beta frontier coefficients. The first element is the intercept, so a
#'   vector of length \code{k} implies \code{k - 1} regressors named
#'   \code{x1}, ..., \code{x(k-1)}.
#' @param sigma_v standard deviation of the symmetric error.
#' @param par_u parameter of the inefficiency distribution: the rate for
#'   \code{ineff = "exponential"}, the scale for \code{ineff = "halfnormal"}.
#' @param ineff distribution of the inefficiency term.
#' @param type either \code{"production"} or \code{"cost"}.
#' @param n_time number of periods per unit.
#'
#' @return A data frame with the dependent variable \code{y}, the regressors,
#'   and a unit identifier \code{id}. The realised inefficiency terms and
#'   efficiency scores are attached as the attributes \code{"u"} and
#'   \code{"efficiency"}.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
#' summary(attr(d, "efficiency"))
#'
#' @export
sim_sf <- function(n,
                   beta = c(1, 0.5, 0.3),
                   sigma_v = 0.2,
                   par_u = 4,
                   ineff = c("exponential", "halfnormal"),
                   type = c("production", "cost"),
                   n_time = 1) {

  ineff <- match.arg(ineff)
  type <- match.arg(type)

  k <- length(beta)
  n_obs <- n * n_time

  X <- matrix(stats::rnorm(n_obs * (k - 1)), nrow = n_obs)
  colnames(X) <- paste0("x", seq_len(k - 1))

  u <- switch(ineff,
              exponential = stats::rexp(n, rate = par_u),
              halfnormal = abs(stats::rnorm(n, sd = par_u)))

  id <- rep(seq_len(n), each = n_time)
  s <- if (type == "production") -1 else 1

  y <- as.numeric(cbind(1, X) %*% beta) + s * u[id] +
    stats::rnorm(n_obs, sd = sigma_v)

  out <- data.frame(y = y, X, id = id)
  attr(out, "u") <- u
  attr(out, "efficiency") <- exp(-u)
  out
}
