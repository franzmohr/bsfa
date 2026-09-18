#' Simulate data from a four-component stochastic frontier model
#'
#' Generates artificial panel data from the model estimated by
#' \code{\link{create_sfmodel4_exp}} and \code{\link{create_sfmodel4_hn}}, for
#' use in examples, tests and prior predictive checks.
#'
#' @param n number of units.
#' @param n_time number of periods per unit. Two is the minimum the model can
#'   work with, and considerably more is needed before the split between the
#'   unit effect and persistent inefficiency is informative.
#' @param beta frontier coefficients. The first element is the intercept, so a
#'   vector of length \code{k} implies \code{k - 1} regressors named
#'   \code{x1}, ..., \code{x(k-1)}.
#' @param sigma_v standard deviation of the symmetric error.
#' @param sigma_mu standard deviation of the unit effect.
#' @param par_eta parameter of the distribution of persistent inefficiency: the
#'   rate for \code{ineff = "exponential"}, the scale for
#'   \code{ineff = "halfnormal"}.
#' @param par_u the same for transient inefficiency.
#' @param ineff distribution of the two one-sided terms.
#' @param type either \code{"production"} or \code{"cost"}.
#'
#' @return A data frame with the dependent variable \code{y}, the regressors
#'   and a unit identifier \code{id}. The realised components are attached as
#'   the attributes \code{"mu"}, \code{"eta"} and \code{"u"}, and the three
#'   efficiency scores as \code{"persistent"}, \code{"transient"} and
#'   \code{"overall"}.
#'
#' @seealso \code{\link{sim_sf}} for the two-component models.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf4(n = 50, n_time = 6, beta = c(1, 0.5, 0.3))
#' summary(attr(d, "persistent"))
#'
#' @export
sim_sf4 <- function(n = 100,
                    n_time = 8,
                    beta = c(1, 0.5, 0.3),
                    sigma_v = 0.2,
                    sigma_mu = 0.3,
                    par_eta = 8,
                    par_u = 6,
                    ineff = c("exponential", "halfnormal"),
                    type = c("production", "cost")) {

  ineff <- match.arg(ineff)
  type <- match.arg(type)

  if (n_time < 2) {
    stop("'n_time' must be at least 2.")
  }

  k <- length(beta)
  n_obs <- n * n_time

  X <- matrix(stats::rnorm(n_obs * (k - 1)), nrow = n_obs)
  colnames(X) <- paste0("x", seq_len(k - 1))

  draw_one_sided <- function(size, par) {
    switch(ineff,
           exponential = stats::rexp(size, rate = par),
           halfnormal = abs(stats::rnorm(size, sd = par)))
  }

  mu <- stats::rnorm(n, sd = sigma_mu)
  eta <- draw_one_sided(n, par_eta)
  u <- draw_one_sided(n_obs, par_u)

  id <- rep(seq_len(n), each = n_time)
  s <- if (type == "production") -1 else 1

  y <- as.numeric(cbind(1, X) %*% beta) + mu[id] + s * eta[id] + s * u +
    stats::rnorm(n_obs, sd = sigma_v)

  out <- data.frame(y = y, X, id = id)
  attr(out, "mu") <- mu
  attr(out, "eta") <- eta
  attr(out, "u") <- u
  attr(out, "persistent") <- exp(-eta)
  attr(out, "transient") <- exp(-u)
  attr(out, "overall") <- exp(-(eta[id] + u))
  out
}
