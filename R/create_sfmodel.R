#' Specify a Bayesian stochastic frontier model
#'
#' Sets up the model object that the \code{add_*} functions operate on. The two
#' functions differ only in the distribution assumed for the one-sided
#' inefficiency term, which determines the class of the object and hence the
#' methods that apply to it.
#'
#' The model is
#' \deqn{y_i = x_i' \beta - u_{g(i)} + v_i, \quad v_i \sim N(0, \sigma_v^2),}
#' for a production frontier, with the sign of \eqn{u} reversed for a cost
#' frontier. The frontier \eqn{x_i'\beta} is the maximum output attainable with
#' the inputs \eqn{x_i}, and \eqn{\exp(-u_i)} the efficiency with which unit
#' \eqn{i} reaches it.
#'
#' \code{create_sfmodel_exp} assumes \eqn{u \sim Exp(\lambda)} and
#' \code{create_sfmodel_hn} assumes \eqn{u \sim N^+(0, \sigma_u^2)}. The choice
#' matters mainly for the prior, since the exponential model admits an exact
#' elicitation from a prior median efficiency; see \code{\link{add_priors}}.
#'
#' If \code{id} is supplied, one inefficiency term is drawn per unit and held
#' fixed over that unit's observations, which is the time-invariant panel model
#' of Pitt and Lee (1981). This attributes all persistent heterogeneity between
#' units to inefficiency, and is therefore a strong assumption whenever units
#' differ systematically for reasons unrelated to efficiency. Leaving \code{id}
#' at \code{NULL} gives every observation its own inefficiency term.
#'
#' @param formula a model formula describing the frontier.
#' @param data a data frame containing the variables in \code{formula}.
#' @param id optional. Either the name of a variable in \code{data} or a vector
#'   of the same length as the data, identifying the units that own the
#'   inefficiency terms. See details.
#' @param type either \code{"production"} or \code{"cost"}.
#' @param iterations number of iterations retained after burn-in, before
#'   thinning.
#' @param burnin number of discarded iterations.
#' @param thin thinning interval. \code{iterations} must be a multiple of it.
#'
#' @return An object of class \code{"sfmodel_exp"} or \code{"sfmodel_hn"}, both
#'   inheriting from \code{"sfmodel"}.
#'
#' @seealso \code{\link{add_priors}}, \code{\link{add_initial_values}},
#'   \code{\link{add_seed}}, \code{\link{add_posterior_coefficients}}
#'
#' @references
#' Pitt, M. M., & Lee, L.-F. (1981). The measurement and sources of technical
#' inefficiency in the Indonesian weaving industry. \emph{Journal of
#' Development Economics}, 9(1), 43--64.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
#'
#' model <- create_sfmodel_exp(y ~ x1 + x2, data = d,
#'                             iterations = 500, burnin = 200)
#' model
#'
#' @export
create_sfmodel_exp <- function(formula,
                               data,
                               id = NULL,
                               type = c("production", "cost"),
                               iterations = 20000,
                               burnin = 2000,
                               thin = 1) {

  object <- sfmodel_skeleton(formula = formula, data = data, id = id,
                             type = match.arg(type), iterations = iterations,
                             burnin = burnin, thin = thin,
                             ineff = "exponential", cl = match.call())
  class(object) <- c("sfmodel_exp", "sfmodel")
  object
}

#' @rdname create_sfmodel_exp
#' @export
create_sfmodel_hn <- function(formula,
                              data,
                              id = NULL,
                              type = c("production", "cost"),
                              iterations = 20000,
                              burnin = 2000,
                              thin = 1) {

  object <- sfmodel_skeleton(formula = formula, data = data, id = id,
                             type = match.arg(type), iterations = iterations,
                             burnin = burnin, thin = thin,
                             ineff = "halfnormal", cl = match.call())
  class(object) <- c("sfmodel_hn", "sfmodel")
  object
}

#' Assemble the common parts of a stochastic frontier model object
#'
#' @param formula a model formula.
#' @param data a data frame.
#' @param id unit identifier, or \code{NULL}.
#' @param type \code{"production"} or \code{"cost"}.
#' @param iterations,burnin,thin MCMC settings.
#' @param ineff \code{"exponential"} or \code{"halfnormal"}.
#' @param cl the originating call.
#'
#' @return A list with the model data and specification.
#'
#' @keywords internal
sfmodel_skeleton <- function(formula, data, id, type, iterations, burnin,
                             thin, ineff, cl) {

  if (iterations < 1 || burnin < 0 || thin < 1) {
    stop("'iterations' and 'thin' must be positive and 'burnin' ",
         "non-negative.")
  }
  if (iterations %% thin != 0) {
    stop("'iterations' must be a multiple of 'thin'.")
  }

  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  mt <- attr(mf, "terms")
  y <- as.numeric(stats::model.response(mf))
  X <- stats::model.matrix(mt, mf)
  n <- length(y)
  k <- ncol(X)

  if (n <= k) {
    stop("Model has at least as many coefficients as observations.")
  }

  # Map observations to the units that own the inefficiency terms.
  if (is.null(id)) {
    g <- seq_len(n)
    unit_labels <- rownames(mf)
    panel <- FALSE
  } else {
    if (length(id) == 1L && is.character(id)) {
      if (!id %in% names(data)) {
        stop("Variable '", id, "' not found in 'data'.")
      }
      id <- data[[id]]
    }
    if (length(id) != nrow(mf)) {
      # Align with the rows that survived na.omit.
      id <- id[as.integer(rownames(mf))]
    }
    fid <- factor(id)
    g <- as.integer(fid)
    unit_labels <- levels(fid)
    panel <- TRUE
  }

  list(data = list(y = y, X = X, g = g, n_units = length(unique(g)),
                   unit_labels = unit_labels),
       model = list(ineff = ineff, type = type, panel = panel,
                    par_u_name = if (ineff == "exponential") "lambda" else
                      "sigma_u"),
       formula = formula,
       terms = mt,
       n = n,
       k = k,
       iterations = iterations,
       burnin = burnin,
       thin = thin,
       priors = NULL,
       initial = NULL,
       posterior = NULL,
       call = cl)
}

#' @export
print.sfmodel <- function(x, ...) {

  cat("Bayesian stochastic frontier model\n\n")
  cat("Frontier:           ", x$model$type, "\n", sep = "")
  cat("Inefficiency:       ", x$model$ineff, "\n", sep = "")
  cat("Observations:       ", x$n, "\n", sep = "")
  cat("Coefficients:       ", x$k, "\n", sep = "")
  cat("Inefficiency terms: ", x$data$n_units,
      if (x$model$panel) " (one per unit, time invariant)" else
        " (one per observation)", "\n", sep = "")
  cat("Iterations:         ", x$iterations, " after ", x$burnin,
      " burn-in, thinning ", x$thin, "\n", sep = "")

  cat("\nSpecification:\n")
  tick <- function(done) if (done) "set" else "not set"
  cat("  priors         ", tick(!is.null(x$priors)), "\n", sep = "")
  cat("  initial values ", tick(!is.null(x$initial)), "\n", sep = "")
  cat("  seed           ", tick(!is.null(x$model$seed)), "\n", sep = "")
  cat("  posterior      ",
      if (is.null(x$posterior)) "not simulated" else
        paste0("simulated",
               if (is.null(x$posterior$loglik)) "" else ", with loglik"),
      "\n", sep = "")

  if (!is.null(x$priors)) {
    cat("\nPrior median efficiency: ", format(x$priors$r_star), "\n", sep = "")
  }

  invisible(x)
}
