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
#' The response is assumed to be on a logarithmic scale, and the regressors
#' usually are as well. That is what makes \eqn{\exp(-u_i)} an efficiency:
#' \eqn{u_i} is then a proportional shortfall, so a unit with \eqn{u_i = 0.1}
#' produces about ten per cent less than its inputs allow. Fitting output in
#' levels does not merely change the units of \eqn{u}, it changes the answer,
#' because the prior on the inefficiency term is expressed in the units of the
#' response: the same data measured in euros and in millions of euros give
#' different efficiency scores. See \code{\link{add_priors}}.
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

  # The identifier is resolved to a vector here and then handed to
  # model.frame(), so that na.omit() drops it in step with the response and the
  # regressors. Aligning it afterwards from rownames() does not work: the
  # rownames of a data frame need not be integers, and coercing them silently
  # yields NA units, which collapses the whole panel into one.
  id_var <- NULL
  if (!is.null(id)) {
    if (length(id) == 1L && is.character(id)) {
      if (!id %in% names(data)) {
        stop("Variable '", id, "' not found in 'data'.")
      }
      id_var <- data[[id]]
    } else {
      id_var <- id
    }
    if (is.data.frame(data) && length(id_var) != nrow(data)) {
      stop("'id' must name a variable in 'data' or have one element per row ",
           "of it; got ", length(id_var), " for ", nrow(data), " rows.")
    }
  }

  # Rows are dropped explicitly rather than by na.action, so that the
  # identifier is filtered by the same logical vector as the model frame. An
  # observation whose unit is unknown is dropped as well, since it cannot be
  # assigned an inefficiency term.
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  keep <- stats::complete.cases(mf)
  if (!is.null(id_var)) {
    keep <- keep & !is.na(id_var)
  }
  if (!any(keep)) {
    stop("No complete observations remain after dropping missing values.")
  }
  mf <- mf[keep, , drop = FALSE]
  if (!is.null(id_var)) {
    id_var <- id_var[keep]
  }
  mt <- attr(mf, "terms")
  y <- as.numeric(stats::model.response(mf))
  X <- stats::model.matrix(mt, mf)
  n <- length(y)
  k <- ncol(X)

  if (n <= k) {
    stop("Model has at least as many coefficients as observations.")
  }

  # Map observations to the units that own the inefficiency terms.
  if (is.null(id_var)) {
    g <- seq_len(n)
    unit_labels <- rownames(mf)
    panel <- FALSE
  } else {
    fid <- factor(id_var)
    g <- as.integer(fid)
    unit_labels <- levels(fid)
    panel <- TRUE
  }

  list(data = list(y = y, X = X, g = g, n_units = length(unique(g)),
                   unit_labels = unit_labels),
       model = list(ineff = ineff, type = type, panel = panel,
                    components = 2L,
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

  cat(if (identical(x$model$components, 4L))
        "Four-component Bayesian stochastic frontier model\n\n" else
        "Bayesian stochastic frontier model\n\n")
  cat("Frontier:           ", x$model$type, "\n", sep = "")
  cat("Inefficiency:       ", x$model$ineff,
      if (identical(x$model$components, 4L))
        ", persistent and transient" else "", "\n", sep = "")
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
    r_star <- x$priors$r_star
    if (length(r_star) > 1L) {
      cat("\nPrior median efficiency, per component:\n")
      print(r_star)
    } else {
      cat("\nPrior median efficiency: ", format(r_star), "\n", sep = "")
    }
  }

  invisible(x)
}
