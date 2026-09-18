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
#' A term wrapped in \code{\link[stats]{offset}} enters the frontier with its
#' coefficient fixed at one, as it does in \code{\link[stats]{lm}}. This is how
#' a known quantity is imposed rather than estimated: a capacity that the
#' frontier cannot exceed, or a scale factor whose elasticity is known to be
#' one. The efficiency scores are measured from the frontier including it.
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
#'   inheriting from \code{"sfmodel"}. Rows with a missing value in the model
#'   frame, or with an unknown unit, are dropped with a message and recorded in
#'   the element \code{na.action}. An \code{\link[stats]{offset}} in the
#'   formula is subtracted from the response and kept as \code{data$offset},
#'   so \code{data$y} is the response the sampler sees rather than the one
#'   supplied.
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
#' @return A list with the model data and specification, and the element
#'   \code{na.action} recording any rows that were dropped.
#'
#' @keywords internal
sfmodel_skeleton <- function(formula, data, id, type, iterations, burnin,
                             thin, ineff, cl) {

  check_count(iterations, "iterations")
  check_count(burnin, "burnin")
  check_count(thin, "thin")
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
  # Dropping rows quietly is how a third of a panel goes missing unnoticed, so
  # the omission is reported and kept on the object in the form na.omit() uses.
  dropped <- which(!keep)
  if (length(dropped) > 0) {
    names(dropped) <- rownames(mf)[dropped]
    class(dropped) <- "omit"
    message("Dropping ", length(dropped), " of ", length(keep),
            " observations with missing values",
            if (!is.null(id_var) && anyNA(id_var))
              " or an unknown unit" else "", ".")
  } else {
    dropped <- NULL
  }
  mf <- mf[keep, , drop = FALSE]
  if (!is.null(id_var)) {
    id_var <- id_var[keep]
  }
  mt <- attr(mf, "terms")
  y <- as.numeric(stats::model.response(mf))
  X <- stats::model.matrix(mt, mf)

  # An offset is a term whose coefficient is fixed at one, so it is not a
  # column of the design matrix and model.matrix() leaves it out. Subtracting
  # it from the response here is what makes it enter the model at all.
  # Everything downstream -- the sampler, the log-likelihood, the residuals --
  # then works on the adjusted response, and fitted() adds it back.
  offs <- stats::model.offset(mf)
  if (!is.null(offs)) {
    offs <- as.numeric(offs)
    y <- y - offs
  }

  n <- length(y)
  k <- ncol(X)

  if (n <= k) {
    stop("Model has at least as many coefficients as observations.")
  }

  # A rank deficient design leaves a direction the data say nothing about.
  # lm() drops the aliased columns and reports NA for them; here the prior
  # would fill the gap instead, giving coefficients with enormous standard
  # deviations whose sum happens to be right, and a chain that barely moves.
  qx <- qr(X)
  if (qx$rank < k) {
    aliased <- colnames(X)[qx$pivot[(qx$rank + 1L):k]]
    one <- length(aliased) == 1L
    stop("The design matrix has ", k, " columns but rank ", qx$rank, ". ",
         if (one) "Column " else "Columns ",
         paste(aliased, collapse = ", "), if (one) " is " else " are ",
         "a linear combination of the others, so the data cannot tell the ",
         "coefficients apart. Drop ", if (one) "it" else "them",
         " from the formula.")
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

  list(data = list(y = y, X = X, offset = offs, g = g,
                   n_units = length(unique(g)), unit_labels = unit_labels),
       model = list(ineff = ineff, type = type, panel = panel,
                    components = 2L,
                    par_u_name = if (ineff == "exponential") "lambda" else
                      "sigma_u"),
       formula = formula,
       terms = mt,
       # Kept so that predict() codes a factor in newdata against the levels
       # the model was fitted with rather than against whatever happens to
       # appear in the new data.
       xlevels = stats::.getXlevels(mt, mf),
       na.action = dropped,
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
  cat("Observations:       ", x$n,
      if (is.null(x$na.action)) "" else
        paste0(" (", length(x$na.action), " dropped)"), "\n", sep = "")
  if (!is.null(x$data$offset)) {
    cat("Offset:             subtracted from the response\n")
  }
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
