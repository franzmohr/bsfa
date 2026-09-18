#' Check that a model carries the posterior draws a function needs
#'
#' The \code{posterior_function} argument of
#' \code{\link{add_posterior_coefficients}} lets an outside sampler supply the
#' draws, so the contents of \code{object$posterior} are not guaranteed to be
#' what the package's own sampler would have written. Without this check the
#' failure surfaces much later, as a non-conformable matrix or a \code{NULL}
#' passed to a plotting function, neither of which says what is wrong.
#'
#' @param object a model object.
#' @param blocks names of the required elements of \code{object$posterior},
#'   each expected to hold a matrix of draws in its \code{coeffs} element.
#'
#' @return Invisibly \code{TRUE}; called for the error it raises.
#'
#' @keywords internal
check_posterior_blocks <- function(object, blocks) {

  if (is.null(object$posterior)) {
    stop("Object does not contain posterior draws. ",
         "See ?add_posterior_coefficients.")
  }

  usable <- vapply(blocks, function(b) {
    x <- object$posterior[[b]]$coeffs
    !is.null(x) && is.numeric(as.matrix(x)) && NROW(x) > 0 && NCOL(x) > 0
  }, logical(1))

  if (!all(usable)) {
    stop("The posterior of this model has no usable draws in ",
         paste0("posterior$", blocks[!usable], "$coeffs", collapse = ", "),
         ". A posterior supplied through 'posterior_function' has to carry ",
         "the same blocks as the package's own sampler.")
  }

  n_draws <- vapply(blocks, function(b) NROW(object$posterior[[b]]$coeffs),
                    numeric(1))
  if (length(unique(n_draws)) > 1L) {
    stop("The posterior blocks hold different numbers of draws: ",
         paste(paste0(blocks, " (", n_draws, ")"), collapse = ", "), ".")
  }

  if ("beta" %in% blocks &&
      NCOL(object$posterior$beta$coeffs) != object$k) {
    stop("posterior$beta$coeffs has ",
         NCOL(object$posterior$beta$coeffs), " column(s) but the model has ",
         object$k, " coefficient(s).")
  }

  invisible(TRUE)
}

#' Drop posterior draws that no longer match the specification
#'
#' Changing a prior or the starting values after the sampler has run leaves
#' draws behind that were produced under the old specification. Keeping them
#' would let \code{summary()} report draws alongside a prior that did not
#' generate them, so they are discarded and the user is told.
#'
#' @param object a model object.
#'
#' @return The model object without its \code{posterior} element.
#'
#' @keywords internal
drop_stale_posterior <- function(object) {
  if (!is.null(object$posterior)) {
    message("Dropping posterior draws: they were simulated under a different ",
            "specification. Re-run add_posterior_coefficients().")
    object$posterior <- NULL
  }
  object
}

#' Names of the scalar parameter blocks of a model
#'
#' @param object a model object.
#'
#' @return A character vector, excluding \code{beta}.
#'
#' @keywords internal
sf_scalar_blocks <- function(object) {
  if (identical(object$model$components, 4L)) {
    c("sigma_v", "sigma_mu", object$model$par_eta_name,
      object$model$par_u_name)
  } else {
    c("sigma_v", object$model$par_u_name)
  }
}

#' Summarise a matrix of efficiency draws
#'
#' Both \code{efficiency} methods report the same statistics and differ only in
#' which draws they hand over and how the rows are labelled, so the table is
#' built in one place.
#'
#' The quantiles are given their shape explicitly rather than by transposing
#' the result of \code{apply}. For more than one probability \code{apply}
#' returns a matrix with one row per probability, but for a single one it
#' returns a plain vector, and transposing that yields a one-row matrix whose
#' columns are the units.
#'
#' @param r a matrix of efficiency draws, one row per draw and one column per
#'   unit or observation.
#' @param unit the labels of the columns of \code{r}.
#' @param probs quantiles of the posterior to report.
#'
#' @return A data frame with one row per column of \code{r}.
#'
#' @keywords internal
efficiency_table <- function(r, unit, probs) {

  qs <- matrix(apply(r, 2, stats::quantile, probs = probs),
               nrow = ncol(r), byrow = TRUE,
               dimnames = list(NULL,
                               paste0(format(100 * probs, trim = TRUE), "%")))

  data.frame(unit = unit,
             mean = colMeans(r),
             sd = apply(r, 2, stats::sd),
             qs,
             row.names = NULL,
             check.names = FALSE,
             stringsAsFactors = FALSE)
}

#' Check that a value is a single number strictly between zero and one
#'
#' Used for the prior median efficiency and for the credible band of every
#' function that reports one. A copy of the test per call site is how one of
#' them came to accept \code{NA}: a comparison with \code{NA} is \code{NA}
#' rather than \code{FALSE}, so a guard written as \code{x <= 0 || x >= 1}
#' raises R's own "missing value where TRUE/FALSE needed" instead of naming the
#' argument at fault.
#'
#' @param x the value to check.
#' @param what the argument name, used in the error message.
#'
#' @return Invisibly \code{TRUE}; called for the error it raises.
#'
#' @keywords internal
check_probability <- function(x, what) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x <= 0 || x >= 1) {
    stop("'", what, "' must be a single number strictly between 0 and 1.")
  }
  invisible(TRUE)
}

#' Check that a value is a single whole number an integer can hold
#'
#' The MCMC settings reach the sampler through \code{as.integer()}, which
#' truncates a fractional value and turns anything beyond the integer range
#' into \code{NA}, neither of them with an error. They are therefore checked
#' here rather than left to fail somewhere inside C++.
#'
#' @param x the value to check.
#' @param what the argument name, used in the error message.
#'
#' @return Invisibly \code{TRUE}; called for the error it raises.
#'
#' @keywords internal
check_count <- function(x, what) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x != trunc(x) || abs(x) > .Machine$integer.max) {
    stop("'", what, "' must be a single whole number no larger than ",
         .Machine$integer.max, ".")
  }
  invisible(TRUE)
}
