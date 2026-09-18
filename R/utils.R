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
