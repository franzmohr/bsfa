#' Add a seed to a stochastic frontier model
#'
#' Replaces the seed that \code{\link{add_posterior_coefficients}} sets before
#' sampling, so that a model object fully determines its own output.
#'
#' \code{\link{add_initial_values}} already draws a seed and stores it as
#' \code{object$model$seed}; this function replaces it. Keeping the seed with
#' the specification rather than calling \code{set.seed} beforehand matters
#' once several models are estimated in one script: the results then depend on
#' each model's own seed rather than on the order in which they happen to be
#' run.
#'
#' What the seed guarantees is reproducibility on one machine with one linear
#' algebra library, which is not quite the same as reproducibility everywhere.
#' The samplers draw their one-sided terms by rejection, so the number of
#' random deviates a sweep consumes depends on the values it computes. Those
#' values pass through the BLAS and LAPACK that R is linked against, and a
#' different library, or the same library running on a different number of
#' threads, can return a result that differs in the last bit. That is enough to
#' change an acceptance somewhere and shift every draw that follows.
#'
#' The consequence is narrow but worth knowing. Two runs of the same seed under
#' the same setup are identical; under a different one they are two different
#' chains from the same posterior, so every summary agrees to Monte Carlo error
#' and none of the draws match. Where bit-identical output matters -- a
#' regression test, a published figure -- pin the thread count as well as the
#' seed, for instance with \code{OPENBLAS_NUM_THREADS}, and record which
#' library was used.
#'
#' @param object an object of class \code{"sfmodel"}.
#' @param seed a single number passed to \code{\link[base]{set.seed}}.
#' @param ... unused, for compatibility with the generic.
#'
#' @return The model object with the element \code{seed} of \code{model}
#'   replaced. Any posterior draws it carried are discarded, since they were
#'   simulated under the seed that has just been replaced.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
#'
#' model <- create_sfmodel_exp(y ~ x1 + x2, data = d,
#'                             iterations = 500, burnin = 200)
#' model <- add_seed(model, 5555)
#'
#' @export
add_seed <- function(object, ...) {
  UseMethod("add_seed")
}

#' @rdname add_seed
#' @export
add_seed.sfmodel <- function(object, seed, ...) {

  if (missing(seed) || length(seed) != 1L || !is.numeric(seed) ||
      !is.finite(seed)) {
    stop("'seed' must be a single number.")
  }
  # set.seed() takes the seed as an integer, so a fractional value would be
  # truncated and one beyond the integer range would become NA, in both cases
  # leaving the model carrying a seed that is not the one it was given.
  if (seed != trunc(seed) || abs(seed) > .Machine$integer.max) {
    stop("'seed' must be a whole number within the range of an integer, ",
         "since set.seed() reads it as one.")
  }

  object$model$seed <- seed
  drop_stale_posterior(object)
}
