#' Add a seed to a stochastic frontier model
#'
#' Records the seed that \code{\link{draw_posterior}} sets before sampling, so
#' that a model object fully determines its own output.
#'
#' Setting the seed here rather than calling \code{set.seed} beforehand keeps
#' the seed with the specification, which matters once several models are
#' estimated in one script: the results then depend on each model's own seed
#' rather than on the order in which they happen to be run.
#'
#' @param object an object of class \code{"sfmodel"}.
#' @param seed a single number passed to \code{\link[base]{set.seed}}.
#' @param ... unused, for compatibility with the generic.
#'
#' @return The model object with the element \code{seed} attached.
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
      is.na(seed)) {
    stop("'seed' must be a single number.")
  }

  object$seed <- seed
  object
}
