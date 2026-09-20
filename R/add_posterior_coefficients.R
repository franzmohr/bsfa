#' Posterior simulation of stochastic frontier coefficients
#'
#' Runs the Gibbs sampler on a model prepared by
#' \code{\link{create_sfmodel_exp}} or \code{\link{create_sfmodel_hn}} and adds
#' the draws to it.
#'
#' Each sweep draws the inefficiency terms from their truncated normal full
#' conditionals, the frontier coefficients from a normal conditioning on
#' \eqn{y + u} rather than \eqn{y}, and the variance parameters from gamma
#' distributions. This is the data augmentation scheme of van den Broeck, Koop,
#' Osiewalski and Steel (1994).
#'
#' Treating \eqn{u} as a latent variable rather than integrating it out means
#' the efficiency scores come out of the sampler as draws, so
#' \code{\link{efficiency}} reports genuine posterior distributions. Maximum
#' likelihood has to fall back on the Jondrow et al. (1982) conditional mean, a
#' point predictor of an unobserved quantity with no comparable measure of
#' uncertainty.
#'
#' The sampler draws with the seed in \code{object$model$seed}, which
#' \code{\link{add_initial_values}} sets and \code{\link{add_seed}} replaces.
#' R's random number generator is set to that seed for the simulation and put
#' back as it was afterwards, so a call of \code{set.seed()} between
#' \code{add_initial_values()} and this function does not change the draws. A
#' model without a seed draws from R's generator as it stands.
#'
#' That reproducibility holds for a given linear algebra library and thread
#' count. The inefficiency terms are drawn by rejection, so a last-bit
#' difference in the BLAS or LAPACK result that a sweep is built from can
#' change how many deviates it consumes and shift the rest of the chain. The
#' draws remain draws from the same posterior; they are simply not the same
#' draws. See \code{\link{add_seed}}.
#'
#' Priors must have been added. Initial values are added with the default
#' method if they are missing, since least squares starting values are a safe
#' choice where no prior is.
#'
#' @param object an object of class \code{"sfmodel_exp"} or
#'   \code{"sfmodel_hn"}, usually the result of a call to
#'   \code{\link{create_sfmodel_exp}} in combination with
#'   \code{\link{add_priors}} and \code{\link{add_initial_values}}.
#' @param posterior_function the function to be applied to the model in
#'   argument \code{object}. If \code{NULL}, the package's own sampler is used.
#'   It is called as \code{posterior_function(object)} and must return the
#'   object with its \code{posterior} element added.
#' @param keep_u whether to store the inefficiency draws, which
#'   \code{\link{efficiency}} needs. The block is one column per unit and one
#'   row per retained draw, so it can be the largest part of the returned
#'   object. A positive whole number stores every that-many-th retained draw
#'   instead, which divides the storage by it and still leaves genuine
#'   posterior draws to summarise. It cannot exceed the number of retained
#'   draws, since nothing would then be stored at all.
#' @param verbose either \code{FALSE}, \code{TRUE} for progress at every ten
#'   per cent of iterations, or an integer reporting interval.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The object in \code{object} with the element \code{posterior} added.
#'   Each of its elements is a list whose element \code{coeffs} holds the draws
#'   after burn-in as an \code{\link[coda]{mcmc}} object with one row per draw
#'   and one column per parameter:
#'   \describe{
#'     \item{\code{beta}}{the frontier coefficients, one column per regressor.}
#'     \item{\code{sigma_v}}{the standard deviation of the symmetric error.}
#'     \item{\code{lambda} or \code{sigma_u}}{the parameter of the inefficiency
#'       distribution, named for the model: the rate for the exponential
#'       specification and the scale for the half-normal one.}
#'     \item{\code{u}}{the inefficiency terms, one column per unit, if
#'       \code{keep_u} is \code{TRUE}.}
#'     \item{\code{inclusion}}{the inclusion indicators, one column per
#'       coefficient under variable selection, if the model was created with
#'       \code{varsel}. Their posterior means are the posterior inclusion
#'       probabilities that \code{\link{summary.sfmodel}} reports. Under
#'       \code{varsel = "bvs"} the draws in \code{beta} are the coefficients
#'       the frontier was built from, so a regressor an indicator switched off
#'       appears there as an exact zero.}
#'   }
#'
#' @references
#' Jondrow, J., Lovell, C. A. K., Materov, I. S., & Schmidt, P. (1982). On the
#' estimation of technical inefficiency in the stochastic frontier production
#' function model. \emph{Journal of Econometrics}, 19(2--3), 233--238.
#'
#' van den Broeck, J., Koop, G., Osiewalski, J., & Steel, M. F. J. (1994).
#' Stochastic frontier models: A Bayesian perspective. \emph{Journal of
#' Econometrics}, 61(2), 273--303.
#'
#' @seealso \code{\link{add_posterior_loglik}}, \code{\link{efficiency}},
#'   \code{\link{summary.sfmodel}}
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
#'
#' model <- create_sfmodel_exp(y ~ x1 + x2, data = d,
#'                             iterations = 500, burnin = 200)
#' model <- add_priors(model)
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#'
#' summary(model)
#'
#' @family posterior simulation
#' @export
add_posterior_coefficients <- function(object, ...) {
  UseMethod("add_posterior_coefficients")
}

#' @rdname add_posterior_coefficients
#' @export
add_posterior_coefficients.sfmodel_exp <- function(object,
                                                   posterior_function = NULL,
                                                   keep_u = TRUE,
                                                   verbose = FALSE, ...) {
  posterior_coefficients_sf(object, posterior_function, keep_u, verbose)
}

#' @rdname add_posterior_coefficients
#' @export
add_posterior_coefficients.sfmodel_tn <- function(object,
                                                  posterior_function = NULL,
                                                  keep_u = TRUE,
                                                  verbose = FALSE, ...) {
  posterior_coefficients_sf(object, posterior_function, keep_u, verbose)
}

#' @rdname add_posterior_coefficients
#' @export
add_posterior_coefficients.sfmodel_hn <- function(object,
                                                  posterior_function = NULL,
                                                  keep_u = TRUE,
                                                  verbose = FALSE, ...) {
  posterior_coefficients_sf(object, posterior_function, keep_u, verbose)
}

#' Run the sampler for a prepared model object
#'
#' @param object a model object with priors attached.
#' @param posterior_function an alternative sampler, or \code{NULL}.
#' @param keep_u whether to store the inefficiency draws.
#' @param verbose progress reporting.
#'
#' @return The model object with \code{posterior} added.
#'
#' @keywords internal
posterior_coefficients_sf <- function(object, posterior_function, keep_u,
                                      verbose) {

  # This allows the method to be used with other compatible classes.
  class_of_object <- class(object)

  if (!is.null(posterior_function)) {
    object <- posterior_function(object)
    class(object) <- class_of_object
    return(object)
  }

  if (is.null(object$priors)) {
    stop("Add priors before simulating the posterior. See ?add_priors.")
  }
  if (is.null(object$initial)) {
    object <- add_initial_values(object)
  }

  n_iter <- object$burnin + object$iterations
  verbose_int <- verbose_interval(verbose, n_iter)

  n_keep <- object$iterations / object$thin
  u_thin <- augmented_thin(keep_u, n_keep)
  n_keep_u <- if (u_thin > 0) n_keep %/% u_thin else 0
  warn_augmented_size(object$data$n_units, n_keep_u)

  sel <- varsel_args(object)
  ds <- determinant_args(object, "scale_u")
  dm <- determinant_args(object, "mean_u")

  out <- .with_model_seed(
    object$model$seed,
    gibbs_sf(y = object$data$y,
             X = object$data$X,
             g = as.integer(object$data$g) - 1L,
             n_units = object$data$n_units,
             b0 = object$priors$b0,
             B0i = object$priors$B0i,
             ssvs_idx = sel$ssvs_idx,
             tau0 = sel$tau0,
             tau1 = sel$tau1,
             prob_prior = sel$prob_prior,
             varsel = sel$varsel,
             a_v = object$priors$shape_v,
             b_v = object$priors$rate_v,
             a_u = object$priors$shape_u,
             b_u = object$priors$rate_u,
             beta_init = object$initial$beta,
             sigma_v2_init = object$initial$sigma_v2,
             par_u_init = object$initial$par_u,
             u_init = object$initial$u,
             Zs = ds$Z,
             Zm = dm$Z,
             g0 = ds$mu,
             G0i = ds$v_i,
             d0 = dm$mu,
             D0i = dm$v_i,
             gamma_init = ds$init,
             delta_init = dm$init,
             ineff = ineff_code(object),
             s = if (object$model$type == "production") -1 else 1,
             draws = as.integer(object$iterations),
             burnin = as.integer(object$burnin),
             thin = as.integer(object$thin),
             u_thin = u_thin,
             verbose = verbose_int))

  colnames(out$beta) <- colnames(object$data$X)

  posterior <- list(
    beta = list(coeffs = .mcmc_draws(object, out$beta)),
    sigma_v = list(coeffs = .mcmc_draws(object, matrix(out$sigma_v,
                                                       dimnames = list(
                                                         NULL, "sigma_v")))))

  par_u <- matrix(out$par_u,
                  dimnames = list(NULL, object$model$par_u_name))
  posterior[[object$model$par_u_name]] <- list(
    coeffs = .mcmc_draws(object, par_u))

  posterior <- add_inclusion_draws(object, posterior, out$inclusion)
  posterior <- add_determinant_draws(object, posterior, out)

  if (!is.null(out$u)) {
    colnames(out$u) <- object$data$unit_labels
    posterior$u <- list(coeffs = mcmc_augmented(object, out$u, u_thin))
  }

  # The acceptance rates of the Metropolis blocks are a property of the run
  # rather than of any parameter, so they sit beside the posterior and not
  # inside it, where the code that walks the blocks would trip over them.
  object$acceptance <- out$acceptance
  object$posterior <- posterior
  class(object) <- class_of_object
  object
}

# Turns a matrix of retained draws into an mcmc object whose iteration index
# counts in the units the model was specified in, so that thinning and burn-in
# survive into anything that reads the draws.
.mcmc_draws <- function(object, x) {
  coda::mcmc(x,
             start = object$burnin + object$thin,
             end = object$burnin + object$iterations,
             thin = object$thin)
}

# Draws a seed for a model that does not have one.
.draw_model_seed <- function() {
  sample.int(.Machine$integer.max, 1L)
}

# Evaluates 'expr' with R's generator set to 'seed', with R's default kinds, and
# puts the generator back afterwards, kinds included. Without a seed 'expr' is
# evaluated with the generator as it stands. 'expr' is a promise and is only
# forced after set.seed().
.with_model_seed <- function(seed, expr) {

  if (is.null(seed)) {
    return(expr)
  }

  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
            add = TRUE)
  }

  set.seed(seed, kind = "default", normal.kind = "default",
           sample.kind = "default")
  expr
}
