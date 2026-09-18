#' Posterior simulation of four-component frontier coefficients
#'
#' Runs the Gibbs sampler on a model prepared by
#' \code{\link{create_sfmodel4_exp}} or \code{\link{create_sfmodel4_hn}} and
#' adds the draws to it.
#'
#' Each sweep draws the transient inefficiency terms, the persistent
#' inefficiency terms and the unit effects from their full conditionals, then
#' the coefficients and the variance parameters. Every one of those conditionals
#' is standard: the one-sided terms are truncated normals, the unit effects are
#' normal, and the variances are gamma. Adding two components to the sampler of
#' \code{\link{add_posterior_coefficients}} therefore changes the loop but not
#' the algorithm.
#'
#' Convergence deserves more attention here than in the two-component models.
#' The unit effect and persistent inefficiency are both constant within a unit
#' and are separated only by the sign restriction, so the two can trade off
#' against each other from sweep to sweep and mix slowly. Run several chains
#' from \code{add_initial_values(method = "prior")} and look at
#' \code{sigma_mu} alongside the persistent inefficiency parameter before
#' trusting the split.
#'
#' @param object an object of class \code{"sfmodel4_exp"} or
#'   \code{"sfmodel4_hn"}.
#' @param posterior_function the function to be applied to the model in
#'   argument \code{object}. If \code{NULL}, the package's own sampler is used.
#' @param keep_u whether to store the draws of the unit effects and of both
#'   one-sided terms, which \code{\link{efficiency}} needs. The transient terms
#'   alone are one column per observation, so on a panel of any size this block
#'   dominates everything else: a hundred thousand observations and twenty
#'   thousand draws come to sixteen gigabytes. A positive whole number stores
#'   every that-many-th retained draw instead, which divides the storage by it
#'   and cannot exceed the number of retained draws. A warning is issued before
#'   anything is allocated if the block would exceed a gigabyte.
#' @param verbose either \code{FALSE}, \code{TRUE} for progress at every ten
#'   per cent of iterations, or an integer reporting interval.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The object in \code{object} with the element \code{posterior} added.
#'   Alongside \code{beta} and \code{sigma_v} it holds \code{sigma_mu}, the two
#'   inefficiency parameters under the names the model gives them, and, if
#'   \code{keep_u} is \code{TRUE}, the draws of the unit effects \code{mu}, the
#'   persistent terms \code{eta} and the transient terms \code{u}.
#'
#' @seealso \code{\link{efficiency}} for the three efficiency concepts the
#'   model distinguishes.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf4(n = 60, n_time = 8, beta = c(1, 0.5, 0.3))
#'
#' model <- create_sfmodel4_exp(y ~ x1 + x2, data = d, id = "id",
#'                              iterations = 500, burnin = 200)
#' model <- add_posterior_coefficients(add_priors(model))
#'
#' summary(model)
#'
#' @family posterior simulation
#' @export
add_posterior_coefficients.sfmodel4_exp <- function(object,
                                                    posterior_function = NULL,
                                                    keep_u = TRUE,
                                                    verbose = FALSE, ...) {
  posterior_coefficients_sf4(object, posterior_function, keep_u, verbose)
}

#' @rdname add_posterior_coefficients.sfmodel4_exp
#' @export
add_posterior_coefficients.sfmodel4_hn <- function(object,
                                                   posterior_function = NULL,
                                                   keep_u = TRUE,
                                                   verbose = FALSE, ...) {
  posterior_coefficients_sf4(object, posterior_function, keep_u, verbose)
}

#' Run the four-component sampler for a prepared model object
#'
#' @param object a model object with priors attached.
#' @param posterior_function an alternative sampler, or \code{NULL}.
#' @param keep_u whether to store the augmented terms.
#' @param verbose progress reporting.
#'
#' @return The model object with \code{posterior} added.
#'
#' @keywords internal
posterior_coefficients_sf4 <- function(object, posterior_function, keep_u,
                                       verbose) {

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
  verbose_int <- if (isTRUE(verbose)) max(1L, floor(n_iter / 10)) else
    if (isFALSE(verbose)) 0L else as.integer(verbose)

  n_keep <- object$iterations / object$thin
  u_thin <- augmented_thin(keep_u, n_keep)
  n_keep_u <- if (u_thin > 0) n_keep %/% u_thin else 0
  warn_augmented_size(2 * object$data$n_units + object$n, n_keep_u)

  out <- .with_model_seed(
    object$model$seed,
    gibbs_sf4(y = object$data$y,
              X = object$data$X,
              g = as.integer(object$data$g) - 1L,
              n_units = object$data$n_units,
              b0 = object$priors$b0,
              B0i = object$priors$B0i,
              a_v = object$priors$shape_v,
              b_v = object$priors$rate_v,
              a_mu = object$priors$shape_mu,
              b_mu = object$priors$rate_mu,
              a_eta = object$priors$shape_eta,
              b_eta = object$priors$rate_eta,
              a_u = object$priors$shape_u,
              b_u = object$priors$rate_u,
              beta_init = object$initial$beta,
              sigma_v2_init = object$initial$sigma_v2,
              sigma_mu2_init = object$initial$sigma_mu2,
              par_eta_init = object$initial$par_eta,
              par_u_init = object$initial$par_u,
              mu_init = object$initial$mu,
              eta_init = object$initial$eta,
              u_init = object$initial$u,
              ineff = if (object$model$ineff == "halfnormal") 0L else 1L,
              s = if (object$model$type == "production") -1 else 1,
              draws = as.integer(object$iterations),
              burnin = as.integer(object$burnin),
              thin = as.integer(object$thin),
              u_thin = u_thin,
              verbose = verbose_int))

  colnames(out$beta) <- colnames(object$data$X)

  named <- function(x, nm) {
    matrix(x, dimnames = list(NULL, nm))
  }

  posterior <- list(
    beta = list(coeffs = .mcmc_draws(object, out$beta)),
    sigma_v = list(coeffs = .mcmc_draws(object, named(out$sigma_v,
                                                      "sigma_v"))),
    sigma_mu = list(coeffs = .mcmc_draws(object, named(out$sigma_mu,
                                                       "sigma_mu"))))

  posterior[[object$model$par_eta_name]] <- list(
    coeffs = .mcmc_draws(object, named(out$par_eta,
                                       object$model$par_eta_name)))
  posterior[[object$model$par_u_name]] <- list(
    coeffs = .mcmc_draws(object, named(out$par_u, object$model$par_u_name)))

  if (!is.null(out$u)) {
    colnames(out$mu) <- object$data$unit_labels
    colnames(out$eta) <- object$data$unit_labels
    colnames(out$u) <- rownames(object$data$X)
    posterior$mu <- list(coeffs = mcmc_augmented(object, out$mu, u_thin))
    posterior$eta <- list(coeffs = mcmc_augmented(object, out$eta, u_thin))
    posterior$u <- list(coeffs = mcmc_augmented(object, out$u, u_thin))
  }

  object$posterior <- posterior
  class(object) <- class_of_object
  object
}

#' @rdname add_posterior_loglik
#' @export
add_posterior_loglik.sfmodel4 <- function(object, ...) {
  stop("The pointwise log-likelihood is not available for the four-component ",
       "model. Integrating out the unit effect and the persistent ",
       "inefficiency term couples the observations of a unit, so the ",
       "likelihood does not factorise over them. The closed skew normal form ",
       "of Colombi et al. (2014) exists but needs normal distribution ",
       "functions of dimension T_i + 1, which this package does not ",
       "implement.")
}
