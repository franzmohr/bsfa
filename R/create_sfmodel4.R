#' Specify a four-component stochastic frontier model
#'
#' Sets up the model object that the \code{add_*} functions operate on. The two
#' functions differ only in the distribution assumed for the one-sided terms,
#' which determines the class of the object and hence the methods that apply to
#' it.
#'
#' The model is
#' \deqn{y_{it} = x_{it}' \beta + \mu_i - \eta_i - u_{it} + v_{it},
#'       \quad v_{it} \sim N(0, \sigma_v^2),}
#' for a production frontier, with the signs of \eqn{\eta} and \eqn{u} reversed
#' for a cost frontier. Its four components are
#' \describe{
#'   \item{\eqn{\mu_i \sim N(0, \sigma_\mu^2)}}{an unrestricted unit effect,
#'     which is latent heterogeneity and \emph{not} inefficiency: technology,
#'     market or business model, anything that makes a unit persistently
#'     different for reasons unrelated to how well it operates.}
#'   \item{\eqn{\eta_i \ge 0}}{persistent inefficiency, constant over a unit's
#'     observations.}
#'   \item{\eqn{v_{it}}}{symmetric noise.}
#'   \item{\eqn{u_{it} \ge 0}}{transient inefficiency, varying within a unit.}
#' }
#'
#' As in the two-component models, the response is assumed to be on a
#' logarithmic scale, which is what makes \eqn{\exp(-\eta_i)} and
#' \eqn{\exp(-u_{it})} efficiencies rather than quantities in the units of the
#' output. See \code{\link{create_sfmodel_exp}}.
#'
#' This is the specification of Kumbhakar, Lien and Hjalmarsson (2014), Colombi
#' et al. (2014) and Tsionas and Kumbhakar (2014). The reason to prefer it over
#' the time-invariant panel model of \code{\link{create_sfmodel_exp}} is that
#' the latter has no \eqn{\mu_i}, so every persistent difference between units
#' is booked as inefficiency; and no \eqn{u_{it}}, so nothing within a unit can
#' vary. Separating the two matters whenever the question is about movement in
#' efficiency over time rather than about a ranking of units.
#'
#' The price is identification. \eqn{\mu_i} and \eqn{\eta_i} are both
#' unit-specific and are told apart only by the sign restriction on \eqn{\eta}
#' and by the distributional assumptions, while the mean of \eqn{\eta} is in
#' addition hard to separate from the intercept of the frontier.
#'
#' What the data determine sharply is the combination \eqn{\mu_i - \eta_i}.
#' How that combination is divided between the two depends on which of them has
#' the wider spread: where the unit effect dominates, it is recovered well and
#' the persistent inefficiency term is not, and where it is small the reverse
#' holds. The transient term is unaffected, since it is identified within a
#' unit. This is a property of the model rather than of the sampler, and it is
#' worth a prior sensitivity check before reading much into the split. The
#' combination itself, and therefore the fitted frontier, is reliable either
#' way.
#'
#' Both one-sided terms follow the same family, half-normal or exponential,
#' with separate parameters.
#'
#' @param formula a model formula describing the frontier.
#' @param data a data frame containing the variables in \code{formula}.
#' @param id either the name of a variable in \code{data} or a vector of the
#'   same length as the data, identifying the units. Required, since the model
#'   cannot separate its components without repeated observations.
#' @param type either \code{"production"} or \code{"cost"}.
#' @param iterations number of iterations retained after burn-in, before
#'   thinning.
#' @param burnin number of discarded iterations.
#' @param thin thinning interval. \code{iterations} must be a multiple of it.
#'
#' @return An object of class \code{"sfmodel4_exp"} or \code{"sfmodel4_hn"},
#'   both inheriting from \code{"sfmodel4"} and \code{"sfmodel"}.
#'
#' @seealso \code{\link{create_sfmodel_exp}} for the two-component models,
#'   \code{\link{add_priors}}, \code{\link{efficiency}}, \code{\link{sim_sf4}}
#'
#' @references
#' Colombi, R., Kumbhakar, S. C., Martini, G., & Vittadini, G. (2014). Closed
#' skew normal distribution and efficiency analysis. \emph{Journal of
#' Productivity Analysis}, 42(2), 123--136.
#'
#' Kumbhakar, S. C., Lien, G., & Hjalmarsson, L. (2014). Technical efficiency
#' in competing panel data models: A study of Norwegian grain farming.
#' \emph{Journal of Productivity Analysis}, 41(2), 321--337.
#'
#' Tsionas, E. G., & Kumbhakar, S. C. (2014). Firm heterogeneity, persistent
#' and transient technical inefficiency: A generalized true random-effects
#' model. \emph{Journal of Applied Econometrics}, 29(1), 110--132.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf4(n = 60, n_time = 8, beta = c(1, 0.5, 0.3))
#'
#' model <- create_sfmodel4_exp(y ~ x1 + x2, data = d, id = "id",
#'                              iterations = 500, burnin = 200)
#' model
#'
#' @export
create_sfmodel4_exp <- function(formula,
                                data,
                                id,
                                type = c("production", "cost"),
                                iterations = 20000,
                                burnin = 2000,
                                thin = 1) {

  object <- sfmodel4_skeleton(formula = formula, data = data, id = id,
                              type = match.arg(type), iterations = iterations,
                              burnin = burnin, thin = thin,
                              ineff = "exponential", cl = match.call())
  class(object) <- c("sfmodel4_exp", "sfmodel4", "sfmodel")
  object
}

#' @rdname create_sfmodel4_exp
#' @export
create_sfmodel4_hn <- function(formula,
                               data,
                               id,
                               type = c("production", "cost"),
                               iterations = 20000,
                               burnin = 2000,
                               thin = 1) {

  object <- sfmodel4_skeleton(formula = formula, data = data, id = id,
                              type = match.arg(type), iterations = iterations,
                              burnin = burnin, thin = thin,
                              ineff = "halfnormal", cl = match.call())
  class(object) <- c("sfmodel4_hn", "sfmodel4", "sfmodel")
  object
}

#' Assemble the common parts of a four-component model object
#'
#' @param formula a model formula.
#' @param data a data frame.
#' @param id unit identifier.
#' @param type \code{"production"} or \code{"cost"}.
#' @param iterations,burnin,thin MCMC settings.
#' @param ineff \code{"exponential"} or \code{"halfnormal"}.
#' @param cl the originating call.
#'
#' @return A list with the model data and specification.
#'
#' @keywords internal
sfmodel4_skeleton <- function(formula, data, id, type, iterations, burnin,
                              thin, ineff, cl) {

  if (missing(id) || is.null(id)) {
    stop("The four-component model needs an 'id': without repeated ",
         "observations per unit its components cannot be separated.")
  }

  object <- sfmodel_skeleton(formula = formula, data = data, id = id,
                             type = type, iterations = iterations,
                             burnin = burnin, thin = thin, ineff = ineff,
                             cl = cl)

  sizes <- tabulate(object$data$g, nbins = object$data$n_units)
  if (min(sizes) < 2) {
    short <- object$data$unit_labels[sizes < 2]
    stop("Every unit needs at least two observations. Units with a single ",
         "observation cannot contribute to the split between persistent and ",
         "transient inefficiency. ", length(short), " of ", length(sizes),
         " have one: ",
         paste(short[seq_len(min(5L, length(short)))], collapse = ", "),
         if (length(short) > 5) ", ..." else "",
         if (is.null(object$na.action)) "" else
           paste0(" (", length(object$na.action),
                  " rows were dropped for missing values first)"), ".")
  }

  object$model$components <- 4L
  object$model$par_eta_name <- if (ineff == "exponential") "lambda_eta" else
    "sigma_eta"
  object$model$par_u_name <- if (ineff == "exponential") "lambda_u" else
    "sigma_u"

  object
}
