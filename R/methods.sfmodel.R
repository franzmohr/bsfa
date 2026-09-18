#' Standard accessors for a stochastic frontier model
#'
#' The methods that R's modelling conventions lead a user to reach for. Without
#' them \code{coef()} falls through to \code{coef.default()}, which looks for an
#' element named \code{coefficients}, finds none, and returns \code{NULL}
#' rather than saying so.
#'
#' A frontier model has two natural notions of a fitted value, and they are
#' kept apart here. \code{fitted()} returns the frontier, the maximum output
#' the inputs allow, which is what the coefficients describe, including any
#' offset. \code{residuals()} returns the composed error, which still contains
#' the one-sided term and is therefore not centred on zero. Use
#' \code{\link{efficiency}} for the inefficiency part of it. The two are
#' related as usual: the residual is the response as supplied minus the fitted
#' frontier, whether or not there is an offset.
#'
#' Every one of these summarises the posterior by its mean. That is a
#' convenience, not the object: the draws themselves are in
#' \code{object$posterior}, and \code{\link{summary.sfmodel}} reports their
#' spread.
#'
#' @param object an object of class \code{"sfmodel"} carrying posterior draws.
#' @param newdata an optional data frame in which to look for the variables of
#'   the frontier. If omitted, the data the model was built on are used.
#' @param ... further arguments, unused.
#'
#' @return \code{coef()} a named vector of posterior means of the frontier
#'   coefficients; \code{vcov()} their posterior covariance matrix;
#'   \code{nobs()} the number of observations; \code{fitted()} and
#'   \code{predict()} the frontier; \code{residuals()} the composed error; and
#'   \code{logLik()} the log-likelihood at the posterior mean, with the
#'   attributes that \code{\link[stats]{AIC}} and \code{\link[stats]{BIC}} need.
#'
#' @seealso \code{\link{summary.sfmodel}} for the posterior summaries,
#'   \code{\link{efficiency}} for the scores, \code{\link{selection_criteria}}
#'   for the criteria in the form the package reports them.
#'
#' @examples
#' set.seed(1234)
#' d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
#'
#' model <- create_sfmodel_exp(y ~ x1 + x2, data = d,
#'                             iterations = 500, burnin = 200)
#' model <- add_posterior_coefficients(add_priors(model))
#'
#' coef(model)
#' nobs(model)
#' head(residuals(model))
#' AIC(model)
#'
#' @name sfmodel-methods
NULL

#' @rdname sfmodel-methods
#' @export
coef.sfmodel <- function(object, ...) {
  check_posterior_blocks(object, "beta")
  colMeans(as.matrix(object$posterior$beta$coeffs))
}

#' @rdname sfmodel-methods
#' @export
vcov.sfmodel <- function(object, ...) {
  check_posterior_blocks(object, "beta")
  stats::cov(as.matrix(object$posterior$beta$coeffs))
}

#' @rdname sfmodel-methods
#' @export
nobs.sfmodel <- function(object, ...) {
  object$n
}

#' The estimated part of the frontier, without any offset
#'
#' @param object a model object carrying posterior draws.
#'
#' @return A numeric vector with one entry per observation.
#'
#' @keywords internal
sf_frontier <- function(object) {
  as.numeric(object$data$X %*% coef.sfmodel(object))
}

#' @rdname sfmodel-methods
#' @export
fitted.sfmodel <- function(object, ...) {
  out <- sf_frontier(object)
  if (!is.null(object$data$offset)) {
    out <- out + object$data$offset
  }
  names(out) <- rownames(object$data$X)
  out
}

#' @rdname sfmodel-methods
#' @export
residuals.sfmodel <- function(object, ...) {
  # Against the frontier without the offset, and the response the sampler saw,
  # which is the one the offset has already been taken out of. Subtracting
  # fitted() from the original response comes to the same thing.
  out <- object$data$y - sf_frontier(object)
  names(out) <- rownames(object$data$X)
  out
}

#' @rdname sfmodel-methods
#' @export
predict.sfmodel <- function(object, newdata = NULL, ...) {

  if (is.null(newdata)) {
    return(fitted.sfmodel(object))
  }

  # The terms object is stripped of the response, so that newdata need not
  # carry the variable being predicted.
  mt <- stats::delete.response(object$terms)
  mf <- stats::model.frame(mt, data = newdata, xlev = object$xlevels)
  X <- stats::model.matrix(mt, mf)

  if (ncol(X) != object$k) {
    stop("'newdata' gives ", ncol(X), " coefficient(s) but the model has ",
         object$k, ".")
  }

  out <- as.numeric(X %*% coef.sfmodel(object))
  offs <- stats::model.offset(mf)
  if (!is.null(offs)) {
    out <- out + as.numeric(offs)
  }
  names(out) <- rownames(X)
  out
}

#' @rdname sfmodel-methods
#' @export
logLik.sfmodel <- function(object, ...) {

  if (object$model$panel || identical(object$model$components, 4L)) {
    stop("The likelihood of this model does not factorise over the ",
         "observations of a unit, so there is no log-likelihood to return. ",
         "See ?add_posterior_loglik.")
  }
  check_posterior_blocks(object, c("beta", sf_scalar_blocks(object)))

  # Evaluated at the posterior mean, which is the point the degrees of freedom
  # of AIC and BIC belong to, and the same point selection_criteria() uses.
  value <- sum(sf_loglik_point(
    object,
    coef.sfmodel(object),
    mean(object$posterior$sigma_v$coeffs),
    mean(object$posterior[[object$model$par_u_name]]$coeffs)))

  structure(value, df = object$k + 2L, nobs = object$n, class = "logLik")
}
