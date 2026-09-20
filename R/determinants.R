#' Build the design matrix of a set of inefficiency determinants
#'
#' The determinants of an inefficiency term are given as a one-sided formula,
#' evaluated against the same data as the frontier. Two things distinguish
#' this from building the frontier design.
#'
#' The first is the intercept. A set of determinants of a \emph{scale} carries
#' none, because the model already has a scale parameter for the term to be
#' measured against: \eqn{\sigma_u} or \eqn{\lambda} is the value the scale
#' takes where every determinant is zero, and a column of ones beside it would
#' be the same quantity written twice, leaving the pair identified only by
#' their priors. A set of determinants of a \emph{mean} keeps its intercept,
#' since the pre-truncation mean of the truncated normal has no parameter of
#' its own: with \code{~ 1} it is the constant mean of Stevenson (1980), and
#' with covariates the specification of Battese and Coelli (1995).
#'
#' The second is the level. An inefficiency term that belongs to a unit can
#' only be explained by something that is a property of that unit, so a
#' determinant of such a term is required to be constant within a unit and is
#' reduced to one row per unit here. A determinant of a term that belongs to
#' an observation is taken as it comes.
#'
#' @param spec a one-sided formula, or \code{NULL} for no determinants.
#' @param data the data frame the frontier was built from.
#' @param what the argument name, used in error messages.
#' @param keep a logical vector over the rows of \code{data}, marking those
#'   that survived the frontier's handling of missing values.
#' @param g the unit index of each surviving row.
#' @param n_units the number of units.
#' @param level \code{"unit"} if the term is one per unit, \code{"obs"} if it
#'   is one per observation.
#' @param intercept whether to keep the intercept.
#'
#' @return \code{NULL}, or a list with the matrix \code{Z}, its column
#'   \code{names}, the \code{formula} as given, and the \code{xlevels} of any
#'   factor among the determinants.
#'
#' @references
#' Battese, G. E., & Coelli, T. J. (1995). A model for technical inefficiency
#' effects in a stochastic frontier production function for panel data.
#' \emph{Empirical Economics}, 20(2), 325--332.
#'
#' Stevenson, R. E. (1980). Likelihood functions for generalized stochastic
#' frontier estimation. \emph{Journal of Econometrics}, 13(1), 57--66.
#'
#' @keywords internal
determinant_matrix <- function(spec, data, what, keep, g, n_units, level,
                               intercept) {

  if (is.null(spec)) {
    return(NULL)
  }
  if (!inherits(spec, "formula")) {
    stop("'", what, "' must be a one-sided formula, such as ", what,
         " = ~ z1 + z2, or NULL for a term whose scale does not depend on ",
         "anything.")
  }
  if (length(spec) != 2L) {
    stop("'", what, "' must be one-sided. It names the variables that ",
         "explain the inefficiency term, and the term itself is latent, so ",
         "there is nothing to put on the left of the tilde.")
  }

  mf <- stats::model.frame(spec, data = data, na.action = stats::na.pass)
  if (nrow(mf) != length(keep)) {
    stop("'", what, "' produced ", nrow(mf), " rows from 'data' but the ",
         "frontier was built from ", length(keep), ". Both are evaluated ",
         "against the same data, so they have to agree.")
  }
  mf <- mf[keep, , drop = FALSE]

  # Missing values are refused rather than dropped. Dropping them here would
  # take out rows the frontier has already kept, and the two designs would no
  # longer describe the same observations.
  if (!all(stats::complete.cases(mf))) {
    stop("'", what, "' has missing values in rows the frontier kept. A ",
         "determinant has to be observed wherever the response is; drop ",
         "those rows from 'data' or leave the variable out.")
  }

  mt <- attr(mf, "terms")
  if (!intercept) {
    attr(mt, "intercept") <- 0L
  }
  Z <- stats::model.matrix(mt, mf)

  if (ncol(Z) == 0L) {
    stop("'", what, "' has no columns. ",
         if (intercept) "Give it at least an intercept, as ~ 1." else
           paste("Its only term was an intercept, which a set of scale",
                 "determinants does not carry, since the scale parameter",
                 "of the distribution already plays that part."))
  }

  # A unit-level term cannot be explained by something that moves within the
  # unit: the term takes one value for the whole unit, so a determinant that
  # differs across its observations does not say which of its values applies.
  if (identical(level, "unit")) {
    varies <- character(0)
    for (j in seq_len(ncol(Z))) {
      rng <- tapply(Z[, j], g, function(v) max(v) - min(v))
      if (any(rng > 0)) {
        varies <- c(varies, colnames(Z)[j])
      }
    }
    if (length(varies) > 0) {
      one <- length(varies) == 1L
      stop("'", what, "' explains a term that is one per unit, but ",
           paste(varies, collapse = ", "), if (one) " varies" else " vary",
           " within a unit. A persistent term takes one value for the whole ",
           "unit, so a determinant of it has to be a property of the unit. ",
           "Use the time average, or put ", if (one) "it" else "them",
           " on the transient term instead.")
    }
    Z <- Z[!duplicated(g), , drop = FALSE][order(unique(g)), , drop = FALSE]
    if (nrow(Z) != n_units) {
      stop("'", what, "' reduced to ", nrow(Z), " rows for ", n_units,
           " units.")
    }
  }

  # The same rank check the frontier gets, and for the same reason: a
  # direction the data say nothing about is filled by the prior instead, and
  # the chain barely moves along it.
  qz <- qr(Z)
  if (qz$rank < ncol(Z)) {
    aliased <- colnames(Z)[qz$pivot[(qz$rank + 1L):ncol(Z)]]
    one <- length(aliased) == 1L
    stop("'", what, "' has ", ncol(Z), " columns but rank ", qz$rank, ". ",
         if (one) "Column " else "Columns ", paste(aliased, collapse = ", "),
         if (one) " is " else " are ",
         "a linear combination of the others. Drop ",
         if (one) "it" else "them", ".")
  }

  list(Z = Z, names = colnames(Z), formula = spec,
       xlevels = stats::.getXlevels(mt, mf))
}


#' Check that a model was created with the determinants a prior is given for
#'
#' The determinants belong to the model and their prior belongs to
#' \code{add_priors}, exactly as the variable selection algorithm and its
#' prior do, so the same pair of mistakes is possible and is reported the same
#' way: a prior for determinants the model does not have, and a model whose
#' determinants have no prior.
#'
#' @param object a model object.
#' @param spec the prior specification supplied, or \code{NULL}.
#' @param slot the name of the determinant set, such as \code{"scale_u"}.
#'
#' @return The merged specification, or \code{NULL} if the model has no such
#'   determinants.
#'
#' @keywords internal
prior_determinants <- function(object, spec, slot) {

  z <- object$model$determinants[[slot]]

  if (is.null(z)) {
    if (!is.null(spec)) {
      stop("'", slot, "' was given, but the model was not created with ",
           "determinants of that term, so there is nothing for it to set a ",
           "prior for. Pass ", slot, " = ~ z1 + z2 to the create_sfmodel ",
           "call instead.")
    }
    return(NULL)
  }

  q <- length(z$names)
  spec <- merge_prior_list(spec, list(mu = 0, v_i = 0.01), slot)

  mu <- spec$mu
  if (length(mu) == 1L) {
    mu <- rep(mu, q)
  }
  if (length(mu) != q || !is.numeric(mu) || anyNA(mu)) {
    stop("'", slot, "$mu' must be one number or one per determinant; the ",
         "model has ", q, " (", paste(z$names, collapse = ", "), ").")
  }

  v_i <- spec$v_i
  if (length(v_i) == 1L) {
    v_i <- diag(v_i, q)
  } else if (is.matrix(v_i)) {
    if (!identical(dim(v_i), c(q, q))) {
      stop("'", slot, "$v_i' must be ", q, " by ", q, ".")
    }
  } else if (length(v_i) == q) {
    v_i <- diag(v_i, q)
  } else {
    stop("'", slot, "$v_i' must be one number, one per determinant, or a ",
         q, " by ", q, " matrix.")
  }
  check_precision(v_i)

  # A determinant coefficient enters through exp(z'gamma), so the prior has to
  # be proper. A flat prior on a coefficient inside an exponential puts most
  # of its mass on scales no data can rule out, and the Metropolis step then
  # wanders off into a region where the inefficiency term is numerically zero
  # or overflows.
  ev <- eigen(v_i, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0) {
    stop("'", slot, "$v_i' must be positive definite. A determinant ",
         "coefficient acts through an exponential, so an improper prior on ",
         "it leaves the scale of the inefficiency term unbounded and the ",
         "sampler stuck where it started.")
  }

  list(mu = mu, v_i = v_i, names = z$names)
}


#' Assemble the determinant sets of a model
#'
#' Called by the constructors once the frontier design is built, so that the
#' determinants are checked against the rows the frontier kept rather than
#' against the data as supplied.
#'
#' @param specs a named list of one-sided formulas or \code{NULL}s.
#' @param data the data frame.
#' @param keep a logical vector marking the surviving rows.
#' @param g the unit index of each surviving row.
#' @param n_units the number of units.
#' @param levels a named character vector, \code{"unit"} or \code{"obs"} per
#'   element of \code{specs}.
#' @param ineff the inefficiency family.
#'
#' @return A named list of determinant sets, with \code{NULL} for those not
#'   given.
#'
#' @keywords internal
build_determinants <- function(specs, data, keep, g, n_units, levels, ineff) {

  # The pre-truncation mean is a feature of the truncated normal and of no
  # other family: a half-normal or an exponential term has no mean to shift,
  # and silently ignoring the argument would leave the user reading a
  # posterior that does not contain what they asked for.
  mean_slots <- grep("^mean_", names(specs), value = TRUE)
  if (!identical(ineff, "truncnormal")) {
    given <- mean_slots[!vapply(specs[mean_slots], is.null, logical(1))]
    if (length(given) > 0) {
      stop("'", given[1], "' shifts the pre-truncation mean of the ",
           "inefficiency term, which only the truncated normal family has. ",
           "This model is ", ineff, ", whose one-sided term is anchored at ",
           "zero. Create the model with create_sfmodel_tn() or ",
           "create_sfmodel4_tn(), or give '",
           sub("^mean_", "scale_", given[1]),
           "' to make the determinants act on the scale instead.")
    }
  }

  out <- list()
  for (nm in names(specs)) {
    out[[nm]] <- determinant_matrix(
      spec = specs[[nm]], data = data, what = nm, keep = keep, g = g,
      n_units = n_units, level = levels[[nm]],
      intercept = startsWith(nm, "mean_"))
  }

  # The truncated normal always has a pre-truncation mean. Where none was
  # asked for it is the constant of Stevenson (1980), which is what makes the
  # family more than a half-normal in the first place.
  if (identical(ineff, "truncnormal")) {
    for (nm in mean_slots) {
      if (is.null(out[[nm]])) {
        lvl <- levels[[nm]]
        rows_n <- if (identical(lvl, "unit")) n_units else sum(keep)
        out[[nm]] <- list(Z = matrix(1, nrow = rows_n, ncol = 1,
                                     dimnames = list(NULL, "(Intercept)")),
                          names = "(Intercept)", formula = ~1,
                          xlevels = list())
      }
    }
  }

  if (all(vapply(out, is.null, logical(1)))) {
    return(NULL)
  }
  out
}


#' The determinant sets a model of a given kind can carry
#'
#' @param components \code{2L} or \code{4L}.
#'
#' @return A named character vector of levels, \code{"unit"} or \code{"obs"}.
#'
#' @keywords internal
determinant_levels <- function(components) {
  if (identical(components, 4L)) {
    # The persistent term is one per unit and the transient one is one per
    # observation, which is the whole point of the four-component model, so
    # their determinants sit at different levels too.
    c(scale_eta = "unit", mean_eta = "unit",
      scale_u = "obs", mean_u = "obs")
  } else {
    # In the two-component model the inefficiency term belongs to the unit,
    # which for a cross-section is the observation.
    c(scale_u = "unit", mean_u = "unit")
  }
}


#' The names of the determinant coefficient blocks a model carries
#'
#' @param object a model object.
#'
#' @return A character vector, possibly empty.
#'
#' @keywords internal
determinant_blocks <- function(object) {
  d <- object$model$determinants
  if (is.null(d)) {
    return(character(0))
  }
  names(d)[!vapply(d, is.null, logical(1))]
}


#' Collect the priors of every determinant set a model carries
#'
#' @param object a model object.
#' @param specs a named list of prior specifications, one per determinant set
#'   the calling \code{add_priors} method offers.
#'
#' @return A list to splice into \code{object$priors}, empty when the model
#'   has no determinants.
#'
#' @keywords internal
prior_determinant_list <- function(object, specs) {

  # A determinant set the method does not offer would otherwise reach the
  # sampler without a prior. This cannot happen through the constructors as
  # they stand, and is checked because the cost of being wrong is a block
  # drawn from whatever the uninitialised prior happened to be.
  unoffered <- setdiff(determinant_blocks(object), names(specs))
  if (length(unoffered) > 0) {
    stop("The model carries determinants in '", unoffered[1], "', for which ",
         "this method sets no prior. This is an internal inconsistency; ",
         "please report it.")
  }

  out <- list()
  for (nm in names(specs)) {
    out[[nm]] <- prior_determinants(object, specs[[nm]], nm)
  }
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0L) {
    return(list())
  }
  list(determinants = out)
}


#' The integer code the samplers take for an inefficiency family
#'
#' @param object a model object.
#'
#' @return \code{0L}, \code{1L} or \code{2L}.
#'
#' @keywords internal
ineff_code <- function(object) {
  switch(object$model$ineff,
         halfnormal = 0L,
         exponential = 1L,
         truncnormal = 2L,
         stop("Unknown inefficiency family '", object$model$ineff, "'."))
}


#' The sampler arguments of one determinant set
#'
#' A set the model does not carry is handed over as a matrix with no columns,
#' which is what the sampler reads as "no determinants": the scale multiplier
#' is then one, the mean shift zero, and no Metropolis block is created.
#'
#' @param object a model object with priors and initial values attached.
#' @param slot the name of the determinant set.
#'
#' @return A list with \code{Z}, \code{mu}, \code{v_i} and \code{init}.
#'
#' @keywords internal
determinant_args <- function(object, slot) {

  z <- object$model$determinants[[slot]]
  if (is.null(z)) {
    return(list(Z = matrix(0, nrow = 0, ncol = 0),
                mu = numeric(0),
                v_i = matrix(0, nrow = 0, ncol = 0),
                init = numeric(0)))
  }
  pr <- object$priors$determinants[[slot]]
  init <- object$initial$determinants[[slot]]
  if (is.null(pr) || is.null(init)) {
    stop("The model carries determinants in '", slot, "' but no ",
         if (is.null(pr)) "prior" else "starting value",
         " for them. Call add_priors() and add_initial_values() again.")
  }
  list(Z = z$Z, mu = pr$mu, v_i = pr$v_i, init = init)
}


#' Attach the draws of the determinant coefficients to a posterior
#'
#' @param object a model object.
#' @param posterior the posterior list built so far.
#' @param out the sampler output.
#'
#' @return The posterior with one block per determinant set the model has.
#'
#' @keywords internal
add_determinant_draws <- function(object, posterior, out) {

  for (slot in determinant_blocks(object)) {
    draws <- out[[slot]]
    if (is.null(draws)) {
      next
    }
    colnames(draws) <- object$model$determinants[[slot]]$names
    posterior[[slot]] <- list(coeffs = .mcmc_draws(object, draws))
  }
  posterior
}
