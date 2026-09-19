#' Check the variable selection argument of a model constructor
#'
#' @param varsel the value supplied, \code{NULL}, \code{"ssvs"} or
#'   \code{"bvs"}.
#'
#' @return \code{NULL}, \code{"ssvs"} or \code{"bvs"}.
#'
#' @keywords internal
check_varsel <- function(varsel) {

  if (is.null(varsel)) {
    return(NULL)
  }
  if (length(varsel) != 1L || !is.character(varsel) || is.na(varsel) ||
      !varsel %in% c("ssvs", "bvs")) {
    stop("'varsel' must be NULL for no variable selection, \"ssvs\" for ",
         "stochastic search variable selection, or \"bvs\" for the Bayesian ",
         "variable selection of Korobilis (2013).")
  }
  varsel
}

#' The integer code the samplers take for a variable selection algorithm
#'
#' @param varsel \code{NULL}, \code{"ssvs"} or \code{"bvs"}.
#'
#' @return \code{0L}, \code{1L} or \code{2L}.
#'
#' @keywords internal
varsel_code <- function(varsel) {
  if (is.null(varsel)) 0L else if (identical(varsel, "ssvs")) 1L else 2L
}

#' Build the prior of the stochastic search variable selection
#'
#' Turns the \code{varsel} argument of \code{\link{add_priors}} into the
#' vectors the sampler draws the inclusion indicators with: one prior standard
#' deviation for a coefficient that is excluded from the frontier, one for a
#' coefficient that is included, and one prior inclusion probability, each per
#' coefficient under selection.
#'
#' Under \code{"ssvs"} the mixture is the one of George, Sun and Ni (2008),
#' which is also what \pkg{bvartools} uses. It is a prior on the coefficient
#' rather than a restriction on it: an excluded coefficient is not set to zero
#' but given a prior so tight around zero that it cannot move away from it, so
#' the sampler stays a Gibbs sampler and the draws of every other block are
#' unaffected.
#'
#' Under \code{"bvs"} the selection of Korobilis (2013) switches the regressor
#' itself off, so there is no mixture and no pair of standard deviations: the
#' coefficient keeps the normal prior in \code{coef}, and the only thing to
#' elicit is the prior inclusion probability. That prior has to be proper,
#' since an excluded coefficient is drawn from it.
#'
#' @param object a model object.
#' @param varsel the \code{varsel} argument of \code{\link{add_priors}}.
#' @param b0 the prior mean of the coefficients, already expanded.
#' @param B0i the prior precision of the coefficients, already expanded.
#'
#' @return \code{NULL} if the model has no variable selection, otherwise a list
#'   with \code{include}, the positions under selection, their \code{names},
#'   \code{tau0}, \code{tau1}, \code{inprior}, and \code{spec}, the merged
#'   specification as given.
#'
#' @references
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR
#' model restrictions. \emph{Journal of Econometrics}, 142(1), 553--580.
#'
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics}, 28(2), 204--230.
#'
#' @keywords internal
prior_varsel <- function(object, varsel, b0, B0i) {

  alg <- object$model$varsel

  if (is.null(alg)) {
    if (!is.null(varsel)) {
      stop("'varsel' was given, but the model was not created with a ",
           "variable selection algorithm, so there is nothing for it to set ",
           "a prior for. Pass varsel = \"ssvs\" or varsel = \"bvs\" to the ",
           "create_sfmodel call instead.")
    }
    return(NULL)
  }
  if (is.null(varsel)) {
    stop("The model was created with varsel = \"", alg, "\", so 'varsel' ",
         "must be given. It needs ",
         if (alg == "ssvs")
           paste("either 'tau', the prior standard deviations of an excluded",
                 "and an included coefficient, or 'semiautomatic', the",
                 "factors to scale their least squares standard errors by")
         else "'inprior', the prior probability that a regressor is included",
         ". See ?add_priors.")
  }

  nms <- colnames(object$data$X)
  ssvs <- identical(alg, "ssvs")
  defaults <- list(inprior = 0.5, include = NULL, exclude_intercept = TRUE)
  if (ssvs) {
    defaults <- c(defaults, list(tau = NULL, semiautomatic = NULL))
  } else {
    # Reported here rather than left to the unknown-element message, which
    # would say that 'tau' is not an element of 'varsel' without saying that
    # it is an element of the other algorithm's.
    mixture <- intersect(names(varsel), c("tau", "semiautomatic"))
    if (length(mixture) > 0) {
      stop("'varsel$", mixture[1], "' belongs to SSVS, which sets a ",
           "coefficient's prior by its inclusion indicator. The Bayesian ",
           "variable selection of Korobilis (2013) switches the regressor ",
           "off instead and leaves the coefficient the normal prior in ",
           "'coef', so it has no pair of standard deviations to set. Give ",
           "'varsel$inprior', or create the model with varsel = \"ssvs\".")
    }
  }
  spec <- merge_prior_list(varsel, defaults, "varsel")

  idx <- varsel_include(spec, nms)
  inprior <- varsel_inprior(spec$inprior, length(idx))
  taus <- if (ssvs) varsel_tau(spec, object, idx) else
    list(tau0 = numeric(0), tau1 = numeric(0))

  if (ssvs) {
    varsel_check_ssvs_prior(idx, nms, b0, B0i)
  } else {
    varsel_check_bvs_prior(idx, nms, B0i, ols_se(object))
  }

  list(algorithm = alg, include = idx, names = nms[idx], tau0 = taus$tau0,
       tau1 = taus$tau1, inprior = inprior, spec = spec)
}

#' Check the coefficient prior that stochastic search variable selection needs
#'
#' @param idx the positions under selection.
#' @param nms the column names of the design matrix.
#' @param b0 the prior mean of the coefficients.
#' @param B0i the prior precision of the coefficients.
#'
#' @return Invisibly \code{TRUE}; called for the errors it raises.
#'
#' @keywords internal
varsel_check_ssvs_prior <- function(idx, nms, b0, B0i) {

  # The restricted half of the mixture stands for the regressor being absent
  # from the frontier, which it only does if the prior is centred on zero. A
  # prior mean elsewhere would make it stand for the coefficient being pinned
  # at that other value, which is not what the posterior inclusion probability
  # would then be read as, so it is refused rather than quietly overridden.
  wrong <- idx[b0[idx] != 0]
  if (length(wrong) > 0) {
    stop("A coefficient under variable selection is given a prior centred on ",
         "zero, since the excluded half of its prior stands for the ",
         "regressor being absent from the frontier. 'coef$mu' is not zero ",
         "for ", paste(nms[wrong], collapse = ", "), ". Leave the prior mean ",
         "of ", if (length(wrong) == 1L) "that coefficient" else
           "those coefficients",
         " at zero, or take ", if (length(wrong) == 1L) "it" else "them",
         " out of 'varsel$include'.")
  }

  # The sampler rewrites the diagonal entry of a selected coefficient in every
  # sweep and leaves the rest of the precision matrix alone, so a prior that
  # ties such a coefficient to another one would be half replaced and half
  # kept. Their prior is therefore required to be independent of the others.
  off <- B0i[idx, , drop = FALSE]
  off[cbind(seq_along(idx), idx)] <- 0
  if (any(off != 0)) {
    stop("'coef$v_i' ties a coefficient under variable selection to another ",
         "coefficient. The prior of a selected coefficient is the two-point ",
         "mixture in 'varsel', which is independent of the other ",
         "coefficients, so it cannot carry an off-diagonal precision as ",
         "well. Set the entries of 'coef$v_i' outside the diagonal to zero ",
         "for ", paste(nms[idx], collapse = ", "), ".")
  }

  invisible(TRUE)
}

#' Check the coefficient prior that Bayesian variable selection needs
#'
#' A coefficient the indicators have switched off is absent from the
#' likelihood, so the sweep draws it from its prior. A prior precision of zero
#' leaves no prior to draw from, and the flat limiting prior that the package
#' otherwise allows is therefore refused here.
#'
#' A proper but very wide prior is worse than it looks, and is the reason this
#' check is scale-aware rather than a fixed bound. The wider the prior, the
#' further an excluded coefficient wanders while it is out, and a coefficient
#' far from anything the data support is one the likelihood will not admit
#' back. Measured against the least squares standard error of the same
#' coefficient, on 400 observations of simulated data: at a prior standard
#' deviation of six times that error the indicator of an irrelevant regressor
#' changed state 342 times in 2000 sweeps, at twenty times 116, at sixty times
#' 48, and at two hundred times not once, leaving an effective sample size of
#' zero and a posterior inclusion probability of exactly zero that reflects the
#' chain being stuck rather than the data. The bound is therefore set at a
#' hundred, and this is a warning rather than an error because a stuck chain is
#' visible in the draws and the user may have a reason.
#'
#' Note that the prior width moves the answer even where the chain mixes: the
#' posterior inclusion probability of the irrelevant regressor above fell from
#' 0.17 to 0.03 over the range that still mixed. That is Bartlett's paradox and
#' is a property of the algorithm, not of this implementation. It is why the
#' stochastic search of George, Sun and Ni (2008), whose excluded coefficients
#' never leave the neighbourhood of zero, is the easier of the two to use with
#' a prior that was not chosen with the selection in mind.
#'
#' @param idx the positions under selection.
#' @param nms the column names of the design matrix.
#' @param B0i the prior precision of the coefficients.
#' @param se the least squares standard errors of the coefficients.
#'
#' @return Invisibly \code{TRUE}; called for the error it raises.
#'
#' @keywords internal
varsel_check_bvs_prior <- function(idx, nms, B0i, se) {

  prec <- diag(B0i)[idx]

  flat <- idx[prec <= 0]
  if (length(flat) > 0) {
    stop("A coefficient under the Bayesian variable selection of Korobilis ",
         "(2013) needs a proper prior, since the sweeps in which it is ",
         "excluded draw it from that prior rather than from the data. ",
         "'coef$v_i' is zero for ", paste(nms[flat], collapse = ", "),
         ". Give ", if (length(flat) == 1L) "it" else "them",
         " a positive prior precision, or use varsel = \"ssvs\", whose ",
         "excluded coefficients are held at zero by the prior itself.")
  }

  ratio <- (1 / sqrt(prec)) / se[idx]
  wide <- ratio > 100
  if (any(wide)) {
    worst <- max(ratio[wide])
    warning("The prior on ", paste(nms[idx][wide], collapse = ", "),
            " is too wide for the Bayesian variable selection to move: its ",
            "standard deviation is ", format(round(worst)), " times the ",
            "least squares standard error of the coefficient, so a sweep ",
            "that excludes the regressor draws a coefficient the likelihood ",
            "will not admit back, and the indicator is likely to stay where ",
            "it started for the whole run. Give 'coef$v_i' a precision of ",
            "the order of the inverse squared scale of the coefficient, or ",
            "use varsel = \"ssvs\", which does not have to travel between ",
            "the prior and the data to change state.", call. = FALSE)
  }

  invisible(TRUE)
}

#' Resolve which coefficients are placed under variable selection
#'
#' @param spec the merged \code{varsel} specification.
#' @param nms the column names of the design matrix.
#'
#' @return An integer vector of positions, in increasing order.
#'
#' @keywords internal
varsel_include <- function(spec, nms) {

  k <- length(nms)
  idx <- seq_len(k)

  if (!is.null(spec$include)) {
    if (is.character(spec$include)) {
      unknown <- setdiff(spec$include, nms)
      if (length(unknown) > 0) {
        stop("'varsel$include' names ", paste(unknown, collapse = ", "),
             ", which ", if (length(unknown) == 1L) "is not a column" else
               "are not columns",
             " of the design matrix. Its columns are ",
             paste(nms, collapse = ", "), ".")
      }
      idx <- match(spec$include, nms)
    } else {
      if (!is.numeric(spec$include) || !all(is.finite(spec$include)) ||
          any(spec$include != trunc(spec$include)) ||
          any(spec$include < 1) || any(spec$include > k)) {
        stop("'varsel$include' must be column names of the design matrix or ",
             "whole numbers between 1 and ", k, ".")
      }
      idx <- as.integer(spec$include)
    }
  }
  idx <- sort(unique(idx))

  if (length(spec$exclude_intercept) != 1L ||
      !is.logical(spec$exclude_intercept) || is.na(spec$exclude_intercept)) {
    stop("'varsel$exclude_intercept' must be TRUE or FALSE.")
  }
  # The intercept is the level of the frontier rather than the effect of a
  # regressor, and selecting it away would move every efficiency score
  # instead of dropping a variable, so it is left out by default.
  if (spec$exclude_intercept) {
    idx <- setdiff(idx, which(nms == "(Intercept)"))
  }

  if (length(idx) == 0L) {
    stop("No coefficient is left under variable selection. ",
         if (isTRUE(spec$exclude_intercept) && identical(nms, "(Intercept)"))
           paste("The frontier has an intercept and nothing else, and the",
                 "intercept is excluded by default.") else
           "'varsel$include' selects none.")
  }
  idx
}

#' Prior standard deviations of the two mixture components
#'
#' Either taken from \code{tau} as given, or, under the semiautomatic approach
#' of George, Sun and Ni (2008), obtained by scaling the least squares standard
#' error of each coefficient by the two factors in \code{semiautomatic}. The
#' least squares fit ignores the one-sided term, so its residual variance
#' carries the variance of the inefficiency as well and the standard errors it
#' gives are somewhat the wider for it. That widens both components of the
#' mixture by the same factor, which is the point of tying them to the scale of
#' the data in the first place.
#'
#' @param spec the merged \code{varsel} specification.
#' @param object a model object.
#' @param idx the positions under selection.
#'
#' @return A list with \code{tau0} and \code{tau1}, each of the length of
#'   \code{idx}.
#'
#' @keywords internal
varsel_tau <- function(spec, object, idx) {

  if (is.null(spec$tau) == is.null(spec$semiautomatic)) {
    stop("Exactly one of 'varsel$tau' and 'varsel$semiautomatic' must be ",
         "given: 'tau' sets the two prior standard deviations directly, ",
         "'semiautomatic' obtains them by scaling the least squares standard ",
         "errors of the coefficients.")
  }

  pair <- function(x, what) {
    x <- as.numeric(x)
    if (length(x) != 2L || !all(is.finite(x)) || any(x <= 0) || x[1] >= x[2]) {
      stop("'varsel$", what, "' must be two positive numbers, the one ",
           "belonging to an excluded coefficient first and the smaller of ",
           "the two. The mixture only separates the two states if the prior ",
           "of an excluded coefficient is the tighter one.")
    }
    x
  }

  if (!is.null(spec$tau)) {
    tau <- pair(spec$tau, "tau")
    return(list(tau0 = rep(tau[1], length(idx)),
                tau1 = rep(tau[2], length(idx))))
  }

  fac <- pair(spec$semiautomatic, "semiautomatic")
  se <- ols_se(object)[idx]
  list(tau0 = fac[1] * se, tau1 = fac[2] * se)
}

#' Prior inclusion probabilities
#'
#' @param inprior the value supplied, one number or one per coefficient.
#' @param n the number of coefficients under selection.
#'
#' @return A numeric vector of length \code{n}.
#'
#' @keywords internal
varsel_inprior <- function(inprior, n) {

  if (!is.numeric(inprior) ||
      !(length(inprior) == 1L || length(inprior) == n)) {
    stop("'varsel$inprior' must be a single number or one per coefficient ",
         "under selection, of which there are ", n, ".")
  }
  for (p in inprior) {
    check_probability(p, "varsel$inprior")
  }
  if (length(inprior) == 1L) rep(inprior, n) else as.numeric(inprior)
}

#' Least squares standard errors of the frontier coefficients
#'
#' @param object a model object.
#'
#' @return A numeric vector with one entry per coefficient.
#'
#' @keywords internal
ols_se <- function(object) {

  X <- object$data$X
  fit <- stats::.lm.fit(X, object$data$y)
  s2 <- sum(as.numeric(fit$residuals)^2) / (object$n - object$k)
  # The design is checked for full rank when the model is created, so the
  # cross-product is positive definite and this inverse exists.
  xtxi <- chol2inv(chol(crossprod(X)))
  sqrt(s2 * diag(xtxi))
}

#' Arguments the samplers take for the variable selection block
#'
#' Both samplers take the same five, and take them empty when the model has no
#' variable selection, so that the block is simply skipped. The two standard
#' deviations are empty under \code{"bvs"} as well, which has no mixture.
#'
#' @param object a model object with priors attached.
#'
#' @return A list with \code{ssvs_idx}, \code{tau0}, \code{tau1},
#'   \code{prob_prior} and \code{varsel}.
#'
#' @keywords internal
varsel_args <- function(object) {

  sel <- object$priors$varsel
  if (is.null(sel)) {
    return(list(ssvs_idx = integer(0), tau0 = numeric(0), tau1 = numeric(0),
                prob_prior = numeric(0), varsel = 0L))
  }
  list(ssvs_idx = as.integer(sel$include) - 1L,
       tau0 = as.numeric(sel$tau0),
       tau1 = as.numeric(sel$tau1),
       prob_prior = as.numeric(sel$inprior),
       varsel = varsel_code(sel$algorithm))
}

#' Add the posterior inclusion probabilities to a summary table
#'
#' The column is placed before the effective sample size rather than after it,
#' so that the statistics of the posterior stay together and the two columns
#' that qualify them rather than describe it come last.
#'
#' Rows that are not under selection take \code{NA}: the error and inefficiency
#' parameters are always in the model, and so is a coefficient the selection
#' was not applied to, such as the intercept. Reporting 1 for them would not
#' distinguish a parameter that was never a candidate from one that the data
#' kept.
#'
#' @param object a model object.
#' @param tab the summary table built by \code{\link{summary.sfmodel}}.
#'
#' @return The table, with a \code{PIP} column if the model has inclusion
#'   draws.
#'
#' @keywords internal
add_pip_column <- function(object, tab) {

  draws <- object$posterior$inclusion$coeffs
  if (is.null(draws)) {
    return(tab)
  }

  pip <- rep(NA_real_, nrow(tab))
  names(pip) <- rownames(tab)
  inc <- colMeans(as.matrix(draws))
  known <- intersect(names(inc), rownames(tab))
  pip[known] <- inc[known]

  ess <- match("ESS", colnames(tab))
  cbind(tab[, -ess, drop = FALSE], PIP = pip,
        ESS = tab[, ess])
}

#' Attach the inclusion draws to a posterior
#'
#' @param object a model object.
#' @param posterior the posterior list being assembled.
#' @param inclusion the matrix of inclusion draws, or \code{NULL}.
#'
#' @return The posterior list, with \code{inclusion} added if there was one.
#'
#' @keywords internal
add_inclusion_draws <- function(object, posterior, inclusion) {

  if (is.null(inclusion)) {
    return(posterior)
  }
  colnames(inclusion) <- object$priors$varsel$names
  # Not named 'lambda', which in the exponential models is already the rate of
  # the inefficiency distribution.
  posterior$inclusion <- list(coeffs = .mcmc_draws(object, inclusion))
  posterior
}
