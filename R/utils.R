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
#' @param obs observation labels, for scores that are one per observation
#'   rather than one per unit, or \code{NULL}.
#'
#' @return A data frame with one row per column of \code{r}.
#'
#' @keywords internal
efficiency_table <- function(r, unit, probs, obs = NULL) {

  qs <- matrix(apply(r, 2, stats::quantile, probs = probs),
               nrow = ncol(r), byrow = TRUE,
               dimnames = list(NULL,
                               paste0(format(100 * probs, trim = TRUE), "%")))

  out <- data.frame(unit = unit,
                    mean = colMeans(r),
                    sd = apply(r, 2, stats::sd),
                    qs,
                    row.names = NULL,
                    check.names = FALSE,
                    stringsAsFactors = FALSE)

  if (is.null(obs)) {
    return(out)
  }
  cbind(out["unit"], obs = obs, out[setdiff(names(out), "unit")],
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

#' Resolve keep_u to a thinning interval for the augmented blocks
#'
#' The augmented terms are the largest thing the sampler returns, and in the
#' four-component model they are one column per observation. Storing every
#' retained draw of them is what makes a long chain on a wide panel run out of
#' memory, so \code{keep_u} also accepts an interval at which to keep them.
#'
#' @param keep_u \code{TRUE}, \code{FALSE}, or a positive whole number.
#' @param n_keep the number of draws the chain retains, which the interval
#'   cannot exceed without leaving nothing to store.
#'
#' @return \code{0L} if the draws are not to be stored, otherwise the interval.
#'
#' @keywords internal
augmented_thin <- function(keep_u, n_keep) {

  if (is.logical(keep_u)) {
    if (length(keep_u) != 1L || is.na(keep_u)) {
      stop("'keep_u' must be TRUE, FALSE or a positive whole number.")
    }
    return(if (keep_u) 1L else 0L)
  }

  if (length(keep_u) != 1L || !is.numeric(keep_u) || !is.finite(keep_u) ||
      keep_u != trunc(keep_u) || keep_u < 0 ||
      keep_u > .Machine$integer.max) {
    stop("'keep_u' must be TRUE, FALSE or a positive whole number.")
  }
  # An interval wider than the chain would store nothing at all, which as a
  # silent outcome is indistinguishable from keep_u = FALSE and leaves
  # efficiency() advising the very argument that was just passed.
  if (keep_u > n_keep) {
    stop("'keep_u' is ", keep_u, ", but the chain retains only ", n_keep,
         " draw", if (n_keep == 1L) "" else "s",
         ", so not one augmented draw would be kept. The interval cannot ",
         "exceed 'iterations' divided by 'thin'.")
  }
  as.integer(keep_u)
}

#' Warn before the augmented draws fill the machine
#'
#' The sampler allocates the whole block up front, so an over-ambitious
#' \code{keep_u} on a wide panel fails on allocation rather than gradually.
#' Saying so beforehand is more use than the allocation error would be.
#'
#' @param columns total number of columns across the stored blocks.
#' @param n_keep_u number of draws of them that will be stored.
#'
#' @return Invisibly \code{NULL}; called for the warning it raises.
#'
#' @keywords internal
warn_augmented_size <- function(columns, n_keep_u) {

  bytes <- 8 * as.numeric(columns) * as.numeric(n_keep_u)
  if (bytes > 1024^3) {
    warning("The augmented draws will take about ",
            format(round(bytes / 1024^3, 1), nsmall = 1), " GB. Pass ",
            "'keep_u' an interval to thin them, or FALSE to drop them if ",
            "the efficiency scores are not needed.", call. = FALSE)
  }
  invisible(NULL)
}

#' Iteration index of a thinned block of augmented draws
#'
#' The augmented blocks may be thinned beyond the chain's own \code{thin}, so
#' their iteration index counts in steps of the two intervals multiplied.
#'
#' @param object a model object.
#' @param x a matrix of stored draws.
#' @param u_thin the interval the block was stored at.
#'
#' @return An \code{\link[coda]{mcmc}} object.
#'
#' @keywords internal
mcmc_augmented <- function(object, x, u_thin) {
  step <- object$thin * u_thin
  coda::mcmc(x,
             start = object$burnin + step,
             end = object$burnin + step * NROW(x),
             thin = step)
}

#' Draw a starting value from a prior that may overflow
#'
#' A draw from a vague prior on a variance is a draw from an inverse gamma
#' whose shape is near zero, and those routinely exceed what a double can hold:
#' under the package's own default of shape and rate 0.01, about one draw in
#' twelve hundred comes back as \code{Inf}, and under 0.001 nearly half of them
#' do. Such a value is not a dispersed starting point but an unusable one. It
#' turns the first sweep's inefficiency terms into \code{NaN}, and the run then
#' fails inside the linear algebra with a message about a matrix inverse rather
#' than about the prior it came from.
#'
#' The draw is therefore repeated until it is representable. That discards only
#' values no chain could have started from in any case, and leaves the draw a
#' draw from the prior everywhere it can be held.
#'
#' @param draw a function of no arguments returning one draw.
#' @param what the name of the quantity, used in the error message.
#'
#' @return A finite positive draw.
#'
#' @keywords internal
draw_finite <- function(draw, what) {

  for (i in seq_len(100L)) {
    x <- draw()
    if (length(x) == 1L && is.finite(x) && x > 0) {
      return(x)
    }
  }

  stop("Could not draw a usable starting value for '", what, "': every one ",
       "of 100 draws from its prior fell outside the range a double can ",
       "hold. The prior is too vague to start a chain from. Use ",
       "method = \"ols\", or give it a shape and rate that put some mass on ",
       "values a variance could take.")
}

#' Log of the standard normal distribution function, plus half its argument
#' squared
#'
#' The quantity \eqn{\log \Phi(x) + x^2/2}. It is what the exponential
#' composed-error density reduces to once the terms that cancel are cancelled
#' by hand, and it has to be computed as one thing rather than as the sum it is
#' written as: in the left tail \eqn{\log \Phi(x)} is close to \eqn{-x^2/2}, so
#' forming the two separately and adding them loses about \eqn{x^2/2} times the
#' machine epsilon. At \eqn{x = -10^4} that is already \eqn{10^{-8}}, and by
#' \eqn{x = -10^6} nothing of the answer survives.
#'
#' Far enough into the tail the asymptotic expansion of the Mills ratio is used
#' instead, in which the cancelling part never appears. Its next omitted term
#' is \eqn{10395/x^{12}}, which at the crossover is of the order of
#' \eqn{10^{-12}}.
#'
#' @param x a numeric vector.
#'
#' @return A numeric vector of the same length.
#'
#' @keywords internal
log_phi_ratio <- function(x) {

  out <- numeric(length(x))
  tail <- !is.na(x) & x < -20

  if (any(!tail)) {
    out[!tail] <- stats::pnorm(x[!tail], log.p = TRUE) + x[!tail]^2 / 2
  }
  if (any(tail)) {
    w <- 1 / x[tail]^2
    out[tail] <- -0.5 * log(2 * pi) - log(-x[tail]) +
      log1p(w * (-1 + w * (3 + w * (-15 + w * (105 + w * -945)))))
  }
  out
}

#' Resolve the verbose argument to a reporting interval
#'
#' @param verbose \code{TRUE}, \code{FALSE}, or a non-negative whole number.
#' @param n_iter the number of iterations the sampler will run.
#'
#' @return The interval in iterations, \code{0L} for no reporting.
#'
#' @keywords internal
verbose_interval <- function(verbose, n_iter) {

  if (is.logical(verbose)) {
    if (length(verbose) != 1L || is.na(verbose)) {
      stop("'verbose' must be TRUE, FALSE or a non-negative whole number.")
    }
    return(if (verbose) max(1L, as.integer(floor(n_iter / 10))) else 0L)
  }

  if (length(verbose) != 1L || !is.numeric(verbose) || !is.finite(verbose) ||
      verbose != trunc(verbose) || verbose < 0 ||
      verbose > .Machine$integer.max) {
    stop("'verbose' must be TRUE, FALSE or a non-negative whole number.")
  }
  as.integer(verbose)
}

#' Effective sample size of a matrix of draws
#'
#' A thin guard around \code{\link[coda]{effectiveSize}}. The draws a sampler
#' returns are autocorrelated, so the number of them overstates how much the
#' posterior summaries are worth: a chain whose effective size is a hundredth
#' of its length carries a Monte Carlo error ten times what an independent
#' sample of the same length would. That ratio is not visible in the summaries
#' themselves, which is why it is reported beside them.
#'
#' A posterior supplied through \code{posterior_function} need not be long
#' enough, or varied enough, for the estimate to exist, so a failure is
#' reported as \code{NA} rather than allowed to stop the summary.
#'
#' @param draws a matrix of draws, one row per draw.
#'
#' @return A numeric vector with one entry per column of \code{draws}.
#'
#' @keywords internal
sf_ess <- function(draws) {

  nms <- colnames(draws)
  if (NROW(draws) < 3L) {
    return(stats::setNames(rep(NA_real_, NCOL(draws)), nms))
  }

  out <- tryCatch(as.numeric(coda::effectiveSize(coda::mcmc(draws))),
                  error = function(e) rep(NA_real_, NCOL(draws)))
  stats::setNames(out, nms)
}
