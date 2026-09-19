sample_with <- function(lambda, n = 300, seed = 170) {
  set.seed(seed)
  x1 <- rnorm(n)
  u <- if (is.na(lambda)) rep(0, n) else rexp(n, lambda)
  data.frame(y = 1 + 0.5 * x1 - u + rnorm(n, sd = 0.25), x1 = x1)
}

test_that("it returns the documented shape", {
  m <- create_sfmodel_exp(y ~ x1, data = sample_with(4), iterations = 400,
                          burnin = 200)
  s <- prior_sensitivity(m, r_star = c(0.6, 0.8))

  expect_s3_class(s, "sfsens")
  expect_named(s, c("efficiency", "coefficients", "spread", "relative",
                    "coef_spread", "seed", "model"))
  expect_equal(names(s$efficiency), c("r_star", "mean", "sd", "min", "max"))
  expect_equal(s$efficiency$r_star, c(0.6, 0.8))
  expect_equal(dim(s$coefficients), c(2L, 2L))
  expect_equal(colnames(s$coefficients), c("(Intercept)", "x1"))
  expect_equal(s$spread, diff(range(s$efficiency$mean)))
})

test_that("it separates a prior driven fit from an estimated one", {
  fit <- function(lambda) {
    prior_sensitivity(create_sfmodel_exp(y ~ x1, data = sample_with(lambda),
                                         iterations = 1500, burnin = 500),
                      r_star = c(0.5, 0.75, 0.9))
  }
  # A one-sided term twice the size of the noise, so that the sample really is
  # informative about the level. At a signal the size of the noise the data
  # genuinely do not pin the level down, and no assertion here should pretend
  # otherwise: an earlier version of this test compared the two spreads by a
  # factor of five and failed about a third of the time.
  empty <- fit(NA)
  plenty <- fit(2)

  # The figure to read is the movement relative to the spread between units.
  expect_gt(empty$relative, 0.5)
  expect_lt(plenty$relative, 0.5)
  expect_gt(empty$relative, plenty$relative)
  expect_equal(empty$relative, empty$spread / mean(empty$efficiency$sd))

  # The mean score rises with the anchor when the prior is doing the work.
  expect_true(all(diff(empty$efficiency$mean) > 0))

  # The frontier holds still either way: it is the distance to it that is in
  # question, not the frontier.
  expect_lt(max(plenty$coef_spread), 0.05)
  expect_output(print(empty), "largely what")
  expect_output(print(plenty), "estimated from the data")
})

test_that("one seed is used for every fit", {
  m <- add_seed(create_sfmodel_exp(y ~ x1, data = sample_with(4),
                                   iterations = 400, burnin = 200), 4242)
  s <- prior_sensitivity(m, r_star = c(0.6, 0.8))
  expect_equal(s$seed, 4242)

  # Repeating the whole thing gives the same table, so a difference between
  # rows is the prior rather than the sampler.
  expect_equal(prior_sensitivity(m, r_star = c(0.6, 0.8))$efficiency,
               s$efficiency)

  # A model without a seed gets one, and uses it throughout.
  set.seed(11)
  a <- prior_sensitivity(create_sfmodel_exp(y ~ x1, data = sample_with(4),
                                            iterations = 400, burnin = 200),
                         r_star = c(0.6, 0.8))
  expect_true(is.numeric(a$seed))
  expect_length(a$seed, 1L)
})

test_that("the rest of the prior specification survives the refit", {
  d <- sample_with(4)
  v <- matrix(c(4, 1, 1, 2), 2)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d, iterations = 200),
                  coef = list(mu = c(3, -1), v_i = v),
                  sigma = list(shape = 2, rate = 0.5),
                  lambda = list(r_star = 0.75, shape = 3))

  rebuilt <- bsfa:::sensitivity_priors(m, 0.9)
  # Only the anchor moves.
  expect_equal(rebuilt$priors$b0, m$priors$b0)
  expect_equal(rebuilt$priors$B0i, m$priors$B0i)
  expect_equal(rebuilt$priors$shape_v, 2)
  expect_equal(rebuilt$priors$rate_v, 0.5)
  expect_equal(rebuilt$priors$shape_u, 3)
  expect_equal(rebuilt$priors$r_star, 0.9)
  # The rate carries the shape, because the anchor is matched on the median of
  # the marginal prior rather than on the rate alone.
  expect_equal(rebuilt$priors$rate_u, -log(0.9) / (2^(1 / 3) - 1))

  # And a model with no priors at all simply gets the defaults.
  bare <- bsfa:::sensitivity_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 200), 0.8)
  expect_equal(bare$priors$r_star, 0.8)
  expect_equal(bare$priors$shape_u, 1)
})

test_that("it works for every model class", {
  d <- sample_with(4)
  expect_s3_class(
    prior_sensitivity(create_sfmodel_hn(y ~ x1, data = d, iterations = 400,
                                        burnin = 200),
                      r_star = c(0.6, 0.8)), "sfsens")

  p <- sim_sf4(n = 25, n_time = 4, beta = c(1, 0.5))
  for (mk in list(create_sfmodel4_exp, create_sfmodel4_hn)) {
    s <- prior_sensitivity(mk(y ~ x1, data = p, id = "id", iterations = 400,
                              burnin = 200), r_star = c(0.8, 0.95))
    expect_s3_class(s, "sfsens")
    expect_equal(nrow(s$efficiency), 2L)
    expect_output(print(s), "both anchors moved together")
  }
})

test_that("the grid is validated", {
  m <- create_sfmodel_exp(y ~ x1, data = sample_with(4), iterations = 200)

  expect_error(prior_sensitivity(m, r_star = 0.75), "at least two values")
  expect_error(prior_sensitivity(m, r_star = c(0.5, 1)), "between 0 and 1")
  expect_error(prior_sensitivity(m, r_star = c(0.5, NA)), "between 0 and 1")
  expect_error(prior_sensitivity(1:10), "stochastic frontier model")

  # Duplicates collapse and the grid is ordered.
  s <- prior_sensitivity(m, r_star = c(0.8, 0.6, 0.8))
  expect_equal(s$efficiency$r_star, c(0.6, 0.8))
})
