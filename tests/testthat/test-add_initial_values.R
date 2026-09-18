test_that("ols starting values match a least squares fit", {
  d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
  m <- add_initial_values(add_priors(
    create_sfmodel_exp(y ~ x1 + x2, data = d)))

  fit <- lm(y ~ x1 + x2, data = d)
  expect_equal(m$initial$beta, unname(coef(fit)))
  expect_equal(m$initial$sigma_v2, summary(fit)$sigma^2)
  expect_equal(m$initial$u, rep(0, 200))
  expect_equal(m$initial$method, "ols")
})

test_that("each model starts its inefficiency parameter at its prior mean", {
  d <- sim_sf(n = 100, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)

  e <- add_initial_values(add_priors(create_sfmodel_exp(y ~ x1, data = d)))
  expect_equal(e$initial$par_u, e$priors$shape_u / e$priors$rate_u)

  # The half-normal prior is on the precision, so the starting value is the
  # square root of the prior mean of sigma_u^2.
  h <- add_initial_values(add_priors(create_sfmodel_hn(y ~ x1, data = d)))
  expect_equal(h$initial$par_u,
               sqrt(h$priors$rate_u / (h$priors$shape_u - 1)))
})

test_that("prior starting values are dispersed across chains", {
  d <- sim_sf(n = 100, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d))

  set.seed(1)
  a <- add_initial_values(m, method = "prior")
  b <- add_initial_values(m, method = "prior")

  expect_false(isTRUE(all.equal(a$initial$beta, b$initial$beta)))
  expect_equal(a$initial$method, "prior")
  expect_true(a$initial$sigma_v2 > 0)
})

test_that("initial values require priors and a known method", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- create_sfmodel_exp(y ~ x1, data = d)

  expect_error(add_initial_values(m), "Add priors before initial values")
  expect_error(add_initial_values(add_priors(m), method = "magic"),
               "should be one of")
})

test_that("a seed is stored and validated", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- create_sfmodel_exp(y ~ x1, data = d)

  expect_equal(add_seed(m, 99)$model$seed, 99)
  expect_error(add_seed(m, "a"), "single number")
  expect_error(add_seed(m, c(1, 2)), "single number")
})

test_that("replacing the seed discards draws it did not produce", {
  set.seed(63)
  d <- sim_sf(n = 80, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 100, burnin = 50)))

  expect_message(reseeded <- add_seed(est, 999), "Dropping posterior draws")
  expect_null(reseeded$posterior)
  expect_equal(reseeded$model$seed, 999)

  # A model without draws has nothing to drop and says nothing.
  expect_silent(add_seed(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 100, burnin = 50)), 999))
})

test_that("a seed set.seed() cannot take is rejected", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- create_sfmodel_exp(y ~ x1, data = d)

  expect_error(add_seed(m, NA), "single number")
  expect_error(add_seed(m, Inf), "single number")
  # Accepted before, then silently coerced to NA when the sampler ran.
  expect_error(add_seed(m, 2^31), "range of an integer")
  expect_error(add_seed(m, 1.5), "whole number")

  expect_equal(add_seed(m, -99)$model$seed, -99)
})

test_that("a prior draw that overflows is redrawn rather than kept", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)

  # Under the default shape and rate of 0.01, about one draw in twelve hundred
  # of 1 / rgamma() comes back as Inf. That start turns the first sweep into
  # NaN and the run then fails inside inv_sympd().
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d, iterations = 50,
                                     burnin = 0))
  set.seed(96)
  starts <- replicate(300, add_initial_values(m, method = "prior")$initial)
  sigma_v2 <- unlist(starts["sigma_v2", ])
  expect_true(all(is.finite(sigma_v2)))
  expect_true(all(sigma_v2 > 0))
  expect_true(all(is.finite(unlist(starts["par_u", ]))))

  # An even vaguer prior overflows about half the time and must still work.
  vague <- add_priors(create_sfmodel_exp(y ~ x1, data = d, iterations = 50,
                                         burnin = 0),
                      sigma = list(shape = 0.001, rate = 0.001))
  set.seed(97)
  v <- replicate(200, add_initial_values(vague, method = "prior")$initial$sigma_v2)
  expect_true(all(is.finite(v)))

  # And the chain it starts is finite throughout.
  fit <- add_posterior_coefficients(
    add_seed(add_initial_values(vague, method = "prior"), 5), keep_u = FALSE)
  expect_true(all(is.finite(as.matrix(fit$posterior$beta$coeffs))))
})

test_that("the four-component model redraws its variances too", {
  d <- sim_sf4(n = 15, n_time = 3, beta = c(1, 0.5))
  vague <- add_priors(create_sfmodel4_exp(y ~ x1, data = d, id = "id"),
                      sigma = list(shape = 0.001, rate = 0.001),
                      sigma_mu = list(shape = 0.001, rate = 0.001))
  set.seed(98)
  s <- replicate(200, {
    i <- add_initial_values(vague, method = "prior")$initial
    c(i$sigma_v2, i$sigma_mu2, i$par_eta, i$par_u)
  })
  expect_true(all(is.finite(s)))
  expect_true(all(s > 0))
})

test_that("an improper prior on the coefficients has no draw", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)

  # v_i = 0 is the flat limiting prior, which add_priors() allows. It has no
  # covariance, so there is nothing to draw a starting value from; this used to
  # surface as a LAPACK message about a singular system.
  flat <- add_priors(create_sfmodel_exp(y ~ x1, data = d),
                     coef = list(v_i = 0))
  expect_error(add_initial_values(flat, method = "prior"), "improper prior")
  expect_error(add_initial_values(flat, method = "prior"), "method = \"ols\"")

  rank_deficient <- add_priors(create_sfmodel_exp(y ~ x1, data = d),
                               coef = list(v_i = tcrossprod(c(1, 1))))
  expect_error(add_initial_values(rank_deficient, method = "prior"),
               "improper prior")

  # Least squares starting values do not need the prior and still work.
  expect_equal(add_initial_values(flat)$initial$beta,
               unname(coef(lm(y ~ x1, data = d))))
})

test_that("prior draws of the coefficients match the prior they came from", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  v_i <- matrix(c(4, 1, 1, 2), 2)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d),
                  coef = list(mu = c(3, -1), v_i = v_i))

  set.seed(99)
  b <- t(replicate(20000, add_initial_values(m, method = "prior")$initial$beta))
  expect_equal(colMeans(b), c(3, -1), tolerance = 0.05)
  expect_equal(cov(b), solve(v_i), tolerance = 0.05)
})
