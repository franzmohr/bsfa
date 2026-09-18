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
