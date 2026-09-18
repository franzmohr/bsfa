test_that("output has the expected structure", {
  set.seed(1)
  d <- sim_sf(n = 80, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
  est <- draw_posterior(add_priors(
    create_sfmodel_exp(y ~ x1 + x2, data = d,
                       iterations = 200, burnin = 100, thin = 2)))

  expect_s3_class(est, "bsfa_exp")
  expect_s3_class(est, "bsfa")
  expect_equal(dim(est$draws$beta), c(100L, 3L))
  expect_equal(colnames(est$draws$beta), c("(Intercept)", "x1", "x2"))
  expect_length(est$draws$sigma_v, 100L)
  expect_length(est$draws$lambda, 100L)
  expect_null(est$draws$sigma_u)
  expect_equal(dim(est$u), c(100L, 80L))
  expect_true(all(est$u >= 0))
  expect_true(all(est$draws$sigma_v > 0))
})

test_that("the half-normal result names its own parameter", {
  set.seed(2)
  d <- sim_sf(n = 80, beta = c(1, 0.5), sigma_v = 0.2, par_u = 0.3,
              ineff = "halfnormal")
  est <- draw_posterior(add_priors(
    create_sfmodel_hn(y ~ x1, data = d, iterations = 200, burnin = 100)))

  expect_s3_class(est, "bsfa_hn")
  expect_length(est$draws$sigma_u, 200L)
  expect_null(est$draws$lambda)
  # The signal-to-noise ratio is reported for the half-normal model only.
  expect_true("lambda" %in% rownames(summary(est)$coefficients))
})

test_that("the exponential model recovers its parameters", {
  set.seed(123)
  beta <- c(1, 0.5, 0.3)
  d <- sim_sf(n = 1500, beta = beta, sigma_v = 0.2, par_u = 4)
  est <- draw_posterior(add_priors(
    create_sfmodel_exp(y ~ x1 + x2, data = d,
                       iterations = 2000, burnin = 1000)))

  pm <- colMeans(est$draws$beta)
  expect_equal(unname(pm[2:3]), beta[2:3], tolerance = 0.05)
  # The intercept is only weakly separated from the mean of the one-sided term,
  # so it is checked against a looser bound than the slopes.
  expect_equal(unname(pm[1]), beta[1], tolerance = 0.15)
  expect_equal(mean(est$draws$sigma_v), 0.2, tolerance = 0.1)
  expect_equal(mean(est$draws$lambda), 4, tolerance = 1.5)
})

test_that("the half-normal model recovers its parameters", {
  set.seed(321)
  d <- sim_sf(n = 1500, beta = c(1, 0.5), sigma_v = 0.2, par_u = 0.3,
              ineff = "halfnormal")
  est <- draw_posterior(add_priors(
    create_sfmodel_hn(y ~ x1, data = d, iterations = 2000, burnin = 1000)))

  expect_equal(unname(colMeans(est$draws$beta)[2]), 0.5, tolerance = 0.05)
  expect_equal(mean(est$draws$sigma_u), 0.3, tolerance = 0.15)
})

test_that("the cost frontier flips the sign of the one-sided term", {
  set.seed(7)
  d <- sim_sf(n = 800, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4,
              type = "cost")

  est <- draw_posterior(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, type = "cost",
                       iterations = 1000, burnin = 500)))
  expect_equal(unname(colMeans(est$draws$beta)[2]), 0.5, tolerance = 0.05)

  # Estimating the same data as a production frontier must inflate the
  # residual variance, because the one-sided term is then pushed the wrong way.
  wrong <- draw_posterior(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, type = "production",
                       iterations = 1000, burnin = 500)))
  expect_gt(mean(wrong$draws$sigma_v), mean(est$draws$sigma_v))
})

test_that("efficiency scores track the simulated truth", {
  set.seed(99)
  d <- sim_sf(n = 400, beta = c(1, 0.5), sigma_v = 0.1, par_u = 4)
  est <- draw_posterior(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 1000, burnin = 500)))

  eff <- efficiency(est)
  expect_equal(nrow(eff), 400L)
  expect_true(all(eff$mean > 0 & eff$mean <= 1))
  expect_true(all(eff[["5%"]] <= eff[["95%"]]))
  expect_gt(cor(eff$mean, attr(d, "efficiency")), 0.8)
})

test_that("the panel model gives one inefficiency term per unit", {
  set.seed(11)
  d <- sim_sf(n = 100, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4,
              n_time = 5)
  est <- draw_posterior(add_priors(
    create_sfmodel_exp(y ~ x1 + x2, data = d, id = "id",
                       iterations = 500, burnin = 250)))

  expect_equal(est$n, 500L)
  expect_equal(est$n_units, 100L)
  expect_equal(ncol(est$u), 100L)
  expect_gt(cor(colMeans(exp(-est$u)), attr(d, "efficiency")), 0.8)
})

test_that("a stored seed makes the run reproducible", {
  set.seed(3)
  d <- sim_sf(n = 100, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_seed(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 200, burnin = 100)), 777)

  expect_equal(draw_posterior(m)$draws$beta, draw_posterior(m)$draws$beta)
})

test_that("starting from the prior reaches the same posterior", {
  set.seed(17)
  d <- sim_sf(n = 600, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d,
                                     iterations = 2000, burnin = 1000))

  a <- draw_posterior(add_initial_values(m, method = "ols"))
  b <- draw_posterior(add_initial_values(m, method = "prior"))

  expect_equal(colMeans(a$draws$beta), colMeans(b$draws$beta),
               tolerance = 0.05)
})

test_that("priors must be added before drawing", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- create_sfmodel_exp(y ~ x1, data = d, iterations = 100, burnin = 50)

  expect_error(draw_posterior(m), "Add priors before drawing")
  # Initial values, unlike priors, are filled in with their default.
  expect_s3_class(draw_posterior(add_priors(m)), "bsfa")
})

test_that("pointwise log-likelihood is stored only without an id", {
  set.seed(5)
  d <- sim_sf(n = 100, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4, n_time = 2)

  est <- draw_posterior(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 200, burnin = 100)),
    keep_ll = TRUE)
  expect_equal(dim(est$log_lik), c(200L, 200L))
  expect_true(all(is.finite(est$log_lik)))

  expect_warning(
    draw_posterior(add_priors(
      create_sfmodel_exp(y ~ x1, data = d, id = "id",
                         iterations = 200, burnin = 100)), keep_ll = TRUE),
    "not available for panel models")
})
