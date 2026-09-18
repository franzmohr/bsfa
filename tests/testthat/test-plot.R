estimated <- function(keep_u = TRUE) {
  set.seed(31)
  d <- sim_sf(n = 120, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
  add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1 + x2, data = d,
                       iterations = 200, burnin = 100)), keep_u = keep_u)
}

test_that("every plot type draws without error and restores par", {
  est <- estimated()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  before <- par(no.readonly = TRUE)
  for (type in c("hist", "trace", "boxplot", "efficiency")) {
    expect_silent(plot(est, type = type))
  }
  expect_equal(par("mfrow"), before$mfrow)
  expect_equal(par("mar"), before$mar)
})

test_that("plot returns its input invisibly", {
  est <- estimated()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_invisible(plot(est))
  expect_identical(plot(est), est)
})

test_that("the efficiency plot thins the caterpillar to the requested units", {
  est <- estimated()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  # plot_sf_efficiency returns the summary it drew, over all units; the
  # thinning applies to the panel rather than to the numbers.
  eff <- plot(est, type = "efficiency", units = 20)
  expect_equal(nrow(eff), 120L)
  expect_silent(plot(est, type = "efficiency", units = NULL))
  expect_silent(plot(est, type = "efficiency", units = Inf))
})

test_that("the parameter matrix carries the right columns", {
  est <- estimated()
  expect_equal(colnames(sf_par_draws(est)),
               c("(Intercept)", "x1", "x2", "sigma_v", "lambda"))

  set.seed(32)
  d <- sim_sf(n = 80, beta = c(1, 0.5), sigma_v = 0.2, par_u = 0.3,
              ineff = "halfnormal")
  hn <- add_posterior_coefficients(add_priors(
    create_sfmodel_hn(y ~ x1, data = d, iterations = 100, burnin = 50)))
  # The half-normal model appends the signal-to-noise ratio.
  expect_equal(colnames(sf_par_draws(hn)),
               c("(Intercept)", "x1", "sigma_v", "sigma_u", "lambda"))
})

test_that("panel layouts fill rows up to max_cols", {
  expect_equal(grid_dim(1, 3), c(1, 1))
  expect_equal(grid_dim(3, 3), c(1, 3))
  expect_equal(grid_dim(5, 3), c(2, 3))
  expect_equal(grid_dim(5, 2), c(3, 2))
})

test_that("plotting is refused without the draws it needs", {
  set.seed(33)
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d,
                                     iterations = 100, burnin = 50))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_error(plot(m), "does not contain posterior draws")
  expect_error(plot(estimated(keep_u = FALSE), type = "efficiency"), "keep_u")
  expect_error(plot(estimated(), type = "efficiency", ci = 2),
               "between 0 and 1")
  expect_error(plot(estimated(), type = "nonsense"), "should be one of")
})
