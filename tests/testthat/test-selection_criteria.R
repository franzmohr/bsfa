fitted_model <- function(ineff = "exponential", n = 150, iterations = 300) {
  set.seed(24)
  d <- sim_sf(n = n, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  ctor <- if (ineff == "exponential") create_sfmodel_exp else create_sfmodel_hn
  m <- ctor(y ~ x1, data = d, iterations = iterations, burnin = 100)
  add_posterior_loglik(add_posterior_coefficients(add_priors(m)))
}

test_that("criteria are returned in the documented shape", {
  crit <- selection_criteria(fitted_model())

  expect_s3_class(crit, "selcrit")
  expect_s3_class(crit, "sfmodel_exp")
  expect_named(crit, c("model", "LL", "AIC", "BIC", "HQ", "WAIC"))

  for (nm in c("LL", "AIC", "BIC", "HQ", "WAIC")) {
    expect_named(crit[[nm]], c("mean", "median", "qlower", "qupper"))
    expect_equal(nrow(crit[[nm]]), 1L)
  }

  # AIC, BIC and HQ are point estimates, so they carry no band.
  for (nm in c("AIC", "BIC", "HQ")) {
    expect_true(is.na(crit[[nm]]$qlower))
    expect_true(is.na(crit[[nm]]$qupper))
  }
  # WAIC gets a normal interval from its standard error.
  expect_false(is.na(crit$WAIC$qlower))
  expect_lt(crit$WAIC$qlower, crit$WAIC$mean)
  expect_gt(crit$WAIC$qupper, crit$WAIC$mean)
})

test_that("LL summarises the draws of the log-likelihood", {
  est <- fitted_model()
  crit <- selection_criteria(est)

  ll_draws <- rowSums(as.matrix(est$posterior$loglik))
  expect_equal(crit$LL$mean, mean(ll_draws))
  expect_equal(crit$LL$median, median(ll_draws))
  expect_equal(crit$LL$qlower, unname(quantile(ll_draws, 0.025)))
})

test_that("the penalties are built on the deviance at the posterior mean", {
  est <- fitted_model()
  crit <- selection_criteria(est)

  beta_bar <- colMeans(as.matrix(est$posterior$beta$coeffs))
  sv_bar <- mean(est$posterior$sigma_v$coeffs)
  lam_bar <- mean(est$posterior$lambda$coeffs)
  e <- as.numeric(est$data$y - est$data$X %*% beta_bar)
  ll_bar <- log(lam_bar) + lam_bar * e + 0.5 * lam_bar^2 * sv_bar^2 +
    pnorm(-e / sv_bar - lam_bar * sv_bar, log.p = TRUE)
  dev <- -2 * sum(ll_bar)

  kappa <- est$k + 2
  tt <- est$n
  expect_equal(crit$AIC$mean, dev + 2 * kappa)
  expect_equal(crit$BIC$mean, dev + log(tt) * kappa)
  expect_equal(crit$HQ$mean, dev + 2 * log(log(tt)) * kappa)

  # The deviance at the point estimate is smaller than the mean deviance over
  # the posterior, which is exactly why the penalties belong to the former.
  expect_lt(dev, mean(-2 * rowSums(as.matrix(est$posterior$loglik))))
})

test_that("WAIC matches its definition", {
  est <- fitted_model()
  ll <- as.matrix(est$posterior$loglik)
  expected <- -2 * sum(log(colMeans(exp(ll))) - apply(ll, 2, var))

  expect_equal(selection_criteria(est)$WAIC$mean, expected)
})

test_that("the log-sum-exp is stable where the naive version overflows", {
  est <- fitted_model()
  # Shifting the log-likelihood by a large constant would overflow exp(), but
  # must move WAIC by exactly -2 * n * shift.
  base <- selection_criteria(est)$WAIC$mean

  shifted <- est
  shifted$posterior$loglik <- est$posterior$loglik + 800
  expect_true(all(is.infinite(exp(as.matrix(shifted$posterior$loglik)))))
  expect_equal(selection_criteria(shifted)$WAIC$mean,
               base - 2 * est$n * 800)
})

test_that("the criteria work for the half-normal model too", {
  crit <- selection_criteria(fitted_model(ineff = "halfnormal"))

  expect_s3_class(crit, "sfmodel_hn")
  expect_true(is.finite(crit$AIC$mean))
  expect_true(is.finite(crit$WAIC$mean))
})

test_that("criteria prefer the distribution the data came from", {
  set.seed(555)
  d <- sim_sf(n = 800, beta = c(1, 0.5), sigma_v = 0.15, par_u = 5)

  fit <- function(ctor) {
    m <- ctor(y ~ x1, data = d, iterations = 2000, burnin = 500)
    add_posterior_loglik(add_posterior_coefficients(add_priors(m)))
  }
  w_exp <- selection_criteria(fit(create_sfmodel_exp))$WAIC$mean
  w_hn <- selection_criteria(fit(create_sfmodel_hn))$WAIC$mean

  expect_lt(w_exp, w_hn)
})

test_that("criteria need a log-likelihood and a valid band", {
  set.seed(6)
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 100, burnin = 50)))

  expect_error(selection_criteria(m), "does not contain draws")
  expect_error(selection_criteria(add_posterior_loglik(m), ci = 0),
               "between 0 and 1")
})

test_that("printing lays the criteria out as a table", {
  crit <- selection_criteria(fitted_model())
  expect_output(print(crit), "WAIC")
  expect_output(print(crit), "exponential")
})
