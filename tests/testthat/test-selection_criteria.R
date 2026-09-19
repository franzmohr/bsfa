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
  # It used to inherit the model's class as well, which bought nothing and sent
  # summary(), plot() and efficiency() to methods that then complained about
  # missing posterior draws.
  expect_equal(class(crit), "selcrit")
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

  expect_s3_class(crit, "selcrit")
  # Which model it came from is recorded in the object rather than in its
  # class, so that it does not answer to the model's own methods.
  expect_equal(crit$model$ineff, "halfnormal")
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

test_that("models that can never have criteria say so", {
  # Pointing the user at add_posterior_loglik() is unhelpful when that function
  # will refuse them too, so each case explains why rather than redirecting.
  set.seed(60)
  d4 <- sim_sf4(n = 30, n_time = 4, beta = c(1, 0.5))
  m4 <- add_posterior_coefficients(add_priors(
    create_sfmodel4_exp(y ~ x1, data = d4, id = "id",
                        iterations = 100, burnin = 50)))
  expect_error(selection_criteria(m4), "not available for the four-component")
  expect_error(selection_criteria(m4), "does not factorise")

  d2 <- sim_sf(n = 40, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4, n_time = 3)
  mp <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d2, id = "id",
                       iterations = 100, burnin = 50)))
  expect_error(selection_criteria(mp), "panel model cannot have one")
})

test_that("printing lays the criteria out as a table", {
  crit <- selection_criteria(fitted_model())
  expect_output(print(crit), "WAIC")
  expect_output(print(crit), "exponential")
})

test_that("a non-finite credible band names the argument at fault", {
  set.seed(64)
  d <- sim_sf(n = 60, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  est <- add_posterior_loglik(add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 200, burnin = 100))))

  # NA compares to NA, not to FALSE, so a guard written as a chain of
  # comparisons used to raise R's own error instead of this one.
  for (bad in list(NA, NaN, Inf, "a", c(0.9, 0.95))) {
    expect_error(selection_criteria(est, ci = bad), "between 0 and 1")
    expect_error(summary(est, ci = bad), "between 0 and 1")
    expect_error(plot(est, type = "efficiency", ci = bad), "between 0 and 1")
  }
})

test_that("a posterior that cannot support the criteria is refused", {
  m <- fitted_model()
  damage <- function(f) f(m)

  # The parameter blocks. The deviance the penalties are added to is taken at
  # the posterior mean of these, so a missing one used to produce a table of
  # NA behind a warning from mean(), and a mis-shaped beta reached the
  # likelihood and failed there with R's own message about conformability.
  expect_error(selection_criteria(damage(function(x) {
    x$posterior$beta <- NULL; x })), "no usable draws")
  expect_error(selection_criteria(damage(function(x) {
    x$posterior$sigma_v <- NULL; x })), "no usable draws")
  expect_error(selection_criteria(damage(function(x) {
    x$posterior$lambda <- NULL; x })), "no usable draws")
  expect_error(selection_criteria(damage(function(x) {
    x$posterior$beta$coeffs <- x$posterior$beta$coeffs[, 1, drop = FALSE]
    x })), "but the model has 2 coefficient(s)", fixed = TRUE)
  expect_error(selection_criteria(damage(function(x) {
    x$posterior$sigma_v$coeffs <- x$posterior$sigma_v$coeffs[1:5, , drop = FALSE]
    x })), "different numbers of draws")

  # The log-likelihood itself. A matrix of the wrong width is the one that
  # matters: every criterion is built from it, so the result came back looking
  # entirely reasonable and was the WAIC of nothing in particular.
  expect_error(selection_criteria(damage(function(x) {
    x$posterior$loglik <- x$posterior$loglik[, 1:5]
    x })), "has 5 column(s) but the model has 150 observation(s)", fixed = TRUE)
  expect_error(selection_criteria(damage(function(x) {
    x$posterior$loglik <- matrix("a", 10, 150); x })), "no usable draws")
  expect_error(selection_criteria(damage(function(x) {
    x$posterior$loglik <- x$posterior$loglik[1:5, , drop = FALSE]
    x })), "holds 5 draw(s) and posterior$beta holds", fixed = TRUE)
})

test_that("thinning and dropped rows are not mistaken for damage", {
  set.seed(24)
  d <- sim_sf(n = 120, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)

  # Thinning shortens the draws and an incomplete row shortens the sample.
  # Both change the shape of the log-likelihood legitimately, so the checks
  # are made against the model rather than against a remembered size.
  thinned <- add_posterior_loglik(add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 500, burnin = 200,
                       thin = 5))))
  expect_equal(nrow(thinned$posterior$loglik), 100L)
  expect_s3_class(selection_criteria(thinned), "selcrit")

  d$y[1:7] <- NA
  dropped <- suppressMessages(add_posterior_loglik(add_posterior_coefficients(
    add_priors(create_sfmodel_exp(y ~ x1, data = d, iterations = 300,
                                  burnin = 100)))))
  expect_equal(ncol(dropped$posterior$loglik), 113L)
  expect_s3_class(selection_criteria(dropped), "selcrit")
})
