sim_with_noise <- function(n = 200, seed = 1234) {
  set.seed(seed)
  d <- sim_sf(n = n, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
  # Two regressors the frontier has no use for.
  d$z1 <- stats::rnorm(n)
  d$z2 <- stats::rnorm(n)
  d
}

ssvs_model <- function(d, ineff = "exponential", ...) {
  create <- if (ineff == "exponential") create_sfmodel_exp else
    create_sfmodel_hn
  create(y ~ x1 + x2 + z1 + z2, data = d, varsel = "ssvs",
         iterations = 400, burnin = 200, ...)
}

test_that("varsel is checked when the model is created", {
  d <- sim_with_noise()
  expect_null(create_sfmodel_exp(y ~ x1, data = d)$model$varsel)
  expect_identical(
    create_sfmodel_exp(y ~ x1, data = d, varsel = "ssvs")$model$varsel,
    "ssvs")
  expect_error(create_sfmodel_exp(y ~ x1, data = d, varsel = "bvs"),
               "must be NULL")
  expect_error(create_sfmodel_exp(y ~ x1, data = d, varsel = TRUE),
               "must be NULL")
  expect_identical(
    create_sfmodel4_hn(y ~ x1, data = sim_sf4(n = 20, n_time = 4),
                       id = "id", varsel = "ssvs")$model$varsel,
    "ssvs")
})

test_that("the varsel prior is required exactly when the model has one", {
  d <- sim_with_noise()
  expect_error(add_priors(create_sfmodel_exp(y ~ x1 + x2, data = d),
                          varsel = list(tau = c(0.1, 10))),
               "not created with varsel")
  expect_error(add_priors(ssvs_model(d)), "'varsel' must be given")
})

test_that("tau and semiautomatic are mutually exclusive and ordered", {
  d <- sim_with_noise()
  m <- ssvs_model(d)
  expect_error(add_priors(m, varsel = list(tau = c(0.1, 10),
                                           semiautomatic = c(0.1, 10))),
               "Exactly one")
  expect_error(add_priors(m, varsel = list(tau = c(10, 0.1))),
               "the smaller of")
  expect_error(add_priors(m, varsel = list(tau = c(0, 10))), "positive")
  expect_error(add_priors(m, varsel = list(tau = 0.1)), "two positive")
  expect_error(add_priors(m, varsel = list(semiautomatic = c(-1, 10))),
               "positive")
})

test_that("tau is expanded and semiautomatic follows the least squares fit", {
  d <- sim_with_noise()
  m <- add_priors(ssvs_model(d), varsel = list(tau = c(0.1, 10)))
  sel <- m$priors$varsel
  # The intercept is left out by default, so four regressors give three.
  expect_identical(sel$names, c("x1", "x2", "z1", "z2"))
  expect_identical(sel$include, 2:5)
  expect_equal(sel$tau0, rep(0.1, 4))
  expect_equal(sel$tau1, rep(10, 4))
  expect_equal(sel$inprior, rep(0.5, 4))

  auto <- add_priors(ssvs_model(d),
                     varsel = list(semiautomatic = c(0.1, 10)))
  se <- ols_se(auto)[2:5]
  expect_equal(auto$priors$varsel$tau0, 0.1 * se)
  expect_equal(auto$priors$varsel$tau1, 10 * se)
  expect_equal(auto$priors$varsel$tau1 / auto$priors$varsel$tau0,
               rep(100, 4))
})

test_that("ols_se matches what lm reports", {
  d <- sim_with_noise()
  m <- ssvs_model(d)
  expect_equal(unname(ols_se(m)),
               unname(coef(summary(stats::lm(y ~ x1 + x2 + z1 + z2,
                                             data = d)))[, "Std. Error"]))
})

test_that("include and exclude_intercept pick the coefficients", {
  d <- sim_with_noise()

  by_name <- add_priors(ssvs_model(d),
                        varsel = list(tau = c(0.1, 10),
                                      include = c("z2", "x1")))
  expect_identical(by_name$priors$varsel$names, c("x1", "z2"))

  by_pos <- add_priors(ssvs_model(d),
                       varsel = list(tau = c(0.1, 10), include = c(5, 2)))
  expect_identical(by_pos$priors$varsel$names, c("x1", "z2"))

  with_const <- add_priors(ssvs_model(d),
                           varsel = list(tau = c(0.1, 10),
                                         exclude_intercept = FALSE))
  expect_identical(with_const$priors$varsel$names,
                   c("(Intercept)", "x1", "x2", "z1", "z2"))

  expect_error(add_priors(ssvs_model(d),
                          varsel = list(tau = c(0.1, 10), include = "x3")),
               "not a column")
  expect_error(add_priors(ssvs_model(d),
                          varsel = list(tau = c(0.1, 10), include = 9)),
               "whole numbers")
  expect_error(add_priors(ssvs_model(d),
                          varsel = list(tau = c(0.1, 10),
                                        exclude_intercept = NA)),
               "TRUE or FALSE")
  expect_error(
    add_priors(create_sfmodel_exp(y ~ 1, data = d, varsel = "ssvs"),
               varsel = list(tau = c(0.1, 10))),
    "No coefficient is left")
})

test_that("inprior takes one value or one per coefficient", {
  d <- sim_with_noise()
  m <- add_priors(ssvs_model(d), varsel = list(tau = c(0.1, 10),
                                               inprior = c(0.1, 0.2, 0.3,
                                                           0.4)))
  expect_equal(m$priors$varsel$inprior, c(0.1, 0.2, 0.3, 0.4))
  expect_error(add_priors(ssvs_model(d),
                          varsel = list(tau = c(0.1, 10),
                                        inprior = c(0.1, 0.2))),
               "one per coefficient")
  expect_error(add_priors(ssvs_model(d),
                          varsel = list(tau = c(0.1, 10), inprior = 1)),
               "strictly between 0 and 1")
})

test_that("an unknown element of varsel is reported", {
  d <- sim_with_noise()
  expect_error(add_priors(ssvs_model(d),
                          varsel = list(tau = c(0.1, 10), taus = 1)),
               "Unknown element")
})

test_that("a selected coefficient needs a prior centred on zero", {
  d <- sim_with_noise()
  expect_error(add_priors(ssvs_model(d), coef = list(mu = 0.5, v_i = 0.01),
                          varsel = list(tau = c(0.1, 10))),
               "centred on zero")
  # A prior mean on a coefficient that is not selected is fine.
  expect_silent(add_priors(ssvs_model(d),
                           coef = list(mu = c(1, 0, 0, 0, 0), v_i = 0.01),
                           varsel = list(tau = c(0.1, 10))))
})

test_that("a selected coefficient cannot carry an off-diagonal precision", {
  d <- sim_with_noise()
  v <- diag(0.01, 5)
  v[2, 3] <- v[3, 2] <- 0.001
  expect_error(add_priors(ssvs_model(d), coef = list(mu = 0, v_i = v),
                          varsel = list(tau = c(0.1, 10))),
               "off-diagonal")
  # The same entry between two coefficients outside the selection passes.
  v2 <- diag(0.01, 5)
  v2[1, 2] <- v2[2, 1] <- 0.001
  expect_silent(add_priors(ssvs_model(d), coef = list(mu = 0, v_i = v2),
                           varsel = list(tau = c(0.1, 10),
                                         include = c("z1", "z2"))))
})

test_that("the stored prior precision is the one of an included coefficient", {
  d <- sim_with_noise()
  m <- add_priors(ssvs_model(d), varsel = list(tau = c(0.1, 10)))
  expect_equal(diag(m$priors$B0i), c(0.01, rep(1 / 100, 4)))
})

test_that("the sampler returns inclusion draws of the right shape", {
  d <- sim_with_noise()
  m <- add_priors(ssvs_model(d), varsel = list(semiautomatic = c(0.1, 10)))
  m <- add_posterior_coefficients(add_seed(add_initial_values(m), 42))

  inc <- m$posterior$inclusion$coeffs
  expect_s3_class(inc, "mcmc")
  expect_identical(colnames(inc), c("x1", "x2", "z1", "z2"))
  expect_identical(nrow(inc), 400L)
  expect_true(all(as.matrix(inc) %in% c(0, 1)))
  expect_identical(coda::thin(inc), 1)
})

test_that("the posterior keeps the regressors that matter and drops the rest", {
  d <- sim_with_noise(n = 400)
  m <- create_sfmodel_exp(y ~ x1 + x2 + z1 + z2, data = d, varsel = "ssvs",
                          iterations = 2000, burnin = 1000)
  m <- add_priors(m, varsel = list(semiautomatic = c(0.1, 10)))
  m <- add_posterior_coefficients(add_seed(add_initial_values(m), 42))

  pip <- colMeans(as.matrix(m$posterior$inclusion$coeffs))
  expect_true(all(pip[c("x1", "x2")] > 0.95))
  expect_true(all(pip[c("z1", "z2")] < 0.5))
  # The coefficients of the regressors that belong are unharmed by the
  # selection, which is the point of leaving the rest of the sampler alone.
  expect_equal(unname(coef(m)[c("x1", "x2")]), c(0.5, 0.3), tolerance = 0.05)
})

test_that("a model without variable selection is left exactly as it was", {
  d <- sim_with_noise()
  fit <- function(...) {
    m <- create_sfmodel_exp(y ~ x1 + x2, data = d, iterations = 300,
                            burnin = 100, ...)
    add_posterior_coefficients(add_seed(add_priors(m), 7))
  }
  m <- fit()
  expect_null(m$priors$varsel)
  expect_null(m$posterior$inclusion)
  expect_null(summary(m)$specifications$varsel)
  expect_false("PIP" %in% colnames(summary(m)$coefficients))
})

test_that("summary reports the inclusion probabilities beside the coefficients", {
  d <- sim_with_noise()
  m <- add_priors(ssvs_model(d), varsel = list(tau = c(0.05, 5)))
  m <- add_posterior_coefficients(add_seed(add_initial_values(m), 42))

  tab <- summary(m)$coefficients
  expect_identical(colnames(tab),
                   c("mean", "sd", "2.5%", "median", "97.5%", "PIP", "ESS"))
  expect_equal(tab[c("x1", "x2", "z1", "z2"), "PIP"],
               colMeans(as.matrix(m$posterior$inclusion$coeffs)))
  # The intercept and the variance parameters were never candidates.
  expect_true(all(is.na(tab[c("(Intercept)", "sigma_v", "lambda"), "PIP"])))
  expect_identical(summary(m)$specifications$varsel, "ssvs")
  expect_output(print(summary(m)), "Variable selection: SSVS")
  expect_output(print(m), "SSVS variable selection")
})

test_that("both inefficiency distributions carry the selection", {
  d <- sim_with_noise()
  for (ineff in c("exponential", "halfnormal")) {
    m <- add_priors(ssvs_model(d, ineff), varsel = list(tau = c(0.05, 5)))
    m <- add_posterior_coefficients(add_seed(add_initial_values(m), 3))
    expect_identical(colnames(m$posterior$inclusion$coeffs),
                     c("x1", "x2", "z1", "z2"))
    # The name of the inefficiency parameter is untouched by the inclusion
    # block, which in the exponential model is also called lambda in the
    # literature and is deliberately not called that here.
    expect_true(m$model$par_u_name %in% names(m$posterior))
  }
})

test_that("the four-component model carries the selection too", {
  set.seed(99)
  d4 <- sim_sf4(n = 40, n_time = 6, beta = c(1, 0.5, 0.3))
  d4$z1 <- stats::rnorm(nrow(d4))

  m <- create_sfmodel4_exp(y ~ x1 + x2 + z1, data = d4, id = "id",
                           varsel = "ssvs", iterations = 600, burnin = 300)
  m <- add_priors(m, varsel = list(semiautomatic = c(0.1, 10)))
  m <- add_posterior_coefficients(add_seed(add_initial_values(m), 5))

  pip <- colMeans(as.matrix(m$posterior$inclusion$coeffs))
  expect_identical(names(pip), c("x1", "x2", "z1"))
  expect_true(all(pip[c("x1", "x2")] > 0.95))
  expect_true(pip[["z1"]] < 0.5)
  expect_true("PIP" %in% colnames(summary(m)$coefficients))
})

test_that("the run is reproducible from the seed", {
  d <- sim_with_noise()
  run <- function() {
    m <- add_priors(ssvs_model(d), varsel = list(tau = c(0.05, 5)))
    add_posterior_coefficients(add_seed(add_initial_values(m), 123))
  }
  a <- run()
  b <- run()
  expect_equal(as.matrix(a$posterior$inclusion$coeffs),
               as.matrix(b$posterior$inclusion$coeffs))
  expect_equal(as.matrix(a$posterior$beta$coeffs),
               as.matrix(b$posterior$beta$coeffs))
})

test_that("an excluded coefficient stays within its tight prior", {
  d <- sim_with_noise(n = 300)
  m <- create_sfmodel_exp(y ~ x1 + x2 + z1 + z2, data = d, varsel = "ssvs",
                          iterations = 1000, burnin = 500)
  # A prior inclusion probability this low keeps the noise regressors out of
  # the frontier in all but a handful of sweeps, and while they are out their
  # coefficients cannot be further from zero than tau0 allows.
  m <- add_priors(m, varsel = list(tau = c(0.001, 10), inprior = 0.01))
  m <- add_posterior_coefficients(add_seed(add_initial_values(m), 42))

  inc <- as.matrix(m$posterior$inclusion$coeffs)
  beta <- as.matrix(m$posterior$beta$coeffs)
  out <- inc[, "z1"] == 0
  expect_true(any(out))
  expect_true(max(abs(beta[out, "z1"])) < 0.01)
})

test_that("prior sensitivity keeps the selection across its fits", {
  d <- sim_with_noise()
  m <- add_priors(ssvs_model(d), varsel = list(semiautomatic = c(0.1, 10)))
  s <- prior_sensitivity(m, r_star = c(0.7, 0.8), keep_u = TRUE)
  expect_identical(colnames(s$coefficients),
                   c("(Intercept)", "x1", "x2", "z1", "z2"))
  expect_identical(nrow(s$coefficients), 2L)
})

test_that("changing the selection prior drops stale posterior draws", {
  d <- sim_with_noise()
  m <- add_priors(ssvs_model(d), varsel = list(tau = c(0.05, 5)))
  m <- add_posterior_coefficients(add_seed(add_initial_values(m), 1))
  expect_message(m2 <- add_priors(m, varsel = list(tau = c(0.1, 10))),
                 "Dropping posterior draws")
  expect_null(m2$posterior)
})
