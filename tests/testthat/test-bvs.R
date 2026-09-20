sim_with_noise <- function(n = 300, seed = 1234) {
  set.seed(seed)
  d <- sim_sf(n = n, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
  d$z1 <- stats::rnorm(n)
  d$z2 <- stats::rnorm(n)
  d
}

bvs_model <- function(d, ineff = "exponential", ...) {
  create <- if (ineff == "exponential") create_sfmodel_exp else
    create_sfmodel_hn
  create(y ~ x1 + x2 + z1 + z2, data = d, varsel = "bvs",
         iterations = 400, burnin = 200, ...)
}

# A prior that is proper and on the scale of the coefficients, which is what
# the Bayesian variable selection needs to move at all.
bvs_coef <- list(mu = 0, v_i = 1)

test_that("bvs is accepted where ssvs is, and nothing else is", {
  d <- sim_with_noise()
  expect_identical(
    create_sfmodel_exp(y ~ x1, data = d, varsel = "bvs")$model$varsel, "bvs")
  expect_identical(
    create_sfmodel4_hn(y ~ x1, data = sim_sf4(n = 20, n_time = 4),
                       id = "id", varsel = "bvs")$model$varsel, "bvs")
  expect_error(create_sfmodel_exp(y ~ x1, data = d, varsel = "lasso"),
               "Korobilis")
})

test_that("the two algorithms reach the sampler with distinct codes", {
  d <- sim_with_noise()
  none <- add_priors(create_sfmodel_exp(y ~ x1 + x2, data = d))
  expect_identical(varsel_args(none)$varsel, 0L)

  ssvs <- add_priors(create_sfmodel_exp(y ~ x1 + x2, data = d,
                                        varsel = "ssvs"),
                     varsel = list(tau = c(0.05, 5)))
  expect_identical(varsel_args(ssvs)$varsel, 1L)

  bvs <- add_priors(bvs_model(d), coef = bvs_coef,
                    varsel = list(inprior = 0.5))
  expect_identical(varsel_args(bvs)$varsel, 2L)
  # BVS has no mixture, so it carries no pair of standard deviations.
  expect_identical(varsel_args(bvs)$tau0, numeric(0))
  expect_identical(bvs$priors$varsel$algorithm, "bvs")
})

test_that("the mixture arguments of ssvs are refused for bvs", {
  d <- sim_with_noise()
  m <- bvs_model(d)
  expect_error(add_priors(m, coef = bvs_coef,
                          varsel = list(tau = c(0.05, 5))),
               "belongs to SSVS")
  expect_error(add_priors(m, coef = bvs_coef,
                          varsel = list(semiautomatic = c(0.1, 10))),
               "belongs to SSVS")
  expect_error(add_priors(m), "'inprior'")
})

test_that("bvs shares include, exclude_intercept and inprior with ssvs", {
  d <- sim_with_noise()
  m <- add_priors(bvs_model(d), coef = bvs_coef,
                  varsel = list(inprior = 0.2, include = c("z1", "x1")))
  expect_identical(m$priors$varsel$names, c("x1", "z1"))
  expect_equal(m$priors$varsel$inprior, c(0.2, 0.2))

  const <- add_priors(bvs_model(d), coef = bvs_coef,
                      varsel = list(inprior = 0.5,
                                    exclude_intercept = FALSE))
  expect_identical(const$priors$varsel$names,
                   c("(Intercept)", "x1", "x2", "z1", "z2"))
  expect_error(add_priors(bvs_model(d), coef = bvs_coef,
                          varsel = list(inprior = 2)),
               "strictly between 0 and 1")
})

test_that("bvs leaves the coefficient prior exactly as it was given", {
  d <- sim_with_noise()
  v <- diag(c(0.5, 1, 2, 3, 4))
  v[2, 3] <- v[3, 2] <- 0.1
  m <- add_priors(bvs_model(d), coef = list(mu = c(0, 1, 0, 0, 0), v_i = v),
                  varsel = list(inprior = 0.5))
  # Unlike SSVS, the selection does not write the precision, so an
  # off-diagonal entry and a non-zero prior mean are both allowed.
  expect_equal(m$priors$B0i, v)
  expect_equal(m$priors$b0, c(0, 1, 0, 0, 0))
})

test_that("bvs needs a proper prior on the coefficients it selects", {
  d <- sim_with_noise()
  expect_error(add_priors(bvs_model(d), coef = list(mu = 0, v_i = 0),
                          varsel = list(inprior = 0.5)),
               "needs a proper prior")
  # A flat prior on a coefficient outside the selection is still allowed.
  expect_silent(add_priors(bvs_model(d),
                           coef = list(mu = 0, v_i = diag(c(0, 1, 1, 1, 1))),
                           varsel = list(inprior = 0.5)))
})

test_that("a prior too wide for the selection to move is reported", {
  d <- sim_with_noise()
  # The package default is vague enough that an excluded coefficient would
  # never be admitted back, which is a stuck chain rather than an answer.
  expect_warning(add_priors(bvs_model(d), varsel = list(inprior = 0.5)),
                 "too wide")
  expect_silent(add_priors(bvs_model(d), coef = bvs_coef,
                           varsel = list(inprior = 0.5)))
  # SSVS is unaffected: its excluded coefficients never leave zero.
  expect_silent(add_priors(create_sfmodel_exp(y ~ x1 + x2 + z1 + z2, data = d,
                                              varsel = "ssvs"),
                           varsel = list(tau = c(0.05, 5))))
})

test_that("an excluded regressor is stored as an exact zero", {
  d <- sim_with_noise(n = 400)
  m <- create_sfmodel_exp(y ~ x1 + x2 + z1 + z2, data = d, varsel = "bvs",
                          iterations = 2000, burnin = 1000)
  m <- add_priors(m, coef = bvs_coef, varsel = list(inprior = 0.5))
  m <- add_posterior_coefficients(add_seed(add_initial_values(m), 42))

  beta <- as.matrix(m$posterior$beta$coeffs)
  inc <- as.matrix(m$posterior$inclusion$coeffs)

  # The indicator and the zero are the same event, which is what makes the
  # stored coefficient the one the frontier was built from.
  for (nm in colnames(inc)) {
    expect_identical(beta[, nm] == 0, inc[, nm] == 0)
  }
  # A coefficient outside the selection is never zeroed.
  expect_true(all(beta[, "(Intercept)"] != 0))

  pip <- colMeans(inc)
  expect_true(all(pip[c("x1", "x2")] > 0.95))
  expect_true(all(pip[c("z1", "z2")] < 0.5))
  expect_equal(unname(coef(m)[c("x1", "x2")]), c(0.5, 0.3), tolerance = 0.05)
})

test_that("the indicators of a relevant regressor do not wander", {
  d <- sim_with_noise(n = 400)
  m <- create_sfmodel_hn(y ~ x1 + x2 + z1, data = d, varsel = "bvs",
                         iterations = 1000, burnin = 500)
  m <- add_priors(m, coef = bvs_coef, varsel = list(inprior = 0.5))
  m <- add_posterior_coefficients(add_seed(add_initial_values(m), 7))

  inc <- as.matrix(m$posterior$inclusion$coeffs)
  expect_true(all(inc[, "x1"] == 1))
  expect_true(all(inc[, "x2"] == 1))
  # The irrelevant one does change state, so the chain is moving rather than
  # merely starting where it stayed.
  expect_true(sum(diff(inc[, "z1"]) != 0) > 0)
})

test_that("the four-component model carries bvs too", {
  set.seed(99)
  d4 <- sim_sf4(n = 40, n_time = 6, beta = c(1, 0.5, 0.3))
  d4$z1 <- stats::rnorm(nrow(d4))

  m <- create_sfmodel4_exp(y ~ x1 + x2 + z1, data = d4, id = "id",
                           varsel = "bvs", iterations = 800, burnin = 400)
  m <- add_priors(m, coef = bvs_coef, varsel = list(inprior = 0.5))
  m <- add_posterior_coefficients(add_seed(add_initial_values(m), 5))

  pip <- colMeans(as.matrix(m$posterior$inclusion$coeffs))
  expect_identical(names(pip), c("x1", "x2", "z1"))
  expect_true(all(pip[c("x1", "x2")] > 0.95))
  expect_true(pip[["z1"]] < 0.5)
  expect_identical(as.matrix(m$posterior$beta$coeffs)[, "z1"] == 0,
                   as.matrix(m$posterior$inclusion$coeffs)[, "z1"] == 0)
})

test_that("summary reports bvs like ssvs", {
  d <- sim_with_noise()
  m <- add_priors(bvs_model(d), coef = bvs_coef,
                  varsel = list(inprior = 0.5))
  m <- add_posterior_coefficients(add_seed(add_initial_values(m), 42))

  tab <- summary(m)$coefficients
  expect_identical(colnames(tab),
                   c("mean", "sd", "2.5%", "median", "97.5%", "PIP", "ESS"))
  expect_equal(tab[c("x1", "x2", "z1", "z2"), "PIP"],
               colMeans(as.matrix(m$posterior$inclusion$coeffs)))
  expect_identical(summary(m)$specifications$varsel, "bvs")
  expect_output(print(summary(m)), "Variable selection: BVS")
  expect_output(print(m), "BVS variable selection")
})

test_that("the run is reproducible from the seed", {
  d <- sim_with_noise()
  run <- function() {
    m <- add_priors(bvs_model(d), coef = bvs_coef,
                    varsel = list(inprior = 0.5))
    add_posterior_coefficients(add_seed(add_initial_values(m), 123))
  }
  a <- run()
  b <- run()
  expect_equal(as.matrix(a$posterior$inclusion$coeffs),
               as.matrix(b$posterior$inclusion$coeffs))
  expect_equal(as.matrix(a$posterior$beta$coeffs),
               as.matrix(b$posterior$beta$coeffs))
})

test_that("the rest of the sampler is untouched by the selection", {
  d <- sim_with_noise(n = 400)

  # With every regressor certain to be included the frontier is the full one
  # in every sweep, so the two runs are sampling one and the same posterior.
  #
  # They cannot be compared draw for draw. The selection draws an indicator
  # per coefficient in every sweep, so the two chains take different paths
  # through the generator and are two independent samples of that posterior
  # rather than one sample repeated. What has to agree is therefore the
  # posterior itself, and the scale on which a gap between two sample means
  # is small or large is the Monte Carlo error of those means.
  #
  # A fixed percentage is not that scale. This test used to allow two per
  # cent, which on its chains was a little over one standard error, so it
  # passed or failed on which way the generator happened to fall. Any change
  # that reshuffled the stream -- including one that altered no arithmetic
  # anyone could see -- could tip it, and one eventually did.
  fit <- function(varsel) {
    m <- create_sfmodel_exp(y ~ x1 + x2, data = d, varsel = varsel,
                            iterations = 4000, burnin = 1000)
    m <- add_priors(m, coef = bvs_coef,
                    varsel = if (is.null(varsel)) NULL else
                      list(inprior = 1 - 1e-12))
    add_posterior_coefficients(add_seed(add_initial_values(m), 21))
  }

  m <- fit("bvs")
  plain <- fit(NULL)

  expect_true(all(as.matrix(m$posterior$inclusion$coeffs) == 1))

  # The standard error of a mean over autocorrelated draws, which is what
  # the effective sample size is for. Six of them, on the difference of two
  # independent means, leaves room for the effective sample size itself
  # being estimated rather than known, and still fails long before a real
  # difference between the two posteriors could hide.
  mcse <- function(x) stats::sd(x) / sqrt(coda::effectiveSize(x))
  within_mcse <- function(a, b) {
    abs(mean(a) - mean(b)) < 6 * sqrt(mcse(a)^2 + mcse(b)^2)
  }

  ba <- as.matrix(m$posterior$beta$coeffs)
  bb <- as.matrix(plain$posterior$beta$coeffs)
  for (j in colnames(ba)) {
    expect_true(within_mcse(ba[, j], bb[, j]),
                label = paste("coefficient", j))
  }

  # The efficiency scores are compared per draw rather than through the
  # single number efficiency() reports, so that the quantity has a spread to
  # be judged against. Its mean is the same number either way, since
  # averaging over draws and over units commutes.
  ea <- rowMeans(exp(-as.matrix(m$posterior$u$coeffs)))
  eb <- rowMeans(exp(-as.matrix(plain$posterior$u$coeffs)))
  expect_true(within_mcse(ea, eb), label = "mean efficiency")
})

test_that("the log-likelihood follows the selected frontier", {
  d <- sim_with_noise()
  m <- add_priors(bvs_model(d), coef = bvs_coef,
                  varsel = list(inprior = 0.5))
  m <- add_posterior_loglik(
    add_posterior_coefficients(add_seed(add_initial_values(m), 9)))
  ll <- as.matrix(m$posterior$loglik)
  expect_identical(ncol(ll), 300L)
  expect_true(all(is.finite(ll)))
})

test_that("prior sensitivity keeps a bvs selection across its fits", {
  d <- sim_with_noise()
  m <- add_priors(bvs_model(d), coef = bvs_coef,
                  varsel = list(inprior = 0.5))
  s <- prior_sensitivity(m, r_star = c(0.7, 0.8), keep_u = TRUE)
  expect_identical(colnames(s$coefficients),
                   c("(Intercept)", "x1", "x2", "z1", "z2"))
  expect_identical(nrow(s$coefficients), 2L)
})
