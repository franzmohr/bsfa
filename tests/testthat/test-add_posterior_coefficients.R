test_that("draws are attached to the model object as mcmc blocks", {
  set.seed(1)
  d <- sim_sf(n = 80, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1 + x2, data = d,
                       iterations = 200, burnin = 100, thin = 2)))

  # The object keeps its class; estimation adds to it rather than replacing it.
  expect_s3_class(est, "sfmodel_exp")
  expect_s3_class(est, "sfmodel")

  expect_s3_class(est$posterior$beta$coeffs, "mcmc")
  expect_equal(dim(est$posterior$beta$coeffs), c(100L, 3L))
  expect_equal(colnames(est$posterior$beta$coeffs),
               c("(Intercept)", "x1", "x2"))
  expect_equal(coda::mcpar(est$posterior$beta$coeffs), c(102, 300, 2))

  expect_equal(dim(est$posterior$sigma_v$coeffs), c(100L, 1L))
  expect_equal(dim(est$posterior$lambda$coeffs), c(100L, 1L))
  expect_null(est$posterior$sigma_u)
  expect_equal(dim(est$posterior$u$coeffs), c(100L, 80L))
  expect_true(all(est$posterior$u$coeffs >= 0))
})

test_that("the half-normal model names its own parameter block", {
  set.seed(2)
  d <- sim_sf(n = 80, beta = c(1, 0.5), sigma_v = 0.2, par_u = 0.3,
              ineff = "halfnormal")
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel_hn(y ~ x1, data = d, iterations = 200, burnin = 100)))

  expect_s3_class(est, "sfmodel_hn")
  expect_equal(dim(est$posterior$sigma_u$coeffs), c(200L, 1L))
  expect_null(est$posterior$lambda)
  # The signal-to-noise ratio is reported for the half-normal model only.
  expect_true("lambda" %in% rownames(summary(est)$coefficients))
})

test_that("keep_u = FALSE omits the inefficiency block", {
  set.seed(4)
  d <- sim_sf(n = 80, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 100, burnin = 50)),
    keep_u = FALSE)

  expect_null(est$posterior$u)
  expect_error(efficiency(est), "keep_u")
})

test_that("the exponential model recovers its parameters", {
  set.seed(123)
  beta <- c(1, 0.5, 0.3)
  d <- sim_sf(n = 1500, beta = beta, sigma_v = 0.2, par_u = 4)
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1 + x2, data = d,
                       iterations = 2000, burnin = 1000)))

  pm <- colMeans(est$posterior$beta$coeffs)
  expect_equal(unname(pm[2:3]), beta[2:3], tolerance = 0.05)
  # The intercept is only weakly separated from the mean of the one-sided term,
  # so it is checked against a looser bound than the slopes.
  expect_equal(unname(pm[1]), beta[1], tolerance = 0.15)
  expect_equal(mean(est$posterior$sigma_v$coeffs), 0.2, tolerance = 0.1)
  expect_equal(mean(est$posterior$lambda$coeffs), 4, tolerance = 1.5)
})

test_that("the half-normal model recovers its parameters", {
  set.seed(321)
  d <- sim_sf(n = 1500, beta = c(1, 0.5), sigma_v = 0.2, par_u = 0.3,
              ineff = "halfnormal")
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel_hn(y ~ x1, data = d, iterations = 2000, burnin = 1000)))

  expect_equal(unname(colMeans(est$posterior$beta$coeffs)[2]), 0.5,
               tolerance = 0.05)
  expect_equal(mean(est$posterior$sigma_u$coeffs), 0.3, tolerance = 0.15)
})

test_that("the cost frontier flips the sign of the one-sided term", {
  set.seed(7)
  d <- sim_sf(n = 800, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4,
              type = "cost")

  est <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, type = "cost",
                       iterations = 1000, burnin = 500)))
  expect_equal(unname(colMeans(est$posterior$beta$coeffs)[2]), 0.5,
               tolerance = 0.05)

  # Estimating the same data as a production frontier must inflate the
  # residual variance, because the one-sided term is then pushed the wrong way.
  wrong <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, type = "production",
                       iterations = 1000, burnin = 500)))
  expect_gt(mean(wrong$posterior$sigma_v$coeffs),
            mean(est$posterior$sigma_v$coeffs))
})

test_that("efficiency scores track the simulated truth", {
  set.seed(99)
  d <- sim_sf(n = 400, beta = c(1, 0.5), sigma_v = 0.1, par_u = 4)
  est <- add_posterior_coefficients(add_priors(
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
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1 + x2, data = d, id = "id",
                       iterations = 500, burnin = 250)))

  expect_equal(est$n, 500L)
  expect_equal(est$data$n_units, 100L)
  expect_equal(ncol(est$posterior$u$coeffs), 100L)
  expect_gt(cor(colMeans(exp(-est$posterior$u$coeffs)),
                attr(d, "efficiency")), 0.8)
})

test_that("the stored seed fixes the draws", {
  set.seed(3)
  d <- sim_sf(n = 100, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_seed(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 200, burnin = 100)), 777)

  a <- add_posterior_coefficients(m)
  b <- add_posterior_coefficients(m)
  expect_equal(a$posterior$beta$coeffs, b$posterior$beta$coeffs)

  # add_seed keeps the seed the model already had out of it, and
  # add_initial_values does not overwrite one that is already there.
  expect_equal(add_initial_values(m)$model$seed, 777)
})

test_that("set.seed after add_initial_values does not change the draws", {
  d <- sim_sf(n = 100, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_initial_values(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 200, burnin = 100)))

  set.seed(1)
  a <- add_posterior_coefficients(m)
  set.seed(2)
  b <- add_posterior_coefficients(m)
  expect_equal(a$posterior$beta$coeffs, b$posterior$beta$coeffs)
})

test_that("the generator is left as it was found", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_seed(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 100, burnin = 50)), 1)

  set.seed(4321)
  before <- .Random.seed
  invisible(add_posterior_coefficients(m))
  expect_equal(.Random.seed, before)
})

test_that("starting from the prior reaches the same posterior", {
  set.seed(17)
  d <- sim_sf(n = 600, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d,
                                     iterations = 2000, burnin = 1000))

  a <- add_posterior_coefficients(add_initial_values(m, method = "ols"))
  b <- add_posterior_coefficients(add_initial_values(m, method = "prior"))

  expect_equal(colMeans(a$posterior$beta$coeffs),
               colMeans(b$posterior$beta$coeffs), tolerance = 0.05)
})

test_that("an alternative sampler can be supplied", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d,
                                     iterations = 100, burnin = 50))

  fake <- function(object) {
    object$posterior <- list(beta = list(coeffs = "supplied elsewhere"))
    object
  }
  est <- add_posterior_coefficients(m, posterior_function = fake)

  expect_equal(est$posterior$beta$coeffs, "supplied elsewhere")
  expect_s3_class(est, "sfmodel_exp")
})

test_that("changing the specification discards draws that no longer match", {
  set.seed(21)
  d <- sim_sf(n = 60, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 100, burnin = 50)))
  expect_false(is.null(est$posterior))

  expect_message(re_prior <- add_priors(est, lambda = list(r_star = 0.5)),
                 "Dropping posterior draws")
  expect_null(re_prior$posterior)
  expect_equal(re_prior$priors$rate_u, -log(0.5))

  expect_message(re_init <- add_initial_values(est), "Dropping posterior")
  expect_null(re_init$posterior)

  # Setting the specification on a model that has no draws says nothing.
  fresh <- create_sfmodel_exp(y ~ x1, data = d, iterations = 100, burnin = 50)
  expect_silent(add_initial_values(add_priors(fresh)))
})

test_that("a posterior from elsewhere is checked before it is used", {
  # posterior_function is a documented extension point, so a posterior that
  # does not carry the blocks the package writes has to be reported as such
  # rather than surfacing later as a non-conformable matrix.
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d,
                                     iterations = 100, burnin = 50))

  only_beta <- add_posterior_coefficients(m, posterior_function = function(x) {
    x$posterior <- list(beta = list(coeffs = matrix(0, 100, 2)))
    x
  })
  expect_error(summary(only_beta), "no usable draws in posterior\\$sigma_v")
  expect_error(add_posterior_loglik(only_beta), "no usable draws")

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_error(plot(only_beta), "no usable draws")

  wrong_width <- add_posterior_coefficients(m, posterior_function = function(x) {
    x$posterior <- list(beta = list(coeffs = matrix(0, 100, 5)),
                        sigma_v = list(coeffs = matrix(1, 100, 1)),
                        lambda = list(coeffs = matrix(1, 100, 1)))
    x
  })
  expect_error(summary(wrong_width), "5 column\\(s\\) but the model has 2")

  ragged <- add_posterior_coefficients(m, posterior_function = function(x) {
    x$posterior <- list(beta = list(coeffs = matrix(0, 100, 2)),
                        sigma_v = list(coeffs = matrix(1, 50, 1)),
                        lambda = list(coeffs = matrix(1, 100, 1)))
    x
  })
  expect_error(summary(ragged), "different numbers of draws")
})

test_that("priors must be added before simulating", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- create_sfmodel_exp(y ~ x1, data = d, iterations = 100, burnin = 50)

  expect_error(add_posterior_coefficients(m), "Add priors before simulating")
  # Initial values, unlike priors, are filled in with their default.
  expect_false(is.null(add_posterior_coefficients(add_priors(m))$initial))
})

test_that("summary requires draws and a valid credible band", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d,
                                     iterations = 100, burnin = 50))

  expect_error(summary(m), "does not contain posterior draws")

  est <- add_posterior_coefficients(m)
  expect_s3_class(summary(est), "summary.sfmodel")
  expect_error(summary(est, ci = 1.5), "between 0 and 1")
  expect_equal(colnames(summary(est, ci = 0.9)$coefficients),
               c("mean", "sd", "5%", "median", "95%", "ESS"))
})

test_that("thinning keeps the last sweep of each block", {
  set.seed(61)
  d <- sim_sf(n = 80, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)

  fit <- function(iterations, burnin, thin) {
    add_posterior_coefficients(add_seed(add_priors(
      create_sfmodel_exp(y ~ x1, data = d, iterations = iterations,
                         burnin = burnin, thin = thin)), 4321))
  }

  # Ten sweeps thinned to one draw must return the tenth sweep, not the first.
  thinned <- fit(iterations = 10, burnin = 0, thin = 10)
  first <- fit(iterations = 1, burnin = 0, thin = 1)
  last <- fit(iterations = 1, burnin = 9, thin = 1)

  expect_equal(nrow(as.matrix(thinned$posterior$beta$coeffs)), 1L)
  expect_equal(as.numeric(thinned$posterior$beta$coeffs),
               as.numeric(last$posterior$beta$coeffs))
  expect_false(isTRUE(all.equal(as.numeric(thinned$posterior$beta$coeffs),
                                as.numeric(first$posterior$beta$coeffs))))
})

test_that("the iteration index names the sweeps the draws came from", {
  set.seed(62)
  d <- sim_sf(n = 80, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)

  est <- add_posterior_coefficients(add_seed(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 20, burnin = 10,
                       thin = 5)), 4321))

  expect_equal(as.numeric(stats::time(est$posterior$beta$coeffs)),
               c(15, 20, 25, 30))

  # The draw labelled 30 is the state after the thirtieth sweep, which is the
  # last one the sampler runs.
  full <- add_posterior_coefficients(add_seed(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 20, burnin = 10)), 4321))
  expect_equal(as.numeric(tail(as.matrix(est$posterior$beta$coeffs), 1)),
               as.numeric(tail(as.matrix(full$posterior$beta$coeffs), 1)))
})

test_that("keep_u takes a thinning interval for the augmented draws", {
  set.seed(65)
  d <- sim_sf(n = 40, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  fit <- function(keep_u) {
    add_posterior_coefficients(add_seed(add_priors(
      create_sfmodel_exp(y ~ x1, data = d, iterations = 200,
                         burnin = 100)), 8080), keep_u = keep_u)
  }

  full <- fit(TRUE)
  thinned <- fit(4)

  expect_equal(nrow(as.matrix(full$posterior$u$coeffs)), 200L)
  expect_equal(nrow(as.matrix(thinned$posterior$u$coeffs)), 50L)

  # The scalar blocks are untouched by it, and the draws that are kept are the
  # very draws the unthinned run kept.
  expect_equal(nrow(as.matrix(thinned$posterior$beta$coeffs)), 200L)
  expect_equal(as.matrix(thinned$posterior$u$coeffs),
               as.matrix(full$posterior$u$coeffs)[seq(4, 200, by = 4), ])

  # Their iteration index counts in the two intervals multiplied.
  expect_equal(as.numeric(stats::time(thinned$posterior$u$coeffs)),
               seq(104, 300, by = 4))
  expect_equal(coda::thin(thinned$posterior$u$coeffs), 4)

  expect_equal(efficiency(thinned)$unit, efficiency(full)$unit)
  expect_null(fit(FALSE)$posterior$u)
})

test_that("keep_u is validated", {
  d <- sim_sf(n = 40, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d, iterations = 20,
                                     burnin = 5))

  # NA used to be taken for TRUE, since Rcpp reads NA_LOGICAL as true.
  for (bad in list(NA, -1, 1.5, "yes", c(TRUE, FALSE))) {
    expect_error(add_posterior_coefficients(m, keep_u = bad), "'keep_u'")
  }
})

test_that("the four-component model thins all three augmented blocks", {
  set.seed(66)
  d <- sim_sf4(n = 15, n_time = 4, beta = c(1, 0.5))
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel4_exp(y ~ x1, data = d, id = "id", iterations = 200,
                        burnin = 100)), keep_u = 5)

  for (b in c("mu", "eta", "u")) {
    expect_equal(nrow(as.matrix(est$posterior[[b]]$coeffs)), 40L)
  }
  expect_equal(nrow(as.matrix(est$posterior$beta$coeffs)), 200L)
  expect_equal(nrow(efficiency(est, "transient")), 60L)
  expect_equal(nrow(efficiency(est, "persistent")), 15L)
})

test_that("a keep_u interval wider than the chain is refused", {
  d <- sim_sf(n = 40, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d, iterations = 100,
                                     burnin = 10))

  # Storing every 101st of 100 draws stores none, which used to come back as
  # a missing block and sent efficiency() on to advise keep_u = TRUE.
  expect_error(add_posterior_coefficients(m, keep_u = 101),
               "retains only 100 draws")
  expect_error(add_posterior_coefficients(m, keep_u = 500), "'keep_u' is 500")

  # The boundary itself is fine and keeps the last draw.
  at_limit <- add_posterior_coefficients(m, keep_u = 100)
  expect_equal(nrow(as.matrix(at_limit$posterior$u$coeffs)), 1L)

  # Thinning counts retained draws, so it is 'iterations / thin' that binds.
  thinned <- add_priors(create_sfmodel_exp(y ~ x1, data = d,
                                           iterations = 100, burnin = 10,
                                           thin = 4))
  expect_error(add_posterior_coefficients(thinned, keep_u = 26),
               "retains only 25 draws")
  expect_equal(nrow(as.matrix(
    add_posterior_coefficients(thinned, keep_u = 25)$posterior$u$coeffs)), 1L)

  m4 <- add_priors(create_sfmodel4_exp(
    y ~ x1, data = sim_sf4(n = 10, n_time = 3, beta = c(1, 0.5)),
    id = "id", iterations = 100, burnin = 10))
  expect_error(add_posterior_coefficients(m4, keep_u = 101),
               "retains only 100 draws")
})

test_that("progress is reported at the requested interval", {
  d <- sim_sf(n = 30, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d, iterations = 20,
                                     burnin = 10))

  expect_output(add_posterior_coefficients(m, verbose = 10),
                "Iteration 10 of 30")
  expect_output(add_posterior_coefficients(m, verbose = TRUE),
                "Iteration 30 of 30")
  expect_silent(add_posterior_coefficients(m, verbose = FALSE))

  m4 <- add_priors(create_sfmodel4_exp(
    y ~ x1, data = sim_sf4(n = 8, n_time = 3, beta = c(1, 0.5)),
    id = "id", iterations = 20, burnin = 10))
  expect_output(add_posterior_coefficients(m4, verbose = 10),
                "Iteration 30 of 30")
})

test_that("the size of the augmented draws is warned about in advance", {
  # The warning has to fire before anything is allocated, so the helper is
  # called directly rather than through a model nobody could hold in memory.
  expect_warning(warn_augmented_size(2e4, 2e4), "about 3.0 GB")
  expect_warning(warn_augmented_size(2e4, 2e4), "keep_u")
  expect_silent(warn_augmented_size(100, 1000))
  expect_silent(warn_augmented_size(0, 0))
})

test_that("verbose is validated", {
  d <- sim_sf(n = 40, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d, iterations = 20,
                                     burnin = 2))

  # "yes" used to become NA_integer_ and quietly disable reporting; 2.7 was
  # truncated to 2 without a word.
  for (bad in list(NA, -1, 1.5, "yes", c(TRUE, FALSE))) {
    expect_error(add_posterior_coefficients(m, verbose = bad), "'verbose'")
  }
  expect_silent(add_posterior_coefficients(m, verbose = 0))
})
