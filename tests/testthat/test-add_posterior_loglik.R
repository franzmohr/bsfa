test_that("the exponential log-likelihood matches numerical integration", {
  set.seed(8)
  d <- sim_sf(n = 60, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  est <- add_posterior_loglik(add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 50, burnin = 50))))

  expect_s3_class(est$posterior$loglik, "mcmc")
  expect_equal(dim(est$posterior$loglik), c(50L, 60L))

  # The composed error of a production frontier is v - u, so its density is the
  # convolution of the normal and the exponential. Integrating it directly is a
  # slow but assumption-free check on the closed form the package uses.
  draw <- 1
  beta <- as.matrix(est$posterior$beta$coeffs)[draw, ]
  sv <- as.numeric(est$posterior$sigma_v$coeffs)[draw]
  lam <- as.numeric(est$posterior$lambda$coeffs)[draw]
  e <- as.numeric(est$data$y - est$data$X %*% beta)

  for (i in c(1L, 17L, 42L)) {
    num <- stats::integrate(
      function(u) dexp(u, rate = lam) * dnorm(e[i] + u, sd = sv),
      lower = 0, upper = Inf)$value
    expect_equal(est$posterior$loglik[draw, i], log(num), tolerance = 1e-6)
  }
})

test_that("the half-normal log-likelihood matches numerical integration", {
  set.seed(9)
  d <- sim_sf(n = 60, beta = c(1, 0.5), sigma_v = 0.2, par_u = 0.3,
              ineff = "halfnormal")
  est <- add_posterior_loglik(add_posterior_coefficients(add_priors(
    create_sfmodel_hn(y ~ x1, data = d, iterations = 50, burnin = 50))))

  draw <- 3
  beta <- as.matrix(est$posterior$beta$coeffs)[draw, ]
  sv <- as.numeric(est$posterior$sigma_v$coeffs)[draw]
  su <- as.numeric(est$posterior$sigma_u$coeffs)[draw]
  e <- as.numeric(est$data$y - est$data$X %*% beta)

  for (i in c(2L, 23L, 55L)) {
    num <- stats::integrate(
      function(u) 2 * dnorm(u, sd = su) * dnorm(e[i] + u, sd = sv),
      lower = 0, upper = Inf)$value
    expect_equal(est$posterior$loglik[draw, i], log(num), tolerance = 1e-6)
  }
})

test_that("the cost frontier reverses the sign of the one-sided term", {
  set.seed(10)
  d <- sim_sf(n = 40, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4,
              type = "cost")
  est <- add_posterior_loglik(add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, type = "cost",
                       iterations = 50, burnin = 50))))

  draw <- 5
  beta <- as.matrix(est$posterior$beta$coeffs)[draw, ]
  sv <- as.numeric(est$posterior$sigma_v$coeffs)[draw]
  lam <- as.numeric(est$posterior$lambda$coeffs)[draw]
  e <- as.numeric(est$data$y - est$data$X %*% beta)

  num <- stats::integrate(
    function(u) dexp(u, rate = lam) * dnorm(e[7] - u, sd = sv),
    lower = 0, upper = Inf)$value
  expect_equal(est$posterior$loglik[draw, 7], log(num), tolerance = 1e-6)
})

test_that("the log-likelihood inherits the iteration index of the draws", {
  set.seed(12)
  d <- sim_sf(n = 40, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  est <- add_posterior_loglik(add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d,
                       iterations = 100, burnin = 50, thin = 2))))

  expect_equal(coda::mcpar(est$posterior$loglik),
               coda::mcpar(est$posterior$beta$coeffs))
})

test_that("log-likelihood requires draws and rejects panel models", {
  d <- sim_sf(n = 40, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4, n_time = 2)

  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d,
                                     iterations = 50, burnin = 50))
  expect_error(add_posterior_loglik(m), "does not contain posterior draws")

  mp <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, id = "id",
                       iterations = 50, burnin = 50)))
  expect_error(add_posterior_loglik(mp), "not available for panel models")
})

test_that("printing reports that the log-likelihood is present", {
  set.seed(13)
  d <- sim_sf(n = 40, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 50, burnin = 50)))

  expect_output(print(est), "posterior      simulated")
  expect_output(print(add_posterior_loglik(est)), "with loglik")
})
