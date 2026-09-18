test_that("constructors set the class that carries the model variant", {
  d <- sim_sf4(n = 30, n_time = 4, beta = c(1, 0.5))

  m_exp <- create_sfmodel4_exp(y ~ x1, data = d, id = "id")
  m_hn <- create_sfmodel4_hn(y ~ x1, data = d, id = "id")

  expect_s3_class(m_exp, "sfmodel4_exp")
  expect_s3_class(m_exp, "sfmodel4")
  expect_s3_class(m_exp, "sfmodel")
  expect_equal(m_exp$model$components, 4L)

  expect_equal(m_exp$model$par_eta_name, "lambda_eta")
  expect_equal(m_exp$model$par_u_name, "lambda_u")
  expect_equal(m_hn$model$par_eta_name, "sigma_eta")
  expect_equal(m_hn$model$par_u_name, "sigma_u")
})

test_that("the model refuses data it cannot separate", {
  d <- sim_sf4(n = 30, n_time = 4, beta = c(1, 0.5))

  expect_error(create_sfmodel4_exp(y ~ x1, data = d), "needs an 'id'")

  # One unit observed once cannot contribute to the split.
  d$id[1] <- 999
  expect_error(create_sfmodel4_exp(y ~ x1, data = d, id = "id"),
               "at least two observations")
})

test_that("both one-sided terms are elicited separately", {
  d <- sim_sf4(n = 30, n_time = 4, beta = c(1, 0.5))

  m <- add_priors(create_sfmodel4_exp(y ~ x1, data = d, id = "id"),
                  lambda_eta = list(r_star = 0.95),
                  lambda_u = list(r_star = 0.8))

  expect_equal(m$priors$rate_eta, -log(0.95))
  expect_equal(m$priors$rate_u, -log(0.8))
  expect_equal(m$priors$r_star, c(persistent = 0.95, transient = 0.8))

  h <- add_priors(create_sfmodel4_hn(y ~ x1, data = d, id = "id"),
                  sigma_eta = list(r_star = 0.95))
  expect_equal(h$priors$rate_eta / (h$priors$shape_eta - 1),
               (-log(0.95) / qnorm(0.75))^2)
})

test_that("four-component prior arguments are validated", {
  d <- sim_sf4(n = 30, n_time = 4, beta = c(1, 0.5))
  m <- create_sfmodel4_exp(y ~ x1, data = d, id = "id")

  expect_error(add_priors(m, lambda_eta = list(rstar = 0.9)),
               "Unknown element")
  expect_error(add_priors(m, lambda_u = list(r_star = 0)), "between 0 and 1")
  expect_error(add_priors(m, sigma_mu = list(shape = 0)), "must be positive")
  expect_error(add_priors(create_sfmodel4_hn(y ~ x1, data = d, id = "id"),
                          sigma_eta = list(shape = 1)), "must exceed 1")
})

test_that("starting values split the least squares residual variance", {
  d <- sim_sf4(n = 30, n_time = 4, beta = c(1, 0.5))
  m <- add_initial_values(add_priors(
    create_sfmodel4_exp(y ~ x1, data = d, id = "id")))

  s2 <- summary(lm(y ~ x1, data = d))$sigma^2
  expect_equal(m$initial$sigma_v2, s2 / 2)
  expect_equal(m$initial$sigma_mu2, s2 / 4)

  expect_equal(m$initial$par_eta, m$priors$shape_eta / m$priors$rate_eta)
  expect_equal(m$initial$par_u, m$priors$shape_u / m$priors$rate_u)

  # One latent term per unit for the first two, one per observation for the
  # transient one.
  expect_length(m$initial$mu, 30L)
  expect_length(m$initial$eta, 30L)
  expect_length(m$initial$u, 120L)
})

test_that("draws are attached as the documented blocks", {
  set.seed(41)
  d <- sim_sf4(n = 40, n_time = 5, beta = c(1, 0.5, 0.3))
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel4_exp(y ~ x1 + x2, data = d, id = "id",
                        iterations = 200, burnin = 100)))

  expect_s3_class(est, "sfmodel4_exp")
  expect_equal(dim(est$posterior$beta$coeffs), c(200L, 3L))
  expect_equal(dim(est$posterior$sigma_mu$coeffs), c(200L, 1L))
  expect_equal(dim(est$posterior$lambda_eta$coeffs), c(200L, 1L))
  expect_equal(dim(est$posterior$lambda_u$coeffs), c(200L, 1L))

  expect_equal(dim(est$posterior$mu$coeffs), c(200L, 40L))
  expect_equal(dim(est$posterior$eta$coeffs), c(200L, 40L))
  expect_equal(dim(est$posterior$u$coeffs), c(200L, 200L))

  expect_true(all(est$posterior$eta$coeffs >= 0))
  expect_true(all(est$posterior$u$coeffs >= 0))
  # The unit effect is unrestricted, which is what separates it from eta.
  expect_true(any(est$posterior$mu$coeffs < 0))
})

test_that("the sampler recovers the frontier and the components", {
  set.seed(4242)
  beta <- c(1, 0.5, 0.3)
  # The unit effect and persistent inefficiency are given comparable spread.
  # Where one dominates the other, the split between them is weakly determined
  # even though the model is correct; see the test below.
  d <- sim_sf4(n = 200, n_time = 10, beta = beta, sigma_v = 0.2,
               sigma_mu = 0.15, par_eta = 4, par_u = 6)

  est <- add_posterior_coefficients(add_priors(
    create_sfmodel4_exp(y ~ x1 + x2, data = d, id = "id",
                        iterations = 3000, burnin = 1000)))

  pm <- colMeans(est$posterior$beta$coeffs)
  expect_equal(unname(pm[2:3]), beta[2:3], tolerance = 0.05)
  expect_equal(mean(est$posterior$sigma_v$coeffs), 0.2, tolerance = 0.1)

  # What the model is for is the split, so that is what is checked: each set of
  # scores has to track the component it is meant to measure.
  expect_gt(cor(efficiency(est, "persistent")$mean, attr(d, "persistent")),
            0.7)
  expect_gt(cor(efficiency(est, "transient")$mean, attr(d, "transient")), 0.5)
  expect_gt(cor(efficiency(est, "overall")$mean, attr(d, "overall")), 0.6)

  # The unit effect must not be absorbed into the efficiency scores.
  expect_gt(cor(colMeans(est$posterior$mu$coeffs), attr(d, "mu")), 0.5)
})

test_that("the identified combination is recovered whatever the split", {
  # mu_i - eta_i is what the data determine; how it is divided between the two
  # depends on which of them dominates. Both settings below must recover the
  # combination sharply, and they differ in which component comes out well.
  fit <- function(sigma_mu, par_eta) {
    set.seed(4242)
    d <- sim_sf4(n = 200, n_time = 10, beta = c(1, 0.5, 0.3), sigma_v = 0.2,
                 sigma_mu = sigma_mu, par_eta = par_eta, par_u = 6)
    est <- add_posterior_coefficients(add_priors(
      create_sfmodel4_exp(y ~ x1 + x2, data = d, id = "id",
                          iterations = 3000, burnin = 1000)))
    mu_hat <- colMeans(est$posterior$mu$coeffs)
    eta_hat <- colMeans(est$posterior$eta$coeffs)
    list(combination = cor(mu_hat - eta_hat, attr(d, "mu") - attr(d, "eta")),
         mu = cor(mu_hat, attr(d, "mu")),
         eta = cor(eta_hat, attr(d, "eta")))
  }

  wide_mu <- fit(sigma_mu = 0.30, par_eta = 8)
  wide_eta <- fit(sigma_mu = 0.05, par_eta = 4)

  expect_gt(wide_mu$combination, 0.95)
  expect_gt(wide_eta$combination, 0.95)

  # A dominant unit effect is pinned down and the persistent term is not; with
  # a small unit effect it is the other way round.
  expect_gt(wide_mu$mu, wide_mu$eta)
  expect_gt(wide_eta$eta, wide_eta$mu)
})

test_that("the unit effect keeps heterogeneity out of the scores", {
  # With a large unit effect and no persistent inefficiency, a two-component
  # panel model books the heterogeneity as inefficiency while the
  # four-component model does not.
  set.seed(77)
  d <- sim_sf4(n = 120, n_time = 8, beta = c(1, 0.5), sigma_v = 0.15,
               sigma_mu = 0.6, par_eta = 200, par_u = 6)

  four <- add_posterior_coefficients(add_priors(
    create_sfmodel4_exp(y ~ x1, data = d, id = "id",
                        iterations = 2000, burnin = 1000)))
  two <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, id = "id",
                       iterations = 2000, burnin = 1000)))

  expect_gt(mean(efficiency(four, "persistent")$mean),
            mean(efficiency(two)$mean))
})

test_that("the three efficiency types have the right shape", {
  set.seed(43)
  d <- sim_sf4(n = 40, n_time = 5, beta = c(1, 0.5))
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel4_exp(y ~ x1, data = d, id = "id",
                        iterations = 200, burnin = 100)))

  expect_equal(nrow(efficiency(est, "persistent")), 40L)
  expect_equal(nrow(efficiency(est, "transient")), 200L)
  expect_equal(nrow(efficiency(est, "overall")), 200L)
  # The default is the overall score.
  expect_equal(efficiency(est), efficiency(est, "overall"))

  # Overall efficiency is the product of the two parts, so it cannot exceed
  # either of them.
  ov <- efficiency(est, "overall")$mean
  tr <- efficiency(est, "transient")$mean
  expect_true(all(ov <= tr + 1e-12))

  expect_error(efficiency(est, "nonsense"), "should be one of")
})

test_that("the cost frontier flips both one-sided terms", {
  set.seed(44)
  d <- sim_sf4(n = 60, n_time = 6, beta = c(1, 0.5), type = "cost")

  est <- add_posterior_coefficients(add_priors(
    create_sfmodel4_exp(y ~ x1, data = d, id = "id", type = "cost",
                        iterations = 1000, burnin = 500)))
  expect_equal(unname(colMeans(est$posterior$beta$coeffs)[2]), 0.5,
               tolerance = 0.05)

  wrong <- add_posterior_coefficients(add_priors(
    create_sfmodel4_exp(y ~ x1, data = d, id = "id", type = "production",
                        iterations = 1000, burnin = 500)))
  expect_gt(mean(wrong$posterior$sigma_v$coeffs),
            mean(est$posterior$sigma_v$coeffs))
})

test_that("summary and plot cover the extra components", {
  set.seed(45)
  d <- sim_sf4(n = 40, n_time = 5, beta = c(1, 0.5))
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel4_hn(y ~ x1, data = d, id = "id",
                       iterations = 200, burnin = 100)))

  tab <- summary(est)$coefficients
  expect_equal(rownames(tab),
               c("(Intercept)", "x1", "sigma_v", "sigma_mu", "sigma_eta",
                 "sigma_u"))
  # No single signal-to-noise ratio is defined with two one-sided terms.
  expect_false("lambda" %in% rownames(tab))

  expect_output(print(summary(est)), "Four-component")
  expect_output(print(est), "persistent and transient")

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  for (type in c("hist", "trace", "boxplot", "efficiency")) {
    expect_silent(plot(est, type = type))
  }
})

test_that("the log-likelihood and the criteria are refused", {
  set.seed(46)
  d <- sim_sf4(n = 30, n_time = 4, beta = c(1, 0.5))
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel4_exp(y ~ x1, data = d, id = "id",
                        iterations = 100, burnin = 50)))

  expect_error(add_posterior_loglik(est), "not available for the four")
  expect_error(selection_criteria(est), "not available for the four-component")
})

test_that("the stored seed fixes the draws", {
  set.seed(47)
  d <- sim_sf4(n = 30, n_time = 4, beta = c(1, 0.5))
  m <- add_seed(add_priors(
    create_sfmodel4_exp(y ~ x1, data = d, id = "id",
                        iterations = 100, burnin = 50)), 909)

  expect_equal(add_posterior_coefficients(m)$posterior$beta$coeffs,
               add_posterior_coefficients(m)$posterior$beta$coeffs)
})

test_that("sim_sf4 returns the components it generated", {
  set.seed(48)
  d <- sim_sf4(n = 25, n_time = 4, beta = c(1, 0.5, 0.3))

  expect_equal(nrow(d), 100L)
  expect_length(attr(d, "mu"), 25L)
  expect_length(attr(d, "eta"), 25L)
  expect_length(attr(d, "u"), 100L)
  expect_equal(attr(d, "overall"),
               attr(d, "persistent")[d$id] * attr(d, "transient"))
  expect_error(sim_sf4(n = 10, n_time = 1), "at least 2")
})
