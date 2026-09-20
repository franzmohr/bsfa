sim_scale <- function(n = 400, gamma = 0.8, sigma_u = 0.3, seed = 42) {
  set.seed(seed)
  z1 <- stats::rnorm(n)
  x1 <- stats::rnorm(n)
  u <- abs(stats::rnorm(n, 0, sigma_u * sqrt(exp(z1 * gamma))))
  y <- 1 + 0.5 * x1 - u + stats::rnorm(n, 0, 0.2)
  data.frame(y = y, x1 = x1, z1 = z1)
}

sim_panel <- function(n_units = 60, n_time = 6, seed = 7) {
  set.seed(seed)
  id <- rep(seq_len(n_units), each = n_time)
  data.frame(y = stats::rnorm(n_units * n_time),
             x1 = stats::rnorm(n_units * n_time),
             zu = stats::rnorm(n_units * n_time),
             ze = rep(stats::rnorm(n_units), each = n_time),
             id = id)
}

fit <- function(object, seed = 1) {
  add_posterior_coefficients(add_seed(add_priors(object), seed))
}

# --- specification -----------------------------------------------------------

test_that("a scale determinant set carries no intercept", {
  d <- sim_scale()
  m <- create_sfmodel_hn(y ~ x1, data = d, scale_u = ~ z1,
                         iterations = 10, burnin = 0)
  expect_identical(m$model$determinants$scale_u$names, "z1")
})

test_that("a set of nothing but an intercept is refused for a scale", {
  d <- sim_scale()
  expect_error(
    create_sfmodel_hn(y ~ x1, data = d, scale_u = ~ 1,
                      iterations = 10, burnin = 0),
    "already plays that part")
})

test_that("the truncated normal defaults to a constant pre-truncation mean", {
  d <- sim_scale()
  m <- create_sfmodel_tn(y ~ x1, data = d, iterations = 10, burnin = 0)
  expect_identical(m$model$determinants$mean_u$names, "(Intercept)")
})

test_that("a mean shift is refused for a family anchored at zero", {
  d <- sim_scale()
  for (f in list(create_sfmodel_hn, create_sfmodel_exp)) {
    expect_error(f(y ~ x1, data = d, mean_u = ~ z1,
                   iterations = 10, burnin = 0),
                 "only the truncated normal family has")
  }
})

test_that("a two-sided or non-formula specification is refused", {
  d <- sim_scale()
  expect_error(create_sfmodel_hn(y ~ x1, data = d, scale_u = y ~ z1,
                                 iterations = 10, burnin = 0),
               "must be one-sided")
  expect_error(create_sfmodel_hn(y ~ x1, data = d, scale_u = "z1",
                                 iterations = 10, burnin = 0),
               "one-sided formula")
})

test_that("a persistent determinant has to be constant within a unit", {
  d <- sim_panel()
  expect_error(
    create_sfmodel4_hn(y ~ x1, data = d, id = "id", scale_eta = ~ zu,
                       iterations = 10, burnin = 0),
    "varies within a unit")
  expect_silent(
    create_sfmodel4_hn(y ~ x1, data = d, id = "id", scale_eta = ~ ze,
                       iterations = 10, burnin = 0))
})

test_that("a persistent set is reduced to one row per unit", {
  d <- sim_panel()
  m <- create_sfmodel4_hn(y ~ x1, data = d, id = "id", scale_eta = ~ ze,
                          scale_u = ~ zu, iterations = 10, burnin = 0)
  expect_identical(nrow(m$model$determinants$scale_eta$Z),
                   m$data$n_units)
  expect_identical(nrow(m$model$determinants$scale_u$Z), m$n)
})

test_that("a rank deficient determinant set is refused", {
  d <- sim_scale()
  d$z2 <- 2 * d$z1
  expect_error(create_sfmodel_hn(y ~ x1, data = d, scale_u = ~ z1 + z2,
                                 iterations = 10, burnin = 0),
               "linear combination")
})

test_that("a determinant missing where the response is not is refused", {
  d <- sim_scale()
  d$z1[3] <- NA
  expect_error(create_sfmodel_hn(y ~ x1, data = d, scale_u = ~ z1,
                                 iterations = 10, burnin = 0),
               "missing values in rows the frontier kept")
})

# --- priors ------------------------------------------------------------------

test_that("a prior is refused for determinants the model does not have", {
  d <- sim_scale()
  m <- create_sfmodel_hn(y ~ x1, data = d, iterations = 10, burnin = 0)
  expect_error(add_priors(m, scale_u = list(mu = 0, v_i = 1)),
               "was given, but the model was not created with")
})

test_that("the determinant prior takes a default and expands scalars", {
  d <- sim_scale()
  m <- add_priors(create_sfmodel_hn(y ~ x1, data = d, scale_u = ~ z1,
                                    iterations = 10, burnin = 0))
  pr <- m$priors$determinants$scale_u
  expect_identical(pr$mu, 0)
  expect_identical(dim(pr$v_i), c(1L, 1L))
})

test_that("an improper determinant prior is refused", {
  d <- sim_scale()
  m <- create_sfmodel_hn(y ~ x1, data = d, scale_u = ~ z1,
                         iterations = 10, burnin = 0)
  expect_error(add_priors(m, scale_u = list(mu = 0, v_i = 0)),
               "must be positive definite")
})

test_that("a determinant prior of the wrong length is refused", {
  d <- sim_scale()
  m <- create_sfmodel_hn(y ~ x1, data = d, scale_u = ~ z1,
                         iterations = 10, burnin = 0)
  expect_error(add_priors(m, scale_u = list(mu = c(0, 0))),
               "one number or one per determinant")
})

# --- sampling ----------------------------------------------------------------

test_that("the scale determinant coefficient is recovered", {
  d <- sim_scale(n = 500, gamma = 0.8)
  m <- fit(create_sfmodel_hn(y ~ x1, data = d, scale_u = ~ z1,
                             iterations = 4000, burnin = 2000))
  band <- stats::quantile(as.matrix(m$posterior$scale_u$coeffs)[, "z1"],
                          c(0.025, 0.975))
  expect_lt(band[1], 0.8)
  expect_gt(band[2], 0.8)
})

test_that("determinants leave the model alone when their coefficient is zero", {
  # Fixing gamma at zero by an overwhelmingly tight prior has to give back the
  # homoskedastic fit, which is the sense in which the parameterisation is a
  # generalisation rather than a different model.
  #
  # The two chains cannot be compared draw for draw even from one seed: the
  # Metropolis block consumes random numbers of its own, so the pinned chain
  # follows a different path through the generator. What has to agree is the
  # posterior they are both drawing from, and the scale on which to judge
  # that is the posterior's own spread rather than a fixed percentage.
  d <- sim_scale(n = 300)
  plain <- fit(create_sfmodel_hn(y ~ x1, data = d,
                                 iterations = 6000, burnin = 2000), seed = 4)
  # Not through fit(), which would add the default priors over the top of
  # the tight one this test is about.
  pinned <- add_posterior_coefficients(add_seed(add_priors(
    create_sfmodel_hn(y ~ x1, data = d, scale_u = ~ z1,
                      iterations = 6000, burnin = 2000),
    scale_u = list(mu = 0, v_i = 1e12)), 4))

  a <- as.matrix(plain$posterior$sigma_u$coeffs)
  b <- as.matrix(pinned$posterior$sigma_u$coeffs)
  expect_lt(abs(mean(a) - mean(b)), 0.3 * stats::sd(a))

  # And the pinned coefficient really did stay where it was put.
  expect_lt(max(abs(as.matrix(pinned$posterior$scale_u$coeffs))), 0.01)
})

test_that("every family accepts scale determinants", {
  d <- sim_scale(n = 250)
  for (f in list(create_sfmodel_hn, create_sfmodel_exp, create_sfmodel_tn)) {
    m <- fit(f(y ~ x1, data = d, scale_u = ~ z1,
               iterations = 600, burnin = 300))
    expect_true("scale_u" %in% names(m$posterior))
    expect_identical(colnames(m$posterior$scale_u$coeffs), "z1")
  }
})

test_that("the pre-truncation mean coefficient is recovered", {
  set.seed(3)
  n <- 500
  z1 <- stats::rnorm(n)
  x1 <- stats::rnorm(n)
  m_u <- 0.2 + 0.45 * z1
  u <- numeric(n)
  for (i in seq_len(n)) {
    repeat {
      v <- stats::rnorm(1, m_u[i], 0.25)
      if (v > 0) {
        u[i] <- v
        break
      }
    }
  }
  d <- data.frame(y = 1 + 0.5 * x1 - u + stats::rnorm(n, 0, 0.2),
                  x1 = x1, z1 = z1)
  m <- fit(create_sfmodel_tn(y ~ x1, data = d, mean_u = ~ z1,
                             iterations = 4000, burnin = 2000), seed = 3)
  band <- stats::quantile(as.matrix(m$posterior$mean_u$coeffs)[, "z1"],
                          c(0.025, 0.975))
  expect_lt(band[1], 0.45)
  expect_gt(band[2], 0.45)
})

test_that("the four-component model takes a set per term", {
  d <- sim_panel(n_units = 50, n_time = 6)
  m <- fit(create_sfmodel4_hn(y ~ x1, data = d, id = "id",
                              scale_eta = ~ ze, scale_u = ~ zu,
                              iterations = 800, burnin = 400), seed = 6)
  expect_identical(colnames(m$posterior$scale_eta$coeffs), "ze")
  expect_identical(colnames(m$posterior$scale_u$coeffs), "zu")
  expect_identical(nrow(as.matrix(m$posterior$scale_eta$coeffs)),
                   nrow(as.matrix(m$posterior$beta$coeffs)))
})

test_that("acceptance rates are reported for every Metropolis block", {
  d <- sim_scale(n = 250)
  m <- fit(create_sfmodel_tn(y ~ x1, data = d, scale_u = ~ z1,
                             iterations = 2000, burnin = 1000))
  expect_setequal(names(m$acceptance), c("scale_u", "mean_u", "par_u"))
  expect_true(all(m$acceptance > 0 & m$acceptance < 1))
})

test_that("a model without determinants reports no acceptance rates", {
  d <- sim_scale(n = 200)
  m <- fit(create_sfmodel_hn(y ~ x1, data = d,
                             iterations = 400, burnin = 200))
  expect_null(m$acceptance)
})

# --- what determinants withdraw ----------------------------------------------

test_that("the pointwise log-likelihood is refused with determinants", {
  d <- sim_scale(n = 200)
  m <- fit(create_sfmodel_hn(y ~ x1, data = d, scale_u = ~ z1,
                             iterations = 400, burnin = 200))
  expect_error(add_posterior_loglik(m), "carries determinants")
})

test_that("the pointwise log-likelihood is refused for the truncated normal", {
  d <- sim_scale(n = 200)
  m <- fit(create_sfmodel_tn(y ~ x1, data = d,
                             iterations = 400, burnin = 200))
  expect_error(add_posterior_loglik(m), "truncated normal family")
})

test_that("summary reports the determinant coefficients", {
  d <- sim_scale(n = 250)
  m <- fit(create_sfmodel_hn(y ~ x1, data = d, scale_u = ~ z1,
                             iterations = 600, burnin = 300))
  expect_true("z1" %in% rownames(summary(m)$coefficients))
})

test_that("efficiency scores are still produced with determinants", {
  d <- sim_scale(n = 250)
  m <- fit(create_sfmodel_hn(y ~ x1, data = d, scale_u = ~ z1,
                             iterations = 600, burnin = 300))
  eff <- efficiency(m)
  expect_identical(nrow(eff), m$n)
  expect_true(all(eff$mean > 0 & eff$mean <= 1))
})
