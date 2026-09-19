test_that("the exponential prior reproduces the target median efficiency", {
  # With shape 1 and rate c = -log(r*) on lambda, the marginal prior of u is
  # c / (u + c)^2, whose median is c. So the prior median of exp(-u) is r*.
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)

  for (r_star in c(0.5, 0.75, 0.9)) {
    m <- add_priors(create_sfmodel_exp(y ~ x1, data = d),
                    lambda = list(r_star = r_star))
    expect_equal(m$priors$rate_u, -log(r_star))
    expect_equal(m$priors$shape_u, 1)

    set.seed(42)
    lambda <- rgamma(2e5, shape = m$priors$shape_u, rate = m$priors$rate_u)
    u <- rexp(length(lambda), rate = lambda)
    expect_equal(median(exp(-u)), r_star, tolerance = 0.01)
  }
})

test_that("the half-normal prior reproduces the target median efficiency", {
  # The prior sits on sigma_u^2 and is inverse gamma, so u is marginally half
  # t on 2 * shape degrees of freedom with scale sqrt(rate / shape), and the
  # rate is set from its median. An earlier version matched the prior mean of
  # sigma_u^2 instead, which left the anchor optimistic by a factor that
  # depended on the shape alone: at the default of 2.5 asking for 0.75 gave a
  # prior median efficiency of 0.787.
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)

  for (r_star in c(0.5, 0.75, 0.9)) {
    for (shape in c(1.5, 2.5, 5)) {
      m <- add_priors(create_sfmodel_hn(y ~ x1, data = d),
                      sigma_u = list(r_star = r_star, shape = shape))
      expect_equal(m$priors$shape_u, shape)
      expect_equal(m$priors$rate_u,
                   shape * (-log(r_star) / qt(0.75, df = 2 * shape))^2)

      # The median of the marginal prior of exp(-u) is r_star itself, which is
      # what the anchor claims and is the whole point of the elicitation.
      set.seed(42)
      s2 <- 1 / rgamma(2e5, shape = m$priors$shape_u, rate = m$priors$rate_u)
      u <- sqrt(s2) * abs(rnorm(length(s2)))
      expect_equal(median(exp(-u)), r_star, tolerance = 0.01)
    }
  }
})

test_that("the two methods elicit the same anchor differently", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)

  e <- add_priors(create_sfmodel_exp(y ~ x1, data = d))
  h <- add_priors(create_sfmodel_hn(y ~ x1, data = d))

  expect_equal(e$priors$r_star, h$priors$r_star)
  expect_false(isTRUE(all.equal(e$priors$rate_u, h$priors$rate_u)))
  # The blocks that do not depend on the distribution must agree.
  expect_equal(e$priors$b0, h$priors$b0)
  expect_equal(e$priors$B0i, h$priors$B0i)
  expect_equal(e$priors$shape_v, h$priors$shape_v)
})

test_that("coefficient priors accept scalars, vectors and matrices", {
  d <- sim_sf(n = 50, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
  m <- create_sfmodel_exp(y ~ x1 + x2, data = d)

  p1 <- add_priors(m, coef = list(mu = 0, v_i = 0.01))
  expect_length(p1$priors$b0, 3L)
  expect_equal(p1$priors$B0i, diag(0.01, 3))

  p2 <- add_priors(m, coef = list(mu = c(1, 0, 0), v_i = diag(1, 3)))
  expect_equal(p2$priors$b0, c(1, 0, 0))
  expect_equal(p2$priors$B0i, diag(1, 3))
})

test_that("partial prior lists keep the remaining defaults", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d),
                  sigma = list(shape = 2))

  expect_equal(m$priors$shape_v, 2)
  expect_equal(m$priors$rate_v, 0.01)
})

test_that("misspelled prior elements are reported rather than ignored", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- create_sfmodel_exp(y ~ x1, data = d)

  expect_error(add_priors(m, coef = list(mean = 0)), "Unknown element")
  expect_error(add_priors(m, lambda = list(rstar = 0.8)), "Unknown element")
  expect_error(add_priors(m, coef = 0.1), "must be a named list")
})

test_that("the improper limiting prior on the error precision is allowed", {
  # Koop's own programs use a flat prior on the error precision, which in the
  # shape and rate parameterisation means zeros. The posterior stays proper.
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_priors(create_sfmodel_exp(y ~ x1, data = d,
                                     iterations = 100, burnin = 50),
                  sigma = list(shape = 0, rate = 0))

  expect_equal(m$priors$shape_v, 0)
  est <- add_posterior_coefficients(m)
  expect_true(all(is.finite(est$posterior$sigma_v$coeffs)))
  expect_true(all(est$posterior$sigma_v$coeffs > 0))
})

test_that("prior arguments are validated", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m_exp <- create_sfmodel_exp(y ~ x1, data = d)
  m_hn <- create_sfmodel_hn(y ~ x1, data = d)

  expect_error(add_priors(m_exp, lambda = list(r_star = 0)), "between 0 and 1")
  expect_error(add_priors(m_exp, lambda = list(r_star = 1)), "between 0 and 1")
  expect_error(add_priors(m_hn, sigma_u = list(shape = 1)), "must exceed 1")
  expect_error(add_priors(m_exp, sigma = list(shape = -1)),
               "must be non-negative")
  expect_error(add_priors(m_exp, coef = list(mu = c(1, 2, 3))),
               "must be of length 1 or 2")
  expect_error(add_priors(m_exp, coef = list(v_i = matrix(1, 3, 3))),
               "2 x 2 matrix")
})

test_that("the shape of an inefficiency prior is validated", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m_exp <- create_sfmodel_exp(y ~ x1, data = d)
  m_hn <- create_sfmodel_hn(y ~ x1, data = d)

  # An exponential rate prior needs a positive shape to be proper. Without the
  # check a negative one reaches the sampler, and the prior mean it implies is
  # a negative starting value for a rate.
  for (bad in list(-1, 0, NA_real_, Inf, "a", c(1, 2))) {
    expect_error(add_priors(m_exp, lambda = list(shape = bad)),
                 "must be positive")
  }
  expect_silent(add_priors(m_exp, lambda = list(shape = 0.5)))

  for (bad in list(1, 0.5, NA_real_, Inf, "a", c(2, 3))) {
    expect_error(add_priors(m_hn, sigma_u = list(shape = bad)),
                 "must exceed 1")
  }

  # The four-component model elicits two of them, and names the one at fault.
  d4 <- sim_sf4(n = 20, n_time = 3, beta = c(1, 0.5))
  m4 <- create_sfmodel4_exp(y ~ x1, data = d4, id = "id")
  expect_error(add_priors(m4, lambda_eta = list(shape = -1)), "lambda_eta")
  expect_error(add_priors(m4, lambda_u = list(shape = -1)), "lambda_u")
})

test_that("a non-finite prior median efficiency is rejected by name", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- create_sfmodel_exp(y ~ x1, data = d)

  for (bad in list(NA_real_, NaN, Inf, "a", c(0.5, 0.6))) {
    expect_error(add_priors(m, lambda = list(r_star = bad)),
                 "between 0 and 1")
  }
  expect_error(add_priors(m, lambda = list(r_star = NA)),
               "lambda$r_star", fixed = TRUE)
})

test_that("the coefficient precision must be a precision matrix", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- create_sfmodel_exp(y ~ x1, data = d)

  expect_error(add_priors(m, coef = list(v_i = matrix(c(1, 2, 3, 4), 2))),
               "must be symmetric")
  expect_error(add_priors(m, coef = list(v_i = diag(-1, 2))),
               "positive semi-definite")
  expect_error(add_priors(m, coef = list(v_i = matrix(c(1, 2, 2, 1), 2))),
               "positive semi-definite")
  expect_error(add_priors(m, coef = list(v_i = diag(NA_real_, 2))),
               "finite numbers")
  expect_error(add_priors(m, coef = list(mu = c(1, NA))), "finite numbers")

  # Zero is the flat limiting prior and stays allowed, as it is for sigma.
  expect_equal(add_priors(m, coef = list(v_i = 0))$priors$B0i, diag(0, 2))
  expect_equal(add_priors(m, coef = list(v_i = diag(0, 2)))$priors$B0i,
               diag(0, 2))
})

test_that("a non-finite prior on the error precision is rejected", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- create_sfmodel_exp(y ~ x1, data = d)

  expect_error(add_priors(m, sigma = list(shape = NA)), "must be non-negative")
  expect_error(add_priors(m, sigma = list(rate = Inf)), "must be non-negative")

  m4 <- create_sfmodel4_exp(y ~ x1, data = sim_sf4(n = 20, n_time = 3,
                                                   beta = c(1, 0.5)),
                            id = "id")
  expect_error(add_priors(m4, sigma_mu = list(shape = NA)), "must be positive")
})
