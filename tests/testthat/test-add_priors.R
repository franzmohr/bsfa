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

test_that("the half-normal prior centres sigma_u on the elicited value", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- add_priors(create_sfmodel_hn(y ~ x1, data = d),
                  sigma_u = list(r_star = 0.75))

  sigma_u <- -log(0.75) / qnorm(0.75)
  expect_equal(m$priors$rate_u / (m$priors$shape_u - 1), sigma_u^2)
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

test_that("prior arguments are validated", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m_exp <- create_sfmodel_exp(y ~ x1, data = d)
  m_hn <- create_sfmodel_hn(y ~ x1, data = d)

  expect_error(add_priors(m_exp, lambda = list(r_star = 0)), "between 0 and 1")
  expect_error(add_priors(m_exp, lambda = list(r_star = 1)), "between 0 and 1")
  expect_error(add_priors(m_hn, sigma_u = list(shape = 1)), "must exceed 1")
  expect_error(add_priors(m_exp, sigma = list(shape = -1)), "must be positive")
  expect_error(add_priors(m_exp, coef = list(mu = c(1, 2, 3))),
               "must be of length 1 or 2")
  expect_error(add_priors(m_exp, coef = list(v_i = matrix(1, 3, 3))),
               "2 x 2 matrix")
})
