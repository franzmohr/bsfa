test_that("constructors set the class that carries the model variant", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)

  m_exp <- create_sfmodel_exp(y ~ x1, data = d)
  m_hn <- create_sfmodel_hn(y ~ x1, data = d)

  expect_s3_class(m_exp, "sfmodel_exp")
  expect_s3_class(m_exp, "sfmodel")
  expect_s3_class(m_hn, "sfmodel_hn")
  expect_s3_class(m_hn, "sfmodel")

  expect_equal(m_exp$model$ineff, "exponential")
  expect_equal(m_hn$model$ineff, "halfnormal")
  expect_equal(m_exp$model$par_u_name, "lambda")
  expect_equal(m_hn$model$par_u_name, "sigma_u")
})

test_that("the model object carries the data and the MCMC settings", {
  d <- sim_sf(n = 50, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
  m <- create_sfmodel_exp(y ~ x1 + x2, data = d,
                          iterations = 400, burnin = 100, thin = 2)

  expect_equal(m$n, 50L)
  expect_equal(m$k, 3L)
  expect_equal(dim(m$data$X), c(50L, 3L))
  expect_equal(colnames(m$data$X), c("(Intercept)", "x1", "x2"))
  expect_equal(m$iterations, 400)
  expect_equal(m$burnin, 100)
  expect_equal(m$thin, 2)
  expect_null(m$priors)
  expect_null(m$initial)
  expect_null(m$seed)
})

test_that("without an id every observation owns an inefficiency term", {
  d <- sim_sf(n = 40, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4, n_time = 3)

  m <- create_sfmodel_exp(y ~ x1, data = d)
  expect_false(m$model$panel)
  expect_equal(m$data$n_units, 120L)

  mp <- create_sfmodel_exp(y ~ x1, data = d, id = "id")
  expect_true(mp$model$panel)
  expect_equal(mp$data$n_units, 40L)
  expect_equal(length(mp$data$g), 120L)
})

test_that("an id vector works as well as an id name", {
  d <- sim_sf(n = 20, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4, n_time = 2)

  by_name <- create_sfmodel_exp(y ~ x1, data = d, id = "id")
  by_value <- create_sfmodel_exp(y ~ x1, data = d, id = d$id)

  expect_equal(by_name$data$g, by_value$data$g)
  expect_equal(by_name$data$unit_labels, by_value$data$unit_labels)
})

test_that("specification errors are caught early", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)

  expect_error(create_sfmodel_exp(y ~ x1, data = d, iterations = 100,
                                  thin = 3), "multiple of")
  expect_error(create_sfmodel_exp(y ~ x1, data = d, iterations = 0),
               "must be positive")
  expect_error(create_sfmodel_exp(y ~ x1, data = d, id = "nope"),
               "not found")
  expect_error(create_sfmodel_exp(y ~ x1, data = d[1:2, ]),
               "as many coefficients as observations")
})

test_that("printing reports which blocks have been set", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- create_sfmodel_exp(y ~ x1, data = d)

  expect_output(print(m), "priors\\s+not set")
  expect_output(print(add_priors(m)), "priors\\s+set")
  expect_output(print(add_seed(m, 1)), "seed\\s+set")
})
