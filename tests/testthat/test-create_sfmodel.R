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

test_that("the identifier survives dropped rows and non-integer rownames", {
  # Aligning the identifier from rownames() used to coerce them to integers,
  # which silently produced NA units and collapsed the panel into one whenever
  # the rownames were not the default integers and a row had been dropped.
  d <- sim_sf(n = 20, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4, n_time = 3)
  d$x1[3] <- NA
  rownames(d) <- paste0("b", seq_len(nrow(d)))

  by_name <- create_sfmodel_exp(y ~ x1, data = d, id = "id")
  by_value <- create_sfmodel_exp(y ~ x1, data = d, id = d$id)

  for (m in list(by_name, by_value)) {
    expect_equal(m$n, 59L)
    expect_equal(m$data$n_units, 20L)
    expect_false(anyNA(m$data$g))
    expect_length(m$data$unit_labels, 20L)
  }
  expect_equal(by_name$data$g, by_value$data$g)
})

test_that("character identifiers are kept as labels", {
  d <- sim_sf(n = 10, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4, n_time = 3)
  d$id <- paste0("bank_", d$id)

  m <- create_sfmodel_exp(y ~ x1, data = d, id = "id")
  expect_equal(m$data$n_units, 10L)
  expect_true(all(grepl("^bank_", m$data$unit_labels)))
})

test_that("an observation whose unit is unknown is dropped", {
  d <- sim_sf(n = 20, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4, n_time = 3)
  d$id[1:3] <- NA

  m <- create_sfmodel_exp(y ~ x1, data = d, id = "id")
  expect_equal(m$n, 57L)
  expect_equal(m$data$n_units, 19L)
  expect_false(anyNA(m$data$g))
})

test_that("the identifier does not enter the design matrix", {
  d <- sim_sf(n = 10, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4, n_time = 2)
  m <- create_sfmodel_exp(y ~ x1, data = d, id = "id")

  expect_equal(colnames(m$data$X), c("(Intercept)", "x1"))
  expect_equal(ncol(m$data$X), 2L)
})

test_that("a mismatched identifier is rejected rather than recycled", {
  d <- sim_sf(n = 10, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4, n_time = 2)

  expect_error(create_sfmodel_exp(y ~ x1, data = d, id = 1:5),
               "one element per row")
  expect_error(create_sfmodel_exp(y ~ x1, data = d[0, ], id = integer(0)),
               "No complete observations remain")
})

test_that("printing reports which blocks have been set", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  m <- create_sfmodel_exp(y ~ x1, data = d)

  expect_output(print(m), "priors\\s+not set")
  expect_output(print(add_priors(m)), "priors\\s+set")
  expect_output(print(add_seed(m, 1)), "seed\\s+set")
})

test_that("the MCMC settings must be whole numbers an integer can hold", {
  d <- sim_sf(n = 50, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  build <- function(...) create_sfmodel_exp(y ~ x1, data = d, ...)

  # Left unchecked these reach as.integer(), which truncates a fractional
  # value and turns one beyond the integer range into NA, neither with an
  # error of its own.
  expect_error(build(iterations = 100.5), "whole number")
  expect_error(build(iterations = 3e9), "whole number")
  expect_error(build(iterations = "100"), "whole number")
  expect_error(build(iterations = NA), "whole number")
  expect_error(build(iterations = c(100, 200)), "whole number")
  expect_error(build(burnin = 1.5), "'burnin'")
  expect_error(build(thin = NaN), "'thin'")

  # The bounds are still reported separately from the type.
  expect_error(build(iterations = 0), "must be positive")
  expect_error(build(burnin = -1), "non-negative")
})

test_that("dropped rows are reported and recorded", {
  d <- sim_sf(n = 40, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  d$x1[1:12] <- NA

  expect_message(m <- create_sfmodel_exp(y ~ x1, data = d),
                 "Dropping 12 of 40 observations")
  expect_equal(m$n, 28L)
  expect_s3_class(m$na.action, "omit")
  expect_equal(unname(m$na.action), 1:12)
  expect_output(print(m), "28 (12 dropped)", fixed = TRUE)

  # Nothing is said, or recorded, when nothing is dropped.
  expect_silent(clean <- create_sfmodel_exp(y ~ x1, data = d[13:40, ]))
  expect_null(clean$na.action)
})

test_that("a unit left with one observation is named", {
  d <- sim_sf4(n = 5, n_time = 4, beta = c(1, 0.5))
  d$bank <- rep(c("DE001", "FR002", "IT003", "ES004", "NL005"), each = 4)
  d$x1[d$bank == "IT003"][1:3] <- NA

  err <- tryCatch(
    suppressMessages(create_sfmodel4_exp(y ~ x1, data = d, id = "bank")),
    error = function(e) conditionMessage(e))
  expect_match(err, "IT003")
  expect_match(err, "3 rows were dropped for missing values first")
})
