estimated <- function() {
  set.seed(130)
  d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
  add_posterior_coefficients(add_seed(add_priors(
    create_sfmodel_exp(y ~ x1 + x2, data = d, iterations = 1000,
                       burnin = 500)), 31))
}

test_that("the coefficient table ends in the effective sample size", {
  est <- estimated()
  tab <- summary(est)$coefficients

  expect_equal(colnames(tab),
               c("mean", "sd", "2.5%", "median", "97.5%", "ESS"))
  expect_equal(rownames(tab),
               c("(Intercept)", "x1", "x2", "sigma_v", "lambda"))
  expect_true(all(tab[, "ESS"] > 0))
  expect_true(all(is.finite(tab[, "ESS"])))
})

test_that("the effective sample size is the one coda computes", {
  est <- estimated()
  pars <- bsfa:::sf_par_draws(est)

  expect_equal(unname(summary(est)$coefficients[, "ESS"]),
               unname(as.numeric(coda::effectiveSize(coda::mcmc(pars)))))

  # It is a count of draws, not a proportion, so it scales with the chain.
  short <- summary(est)$coefficients[, "ESS"]
  set.seed(130)
  d <- sim_sf(n = 200, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
  long <- summary(add_posterior_coefficients(add_seed(add_priors(
    create_sfmodel_exp(y ~ x1 + x2, data = d, iterations = 4000,
                       burnin = 500)), 31)))$coefficients[, "ESS"]
  expect_true(mean(long / short) > 2)
})

test_that("it is printed as a count rather than to four decimals", {
  est <- estimated()
  out <- utils::capture.output(print(summary(est)))
  line <- grep("^sigma_v", out, value = TRUE)

  expect_length(line, 1L)
  # The last field on the row is the ESS and must carry no decimal point.
  last <- utils::tail(strsplit(trimws(line), " +")[[1]], 1)
  expect_false(grepl(".", last, fixed = TRUE))
  expect_gt(as.numeric(last), 0)
  expect_true(any(grepl("ESS", out, fixed = TRUE)))
})

test_that("the four-component model reports it for every block", {
  set.seed(131)
  d <- sim_sf4(n = 30, n_time = 5, beta = c(1, 0.5))
  est <- add_posterior_coefficients(add_seed(add_priors(
    create_sfmodel4_exp(y ~ x1, data = d, id = "id", iterations = 1000,
                        burnin = 500)), 32), keep_u = FALSE)
  tab <- summary(est)$coefficients

  expect_equal(rownames(tab),
               c("(Intercept)", "x1", "sigma_v", "sigma_mu", "lambda_eta",
                 "lambda_u"))
  expect_true("ESS" %in% colnames(tab))
  expect_true(all(tab[, "ESS"] > 0))
})

test_that("a posterior too short to measure reports NA rather than failing", {
  set.seed(132)
  d <- sim_sf(n = 60, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  const <- function(nm, rows) {
    list(coeffs = matrix(1, rows, length(nm), dimnames = list(NULL, nm)))
  }
  elsewhere <- function(rows) function(object) {
    object$posterior <- list(beta = const(c("(Intercept)", "x1"), rows),
                             sigma_v = const("sigma_v", rows),
                             lambda = const("lambda", rows))
    object
  }
  for (rows in c(1L, 2L)) {
    est <- add_posterior_coefficients(
      create_sfmodel_exp(y ~ x1, data = d, iterations = 50),
      posterior_function = elsewhere(rows))
    tab <- summary(est)$coefficients
    expect_true(all(is.na(tab[, "ESS"])))
    expect_true(all(is.finite(tab[, "mean"])))
  }
})
