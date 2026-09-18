estimated <- function() {
  set.seed(77)
  d <- sim_sf(n = 60, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, iterations = 200, burnin = 100)))
}

test_that("any number of quantiles can be requested", {
  est <- estimated()

  for (probs in list(0.5, c(0.1, 0.9), c(0.05, 0.25, 0.5, 0.75, 0.95))) {
    eff <- efficiency(est, probs = probs)
    expect_equal(nrow(eff), 60L)
    expect_equal(names(eff),
                 c("unit", "mean", "sd",
                   paste0(format(100 * probs, trim = TRUE), "%")))
  }
})

test_that("a single quantile is the one the full set reports", {
  est <- estimated()

  expect_equal(efficiency(est, probs = 0.5)[["50%"]],
               efficiency(est, probs = c(0.05, 0.5, 0.95))[["50%"]])
})

test_that("the quantiles are those of the efficiency draws", {
  est <- estimated()
  r <- exp(-as.matrix(est$posterior$u$coeffs))
  eff <- efficiency(est, probs = c(0.1, 0.9))

  expect_equal(eff$mean, unname(colMeans(r)))
  expect_equal(eff[["10%"]],
               unname(apply(r, 2, stats::quantile, probs = 0.1)))
  expect_equal(eff[["90%"]],
               unname(apply(r, 2, stats::quantile, probs = 0.9)))
})

test_that("the four-component method takes a single quantile too", {
  set.seed(78)
  d <- sim_sf4(n = 20, n_time = 4, beta = c(1, 0.5))
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel4_exp(y ~ x1, data = d, id = "id",
                        iterations = 200, burnin = 100)))

  for (type in c("persistent", "transient", "overall")) {
    eff <- efficiency(est, type = type, probs = 0.5)
    expect_equal(names(eff),
                 if (type == "persistent") c("unit", "mean", "sd", "50%")
                 else c("unit", "obs", "mean", "sd", "50%"))
    expect_equal(nrow(eff), if (type == "persistent") 20L else 80L)
  }
})

test_that("observation-level scores carry the row they came from", {
  set.seed(79)
  d <- sim_sf4(n = 3, n_time = 4, beta = c(1, 0.5))
  d$bank <- rep(c("DE001", "FR002", "IT003"), each = 4)
  d$quarter <- rep(2020:2023, times = 3)
  rownames(d) <- paste0(d$bank, "-", d$quarter)

  est <- add_posterior_coefficients(add_priors(
    create_sfmodel4_exp(y ~ x1, data = d, id = "bank",
                        iterations = 200, burnin = 100)))

  # Labelled only by the unit, a transient score cannot be told from the
  # eleven others that unit has.
  for (type in c("transient", "overall")) {
    eff <- efficiency(est, type = type)
    expect_equal(names(eff)[1:2], c("unit", "obs"))
    expect_equal(eff$obs, rownames(d))
    expect_equal(eff$unit, as.character(d$bank))
    expect_equal(nrow(eff), 12L)
  }

  # One score per unit needs no such column.
  persistent <- efficiency(est, type = "persistent")
  expect_false("obs" %in% names(persistent))
  expect_equal(persistent$unit, c("DE001", "FR002", "IT003"))

  # The scores themselves are unchanged by the extra column.
  eff <- efficiency(est, type = "transient")
  expect_equal(eff$mean,
               unname(colMeans(exp(-as.matrix(est$posterior$u$coeffs)))))
})
