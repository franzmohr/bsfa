mk <- function(lambda, sign = -1, n = 400, seed = 160) {
  set.seed(seed)
  x1 <- rnorm(n)
  u <- if (is.na(lambda)) rep(0, n) else rexp(n, lambda)
  data.frame(y = 1 + 0.5 * x1 + sign * u + rnorm(n, sd = 0.3), x1 = x1)
}

test_that("the statistic is the one Coelli defines", {
  d <- mk(4)
  m <- create_sfmodel_exp(y ~ x1, data = d)
  t <- skewness_test(m)

  r <- as.numeric(stats::.lm.fit(m$data$X, m$data$y)$residuals)
  m2 <- mean(r^2); m3 <- mean(r^3); n <- length(r)
  expect_equal(t$skewness, m3 / m2^1.5)
  expect_equal(t$statistic, m3 / sqrt(6 * m2^3 / n))
  # One sided, towards the tail a production frontier implies.
  expect_equal(t$p.value, stats::pnorm(t$statistic))
  expect_equal(t$n, 400L)
  expect_s3_class(t, "sfskew")
})

test_that("it tells the three cases apart", {
  # Plenty of inefficiency: skewed as the frontier implies.
  clear <- skewness_test(create_sfmodel_exp(y ~ x1, data = mk(2)))
  expect_lt(clear$skewness, 0)
  expect_lt(clear$p.value, 0.01)
  expect_equal(clear$verdict, "informative")

  # Skewed the other way: the wrong skew case.
  wrong <- skewness_test(create_sfmodel_exp(y ~ x1, data = mk(2, sign = 1)))
  expect_gt(wrong$skewness, 0)
  expect_equal(wrong$verdict, "wrong")

  # Symmetric: no evidence either way.
  none <- skewness_test(create_sfmodel_exp(y ~ x1,
                                           data = mk(NA, seed = 161)))
  expect_true(none$verdict %in% c("weak", "wrong"))
})

test_that("the frontier fixes which tail counts as evidence", {
  d <- mk(2, sign = 1)          # inefficiency added, so a cost frontier
  prod <- skewness_test(create_sfmodel_exp(y ~ x1, data = d))
  cost <- skewness_test(create_sfmodel_exp(y ~ x1, data = d, type = "cost"))

  expect_equal(prod$skewness, cost$skewness)
  expect_equal(prod$statistic, cost$statistic)
  expect_equal(prod$expected, "negative")
  expect_equal(cost$expected, "positive")
  # The same residuals are evidence for one orientation and against the other.
  expect_equal(cost$p.value, 1 - prod$p.value)
  expect_equal(cost$verdict, "informative")
  expect_equal(prod$verdict, "wrong")
})

test_that("it needs no posterior and works for every model class", {
  d <- mk(3)
  expect_s3_class(skewness_test(create_sfmodel_exp(y ~ x1, data = d)),
                  "sfskew")
  expect_s3_class(skewness_test(create_sfmodel_hn(y ~ x1, data = d)),
                  "sfskew")

  p <- sim_sf4(n = 40, n_time = 5, beta = c(1, 0.5))
  for (m in list(create_sfmodel_exp(y ~ x1, data = p, id = "id"),
                 create_sfmodel4_exp(y ~ x1, data = p, id = "id"),
                 create_sfmodel4_hn(y ~ x1, data = p, id = "id"))) {
    expect_s3_class(skewness_test(m), "sfskew")
  }

  expect_error(skewness_test(1:10), "stochastic frontier model")
})

test_that("it reads the response the sampler sees, offset and all", {
  d <- mk(3)
  set.seed(162)
  d$z <- runif(nrow(d), 1, 3)
  d$y <- d$y + 2 * log(d$z)

  with_offset <- skewness_test(
    create_sfmodel_exp(y ~ x1 + offset(2 * log(z)), data = d))
  by_hand <- skewness_test(
    create_sfmodel_exp(I(y - 2 * log(z)) ~ x1, data = d))
  ignored <- skewness_test(create_sfmodel_exp(y ~ x1, data = d))

  expect_equal(with_offset$skewness, by_hand$skewness)
  # Leaving the offset out buries the skew in variation it does not explain.
  expect_lt(abs(ignored$skewness), abs(with_offset$skewness))
})

test_that("printing says which case it is", {
  expect_output(print(skewness_test(create_sfmodel_exp(y ~ x1, data = mk(2)))),
                "skewed as the frontier implies")
  expect_output(print(skewness_test(create_sfmodel_exp(y ~ x1,
                                                       data = mk(2, sign = 1)))),
                "skewed the wrong way")
  expect_output(print(skewness_test(create_sfmodel_exp(y ~ x1, data = mk(2)))),
                "M3T statistic")
})
