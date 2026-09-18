sim <- function(n = 300, seed = 140) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), z = runif(n, 1, 3))
  d$u <- rexp(n, 5)
  d$y <- 1 + 0.5 * d$x1 + 2 * log(d$z) - d$u + rnorm(n, sd = 0.15)
  d
}

test_that("an offset enters the model rather than being discarded", {
  d <- sim()
  fit <- function(f) add_posterior_coefficients(add_seed(add_priors(
    create_sfmodel_exp(f, data = d, iterations = 1000, burnin = 500)), 41),
    keep_u = FALSE)

  # Fixing the coefficient at one must come to the same thing as taking the
  # term out of the response beforehand. It used to come to the same thing as
  # leaving the variable out altogether.
  expect_equal(coef(fit(y ~ x1 + offset(2 * log(z)))),
               coef(fit(I(y - 2 * log(z)) ~ x1)))
  expect_false(isTRUE(all.equal(coef(fit(y ~ x1 + offset(2 * log(z)))),
                                coef(fit(y ~ x1)))))
})

test_that("the offset is stored and the response adjusted", {
  d <- sim()
  m <- create_sfmodel_exp(y ~ x1 + offset(2 * log(z)), data = d,
                          iterations = 100)

  expect_equal(m$data$offset, 2 * log(d$z))
  expect_equal(m$data$y, d$y - 2 * log(d$z))
  expect_equal(colnames(m$data$X), c("(Intercept)", "x1"))
  expect_output(print(m), "Offset")

  # Without one, nothing is stored and the response is untouched.
  plain <- create_sfmodel_exp(y ~ x1, data = d, iterations = 100)
  expect_null(plain$data$offset)
  expect_equal(plain$data$y, d$y)
  expect_false(any(grepl("Offset", utils::capture.output(print(plain)))))
})

test_that("fitted, residuals and predict agree with the offset in place", {
  d <- sim()
  est <- add_posterior_coefficients(add_seed(add_priors(
    create_sfmodel_exp(y ~ x1 + offset(2 * log(z)), data = d,
                       iterations = 1000, burnin = 500)), 42), keep_u = FALSE)

  expect_equal(unname(fitted(est)),
               unname(bsfa:::sf_frontier(est)) + 2 * log(d$z))
  # The residual is the response as supplied minus the fitted frontier.
  expect_equal(unname(residuals(est)), d$y - unname(fitted(est)))
  expect_equal(unname(predict(est, d)), unname(fitted(est)))
  expect_equal(unname(predict(est, d[1:5, ])), unname(fitted(est))[1:5])

  # A different offset in newdata must move the prediction by that difference.
  nd <- d[1:5, ]
  nd$z <- nd$z * exp(1)
  expect_equal(unname(predict(est, nd)),
               unname(predict(est, d[1:5, ])) + 2)
})

test_that("an offset works for the panel and four-component models", {
  set.seed(143)
  nu <- 30; tt <- 5; n <- nu * tt
  d <- data.frame(id = rep(seq_len(nu), each = tt), x1 = rnorm(n),
                  z = runif(n, 1, 3))
  d$y <- 1 + 0.5 * d$x1 + log(d$z) - rexp(nu, 5)[d$id] + rnorm(n, sd = 0.15)

  for (mk in list(create_sfmodel_exp, create_sfmodel4_exp)) {
    m <- mk(y ~ x1 + offset(log(z)), data = d, id = "id", iterations = 100)
    expect_equal(m$data$offset, log(d$z))
    expect_equal(m$data$y, d$y - log(d$z))
  }
})

test_that("a row whose offset is missing is dropped with the rest", {
  d <- sim(n = 60, seed = 144)
  d$z[1:5] <- NA

  expect_message(m <- create_sfmodel_exp(y ~ x1 + offset(2 * log(z)),
                                         data = d, iterations = 100),
                 "Dropping 5 of 60")
  expect_equal(m$n, 55L)
  expect_length(m$data$offset, 55L)
  expect_false(anyNA(m$data$offset))
})

test_that("a rank deficient design is refused and the column named", {
  d <- sim(n = 100, seed = 145)
  d$dup <- d$x1

  expect_error(create_sfmodel_exp(y ~ x1 + dup, data = d, iterations = 100),
               "rank 2")
  expect_error(create_sfmodel_exp(y ~ x1 + dup, data = d, iterations = 100),
               "Column dup")
  expect_error(create_sfmodel_exp(y ~ x1 + I(2 * x1), data = d,
                                  iterations = 100),
               "linear combination")

  # Designs that are merely correlated, or that carry a factor, are fine.
  d$near <- d$x1 + rnorm(100, sd = 1e-3)
  d$g <- factor(rep(c("a", "b", "c"), length.out = 100))
  expect_silent(create_sfmodel_exp(y ~ x1 + near, data = d, iterations = 100))
  expect_silent(create_sfmodel_exp(y ~ x1 + g, data = d, iterations = 100))
})
