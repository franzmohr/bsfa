estimated <- function(...) {
  set.seed(120)
  d <- sim_sf(n = 250, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)
  add_posterior_coefficients(add_seed(add_priors(
    create_sfmodel_exp(y ~ x1 + x2, data = d, iterations = 1000,
                       burnin = 500, ...)), 21))
}

test_that("coef and vcov summarise the frontier draws", {
  est <- estimated()
  b <- as.matrix(est$posterior$beta$coeffs)

  # coef.default() used to return NULL here, silently.
  expect_equal(coef(est), colMeans(b))
  expect_named(coef(est), c("(Intercept)", "x1", "x2"))
  expect_equal(vcov(est), cov(b))
  expect_equal(dim(vcov(est)), c(3L, 3L))
  expect_equal(nobs(est), 250L)
})

test_that("fitted is the frontier and residuals the composed error", {
  est <- estimated()

  expect_equal(unname(fitted(est)),
               as.numeric(est$data$X %*% coef(est)))
  expect_equal(unname(residuals(est)), est$data$y - unname(fitted(est)))
  expect_named(fitted(est), rownames(est$data$X))

  # The composed error still holds the one-sided term, so for a production
  # frontier it is centred below zero rather than on it.
  expect_lt(mean(residuals(est)), 0)
})

test_that("predict uses the levels the model was fitted with", {
  set.seed(121)
  d <- sim_sf(n = 200, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  d$g <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1 + g, data = d, iterations = 500, burnin = 200)))

  expect_equal(predict(est), fitted(est))

  full <- data.frame(x1 = c(0, 1),
                     g = factor(c("c", "a"), levels = c("a", "b", "c")))
  partial <- data.frame(x1 = c(0, 1), g = factor(c("c", "a")))
  # The second frame's factor has only two levels of its own, so without the
  # stored levels it would be coded against a different contrast basis.
  expect_equal(predict(est, full), predict(est, partial))
  expect_length(predict(est, full), 2L)

  expect_error(predict(est, data.frame(x1 = 1)), "not found")
})

test_that("logLik carries the degrees of freedom AIC and BIC need", {
  est <- estimated()
  ll <- logLik(est)

  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "df"), 5L)         # 3 coefficients, sigma_v, lambda
  expect_equal(attr(ll, "nobs"), 250L)
  expect_equal(as.numeric(ll),
               sum(bsfa:::sf_loglik_point(
                 est, coef(est), mean(est$posterior$sigma_v$coeffs),
                 mean(est$posterior$lambda$coeffs))))

  # AIC() and BIC() must agree with the package's own criteria, which are
  # built on the same deviance at the same point.
  crit <- selection_criteria(add_posterior_loglik(est))
  expect_equal(AIC(est), crit$AIC$mean)
  expect_equal(BIC(est), crit$BIC$mean)
})

test_that("the methods refuse a model that cannot support them", {
  d <- sim_sf(n = 60, n_time = 3, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  panel <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1, data = d, id = "id", iterations = 300,
                       burnin = 100)))
  expect_error(logLik(panel), "does not factorise")

  d4 <- sim_sf4(n = 20, n_time = 4, beta = c(1, 0.5))
  four <- add_posterior_coefficients(add_priors(
    create_sfmodel4_exp(y ~ x1, data = d4, id = "id", iterations = 300,
                        burnin = 100)))
  expect_error(logLik(four), "does not factorise")
  # The rest apply to any model that carries draws.
  expect_length(coef(four), 2L)
  expect_equal(nobs(four), 80L)

  bare <- create_sfmodel_exp(y ~ x1, data = d, iterations = 300)
  expect_error(coef(bare), "does not contain posterior draws")
  expect_error(vcov(bare), "does not contain posterior draws")
})

test_that("predict rejects a factor level the model never saw", {
  set.seed(123)
  d <- sim_sf(n = 120, beta = c(1, 0.5), sigma_v = 0.2, par_u = 4)
  d$g <- factor(sample(c("a", "b"), nrow(d), TRUE))
  est <- add_posterior_coefficients(add_priors(
    create_sfmodel_exp(y ~ x1 + g, data = d, iterations = 400, burnin = 100)))

  # Carrying the fitted levels is what makes this an error naming the variable
  # rather than a silent recoding against a different contrast basis.
  nd <- data.frame(x1 = 1, g = factor("z", levels = c("a", "b", "z")))
  expect_error(predict(est, nd), "factor g has new level")

  expect_equal(unname(predict(est, d[1:3, ])), unname(fitted(est))[1:3])
})
