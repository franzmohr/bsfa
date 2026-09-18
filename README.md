# bsfa

Bayesian estimation of stochastic frontier models in R.

The package implements Gibbs samplers with data augmentation for composed-error
frontier models in the tradition of van den Broeck, Koop, Osiewalski and Steel
(1994). Because the one-sided inefficiency term is sampled as a latent variable
rather than integrated out, efficiency scores come out of the sampler as draws.
`efficiency()` therefore reports genuine posterior distributions, where maximum
likelihood has to fall back on the Jondrow et al. (1982) conditional mean — a
point predictor of an unobserved quantity, with no comparable measure of
uncertainty attached.

At the time of writing there is no maintained general-purpose Bayesian
stochastic frontier package on CRAN. `sfaR`, `sfa`, `frontier`, `semsfa` and
`ssfa` are all maximum likelihood.

## Installation

```r
# install.packages("remotes")
remotes::install_github("franzmohr/bsfa")
```

## Usage

A model is built up step by step, in the manner of
[bvartools](https://github.com/franzmohr/bvartools): a constructor fixes the
data and the specification, and the `add_*` functions attach the remaining
blocks and then the posterior draws.

```r
library(bsfa)

set.seed(1234)
d <- sim_sf(n = 500, beta = c(1, 0.5, 0.3), sigma_v = 0.2, par_u = 4)

model <- create_sfmodel_exp(y ~ x1 + x2, data = d,
                            iterations = 5000, burnin = 2000)
model <- add_priors(model)
model <- add_initial_values(model)
model <- add_seed(model, 1234)
model <- add_posterior_coefficients(model)
model <- add_posterior_loglik(model)

summary(model)
plot(model, type = "efficiency")
selection_criteria(model)
head(efficiency(model))
```

The draws are attached to the same object rather than returned as a separate
result, and are stored as `coda` objects, one block per parameter. Printing the
model reports which blocks have been set, which is worth a glance before
committing to a long run.

The distribution of the inefficiency term is carried by the model's class rather
than by an argument, so each variant gets its own methods:

| Constructor | Class | Inefficiency |
|---|---|---|
| `create_sfmodel_exp()` | `sfmodel_exp` | exponential, rate `lambda` |
| `create_sfmodel_hn()` | `sfmodel_hn` | half-normal, scale `sigma_u` |

That split is not cosmetic. The prior on the inefficiency term is elicited from
a prior median efficiency `r_star`, and the two models map that anchor
differently:

```r
model <- add_priors(model, lambda = list(r_star = 0.85))   # sfmodel_exp
model <- add_priors(model, sigma_u = list(r_star = 0.85))  # sfmodel_hn
```

For the exponential model the mapping is exact. With a gamma prior of shape 1
and rate `c = -log(r_star)` on the rate parameter, the implied marginal prior
density of `u` is `c / (u + c)^2`, whose median is `c`, so the prior median of
`exp(-u)` is exactly `r_star`. For the half-normal model the same target is
matched only in expectation.

## What is implemented

| | |
|---|---|
| Inefficiency | exponential, half-normal |
| Frontier | production, cost |
| Data | cross-section, time-invariant panel (Pitt and Lee, 1981) |
| Output | posterior draws of coefficients, variances and unit-level efficiency |
| Inference | `summary()`, `plot()`, `efficiency()`, `selection_criteria()` (LL, AIC, BIC, HQ, WAIC) |

## Roadmap

The time-invariant panel model attributes *all* persistent heterogeneity between
units to inefficiency, which is a strong assumption whenever units differ for
reasons unrelated to efficiency. The planned extensions, roughly in order:

1. **Time-varying inefficiency** — `u_it` rather than `u_i`, so that the model
   can speak to cyclical rather than only structural variation in efficiency.
2. **Four-component model** (Kumbhakar, Lien and Hjalmarsson, 2014) — separating
   a unit effect, persistent inefficiency, noise and transient inefficiency. All
   full conditionals remain standard, so this is a hierarchical extension of the
   present sampler rather than a new algorithm.
3. **Inefficiency determinants** — a Battese and Coelli (1995) style regression
   in the mean or scale of the inefficiency distribution, estimated jointly with
   the frontier rather than in an inconsistent second step.
4. **Heteroskedasticity** in both error components, which matters whenever unit
   size is widely dispersed.
5. **Latent class / regime switching**, with class membership driven by
   covariates. Gibbs handles the discrete membership indicator naturally, where
   Hamiltonian Monte Carlo would not.

Each of these is a block of the specification rather than an argument to an
estimator, which is why the package is organised around a model object.

`selection_criteria()` does not provide LOOIC. The Pareto smoothed importance
sampling it needs is a piece of machinery in its own right and deserves its own
testing pass, so it is left out rather than added hastily. Note also that the
pointwise log-likelihood, and therefore every criterion, is unavailable for
panel models: the units share one inefficiency term, so integrating it out
couples the observations that belong to the same unit.

## Validation

The exponential, time-invariant panel specification reproduces the sampler
distributed by Justin Tobias for exercise 14.13 of Koop, Poirier and Tobias
(2007), *Bayesian Econometric Methods*. The `koop-exercise` vignette translates
both that program and its data-generating script into `bsfa` and compares the
results, which is worth having because the data-generating process is known.

The correspondence is exact. The mean of the truncated normal full conditional
for `z_i` carries the `-1/(T h mu_z)` shift contributed by the exponential
prior, which in the half-normal case is replaced by an additional `1/sigma_u^2`
term in the precision; that one line is the whole difference between the two
specifications inside the sampler.

The closed-form log-likelihoods are checked against direct numerical
integration of the composed-error density in the test suite, for both
distributions and for both orientations of the frontier.

## References

Aigner, D., Lovell, C. A. K., & Schmidt, P. (1977). Formulation and estimation
of stochastic frontier production function models. *Journal of Econometrics*,
6(1), 21–37.

Battese, G. E., & Coelli, T. J. (1995). A model for technical inefficiency
effects in a stochastic frontier production function for panel data.
*Empirical Economics*, 20(2), 325–332.

Jondrow, J., Lovell, C. A. K., Materov, I. S., & Schmidt, P. (1982). On the
estimation of technical inefficiency in the stochastic frontier production
function model. *Journal of Econometrics*, 19(2–3), 233–238.

Koop, G. (2003). *Bayesian Econometrics*. Chichester: Wiley.

Kumbhakar, S. C., Lien, G., & Hjalmarsson, L. (2014). Technical efficiency in
competing panel data models: A study of Norwegian grain farming. *Journal of
Productivity Analysis*, 41(2), 321–337.

Pitt, M. M., & Lee, L.-F. (1981). The measurement and sources of technical
inefficiency in the Indonesian weaving industry. *Journal of Development
Economics*, 9(1), 43–64.

Robert, C. P. (1995). Simulation of truncated normal variables. *Statistics and
Computing*, 5(2), 121–125.

van den Broeck, J., Koop, G., Osiewalski, J., & Steel, M. F. J. (1994).
Stochastic frontier models: A Bayesian perspective. *Journal of Econometrics*,
61(2), 273–303.

## License

GPL (>= 2)
