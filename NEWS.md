# bsfa 0.1.0

First release.

## Models

* `create_sfmodel_exp()` and `create_sfmodel_hn()` set up a stochastic frontier
  model with exponential or half-normal inefficiency, as a cross-section or,
  with `id`, as the time-invariant panel model of Pitt and Lee (1981). The
  distribution is carried by the object's class, so each variant has its own
  methods.
* `create_sfmodel4_exp()` and `create_sfmodel4_hn()` set up the four-component
  model of Kumbhakar, Lien and Hjalmarsson (2014), which separates a unit
  effect, persistent inefficiency, noise and transient inefficiency.
* Production and cost frontiers throughout.

## Building a model

* `add_priors()` attaches the prior specification, eliciting the prior on the
  inefficiency term from a prior median efficiency. The mapping is exact for
  the exponential model and matched in expectation for the half-normal one.
* `add_initial_values()` and `add_seed()` attach the remaining blocks. The
  former also draws a seed for the simulation if the model does not have one,
  so that a model object determines its own output.
* `add_posterior_coefficients()` runs the Gibbs sampler and adds the draws to
  the model as `coda` objects, one block per parameter. `keep_u` takes an
  interval at which to store the augmented draws, which on a wide panel are
  otherwise by far the largest part of the object.
* `add_posterior_loglik()` adds the pointwise log-likelihood, with the
  inefficiency term integrated out. It is available for cross-sections only,
  since the marginal likelihood of a panel unit does not factorise over its
  observations.

## Reading the results

* `summary()`, `plot()` and `efficiency()` report the results. `plot()` takes a
  `type` of `"hist"`, `"trace"`, `"boxplot"` or `"efficiency"`.
* The summary table ends in `ESS`, the effective sample size of each block.
  The draws are autocorrelated, and in the four-component model severely so,
  where the unit effect and persistent inefficiency trade off against each
  other from sweep to sweep; a credible band is worth reading only against the
  number of independent draws behind it.
* `efficiency()` summarises the posterior distribution of the scores rather
  than a point predictor of them, which is what sampling the inefficiency term
  as a latent variable buys. For the four-component model it takes a `type` of
  `"overall"`, `"persistent"` or `"transient"`; the two observation-level types
  carry an `obs` column so that the scores can be joined back to the data.
* `selection_criteria()` returns LL, AIC, BIC, HQ and WAIC from the
  log-likelihood draws. Like the log-likelihood itself, it is unavailable for
  panel and four-component models.

## Simulation

* `sim_sf()` and `sim_sf4()` generate artificial data from the two- and
  four-component models, for examples, tests and prior predictive checks.

## Notes

* Rows with a missing value, or with an unknown unit, are dropped with a
  message and recorded in the `na.action` element of the model object.
* The response is assumed to be on a logarithmic scale. The prior on the
  inefficiency term is expressed in the units of the response, so fitting
  output in levels does not merely change the units of the efficiency scores,
  it changes their values.
