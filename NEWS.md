# bsfa 0.0.0.9000

First development version.

* `create_sfmodel_exp()` and `create_sfmodel_hn()` set up a stochastic frontier
  model with exponential or half-normal inefficiency. The distribution is
  carried by the object's class, so each variant has its own methods.
* `add_priors()` attaches the prior specification, eliciting the prior on the
  inefficiency term from a prior median efficiency. The mapping is exact for
  the exponential model and matched in expectation for the half-normal one.
* `add_initial_values()` and `add_seed()` attach the remaining blocks. The
  former also draws a seed for the simulation if the model does not have one.
* `add_posterior_coefficients()` runs the Gibbs sampler and adds the draws to
  the model as `coda` objects, one block per parameter.
* `add_posterior_loglik()` adds the pointwise log-likelihood, with the
  inefficiency term integrated out.
* `selection_criteria()` returns LL, AIC, BIC, HQ and WAIC from those draws.
* `summary()`, `plot()` and `efficiency()` report the results. `plot()` takes
  a `type` of `"hist"`, `"trace"`, `"boxplot"` or `"efficiency"`.
* `sim_sf()` generates artificial data from the model.
