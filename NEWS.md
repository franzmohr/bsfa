# bsfa 0.0.0.9000

First development version.

* `create_sfmodel_exp()` and `create_sfmodel_hn()` set up a stochastic frontier
  model with exponential or half-normal inefficiency. The distribution is
  carried by the object's class, so each variant has its own methods.
* `add_priors()` attaches the prior specification, eliciting the prior on the
  inefficiency term from a prior median efficiency. The mapping is exact for
  the exponential model and matched in expectation for the half-normal one.
* `add_initial_values()` and `add_seed()` attach the remaining blocks.
* `draw_posterior()` runs the Gibbs sampler and returns an object of class
  `bsfa_exp` or `bsfa_hn`.
* `efficiency()` returns posterior summaries of the unit-level efficiency
  scores.
* `sim_sf()` generates artificial data from the model.
