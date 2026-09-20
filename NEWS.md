# bsfa (development version)

## Determinants of inefficiency

* Every model constructor gains `scale_u`, a one-sided formula naming
  covariates that the scale of the inefficiency term depends on. The scale is
  multiplied by `exp(z'gamma)`, so a coefficient of zero is the model without
  determinants and the parameter the scale is measured against keeps the
  prior `add_priors()` gives it. The set carries no intercept, since that
  parameter already plays the part of one. This is the heteroskedastic
  specification of Caudill, Ford and Gropper (1995), and it is what the
  applied four-component literature uses to let covariates explain
  inefficiency. A model created without determinants draws exactly what it
  drew before they existed, down to the last bit, so a script pinned to a
  seed under 0.2.0 reproduces its results.
* `create_sfmodel4_exp()` and `create_sfmodel4_hn()` take `scale_eta` in
  addition, for the persistent term. The two sets sit at different levels:
  the persistent term is one per unit, so its determinants have to be
  constant within a unit and are refused with a message naming the offending
  variable if they are not, while the transient term's may vary within one.
* A new inefficiency family, the truncated normal, in
  `create_sfmodel_tn()` and `create_sfmodel4_tn()`. Its one-sided term is
  `N+(z'delta, sigma^2)`, whose pre-truncation mean the argument `mean_u`,
  and `mean_eta` in the four-component model, gives determinants. With the
  default `~ 1` this is the constant mean of Stevenson (1980); with
  covariates it is Battese and Coelli (1995); together with `scale_u` it is
  Wang (2002). The half-normal and the exponential are anchored at zero and
  refuse `mean_u` rather than ignoring it.
* `add_priors()` gains the matching arguments. A determinant coefficient
  takes a normal prior, defaulting to mean zero and precision 0.01, and the
  prior has to be proper: the coefficient acts through an exponential, so a
  flat prior leaves the scale of the inefficiency term unbounded.
* The draws are added to the posterior as the blocks `scale_u`, `mean_u`,
  `scale_eta` and `mean_eta`, and `summary()` reports them beside the
  frontier coefficients and names which term each set explains.
* Neither the scale nor the mean coefficients have a conjugate full
  conditional, so each is drawn by an adaptive random walk Metropolis step
  inside the Gibbs sweep, as is the baseline scale of a truncated normal,
  whose normalising constant carries it. The proposal is tuned during
  burn-in towards an acceptance rate of 0.234 and then fixed, so the
  retained draws come from a kernel with the right invariant distribution.
  The acceptance rates are returned on the model as `acceptance` and printed
  by `summary()`, with a warning where a block barely moved.
* The pointwise log-likelihood, and with it `selection_criteria()`, is
  refused for a model carrying determinants and for the truncated normal
  family: the closed form the one-sided term is integrated out with assumes
  a single scale and a half-normal shape. The refusal explains itself rather
  than returning a number that is not the log-likelihood of the model.

## Variable selection

* `create_sfmodel_exp()`, `create_sfmodel_hn()`, `create_sfmodel4_exp()` and
  `create_sfmodel4_hn()` gain the argument `varsel`, which places the frontier
  coefficients under one of the two variable selection algorithms of the
  bvartools package. Either way the sampler draws an inclusion indicator per
  selected coefficient in every sweep.
* `varsel = "ssvs"` is the stochastic search variable selection of George, Sun
  and Ni (2008). Each selected coefficient carries a mixture of two normal
  priors centred on zero, a tight one standing for the regressor being absent
  from the frontier and a loose one for its being present, and the indicator
  says which is in force. The regressor never leaves the design.
* `varsel = "bvs"` is the Bayesian variable selection of Korobilis (2013). The
  frontier is `X %*% diag(lambda) %*% beta`, so an excluded regressor leaves
  the likelihood outright and the coefficient stored for it is an exact zero.
  Its prior is the ordinary normal in `coef`, which has to be proper, since the
  sweeps that exclude the regressor draw the coefficient from it. A prior
  precision of zero is refused, and one wide enough that the likelihood would
  never admit the coefficient back is reported: measured against the least
  squares standard error of the same coefficient, a prior standard deviation of
  two hundred times it left the indicator of an irrelevant regressor stuck for
  a whole run of 2000 sweeps.
* `add_priors()` gains the matching argument `varsel`. Both algorithms take
  `inprior`, the prior inclusion probability, and `include` and
  `exclude_intercept` to choose the candidates; SSVS takes in addition either
  `tau`, the two prior standard deviations, or `semiautomatic`, the two factors
  to scale the least squares standard error of each coefficient by. It is
  required for a model created with a `varsel` algorithm and not allowed for
  any other.
* The draws of the indicators are added to the posterior as the block
  `inclusion`, one column per selected coefficient, and `summary()` reports
  their means as `PIP`, the posterior probability that the regressor belongs in
  the frontier, beside the coefficient it belongs to.
* The selection applies to the frontier only. The inefficiency term, the error
  and the efficiency scores are drawn exactly as they are without it, except
  that they are now averaged over the frontiers the selection admits. AIC, BIC
  and HQ still charge the model for every coefficient the frontier was written
  with, because the number the selection keeps is not the same in every draw and
  there is no whole number to subtract; WAIC reads the effective number of
  parameters off the draws and accounts for it.
* A new `variable-selection` vignette covers both algorithms: what an inclusion
  probability means and how it should not be read, what each of them asks of the
  coefficient prior, why a prior too wide for BVS leaves the indicators stuck,
  and which of the two to reach for.

## Documentation

* `add_seed()` and `add_posterior_coefficients()` now state what the seed does
  and does not guarantee. Draws are reproducible for a given linear algebra
  library and thread count, but not across them: the samplers draw their
  one-sided terms by rejection, so a last-bit difference in a BLAS or LAPACK
  result can change how many random deviates a sweep consumes and shift every
  draw after it. The chains remain draws from the same posterior. Pin the
  thread count alongside the seed where bit-identical output is needed.
* `citation("bsfa")` now reports the version that is installed rather than a
  version written into the file by hand, which between releases was the last
  released one and so attributed results to a package the reader could
  install and not reproduce.

# bsfa 0.2.0

## Priors

* `add_priors()` now matches `r_star` on the median of the marginal prior of
  the inefficiency term, after the parameter of its distribution has been
  integrated out, so that the prior median efficiency is `r_star` exactly for
  both models and at any shape.

  This changes the priors that were released in 0.1.0, and so the posteriors
  drawn under them. The half-normal anchor was previously matched through the
  prior mean of the squared scale, which left it optimistic by a factor that
  depended only on the shape: at the default shape of 2.5 the implied prior
  median efficiency was `r_star^0.834`, so that asking for 0.75 gave 0.787 and
  asking for 0.5 gave 0.561. The exponential anchor was matched exactly at the
  default shape of 1, and only there; at shape 2 an `r_star` of 0.75 implied a
  prior median efficiency of 0.888.

  An exponential prior at the default shape is unaffected, which is the
  elicitation of van den Broeck, Koop, Osiewalski and Steel (1994) and was
  already exact. Every half-normal prior moves, as does any exponential prior
  with a shape other than 1. The efficiency scores of a sample that is
  informative about the level will barely notice; those of one that is not
  will follow the anchor, which is what `prior_sensitivity()` measures.

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

## Diagnostics

* `skewness_test()` tests the least squares residuals for skewness in the
  direction the frontier implies, which is the standard check for whether a
  sample supports a one-sided term at all. It matters more here than for a
  maximum likelihood fit: that one collapses onto least squares and says so,
  whereas a posterior cannot, and returns plausible looking efficiency scores
  read off the prior instead. Its standard error comes from a bootstrap over
  the units rather than from the closed form, which assumes the residuals are
  independent and so rejects too often on a panel.

* `prior_sensitivity()` refits the model across a grid of prior median
  efficiencies and reports how far the scores follow the anchor, which is the
  same question asked of the posterior rather than of the residuals. Every fit
  uses one seed, so that the table shows the prior moving the answer rather
  than the sampler.

## Simulation

* `sim_sf()` and `sim_sf4()` generate artificial data from the two- and
  four-component models, for examples, tests and prior predictive checks.

## Notes

* A term wrapped in `offset()` enters the frontier with its coefficient fixed
  at one, as it does in `lm()`, which is how a known capacity or a known
  elasticity is imposed rather than estimated.
* A rank deficient design matrix is refused, naming the column that is a
  linear combination of the others, rather than left for the prior to fill in.
* The response has to be a quantity. A factor would otherwise be fitted as the
  integers coding its levels and a character vector as a column of `NA`.
* An `id` is checked against the model frame rather than against `data`, so a
  mismatched one is refused whether `data` is a data frame or a list.
* Rows with a missing value, or with an unknown unit, are dropped with a
  message and recorded in the `na.action` element of the model object.
* The response is assumed to be on a logarithmic scale. The prior on the
  inefficiency term is expressed in the units of the response, so fitting
  output in levels does not merely change the units of the efficiency scores,
  it changes their values.
