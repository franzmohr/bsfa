# Contributing to bsfa

Bug reports, questions and pull requests are all welcome. A few things about
this package are worth knowing before you start, because they decide *where* a
change belongs.

## How the package is put together

A model is not estimated in one call. A constructor builds a model object, the
`add_*` functions attach one block of the specification each, and the draws are
attached to the same object rather than returned as a separate result:

```r
model <- create_sfmodel_exp(y ~ x1 + x2, data = d,
                            iterations = 5000, burnin = 2000)
model <- add_priors(model)
model <- add_initial_values(model)
model <- add_seed(model, 1234)
model <- add_posterior_coefficients(model)
model <- add_posterior_loglik(model)
```

This mirrors [bvartools](https://github.com/franzmohr/bvartools), so that anyone
who knows one package can move to the other. A change that adds a step should
add it as an `add_*` generic rather than as an argument to an existing function.

**Model variants are classes, not arguments.** The exponential and half-normal
specifications are `sfmodel_exp` and `sfmodel_hn`, and the four-component model
adds `sfmodel4_exp` and `sfmodel4_hn`, all inheriting from `sfmodel`. Each gets
its own method wherever the models genuinely differ — most visibly in
`add_priors()`, since the prior median efficiency maps into a rate differently
for each distribution. A string argument that switches behaviour inside one
function is the thing this structure exists to avoid.

Where a difference is mathematical rather than interfacial, it can stay in one
place: `sf_loglik_point()` holds both composed-error densities side by side,
where they can be compared, and the dispatch lives at the level of the exported
functions.

## Where the sampler lives

`src/gibbs_sf.cpp` holds the two-component sampler and `src/gibbs_sf4.cpp` the
four-component one. Both are plain Gibbs loops with data augmentation: the
one-sided terms are drawn from truncated normals, the coefficients from a
normal, and the variances from gammas.

If you change a full conditional, derive it in a comment above the code and add
a test that pins the result down against something independent — the existing
tests check the closed-form log-likelihoods against direct numerical
integration of the composed-error density, and the Koop vignette checks the
exponential panel sampler against the reference MATLAB program distributed with
*Bayesian Econometric Methods*. A conditional that is subtly wrong will still
produce plausible-looking draws, so "it ran and the numbers seemed fine" is not
evidence.

## Setting up

Building the package needs a C++ toolchain — Rtools on Windows, Xcode command
line tools on macOS — because of the RcppArmadillo code in `src/`.

```r
# install.packages("remotes")
remotes::install_deps("path/to/bsfa", dependencies = TRUE)
```

## Making a change

* **Documentation is generated.** `man/*.Rd` and `NAMESPACE` come from the
  roxygen blocks in `R/`. Edit the `R/` file and run `devtools::document()`;
  never edit an `.Rd` by hand.
* **`src/RcppExports.*` are generated too**, by `Rcpp::compileAttributes()`.
  Run it before `document()` whenever a `// [[Rcpp::export]]` signature changes.
* **Tests.** `tests/testthat/`, testthat edition 3. New behaviour comes with a
  test; a bug fix comes with a test that failed before it. Recovery tests on
  simulated data are welcome but should assert what the model actually
  identifies. The four-component tests are the example worth copying: the split
  between the unit effect and persistent inefficiency is weakly determined, so
  the test asserts that `mu_i - eta_i` is recovered and that *which* of the two
  comes out well depends on their relative spread, rather than putting an
  arbitrary threshold on either.
* **Check before opening a pull request.** `R CMD check` should be clean:

  ```r
  devtools::check()
  ```

* **Style.** Follow the surrounding code rather than a style guide: two-space
  indent, `<-` for assignment, `stats::`/`coda::` prefixes on imported
  functions, and comments that say why rather than what.
* **NEWS.md.** A user-visible change gets an entry under the development
  version.

## Adding a specification

The roadmap in `README.md` lists what is planned. Most of it is a block of the
specification rather than a new algorithm, which is the test of whether it fits:
inefficiency determinants, heteroskedastic error components and latent class
membership all leave the full conditionals standard, so they extend the existing
loop. If a proposal needs a Metropolis step or a numerical integral, say so in
the issue, because that changes what the sampler has to guarantee.

Two things are deliberately absent and would each be a substantial contribution
in its own right:

* **LOOIC**, which needs Pareto smoothed importance sampling.
* **The closed skew normal likelihood** of Colombi et al. (2014), which would
  make information criteria available for the four-component model. It needs
  normal distribution functions of dimension `T_i + 1`.

## Reporting a bug

Include a reproducible example built on `sim_sf()` or `sim_sf4()`, plus the
output of `sessionInfo()`. Give the whole chain, priors included: which code
paths run depends on the model class, on whether `id` was supplied and on
whether the frontier is a production or a cost frontier.

If the report is about draws being wrong rather than an error, say which
parameter block and at what sample size. Frontier intercepts and the parameters
of the one-sided distributions are the weakly determined part of every model
here, and a posterior that looks off may be the shape of the problem rather than
a defect — the documentation of each model says which quantities that applies
to.

## Licensing

Contributions are accepted under the package's license, GPL (>= 2).
