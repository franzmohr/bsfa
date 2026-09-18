<!--
Thanks for the pull request. .github/CONTRIBUTING.md has the details; the
checklist below is the short version.
-->

## What this changes

<!-- One or two sentences, and the issue number if there is one. -->

## Checklist

- [ ] `devtools::check()` is clean (no errors, warnings or new notes)
- [ ] `devtools::test()` passes, and new behaviour or a fixed bug comes with a test
- [ ] `devtools::document()` was run if any roxygen block changed — `man/` and
      `NAMESPACE` are generated, not edited
- [ ] `Rcpp::compileAttributes()` was run if any `// [[Rcpp::export]]` changed
- [ ] `NEWS.md` has an entry if the change is user-visible
- [ ] A changed full conditional is derived in a comment and pinned down by a
      test against something independent — numerical integration, a reference
      implementation, or a known-truth simulation
- [ ] A new step in the workflow is an `add_*` generic, and a new model variant
      is a class with its own methods, rather than an argument that switches
      behaviour inside an existing function
