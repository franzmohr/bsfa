# Runs inside the container built from tools/Dockerfile and repeats the check
# step of .github/workflows/R-CMD-check.yaml, down to the arguments r-lib's
# check-r-package action passes and the warnings it treats as failures.
rcmdcheck::rcmdcheck(
  "/pkg",
  args = c("--no-manual", "--as-cran"),
  build_args = c("--no-manual", "--compact-vignettes=gs+qpdf"),
  error_on = "warning"
)
