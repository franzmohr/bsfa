# Image for running the checks of .github/workflows/R-CMD-check.yaml locally.
#
# It stands in for the ubuntu-latest / R release job of that workflow: a clean
# Linux system with the current R release, the system tools the check needs and
# a package library built from DESCRIPTION.
#
# One NOTE is expected here and nowhere else:
#
#   checking compilation flags used ... NOTE
#     Compilation used the following non-portable flag(s):
#       '-Wdate-time' '-Werror=format-security' '-Wformat'
#
# Those come from the Debian R build this image is based on, which puts them
# in the system Makeconf; src/Makevars sets nothing but PKG_LIBS. Every
# package built here collects the NOTE, and CRAN and the GitHub runners do
# not report it. tools/check.R fails on warnings rather than notes, so it
# does not block a push.
FROM rocker/r-ver:latest

# pandoc builds the vignettes and is read by the top-level-files check,
# ghostscript and qpdf are what --compact-vignettes=gs+qpdf calls.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ghostscript \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        pandoc \
        qpdf \
    && rm -rf /var/lib/apt/lists/*

# Dependencies are taken from DESCRIPTION, so adding or dropping one is enough
# to make the next build install it; nothing here has to be kept in sync by hand.
COPY DESCRIPTION /deps/DESCRIPTION
RUN Rscript -e 'install.packages(c("remotes", "rcmdcheck"))' \
    && Rscript -e 'remotes::install_deps("/deps", dependencies = TRUE)'
