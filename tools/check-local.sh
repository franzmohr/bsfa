#!/usr/bin/env bash
#
# Runs the GitHub check suite locally, in Docker, on the current working tree.
# Nothing should be pushed before this has passed.
#
#   bash tools/check-local.sh
#
set -euo pipefail

pkg_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
image=${BSFA_CHECK_IMAGE:-bsfa-check}

docker build -f "$pkg_dir/tools/Dockerfile" -t "$image" "$pkg_dir"

# The working tree is piped into the container rather than mounted: the check
# then never sees the objects and DLLs a Windows build leaves in src/, and it
# cannot write anything back into the source tree.
tar -cf - -C "$pkg_dir" \
    --exclude=./.git \
    --exclude=./.Rproj.user \
    --exclude=./.Rhistory \
    --exclude=./docs \
    --exclude='*.o' \
    --exclude='*.so' \
    --exclude='*.dll' \
    --exclude='*.dylib' \
    . |
  docker run --rm -i \
    -e _R_CHECK_CRAN_INCOMING_=false \
    "$image" \
    bash -c 'mkdir -p /pkg && tar -xf - -C /pkg && Rscript /pkg/tools/check.R'
