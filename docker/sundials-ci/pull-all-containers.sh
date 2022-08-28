#!/usr/bin/bash

set -x
docker pull ghcr.io/llnl/sundials-ci-spack-e4s-base:latest
docker pull ghcr.io/llnl/sundials-ci-spack-e4s-base:e4s-22.05
docker pull ghcr.io/llnl/sundials-ci-spack-e4s:latest
docker pull ghcr.io/llnl/sundials-ci-spack-e4s:e4s-22.05
docker pull ghcr.io/llnl/sundials-ci-spack-develop:spack-develop
