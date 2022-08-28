#!/usr/bin/bash

set -x
docker push ghcr.io/llnl/sundials-ci-spack-e4s-base:latest
docker push ghcr.io/llnl/sundials-ci-spack-e4s-base:e4s-22.05
docker push ghcr.io/llnl/sundials-ci-spack-e4s:latest
docker push ghcr.io/llnl/sundials-ci-spack-e4s:e4s-22.05
docker push ghcr.io/llnl/sundials-ci-spack-develop:spack-develop
