#!/usr/bin/bash

set -x
docker push ghcr.io/llnl/sundials-ci-spack-latest:spack-latest
docker push ghcr.io/llnl/sundials-ci-spack-develop:spack-develop
