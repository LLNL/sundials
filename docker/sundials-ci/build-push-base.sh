#!/usr/bin/bash

export DOCKER_BUILDKIT=1
export BUILDKIT_PROGRESS=plain
set -x
docker build -t ghcr.io/llnl/sundials-ci-e4s-base:22.05:22.05 e4s-base/
docker push ghcr.io/llnl/sundials-ci-e4s-base:22.05:22.05
