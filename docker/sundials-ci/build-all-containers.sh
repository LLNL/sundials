#!/usr/bin/bash

export DOCKER_BUILDKIT=1
export BUILDKIT_PROGRESS=plain
set -x
docker build -t ghcr.io/llnl/sundials-ci-spack-e4s-base:latest e4s-base/
docker build -t ghcr.io/llnl/sundials-ci-spack-e4s-base:e4s-22.05 e4s-base/
docker build -t ghcr.io/llnl/sundials-ci-spack-e4s:latest --build-arg spack_yaml="spack.yaml" spack-e4s/
docker build -t ghcr.io/llnl/sundials-ci-spack-develop:spack-develop --build-arg spack_yaml="spack.yaml" spack-develop/
