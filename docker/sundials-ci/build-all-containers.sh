#!/usr/bin/bash

export DOCKER_BUILDKIT=1
export BUILDKIT_PROGRESS=plain
set -x
docker build -t ghcr.io/llnl/sundials-ci-spack-latest:spack-latest --build-arg spack_yaml="spack.yaml" spack-latest/
docker build -t ghcr.io/llnl/sundials-ci-spack-develop:spack-develop --build-arg spack_yaml="spack.yaml" spack-develop/
