#!/usr/bin/bash

images=('int32-double' 'int32-extended' 'int32-single' 'int64-double' 'int64-extended' 'int64-single')
set -x
for image in ${images[@]}; do
  docker push ghcr.io/llnl/sundials-ci-${image}:latest
  docker push ghcr.io/llnl/sundials-ci-${image}:e4s-22.05
  docker push ghcr.io/llnl/sundials-ci-${image}:spack-develop
done
