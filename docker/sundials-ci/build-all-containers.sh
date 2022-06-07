#!/usr/bin/bash

export DOCKER_BUILDKIT=1
export BUILDKIT_PROGRESS=plain

images=('int32-double' 'int32-extended' 'int32-single' 'int64-double' 'int64-extended' 'int64-single')
for image in ${images[@]}; do
  cd $image
  docker build -t balos1/sundials-ci:${image}-latest .
  cd -
done