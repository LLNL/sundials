#!/usr/bin/bash

images=('int32-double' 'int32-extended' 'int32-single' 'int64-double' 'int64-extended' 'int64-single')
set -x
for image in ${images[@]}; do
  docker push balos1/sundials-ci-${image}:latest
done