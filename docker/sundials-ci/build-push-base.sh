#!/usr/bin/bash

export DOCKER_BUILDKIT=1
export BUILDKIT_PROGRESS=plain

docker build -t balos1/sundials-ci-e4s-base:22.05:22.05 e4s-base/
docker push balos1/sundials-ci-e4s-base:22.05:22.05