#!/bin/bash

python generate.py sundials/generate.yaml
python generate.py arkode/generate.yaml
python generate.py cvodes/generate.yaml
python generate.py idas/generate.yaml
python generate.py kinsol/generate.yaml
