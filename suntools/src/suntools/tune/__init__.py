#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------

"""Tune SUNDIALS applications by sweeping ``SetOptions`` key/value pairs.

The :mod:`suntools.tune` package exposes validated configuration models.  The
command-line runner and optional DeepHyper, GPTune, and Ytopt integrations are
available from their respective submodules.
"""

from suntools.tune.models import (
    BackendConfig,
    ConstraintConfig,
    ExecutableConfig,
    MetricConfig,
    ObjectiveConfig,
    ParameterSpec,
    SearchConfig,
    TuneConfig,
)

__all__ = [
    "BackendConfig",
    "ConstraintConfig",
    "ExecutableConfig",
    "MetricConfig",
    "ObjectiveConfig",
    "ParameterSpec",
    "SearchConfig",
    "TuneConfig",
]
