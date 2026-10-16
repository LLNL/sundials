# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
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
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Option for testing deprecated options
# ---------------------------------------------------------------

if(DEFINED ENV{SUNDIALS_ENABLE_UNSET_DEPRECATED})
  set(_sundials_enable_unset_deprecated_default
      "$ENV{SUNDIALS_ENABLE_UNSET_DEPRECATED}")
else()
  set(_sundials_enable_unset_deprecated_default OFF)
endif()

sundials_option(
  SUNDIALS_ENABLE_UNSET_DEPRECATED BOOL "Unset deprecated CMake options"
  ${_sundials_enable_unset_deprecated_default} ADVANCED)

if(SUNDIALS_ENABLE_UNSET_DEPRECATED)
  message(WARNING "Unsetting deprecated SUNDIALS options.")
endif()

# ---------------------------------------------------------------
# Deprecated options not handled by DEPRECATED_NAMES
# ---------------------------------------------------------------
