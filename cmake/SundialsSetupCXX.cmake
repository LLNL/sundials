# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
#                Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2024, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# C++-related tests for SUNDIALS CMake-based configuration.
# ---------------------------------------------------------------

enable_language(CXX)
set(CXX_FOUND TRUE)

# ---------------------------------------------------------------
# Option to specify the C++ standard SUNDIALS will use. Defined
# here so it is set in the same configuration pass as the C++
# compiler and related options.
# ---------------------------------------------------------------

# Do not allow decaying to previous standards -- generates error if the standard
# is not supported
sundials_option(CMAKE_CXX_STANDARD_REQUIRED BOOL
  "Require C++ standard version" ON)

if(ENABLE_SYCL)
  set(DOCSTR "The C++ standard to use if C++ is enabled (17, 20)")
  sundials_option(CMAKE_CXX_STANDARD STRING "${DOCSTR}" "17"
                  OPTIONS "17;20")
else()
  set(DOCSTR "The C++ standard to use if C++ is enabled (14, 17, 20)")
  sundials_option(CMAKE_CXX_STANDARD STRING "${DOCSTR}" "14"
                  OPTIONS "14;17;20")
endif()
message(STATUS "CXX standard set to ${CMAKE_CXX_STANDARD}")

set(DOCSTR "Enable C++ compiler specific extensions")
sundials_option(CMAKE_CXX_EXTENSIONS BOOL "${DOCSTR}" ON)
message(STATUS "C++ extensions set to ${CMAKE_CXX_EXTENSIONS}")

# SYCL requries C++17
if(ENABLE_SYCL AND (CMAKE_CXX_STANDARD LESS "17"))
  message(SEND_ERROR "CMAKE_CXX_STANDARD must be >= 17 because ENABLE_SYCL=ON")
endif()
