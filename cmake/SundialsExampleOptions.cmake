# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Options for SUNDIALS examples.
# ---------------------------------------------------------------

# -----------------------------------------------------------------------------
# Options for C/C++ examples
# -----------------------------------------------------------------------------

sundials_option(EXAMPLES_ENABLE_C BOOL "Build SUNDIALS C examples" ON)

# Some TPLs only have C++ examples. Default the C++ examples to ON if any of
# these are enabled on the initial configuration pass.
if (ENABLE_TRILINOS OR ENABLE_SUPERLUDIST OR ENABLE_XBRAID OR ENABLE_HIP OR
    ENABLE_MAGMA OR ENABLE_SYCL OR ENABLE_ONEMKL OR ENABLE_RAJA)
  sundials_option(EXAMPLES_ENABLE_CXX BOOL "Build SUNDIALS C++ examples" ON)
else()
  sundials_option(EXAMPLES_ENABLE_CXX BOOL "Build SUNDIALS C++ examples" OFF)
endif()

# -----------------------------------------------------------------------------
# Options for Fortran Examples
# -----------------------------------------------------------------------------

# F2003 examples (on by default) are an option only if the
# Fortran 2003 interface is enabled.
set(DOCSTR "Build SUNDIALS Fortran 2003 examples")
if(BUILD_FORTRAN_MODULE_INTERFACE)

  sundials_option(EXAMPLES_ENABLE_F2003 BOOL "${DOCSTR}" ON)

  # Fortran 2003 examples only support double precision
  if(EXAMPLES_ENABLE_F2003 AND (NOT (SUNDIALS_PRECISION MATCHES "DOUBLE")))
    print_warning("F2003 examples are not compatible with ${SUNDIALS_PRECISION} precision. "
                  "Setting EXAMPLES_ENABLE_F2003 to OFF.")
    force_variable(EXAMPLES_ENABLE_F2003 BOOL "${DOCSTR}" OFF)
  endif()

  # Fortran 2003 examples only support 64-bit indices
  if(EXAMPLES_ENABLE_F2003 AND (NOT (SUNDIALS_INDEX_SIZE MATCHES "64")))
    print_warning("F2003 examples are not compatible with ${SUNDIALS_INDEX_SIZE}-bit indices. "
                  "Setting EXAMPLES_ENABLE_F2003 to OFF.")
    force_variable(EXAMPLES_ENABLE_F2003 BOOL "${DOCSTR}" OFF)
  endif()

else()

  # set back to OFF (in case it was ON)
  if(EXAMPLES_ENABLE_F2003)
    print_warning("EXAMPLES_ENABLE_F2003 is ON but BUILD_FORTRAN_MODULE_INTERFACE is OFF. "
                  "Setting EXAMPLES_ENABLE_F2003 to OFF.")
    force_variable(EXAMPLES_ENABLE_F2003 BOOL "${DOCSTR}" OFF)
  endif()

endif()

# -----------------------------------------------------------------------------
# Options for CUDA Examples
# -----------------------------------------------------------------------------

sundials_option(EXAMPLES_ENABLE_CUDA BOOL "Build SUNDIALS CUDA examples" ON
                DEPENDS_ON ENABLE_CUDA)

# -----------------------------------------------------------------------------
# Options for installing examples
# -----------------------------------------------------------------------------

# Enable installing examples by default
sundials_option(EXAMPLES_INSTALL BOOL "Install SUNDIALS examples" ON)

sundials_option(EXAMPLES_INSTALL_PATH PATH "Output directory for installing example files" "${CMAKE_INSTALL_PREFIX}/examples")

# If examples are to be exported, check where we should install them.
if(EXAMPLES_INSTALL AND NOT EXAMPLES_INSTALL_PATH)
  print_warning("The example installation path is empty. "
                "Example installation path was reset to its default value")
  set(EXAMPLES_INSTALL_PATH "${CMAKE_INSTALL_PREFIX}/examples" CACHE STRING
      "Output directory for installing example files" FORCE)
endif()

# -----------------------------------------------------------------------------
# Internal variables.
# -----------------------------------------------------------------------------

if(EXAMPLES_ENABLE_C OR
   EXAMPLES_ENABLE_CXX OR
   EXAMPLES_ENABLE_CUDA OR
   EXAMPLES_ENABLE_F2003)
  set(_BUILD_EXAMPLES TRUE CACHE INTERNAL "")
else()
  set(_BUILD_EXAMPLES FALSE CACHE INTERNAL "")
endif()
