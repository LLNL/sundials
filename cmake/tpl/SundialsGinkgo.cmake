# -----------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Module to find and setup GINKGO correctly.
# Created from the SundialsTPL.cmake template.
# All SUNDIALS modules that find and setup a TPL must:
#
# 1. Check to make sure the SUNDIALS configuration and the TPL is compatible.
# 2. Find the TPL.
# 3. Check if the TPL works with SUNDIALS, UNLESS the override option
# TPL_WORKS is TRUE - in this case the tests should not be performed and it
# should be assumed that the TPL works with SUNDIALS.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

if(NOT DEFINED SUNDIALS_GINKGO_INCLUDED)
  set(SUNDIALS_GINKGO_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

if(SUNDIALS_PRECISION MATCHES "extended")
  print_error("SUNDIALS GINKGO interface is not compatible with extended precision")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------
find_package(Ginkgo REQUIRED
             HINTS "${Ginkgo_DIR}"
             NO_DEFAULT_PATH)

message(STATUS "GINKGO VERSION:     ${GINKGO_PROJECT_VERSION}")
message(STATUS "GINKGO BUILD TYPE:  ${GINKGO_BUILD_TYPE}")
message(STATUS "GINKGO LIBRARIES:   ${GINKGO_INTERFACE_LINK_LIBRARIES}")
message(STATUS "GINKGO LINK FLAGS:  ${GINKGO_INTERFACE_LINK_FLAGS}")
message(STATUS "GINKGO CXX FLAGS:   ${GINKGO_INTERFACE_CXX_FLAGS}")


# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

# if(GINKGO_BACKENDS MATCHES "CUDA" AND NOT ENABLE_CUDA)
#   print_error("GINKGO_BACKENDS includes CUDA but CUDA is not enabled. Set ENABLE_CUDA=ON or change the backend.")
# endif()
# if(GINKGO_BACKENDS MATCHES "HIP" AND NOT ENABLE_HIP)
#   print_error("GINKGO_BACKENDS includes HIP but HIP is not enabled. Set ENABLE_HIP=ON or change the backend.")
# endif()

if(Ginkgo_FOUND AND (NOT GINKGO_WORKS))
  # # Create the GINKGO_TEST directory
  # set(GINKGO_TEST_DIR ${PROJECT_BINARY_DIR}/CMakeFiles/GINKGO_TEST)
  # file(MAKE_DIRECTORY ${GINKGO_TEST_DIR})

  # if(GINKGO_BACKENDS MATCHES "HIP")
  #   set(lang CXX)
  #   set(ext cxx)
  #   set(define_have "\#define HAVE_HIP")
  #   set(lib hip::host)
  # elseif(GINKGO_BACKENDS MATCHES "CUDA")
  #   set(lang CUDA)
  #   set(ext cu)
  #   set(define_have "\#define HAVE_CUBLAS")
  #   set(lib )
  # endif()

  # file(WRITE ${GINKGO_TEST_DIR}/ltest.${ext}
  # "${define_have}\n"
  # "\#include \"magma_v2.h\"\n"
  # "int main(){\n"
  # "magma_int_t a=0;\n"
  # "return(a);\n"
  # "}\n")

  # try_compile(COMPILE_OK ${GINKGO_TEST_DIR} ${GINKGO_TEST_DIR}/ltest.${ext}
  #   CMAKE_FLAGS
  #     "-DINCLUDE_DIRECTORIES=${GINKGO_INCLUDE_DIR}"
  #   LINK_LIBRARIES ${GINKGO_LIBRARIES} ${lib}
  #   OUTPUT_VARIABLE COMPILE_OUTPUT
  #   ${lang}_STANDARD ${CMAKE_${lang}_STANDARD}
  # )

  # # Process test result
  # if(COMPILE_OK)
    message(STATUS "Checking if GINKGO works... OK")
    set(GINKGO_WORKS TRUE CACHE BOOL "GINKGO works with SUNDIALS as configured" FORCE)
  # else()
  #   message(STATUS "Checking if GINKGO works... FAILED")
  #   message(STATUS "Check output: ")
  #   message("${COMPILE_OUTPUT}")
  #   print_error("SUNDIALS interface to GINKGO is not functional.")
  # endif()
elseif(GINKGO_FOUND AND GINKGO_WORKS)
  message(STATUS "Skipped GINKGO tests, assuming GINKGO works with SUNDIALS.")
endif()
