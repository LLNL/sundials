# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# CMakeLists.txt file for ARKODE C++ standard library vector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases. The examples use the C++17
# (and later) NVectorBase and the tests do not compare against answer files.
set(ARKODE_examples)

# Header files to install
set(ARKODE_headers)

# std::vector examples require the C++17 NVectorBase
if(CMAKE_CXX_STANDARD GREATER_EQUAL 17)
  list(APPEND ARKODE_examples "harmonic_oscillator.cpp\;\;develop"
       "harmonic_oscillator_serial.cpp\;\;develop")
  list(APPEND ARKODE_headers "nvector_stdvector.hpp")
else()
  message(
    STATUS
      "Skipping ARKODE CXX_std std::vector examples (requires CMAKE_CXX_STANDARD >= 17)"
  )
endif()

# std::mdspan examples require C++23 <mdspan> with multi-index operator[]
# (libstdc++ and libc++ both gate <mdspan> on C++23). The check is compiled with
# the project's C++ standard, and its cached result is keyed on that standard so
# changing CMAKE_CXX_STANDARD in place re-runs it.
if(CMAKE_CXX_STANDARD GREATER_EQUAL 23)
  include(CheckCXXSourceCompiles)
  set(_mdspan_probe
      "
#include <mdspan>
#if !defined(__cpp_lib_mdspan)
#error \"std::mdspan is not available\"
#endif
int main()
{
  double data[4] = {0.0, 0.0, 0.0, 0.0};
  std::mdspan<double, std::dextents<int, 2>> view(data, 2, 2);
  return view[1, 1] == 0.0 ? 0 : 1;
}
")
  set(_mdspan_var SUNDIALS_EXAMPLES_HAVE_STD_MDSPAN_CXX${CMAKE_CXX_STANDARD})
  check_cxx_source_compiles("${_mdspan_probe}" ${_mdspan_var})
  set(_have_std_mdspan ${${_mdspan_var}})
  unset(_mdspan_probe)
  unset(_mdspan_var)
else()
  set(_have_std_mdspan FALSE)
endif()

if(_have_std_mdspan)
  list(APPEND ARKODE_examples "coupled_oscillators.cpp\;\;develop"
       "coupled_oscillators_serial.cpp\;\;develop")
  list(APPEND ARKODE_headers "nvector_stdmdspan.hpp")
else()
  message(
    STATUS
      "Skipping ARKODE CXX_std std::mdspan examples (requires std::mdspan with multi-index operator[], i.e., C++23)"
  )
endif()
unset(_have_std_mdspan)

# Add the build and install targets for each example
foreach(example_tuple ${ARKODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  # add example executable
  if(NOT TARGET ${example_target})
    sundials_add_executable(${example_target} ${example})

    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

    # directories to include
    target_include_directories(${example_target}
                               PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities)

    # libraries to link against
    target_link_libraries(${example_target} sundials_arkode sundials_nvecserial
                          ${EXE_EXTRA_LINK_LIBS})
  endif()

  # add example to regression tests, there are no answer files to compare
  sundials_add_test(${example_target} ${example_target} NODIFF EXAMPLE_TYPE
                    ${example_type})

endforeach(example_tuple ${ARKODE_examples})

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL AND ARKODE_examples)

  sundials_install_examples(
    arkode ARKODE_examples
    CMAKE_TEMPLATE cmakelists_CXX_ex.in
    MAKE_TEMPLATE makefile_serial_CXX_ex.in
    SUNDIALS_TARGETS nvecserial arkode
    DESTINATION arkode/CXX_std
    EXTRA_FILES ${ARKODE_headers}
                "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
    TEST_INSTALL CXX_std)

endif()
