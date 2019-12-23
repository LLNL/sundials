# ---------------------------------------------------------------
# Programmer:  Cody Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Based on the FindPackageMultipass module by Jed Brown.
# ---------------------------------------------------------------
# Copyright Jed Brown
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in
#   the documentation and/or other materials provided with the
#   distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ---------------------------------------------------------------
# PackageMultipass - this module defines two macros
#
# FIND_PACKAGE_MULTIPASS (Name CURRENT
#  STATES VAR0 VAR1 ...
#  DEPENDENTS DEP0 DEP1 ...)
#
#  This function creates a cache entry <UPPERCASED-Name>_CURRENT which
#  the user can set to "NO" to trigger a reconfiguration of the package.
#  The first time this function is called, the values of
#  <UPPERCASED-Name>_VAR0, ... are saved.  If <UPPERCASED-Name>_CURRENT
#  is false or if any STATE has changed since the last time
#  FIND_PACKAGE_MULTIPASS() was called, then CURRENT will be set to "NO",
#  otherwise CURRENT will be "YES".  IF not CURRENT, then
#  <UPPERCASED-Name>_DEP0, ... will be FORCED to NOTFOUND.
#  Example:
#    find_path (FOO_DIR include/foo.h)
#    FIND_PACKAGE_MULTIPASS (Foo foo_current
#      STATES DIR
#      DEPENDENTS INCLUDES LIBRARIES)
#    if (NOT foo_current)
#      # Make temporary files, run programs, etc, to determine FOO_INCLUDES and FOO_LIBRARIES
#    endif (NOT foo_current)
#
# MULTIPASS_SOURCE_RUNS (Name INCLUDES LIBRARIES SOURCE RUNS LANGUAGE)
#  Always runs the given test, use this when you need to re-run tests
#  because parent variables have made old cache entries stale. The LANGUAGE
#  variable is either C or CXX indicating which compiler the test should
#  use.
# MULTIPASS_C_SOURCE_RUNS (Name INCLUDES LIBRARIES SOURCE RUNS)
#  DEPRECATED! This is only included for backwards compatability. Use
#  the more general MULTIPASS_SOURCE_RUNS instead.
#  Always runs the given test, use this when you need to re-run tests
#  because parent variables have made old cache entries stale.
# ---------------------------------------------------------------

macro (FIND_PACKAGE_MULTIPASS _name _current)
  string (TOUPPER ${_name} _NAME)
  set (_args ${ARGV})
  list (REMOVE_AT _args 0 1)

  set (_states_current "YES")
  list (GET _args 0 _cmd)
  if (_cmd STREQUAL "STATES")
    list (REMOVE_AT _args 0)
    list (GET _args 0 _state)
    while (_state AND NOT _state STREQUAL "DEPENDENTS")
      # The name of the stored value for the given state
      set (_stored_var PACKAGE_MULTIPASS_${_NAME}_${_state})
      if (NOT "${${_stored_var}}" STREQUAL "${${_NAME}_${_state}}")
        set (_states_current "NO")
      endif (NOT "${${_stored_var}}" STREQUAL "${${_NAME}_${_state}}")
      set (${_stored_var} "${${_NAME}_${_state}}" CACHE INTERNAL "Stored state for ${_name}." FORCE)
      list (REMOVE_AT _args 0)
      list (GET _args 0 _state)
    endwhile (_state AND NOT _state STREQUAL "DEPENDENTS")
  endif (_cmd STREQUAL "STATES")

  set (_stored ${_NAME}_CURRENT)
  if (NOT ${_stored})
    set (${_stored} "YES" CACHE BOOL "Is the configuration for ${_name} current?  Set to \"NO\" to reconfigure." FORCE)
    set (_states_current "NO")
  endif (NOT ${_stored})

  set (${_current} ${_states_current})
  if (NOT ${_current} AND PACKAGE_MULTIPASS_${_name}_CALLED)
    message (STATUS "Clearing ${_name} dependent variables")
    # Clear all the dependent variables so that the module can reset them
    list (GET _args 0 _cmd)
    if (_cmd STREQUAL "DEPENDENTS")
      list (REMOVE_AT _args 0)
      foreach (dep ${_args})
        set (${_NAME}_${dep} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
      endforeach (dep)
    endif (_cmd STREQUAL "DEPENDENTS")
    set (${_NAME}_FOUND "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
  endif ()
  set (PACKAGE_MULTIPASS_${name}_CALLED YES CACHE INTERNAL "Private" FORCE)
endmacro (FIND_PACKAGE_MULTIPASS)


macro (MULTIPASS_SOURCE_RUNS includes libraries source runs language)
  # Cody Balos on 7/2019: CHECK_<lang>_SOURCE_RUNS does not allow for the
  # compiler to be set which causes problems if MPI is not the C compiler.
  # Therefore, this has been changed to use try_run instead. Everything
  # else was kept so FindPETSC did not need to be modified.

  # This is a ridiculous hack. CHECK_${language}_SOURCE_*
  # thinks that if the *name* of the return variable doesn't change,
  # then the test does not need to be re-run.  We keep an internal count
  # which we increment to guarantee that every test name is unique. If we've
  # gotten here, then the configuration has changed enough that the
  # test *needs* to be rerun.
  if (NOT MULTIPASS_TEST_COUNT)
    set (MULTIPASS_TEST_COUNT 00)
  endif (NOT MULTIPASS_TEST_COUNT)
  math (EXPR _tmp "${MULTIPASS_TEST_COUNT} + 1") # Why can't I add to a cache variable?
  set (MULTIPASS_TEST_COUNT ${_tmp} CACHE INTERNAL "Unique test ID")
  set (testname MULTIPASS_TEST_${MULTIPASS_TEST_COUNT}_${runs})
  set (testdir ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp)
  set (CMAKE_REQUIRED_INCLUDES ${includes})
  set (CMAKE_REQUIRED_LIBRARIES ${libraries})
  # if MPI is available, use it for the test
  if (MPI_${language}_COMPILER)
    set (REQUIRED_COMPILER ${MPI_${language}_COMPILER})
  else ()
    set (REQUIRED_COMPILER ${CMAKE_${language}_COMPILER})
  endif ()
  if(${language} STREQUAL "C")
    file(WRITE ${testdir}/src.c "${source}")
    try_run(${testname} _compiles ${testdir} ${testdir}/src.c
            CMAKE_FLAGS -DCMAKE_C_COMPILER=${REQUIRED_COMPILER} -DINCLUDE_DIRECTORIES=${CMAKE_REQUIRED_INCLUDES}
            LINK_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
  elseif(${language} STREQUAL "CXX")
    file(WRITE ${testdir}/src.cxx "${source}")
    try_run(${testname} _compiles ${testdir} ${testdir}/src.cxx
            CMAKE_FLAGS -DCMAKE_CXX_COMPILER=${REQUIRED_COMPILER} -DINCLUDE_DIRECTORIES=${CMAKE_REQUIRED_INCLUDES}
            LINK_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
  endif()
  # ${testname} is the exit code returned by try_run,
  # so 0 is success and anything else is a failure.
  if (${testname})
    set (${runs} FALSE)
  else ()
    set (${runs} TRUE)
  endif ()
  unset (_compiles)
endmacro (MULTIPASS_SOURCE_RUNS)


macro (MULTIPASS_C_SOURCE_RUNS includes libraries source runs)
  multipass_source_runs("${includes}" "${libraries}" "${source}" ${runs} "C")
endmacro (MULTIPASS_C_SOURCE_RUNS)
