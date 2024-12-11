# ------------------------------------------------------------------------------
# Programmer(s): Steven Smith and David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2024, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------

# ~~~
# sundials_add_test(<test name> <executable>
#                   [NODIFF]
#                   [MPI_NPROCS num_processes]
#                   [FLOAT_PRECISION num_digits]
#                   [INTEGER_PRECISION percent_difference]
#                   [ANSWER_DIR path]
#                   [ANSWER_FIEL file]
#                   [EXAMPLE_TYPE type]
#                   [TEST_ARGS arg1 arg2 ...]
#                   [LABELS label1 label2 ...])
# ~~~
#
# CMake macro to add a SUNDIALS regression test. Keyword input arguments can be
# added after <executable> to set regression test options.
#
# The option NODIFF disables comparison of the test output against the answer
# file
#
# The option MPI_NPROCS sets the number of mpi tasks to use in parallel tests
#
# The option FLOAT_PRECISION set the precision (number of digits) for floating
# point failure comparisons. To use the default value, either don't provide the
# keyword, or provide the value "default".
#
# The option INTEGER_PRECISION sets the integer percentage difference for
# failure comparison.
#
# The option ANSWER_DIR sets the path to the directory containing the test
# answer file
#
# The option ANSWER_FILE set the name of test answer file
#
# The option EXAMPLE_TYPE set the example type i.e., release or develop examples
#
# The option TEST_ARGS are command line arguments to pass to the executable
#
# The options LABELS are labels added to the test properties to easily run (or
# exclude) groups of test with ctest -L <label> (or ctest -LE <label>)
#
# When SUNDIALS_TEST_DEVTESTS is OFF (default) the executable is run and success
# or failure is determined by the executable return value (zero or non-zero
# respectively).
#
# When SUNDIALS_TEST_DEVTESTS is ON the executable is run and its output is
# compared with the corresponding .out file. If the output differs significantly
# then the test fails. The default level of significance is 4 decimal places for
# floating point values and 10% for integer values.
#
# The level of precision can be adjusted for an individual test with the
# FLOAT_PRECISION AND INTEGER_PRECISION keyword inputs to the macro or globally
# for all tests with the cache variables SUNDIALS_TEST_FLOAT_PRECISION and
# SUNDIALS_TEST_INTEGER_PRECISION.
#
# -D SUNDIALS_TEST_FLOAT_PRECISION=<number of digits>
#
# -D SUNDIALS_TEST_INTEGER_PRECISION=<% difference>
#
# By default testing output is written to builddir/Testing/output and the .out
# answer file directory is set using the ANSWER_DIR keyword input to
# sourcedir/examples/package/testdir. These can be changed by setting the cache
# variables SUNDIALS_TEST_OUTPUT_DIR and SUNDIALS_TEST_ANSWER_DIR.
#
# -D SUNDIALS_TEST_OUTPUT_DIR=<path to output directory>
#
# -D SUNDIALS_TEST_ANSWER_DIR=<path to answer directory>
#
# By default the caliper output is written to builddir/Caliper. This can be
# changed by setting the cache variable SUNDIALS_CALIPER_OUTPUT_DIR.
#
# -D SUNDIALS_CALIPER_OUTPUT_DIR=<path to caliper output directory>

macro(SUNDIALS_ADD_TEST NAME EXECUTABLE)

  set(options "NODIFF")
  set(oneValueArgs "MPI_NPROCS" "FLOAT_PRECISION" "INTEGER_PRECISION"
                   "ANSWER_DIR" "ANSWER_FILE" "EXAMPLE_TYPE")
  set(multiValueArgs "LABELS" "TEST_ARGS" "EXTRA_ARGS")

  # parse inputs and create variables SUNDIALS_ADD_TEST_<keyword>
  cmake_parse_arguments(SUNDIALS_ADD_TEST "${options}" "${oneValueArgs}"
                        "${multiValueArgs}" ${ARGN})

  # check if any options for setting command line args have been provided by
  # comparing them against an empty string -- this is needed to avoid missing
  # single args that are a false constant in CMake e.g., 0, FALSE, OFF, etc.
  if("${SUNDIALS_ADD_TEST_TEST_ARGS}" STREQUAL "")
    set(_have_test_args FALSE)
  else()
    set(_have_test_args TRUE)
  endif()

  if("${SUNDIALS_ADD_TEST_EXTRA_ARGS}" STREQUAL "")
    set(_have_extra_test_args FALSE)
  else()
    set(_have_extra_test_args TRUE)
  endif()

  # check that the test is not excluded
  string(TOLOWER "exclude-${SUNDIALS_PRECISION}" _exclude_precision)
  if(("${SUNDIALS_ADD_TEST_EXAMPLE_TYPE}" STREQUAL "exclude")
     OR ("${SUNDIALS_ADD_TEST_EXAMPLE_TYPE}" STREQUAL _exclude_precision))

    message(
      STATUS
        "Skipped test ${NAME} because it had type ${SUNDIALS_ADD_TEST_EXAMPLE_TYPE}"
    )

  else()

    if(SUNDIALS_TEST_DEVTESTS)

      # run all tests (standard and develop) with the test runner

      # command line arguments for the test runner script
      set(TEST_ARGS "--verbose" "--testname=${NAME}"
                    "--executablename=$<TARGET_FILE:${EXECUTABLE}>")

      if(SUNDIALS_TEST_PROFILE)
        list(APPEND TEST_ARGS "--profile")
        if(SUNDIALS_CALIPER_OUTPUT_DIR)
          list(APPEND TEST_ARGS
               "--calidir=${SUNDIALS_CALIPER_OUTPUT_DIR}/Example/${JOB_ID}")
        else()
          list(APPEND TEST_ARGS "--calidir=${TEST_OUTPUT_DIR}/Caliper/Example")
        endif()
      endif()

      # check for a non-default output directory
      if(SUNDIALS_TEST_OUTPUT_DIR)
        list(APPEND TEST_ARGS "--outputdir=${SUNDIALS_TEST_OUTPUT_DIR}")
      else()
        list(APPEND TEST_ARGS "--outputdir=${TEST_OUTPUT_DIR}")
      endif()

      # set a non-default answer directory (default is test/answers)
      if(SUNDIALS_TEST_ANSWER_DIR)
        list(APPEND TEST_ARGS "--answerdir=${SUNDIALS_TEST_ANSWER_DIR}")
      elseif(SUNDIALS_ADD_TEST_ANSWER_DIR)
        list(APPEND TEST_ARGS "--answerdir=${SUNDIALS_ADD_TEST_ANSWER_DIR}")
      endif()

      # set the test answer file name (default is test_name_test_agrs)
      if(SUNDIALS_ADD_TEST_ANSWER_FILE)
        list(APPEND TEST_ARGS "--answerfile=${SUNDIALS_ADD_TEST_ANSWER_FILE}")
      endif()

      # check if a diff is needed and if non-default precisions were provided
      if(SUNDIALS_ADD_TEST_NODIFF OR SUNDIALS_TEST_NODIFF)
        # do not diff the output and answer files
        list(APPEND TEST_ARGS "--nodiff")
      else()
        # set a non-default floating point precision (number of digits, default
        # 4)
        if(SUNDIALS_ADD_TEST_FLOAT_PRECISION
           AND (NOT SUNDIALS_ADD_TEST_FLOAT_PRECISION MATCHES "DEFAULT|default"
               ))
          list(APPEND TEST_ARGS
               "--floatprecision=${SUNDIALS_ADD_TEST_FLOAT_PRECISION}")
        elseif(SUNDIALS_TEST_FLOAT_PRECISION GREATER_EQUAL "0")
          list(APPEND TEST_ARGS
               "--floatprecision=${SUNDIALS_TEST_FLOAT_PRECISION}")
        endif()
        # set a non-default integer precision (percent difference, default 10%)
        if(SUNDIALS_ADD_TEST_INTEGER_PRECISION
           AND (NOT SUNDIALS_ADD_TEST_INTEGER_PRECISION MATCHES
                "DEFAULT|default"))
          list(APPEND TEST_ARGS
               "--integerpercentage=${SUNDIALS_ADD_TEST_INTEGER_PRECISION}")
        elseif(SUNDIALS_TEST_INTEGER_PRECISION GREATER_EQUAL "0")
          list(APPEND TEST_ARGS
               "--integerpercentage=${SUNDIALS_TEST_INTEGER_PRECISION}")
        endif()
      endif()

      # check if this test is run with MPI and set the MPI run command
      if((SUNDIALS_ADD_TEST_MPI_NPROCS) AND ((MPIEXEC_EXECUTABLE)
                                             OR (SUNDIALS_TEST_MPIRUN_COMMAND)))
        if(SUNDIALS_TEST_MPIRUN_COMMAND)
          if(MPIEXEC_PREFLAGS)
            set(RUN_COMMAND
                "${SUNDIALS_TEST_MPIRUN_COMMAND} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${SUNDIALS_ADD_TEST_MPI_NPROCS}"
            )
          else()
            set(RUN_COMMAND
                "${SUNDIALS_TEST_MPIRUN_COMMAND} ${MPIEXEC_NUMPROC_FLAG} ${SUNDIALS_ADD_TEST_MPI_NPROCS}"
            )
          endif()
        elseif(MPIEXEC_EXECUTABLE)
          if(MPIEXEC_PREFLAGS)
            set(RUN_COMMAND
                "${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${SUNDIALS_ADD_TEST_MPI_NPROCS}"
            )
          else()
            set(RUN_COMMAND
                "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SUNDIALS_ADD_TEST_MPI_NPROCS}"
            )
          endif()
        endif()

        # remove trailing white space (empty MPIEXEC_PREFLAGS) as it can cause
        # erroneous test failures with some MPI implementations
        string(STRIP "${RUN_COMMAND}" RUN_COMMAND)
        list(APPEND TEST_ARGS "--runcommand=\"${RUN_COMMAND}\"")
      endif()

      # set the test input args
      if(_have_test_args)
        string(REPLACE ";" " " _user_args "${SUNDIALS_ADD_TEST_TEST_ARGS}")
        set(_run_args "${_user_args}")
        unset(_user_args)
      endif()
      if(_have_extra_test_args)
        string(REPLACE ";" " " _extra_args "${SUNDIALS_ADD_TEST_EXTRA_ARGS}")
        set(_run_args "${_run_args} ${_extra_args}")
        unset(_extra_args)
      endif()
      if(MPIEXEC_POSTFLAGS)
        set(_run_args "${MPIEXEC_POSTFLAGS} ${_run_args}")
      endif()
      if(_have_test_args
         OR _have_extra_test_args
         OR MPIEXEC_POSTFLAGS)
        string(STRIP "${_run_args}" _run_args)
        list(APPEND TEST_ARGS "--runargs=\"${_run_args}\"")
        unset(_run_args)
      endif()

      # create test case with the corresponding test runner command and
      # arguments all tests are added during development and only unlabeled
      # tests when released
      add_test(NAME ${NAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTRUNNER}
                                    ${TEST_ARGS})

    elseif(NOT SUNDIALS_ADD_TEST_EXAMPLE_TYPE)

      # if a test type was not set then it is a standard test that returns
      # pass/fail

      # convert string to list
      if(_have_test_args)
        string(REPLACE " " ";" TEST_ARGS "${SUNDIALS_ADD_TEST_TEST_ARGS}")
      endif()

      # check if this test is run with MPI and add the test run command
      if((SUNDIALS_ADD_TEST_MPI_NPROCS) AND ((MPIEXEC_EXECUTABLE)
                                             OR (SUNDIALS_TEST_MPIRUN_COMMAND)))
        if(MPIEXEC_PREFLAGS)
          string(REPLACE " " ";" PREFLAGS "${MPIEXEC_PREFLAGS}")
        endif()
        if(MPIEXEC_POSTFLAGS)
          string(REPLACE " " ";" POSTLAGS "${MPIEXEC_POSTFLAGS}")
        endif()
        if(SUNDIALS_TEST_MPIRUN_COMMAND)
          string(REPLACE " " ";" MPI_EXEC_ARGS
                         "${SUNDIALS_TEST_MPIRUN_COMMAND}")
          add_test(
            NAME ${NAME}
            COMMAND
              ${MPI_EXEC_ARGS} ${PREFLAGS} ${MPIEXEC_NUMPROC_FLAG}
              ${SUNDIALS_ADD_TEST_MPI_NPROCS} $<TARGET_FILE:${EXECUTABLE}>
              ${POSTFLAGS} ${TEST_ARGS})
        else()
          add_test(
            NAME ${NAME}
            COMMAND
              ${MPIEXEC_EXECUTABLE} ${PREFLAGS} ${MPIEXEC_NUMPROC_FLAG}
              ${SUNDIALS_ADD_TEST_MPI_NPROCS} ${POSTFLAGS}
              $<TARGET_FILE:${EXECUTABLE}> ${POSTFLAGS} ${TEST_ARGS})
        endif()
      else()
        add_test(NAME ${NAME} COMMAND $<TARGET_FILE:${EXECUTABLE}> ${TEST_ARGS})
      endif()

    endif()

    if(SUNDIALS_TEST_DEVTESTS OR NOT SUNDIALS_ADD_TEST_EXAMPLE_TYPE)
      # set any labels (must quote SUNDIALS_ADD_TEST_LABELS)
      if(SUNDIALS_ADD_TEST_LABELS)
        set_tests_properties(${NAME} PROPERTIES LABELS
                                                "${SUNDIALS_ADD_TEST_LABELS}")
      endif()
    endif()

  endif()

endmacro()
