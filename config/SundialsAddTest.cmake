# ---------------------------------------------------------------
# Author:  Steven Smith @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2013, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
#
# SUNDIALS_ADD_TEST(<test name> <executable> )
# 
# Add Sundials regression test.

# Executable is run and output is compared with output in the
# test/answers directory.  If output differs significantly then test
# fails.  Default signicance is 4 decimal points for floating values
# and 10% for integer values.

IF(EXAMPLES_ENABLE)

  find_package(PythonInterp)
  IF(${PYTHON_VERSION_MAJOR} LESS 3)
      IF(${PYTHON_VERSION_MINOR} LESS 7)
	message( WARNING "***************************************************************************\nWARNING\nPython version must be 2.7.x or greater in order to run regression tests.\nExamples will build but 'make test' will fail.\n***************************************************************************")
      ENDIF()
  ENDIF()

  FIND_PROGRAM(TESTRUNNER testRunner PATHS test)
ENDIF(EXAMPLES_ENABLE)

macro(SUNDIALS_ADD_TEST NAME EXECUTABLE)

  set(options "")
  set(oneValueArgs "MPI_NPROCS" "ANSWER_FILE" "FLOAT_PRECISION" "INTEGER_PERCENTAGE")
  set(multiValueArgs "TEST_ARGS")

  CMAKE_PARSE_ARGUMENTS(SUNDIALS_ADD_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # SGS add check to make sure parallel is integer
  # SGS add check for float and integer precision

  set(TEST_ARGS  "-v" "--testname=${NAME}" 
    "--executablename=$<TARGET_FILE:${EXECUTABLE}>"
    "--answersdir=${CMAKE_SOURCE_DIR}/test/answers"
    "--outputdir=${CMAKE_BINARY_DIR}/Testing/output"
    )

  IF("${SUNDIALS_ADD_TEST_MPI_NPROCS}" STREQUAL "")
  ELSE()

    IF(MPI_ENABLE)
      IF(MPI_RUN_COMMAND MATCHES "srun")
	set(RUN_COMMAND "srun -N1 -n${SUNDIALS_ADD_TEST_MPI_NPROCS} -ppdebug")
      ELSE(MPI_RUN_COMMAND MATCHES "srun")
	set(RUN_COMMAND "mpirun -n ${SUNDIALS_ADD_TEST_MPI_NPROCS}")
      ENDIF(MPI_RUN_COMMAND MATCHES "srun")
      
      LIST(APPEND TEST_ARGS "--runcommand=\"${RUN_COMMAND}\"")

    ENDIF(MPI_ENABLE)

  ENDIF()
  
  IF("${SUNDIALS_ADD_TEST_TEST_ARGS}" STREQUAL "")
  ELSE()
    string (REPLACE ";" " " USER_ARGS "${SUNDIALS_ADD_TEST_TEST_ARGS}")
    LIST(APPEND TEST_ARGS "--runargs=\"${USER_ARGS}\"")
  ENDIF()

  IF("${SUNDIALS_ADD_TEST_ANSWER_FILE}" STREQUAL "")
  ELSE()
    LIST(APPEND TEST_ARGS "--answerfile=${SUNDIALS_ADD_TEST_ANSWER_FILE}")
  ENDIF()

  IF("${SUNDIALS_ADD_TEST_FLOAT_PRECISION}" STREQUAL "")
  ELSE()
    LIST(APPEND TEST_ARGS "--floatprecision=${SUNDIALS_ADD_TEST_FLOAT_PRECISION}")
  ENDIF()

  IF("${SUNDIALS_ADD_TEST_INTEGER_PERCENTAGE}" STREQUAL "")
  ELSE()
    LIST(APPEND TEST_ARGS "--integerpercentage=${SUNDIALS_ADD_TEST_INTEGER_PERCENTAGE}")
  ENDIF()

  ADD_TEST(NAME ${NAME}
    COMMAND ${PYTHON_EXECUTABLE} ${TESTRUNNER} ${TEST_ARGS})

endmacro()
