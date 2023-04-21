# Take in an executable, determine how to construct the arguments

# ---------------------------------------------------------------
# CMake macro for adding benchmarks to `make benchmark`.
# ---------------------------------------------------------------

macro(sundials_add_benchmark NAME EXECUTABLE)
  # Define single value parameters the macro takes in to set up the test runner
  #
  # NUM_NODES         = number of nodes (GPU count or CPU count) to run on or number of resource sets
  # BENCHMARK_ARGS    = arguments to pass to the executable
  # IDENTIFIER        = suffix to append to end of benchmark name
  set(oneValueArgs NUM_NODES BENCHMARK_ARGS IDENTIFIER)

  # TEST_RUNNER_ARGS  = command line arguments to pass to the test executable
  set(multiValueArgs TEST_RUNNER_ARGS )

  # ENABLE_GPU        = indicate this benchmark should be run with GPUs
  set(options ENABLE_GPU)

  cmake_parse_arguments(sundials_add_benchmark
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # command line arguments for the test runner script
  set(TEST_RUNNER_ARGS
    "--verbose"
    "--executablename=$<TARGET_FILE:${EXECUTABLE}>"
    "--outputdir=${SUNDIALS_BENCHMARK_OUTPUT_DIR}"
    "--nodiff"
    )

  # set the target name
  if(sundials_add_benchmark_IDENTIFIER)
    set(TARGET_NAME ${NAME}_${sundials_add_benchmark_IDENTIFIER})
  else()
    if(sundials_add_benchmark_BENCHMARK_ARGS)
      string(REPLACE " " "_" TEST_SUFFIX "${sundials_add_benchmark_BENCHMARK_ARGS}")
      set(TARGET_NAME ${NAME}_${TEST_SUFFIX})
    else()
      set(TARGET_NAME ${NAME}_run)
    endif()
  endif()

  # incorporate scheduler arguments into test_runner
  if(SUNDIALS_SCHEDULER_COMMAND STREQUAL "flux run")
    set(SCHEDULER_STRING " -nnodes=${sundials_add_benchmark_NUM_NODES}")
  elseif(SUNDIALS_SCHEDULER_COMMAND STREQUAL "jsrun" AND ${sundials_add_benchmark_ENABLE_GPU})
    set(SCHEDULER_STRING " --smpiargs=\\\"-gpu\\\" -n${sundials_add_benchmark_NUM_NODES} -a1 -c1 -g1")
  elseif(SUNDIALS_SCHEDULER_COMMAND STREQUAL "jsrun")
    set(SCHEDULER_STRING " -n${sundials_add_benchmark_NUM_NODES} -a1 -c1")
  elseif(SUNDIALS_SCHEDULER_COMMAND STREQUAL "srun")
    set(SCHEDULER_STRING " -n${sundials_add_benchmark_NUM_NODES} --cpus-per-task=1 --ntasks-per-node=1")
  endif()
  string(REPLACE " " ";" SCHEDULER_ARGS "${SCHEDULER_STRING}")

  # taken from SundialsAddTest
  string(STRIP "${RUN_COMMAND}" RUN_COMMAND)
  set(RUN_COMMAND ${SUNDIALS_SCHEDULER_COMMAND} ${SCHEDULER_ARGS})
  list(APPEND TEST_RUNNER_ARGS "--runcommand=\"${RUN_COMMAND}\"")

  # enable test runner to set up caliper output paths and configs
  if(SUNDIALS_TEST_PROFILE)
    list(APPEND TEST_RUNNER_ARGS "--profile")
  endif()

  list(APPEND TEST_RUNNER_ARGS "--runargs=${sundials_add_benchmark_BENCHMARK_ARGS}" "--testname=${TARGET_NAME}")
  add_custom_target(${TARGET_NAME}
    COMMENT "Running ${TARGET_NAME}"
    COMMAND ${PYTHON_EXECUTABLE} ${TESTRUNNER} ${TEST_RUNNER_ARGS}
    )
  add_dependencies(benchmark ${TARGET_NAME})

endmacro()
