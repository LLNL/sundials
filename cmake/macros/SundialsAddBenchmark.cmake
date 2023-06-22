# Take in an executable, determine how to construct the arguments

# ---------------------------------------------------------------
# CMake macro for adding benchmarks to `make benchmark`.
# ---------------------------------------------------------------

macro(sundials_add_benchmark NAME EXECUTABLE BASE_BENCHMARK_NAME)

  # Define single value parameters the macro takes in to set up the test runner
  #
  # NUM_CORES         = number of cores (GPU count or CPU count) to run on/number of resource sets
  # BENCHMARK_ARGS    = arguments to pass to the executable
  # IDENTIFIER        = suffix to append to end of benchmark name
  set(oneValueArgs NUM_CORES BENCHMARK_ARGS IDENTIFIER)

  # TEST_RUNNER_ARGS  = command line arguments to pass to the test executable
  set(multiValueArgs TEST_RUNNER_ARGS )

  # ENABLE_GPU        = indicate this benchmark should be run with GPUs
  set(options ENABLE_GPU)

  cmake_parse_arguments(sundials_add_benchmark
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

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

  # Create default benchmark output directory if custom directory not defined
  if(SUNDIALS_TEST_OUTPUT_DIR)
    set(SUNDIALS_BENCHMARK_OUTPUT_DIR ${SUNDIALS_TEST_OUTPUT_DIR}/Benchmarking/${BASE_BENCHMARK_NAME})
  else()
    set(SUNDIALS_BENCHMARK_OUTPUT_DIR ${PROJECT_BINARY_DIR}/Benchmarking/${BASE_BENCHMARK_NAME})
  endif()

  # make the output directory if it doesn't exist
  if(NOT EXISTS ${SUNDIALS_BENCHMARK_OUTPUT_DIR}/${TARGET_NAME})
    file(MAKE_DIRECTORY ${SUNDIALS_BENCHMARK_OUTPUT_DIR}/${TARGET_NAME})
  endif()

  # command line arguments for the test runner script
  set(TEST_RUNNER_ARGS
    "--verbose"
    "--executablename=$<TARGET_FILE:${EXECUTABLE}>"
    "--outputdir=${SUNDIALS_BENCHMARK_OUTPUT_DIR}/${TARGET_NAME}"
    "--nodiff"
    )
  
  # incorporate scheduler arguments into test_runner
  if(SUNDIALS_SCHEDULER_COMMAND STREQUAL "flux run")
    set(SCHEDULER_STRING " -n${sundials_add_benchmark_NUM_CORES}")
  elseif(SUNDIALS_SCHEDULER_COMMAND STREQUAL "jsrun" AND ${sundials_add_benchmark_ENABLE_GPU})
    set(SCHEDULER_STRING " --smpiargs=\\\"-gpu\\\" -n${sundials_add_benchmark_NUM_CORES} -a1 -c1 -g1")
  elseif(SUNDIALS_SCHEDULER_COMMAND STREQUAL "jsrun")
    set(SCHEDULER_STRING " -n${sundials_add_benchmark_NUM_CORES} -a1 -c1")
  elseif(SUNDIALS_SCHEDULER_COMMAND STREQUAL "srun")
    set(SCHEDULER_STRING " -n${sundials_add_benchmark_NUM_CORES} --cpus-per-task=1 --ntasks-per-node=1")
  endif()
  string(REPLACE " " ";" SCHEDULER_ARGS "${SCHEDULER_STRING}")
  string(REPLACE " " ";" SCHEDULER_COMMAND_ARGS "${SUNDIALS_SCHEDULER_COMMAND}")

  # taken from SundialsAddTest
  string(STRIP "${RUN_COMMAND}" RUN_COMMAND)
  set(RUN_COMMAND ${SCHEDULER_COMMAND_ARGS} ${SCHEDULER_ARGS})
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
