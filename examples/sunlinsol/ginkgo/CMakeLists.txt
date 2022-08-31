
# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using SUNDIALS Ginkgo dense matrix
set(sunmatrix_ginkgo_examples
  "test_sunlinsol_ginkgo.cpp\;bicg 100 10 100 0\;"
  "test_sunlinsol_ginkgo.cpp\;bicgstab 100 10 100 0\;"
  "test_sunlinsol_ginkgo.cpp\;cgs 100 10 100 0\;"
  "test_sunlinsol_ginkgo.cpp\;gmres 100 10 100 0\;"
  "test_sunlinsol_ginkgo.cpp\;idr 100 10 100 0\;"
  )

# Add source directory to include directories
include_directories(.. ${PROJECT_SOURCE_DIR}/examples/utilities)

# Add the build and install targets for each example
foreach(example_tuple ${sunmatrix_ginkgo_examples})
  foreach(matrix_type CSR DENSE)
    foreach(backend ${SUNDIALS_GINKGO_BACKENDS})

      # parse the example tuple
      list(GET example_tuple 0 example)
      list(GET example_tuple 1 example_args)
      list(GET example_tuple 2 example_type)

      if(backend MATCHES "CUDA")
        set_source_files_properties(${example} PROPERTIES LANGUAGE CUDA)
        set(vector nveccuda)
      elseif(backend MATCHES "HIP")
        set_source_files_properties(${example} PROPERTIES LANGUAGE CXX)
        set(vector nvechip)
      elseif(backend MATCHES "OMP")
        set(vector nvecserial)
      elseif(backend MATCHES "REF")
        set(vector nvecserial)
      endif()

      # extract the file name without extension
      get_filename_component(example_target ${example} NAME_WE)
      set(example_target "${example_target}.${matrix_type}.${backend}")

      if (NOT TARGET ${example_target})
        # example source files
        add_executable(${example_target} ${example})

        # folder for IDEs
        set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

        # which backend to use
        target_compile_definitions(${example_target} PRIVATE USE_${backend} USE_${matrix_type})

        # libraries to link against
        target_link_libraries(${example_target}
          PRIVATE
            test_sunlinsol_obj
            sundials_${vector}
            Ginkgo::ginkgo
            ${EXTRA_LINK_LIBS}
        )
      endif()

      # check if example args are provided and set the test name
      if("${example_args}" STREQUAL "")
        set(test_name ${example_target})
      else()
        string(REGEX REPLACE " " "_" test_name ${example_target}_${example_args})
      endif()

      # add example to regression tests
      sundials_add_test(${test_name} ${example_target}
        TEST_ARGS ${example_args}
        EXAMPLE_TYPE ${example_type}
        NODIFF)

    endforeach(backend)
  endforeach(matrix_type)
endforeach(example_tuple)

# if(EXAMPLES_INSTALL)

#   sundials_install_examples(sunmatrix_ginkgo sunmatrix_ginkgo_examples
#     EXAMPLES_DEPENDENCIES
#       ${SUNDIALS_SOURCE_DIR}/examples/sunmatrix/test_sunmatrix.c
#       ${SUNDIALS_SOURCE_DIR}/examples/sunmatrix/test_sunmatrix.h
#     CMAKE_TEMPLATE
#       cmakelists_${cuda_or_hip}_ex.in
#     SUNDIALS_TARGETS
#       ${vector}
#       sundials_sunmatrixginkgo
#     EXTRA_INCLUDES
#       ${MAGMA_INCLUDE_DIR}
#     DESTINATION
#       sunmatrix/magmadense
#     )

# endif()
