
install(CODE "MESSAGE(\"\nInstall NVECTOR_KOKKOS\n\")")

# Add variable nveckokkos_SOURCES with the sources for the NVECKOKKOS lib
set(nveckokkos_SOURCES
  nvector_kokkos.cpp
)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECKOKKOS library
set(shared_SOURCES
  sundials_nvector.c
  sundials_math.c
)
add_prefix(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nvec_HEADERS with the exported NVECKOKKOS header files
set(nveckokkos_HEADERS nvector_kokkos.h)
add_prefix(${sundials_SOURCE_DIR}/include/nvector/ nveckokkos_HEADERS)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
add_definitions(-DBUILD_SUNDIALS_LIBRARY)

#checkpoint
if(BUILD_STATIC_LIBS)
  add_library(sundials_nveckokkos_static STATIC ${nveckokkos_SOURCES} ${shared_SOURCES})
  target_compile_features(sundials_nveckokkos_static PUBLIC cxx_std_11)
  target_link_libraries(sundials_nveckokkos_static PUBLIC Kokkos::kokkos)
  target_include_directories(sundials_nveckokkos_static PRIVATE . ${sundials_SOURCE_DIR}/src/sundials)
  set_target_properties(sundials_nveckokkos_static PROPERTIES
                        OUTPUT_NAME sundials_nveckokkos
                        CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_nveckokkos_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_STATIC_LIBS)


if(BUILD_SHARED_LIBS)
  add_library(sundials_nveckokkos_shared SHARED ${nveckokkos_SOURCES} ${shared_SOURCES})
  #TODO Do I need to set cxx_std_11 flag
  # target_compile_features(sundials_nveckokkos_shared PUBLIC cxx_std_11)
  target_link_libraries(sundials_nveckokkos_shared PUBLIC Kokkos::kokkos)
  target_include_directories(sundials_nveckokkos_shared PRIVATE . ${sundials_SOURCE_DIR}/src/sundials)
  set_target_properties(sundials_nveckokkos_shared PROPERTIES
                        OUTPUT_NAME sundials_nveckokkos
                        CLEAN_DIRECT_OUTPUT 1)
  install(TARGETS sundials_nveckokkos_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(BUILD_SHARED_LIBS)

# Install the NVECKOKKOS header files
install(FILES ${nveckokkos_HEADERS} DESTINATION include/nvector)

message(STATUS "Added NVECTOR_KOKKOS module")
