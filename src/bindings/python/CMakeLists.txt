

find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# Add the source file
set(sundials_SOURCES pysundials.cpp sundials_core.cpp arkode.cpp nvector/serial.cpp)

# Create the library
nanobind_add_module(pysundials ${sundials_SOURCES})

target_include_directories(pysundials PRIVATE ${SUNDIALS_SOURCE_DIR}/src ${SUNDIALS_SOURCE_DIR}/src/arkode)

# Link against sundials libraries
target_link_libraries(pysundials PRIVATE sundials_arkode sundials_nvecserial sundials_core)
