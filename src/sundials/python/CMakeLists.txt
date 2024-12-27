

find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# Add the source file
set(sundials_SOURCES bind_sundials_core.cpp)

# Create the library
nanobind_add_module(pysundials_core ${sundials_SOURCES})

# Link against sundials libraries
target_link_libraries(pysundials_core PRIVATE sundials_core)
