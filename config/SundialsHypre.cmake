# ---------------------------------------------------------------
# Programmer:  Jean Sexton @ LLNL
# ---------------------------------------------------------------
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# Hypre tests for SUNDIALS CMake-based configuration.
# 

set(HYPRE_FOUND FALSE)

include(FindHypre)

if(HYPRE_LIBRARY)
  message(STATUS "Looking for HYPRE LIBRARY...")
  # Create the HYPRETest directory
  set(HYPRETest_DIR ${PROJECT_BINARY_DIR}/HYPRETest)
  file(MAKE_DIRECTORY ${HYPRETest_DIR})
  # Create a CMakeLists.txt file 
  file(WRITE ${HYPRETest_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 2.2)\n"
    "PROJECT(ltest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_C_COMPILER ${MPI_MPICC})\n"
    "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "INCLUDE_DIRECTORIES(${HYPRE_INCLUDE_DIR})\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_LINK_LIBRARIES(ltest ${HYPRE_LIBRARY})\n")    
# Create a C source file which calls a KLU function
# SGS TODO what is a simple KLU method to invoke?
  file(WRITE ${HYPRETest_DIR}/ltest.c
    "\#include \"HYPRE_parcsr_ls.h\"\n"
    "int main(){\n"
    "HYPRE_ParVector par_b;\n"
    "HYPRE_IJVector b;\n"
    "HYPRE_IJVectorAssemble(b);\n"
    "HYPRE_IJVectorGetObject(b, (void **) &par_b);\n"
    "return(0);\n"
    "}\n")
  # Attempt to link the "ltest" executable
  try_compile(LTEST_OK ${HYPRETest_DIR} ${HYPRETest_DIR} ltest OUTPUT_VARIABLE MY_OUTPUT)
      
  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${HYPRETest_DIR}/CMakeFiles)
  # Process test result
#PRINT_WARNING("LTEST_OK" "${LTEST_OK}")
  if(LTEST_OK)
#PRINT_WARNING("x SundialsKLU.cmake KLU_LIBRARY" "${KLU_LIBRARY}")
    message(STATUS "Checking if HYPRE works... OK")
    set(HYPRE_FOUND TRUE)
    #print_warning("KLU_FOUND" "${KLU_FOUND}")
  else(LTEST_OK)
    message(STATUS "Checking if HYPRE works... FAILED")
  endif(LTEST_OK)
else(HYPRE_LIBRARY)
  message(STATUS "Looking for HYPRE LIBRARY... FAILED")
endif(HYPRE_LIBRARY)


#    "HYPRE_IJVector b;\n"
#    "HYPRE_IJVectorAssemble(b);\n"
#    "HYPRE_IJVectorGetObject(b, (void **) &par_b);\n"
