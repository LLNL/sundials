# ---------------------------------------------------------------
# $Revision: 1.0 $
# $Date: 2013-02-14 $
# ---------------------------------------------------------------
# Programmer:  Eddy Banks @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2013, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# SuperLU tests for SUNDIALS CMake-based configuration.
# 

# make sure valid thread type - if not, then warn and return
STRING(TOUPPER "${SUPERLU_THREAD_TYPE}" SUPERLU_THREAD_TYPE_UPPER)
If(SUPERLU_THREAD_TYPE AND NOT SUPERLU_THREAD_TYPE_UPPER STREQUAL "OPENMP" AND NOT SUPERLU_THREAD_TYPE_UPPER STREQUAL "PTHREAD")
    PRINT_WARNING("Unknown thread type: ${SUPERLU_THREAD_TYPE}" "Please enter Pthread or OpenMP")
ENDIF(SUPERLU_THREAD_TYPE AND NOT SUPERLU_THREAD_TYPE_UPPER STREQUAL "OPENMP" AND NOT SUPERLU_THREAD_TYPE_UPPER STREQUAL "PTHREAD")

# find the SuperLU include directory path
IF(SUPERLU_THREAD_TYPE)
    # if have user input for thread type - set postfix of library name
    set(POST ${SUPERLU_THREAD_TYPE_UPPER})
    set(temp_SUPERLU_INCLUDE_DIR ${SUPERLU_INCLUDE_DIR})
    unset(SUPERLU_INCLUDE_DIR CACHE)  
    find_path(SUPERLU_INCLUDE_DIR pdsp_defs.h ${temp_SUPERLU_INCLUDE_DIR})
ELSE(SUPERLU_THREAD_TYPE)
    FIND_PATH(SUPERLU_INCLUDE_DIR slu_defs.h ${temp_SUPERLU_INCLUDE_DIR})
ENDIF(SUPERLU_THREAD_TYPE)

IF(MSVC)
  SET(PRE lib)
ENDIF(MSVC)

if (SUPERLU_LIBRARIES)
    #print_warning("FindSUPERLU.cmake SUPERLU_LIBRARIES" "${SUPERLU_LIBRARIES}")
    get_filename_component(SUPERLU_LIBRARY_DIR ${SUPERLU_LIBRARIES} PATH)
    #print_warning("FindSUPERLU.cmake SUPERLU_LIBRARY_DIR" "${SUPERLU_LIBRARY_DIR}")
    
else (SUPERLU_LIBRARIES)

    set(SUPERLU_LIBRARY_NAME superlu_mt_${POST})
    
    # find library path using potential names for static and/or shared libs
    set(temp_SUPERLU_LIBRARY_DIR ${SUPERLU_LIBRARY_DIR})
    unset(SUPERLU_LIBRARY_DIR CACHE)  
    find_path(SUPERLU_LIBRARY_DIR
        NAMES lib${SUPERLU_LIBRARY_NAME}.so lib${SUPERLU_LIBRARY_NAME}.a
        PATHS ${temp_SUPERLU_LIBRARY_DIR}
        )
    
    #unset(SUPERLU_LIBRARIES CACHE)
    FIND_LIBRARY( SUPERLU_LIBRARIES ${PRE}superlu_mt_${POST} ${SUPERLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(SUPERLU_LIBRARIES)
    
    #print_warning("FindSUPERLU.cmake SUPERLU_LIBRARIES" "${SUPERLU_LIBRARIES}")
endif (SUPERLU_LIBRARIES)
#print_warning("FindSUPERLU.cmake SUPERLU_LIBRARIES" "${SUPERLU_LIBRARIES}")

# add threading library (pthread or gomp)
unset(SUPERLU_THREAD_LIBRARY CACHE)
If(SUPERLU_THREAD_TYPE_UPPER STREQUAL "PTHREAD")
  # add pthread to libraries
  find_library(SUPERLU_THREAD_LIBRARY
      NAMES pthread
      PATHS /usr/lib /usr/local/lib
      "$ENV{ProgramFiles}/SUPERLU/Lib"
      )
  # add to SUPERLU_LIBRARIES
ELSE(SUPERLU_THREAD_TYPE_UPPER STREQUAL "PTHREAD")
  # add pthread to libraries
  find_library(SUPERLU_THREAD_LIBRARY
      NAMES gomp
      PATHS /usr/lib /usr/local/lib
      "$ENV{ProgramFiles}/SUPERLU/Lib"
      )
ENDIF(SUPERLU_THREAD_TYPE_UPPER STREQUAL "PTHREAD")
mark_as_advanced(SUPERLU_THREAD_LIBRARY)

# add to SUPERLU_LIBRARIES
set(SUPERLU_LIBRARIES ${SUPERLU_LIBRARIES} ${SUPERLU_THREAD_LIBRARY})


# If LAPACK/BLAS not enabled - find BLAS with SuperLU
if(NOT LAPACK_ENABLE)
    set(SUPERLU_BLAS_LIBRARY_NAME blas_${POST})
    
    #unset(SUPERLU_BLAS_LIBRARIES CACHE)
    FIND_LIBRARY( SUPERLU_BLAS_LIBRARIES ${PRE}blas_${POST} ${SUPERLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    set(SUPERLU_LIBRARIES ${SUPERLU_LIBRARIES} ${SUPERLU_BLAS_LIBRARIES})
    mark_as_advanced(SUPERLU_BLAS_LIBRARIES)
endif(NOT LAPACK_ENABLE)