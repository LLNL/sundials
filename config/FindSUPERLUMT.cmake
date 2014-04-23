# ---------------------------------------------------------------
# $Revision$
# $Date$
# ---------------------------------------------------------------
# Programmer:  Eddy Banks @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2013, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# SUPERLUMT tests for SUNDIALS CMake-based configuration.
# 

# make sure valid thread type - if not, then warn and return
STRING(TOUPPER "${SUPERLUMT_THREAD_TYPE}" SUPERLUMT_THREAD_TYPE_UPPER)
If(SUPERLUMT_THREAD_TYPE AND NOT SUPERLUMT_THREAD_TYPE_UPPER STREQUAL "OPENMP" AND NOT SUPERLUMT_THREAD_TYPE_UPPER STREQUAL "PTHREAD")
    PRINT_WARNING("Unknown thread type: ${SUPERLUMT_THREAD_TYPE}" "Please enter Pthread or OpenMP")
ENDIF(SUPERLUMT_THREAD_TYPE AND NOT SUPERLUMT_THREAD_TYPE_UPPER STREQUAL "OPENMP" AND NOT SUPERLUMT_THREAD_TYPE_UPPER STREQUAL "PTHREAD")

# find the SUPERLUMT include directory path
IF(SUPERLUMT_THREAD_TYPE)
    # if have user input for thread type - set postfix of library name
    set(POST ${SUPERLUMT_THREAD_TYPE_UPPER})
    set(temp_SUPERLUMT_INCLUDE_DIR ${SUPERLUMT_INCLUDE_DIR})
    unset(SUPERLUMT_INCLUDE_DIR CACHE)  
    find_path(SUPERLUMT_INCLUDE_DIR pdsp_defs.h ${temp_SUPERLUMT_INCLUDE_DIR})
ELSE(SUPERLUMT_THREAD_TYPE)
    FIND_PATH(SUPERLUMT_INCLUDE_DIR slu_defs.h ${temp_SUPERLUMT_INCLUDE_DIR})
ENDIF(SUPERLUMT_THREAD_TYPE)

IF(MSVC)
  SET(PRE lib)
ENDIF(MSVC)

if (SUPERLUMT_LIBRARIES)
    #print_warning("FindSUPERLUMT.cmake SUPERLUMT_LIBRARIES" "${SUPERLUMT_LIBRARIES}")
    get_filename_component(SUPERLUMT_LIBRARY_DIR ${SUPERLUMT_LIBRARIES} PATH)
    #print_warning("FindSUPERLUMT.cmake SUPERLUMT_LIBRARY_DIR" "${SUPERLUMT_LIBRARY_DIR}")
    
else (SUPERLUMT_LIBRARIES)

    set(SUPERLUMT_LIBRARY_NAME superlu_mt_${POST})
    
    # find library path using potential names for static and/or shared libs
    set(temp_SUPERLUMT_LIBRARY_DIR ${SUPERLUMT_LIBRARY_DIR})
    unset(SUPERLUMT_LIBRARY_DIR CACHE)  
    find_path(SUPERLUMT_LIBRARY_DIR
        NAMES lib${SUPERLUMT_LIBRARY_NAME}.so lib${SUPERLUMT_LIBRARY_NAME}.a
        PATHS ${temp_SUPERLUMT_LIBRARY_DIR}
        )
    
    #unset(SUPERLUMT_LIBRARIES CACHE)
    FIND_LIBRARY( SUPERLUMT_LIBRARIES ${PRE}superlu_mt_${POST} ${SUPERLUMT_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(SUPERLUMT_LIBRARIES)
    
    #print_warning("FindSUPERLUMT.cmake SUPERLUMT_LIBRARIES" "${SUPERLUMT_LIBRARIES}")
endif (SUPERLUMT_LIBRARIES)
#print_warning("FindSUPERLUMT.cmake SUPERLUMT_LIBRARIES" "${SUPERLUMT_LIBRARIES}")

# add threading library (pthread or gomp)
unset(SUPERLUMT_THREAD_LIBRARY CACHE)
If(SUPERLUMT_THREAD_TYPE_UPPER STREQUAL "PTHREAD")
  # add pthread to libraries
  find_library(SUPERLUMT_THREAD_LIBRARY
      NAMES pthread
      PATHS /usr/lib /usr/local/lib
      "$ENV{ProgramFiles}/SUPERLUMT/Lib"
      )
  # add to SUPERLUMT_LIBRARIES
ELSE(SUPERLUMT_THREAD_TYPE_UPPER STREQUAL "PTHREAD")
  # add pthread to libraries
  find_library(SUPERLUMT_THREAD_LIBRARY
      NAMES gomp
      PATHS /usr/lib /usr/local/lib
      "$ENV{ProgramFiles}/SUPERLUMT/Lib"
      )
ENDIF(SUPERLUMT_THREAD_TYPE_UPPER STREQUAL "PTHREAD")
mark_as_advanced(SUPERLUMT_THREAD_LIBRARY)

# add to SUPERLUMT_LIBRARIES
set(SUPERLUMT_LIBRARIES ${SUPERLUMT_LIBRARIES} ${SUPERLUMT_THREAD_LIBRARY})


# If LAPACK/BLAS not enabled - find BLAS with SUPERLUMT
if(NOT LAPACK_ENABLE)
    set(SUPERLUMT_BLAS_LIBRARY_NAME blas_${POST})
    
    #unset(SUPERLUMT_BLAS_LIBRARIES CACHE)
    FIND_LIBRARY( SUPERLUMT_BLAS_LIBRARIES ${PRE}blas_${POST} ${SUPERLUMT_LIBRARY_DIR} NO_DEFAULT_PATH)
    set(SUPERLUMT_LIBRARIES ${SUPERLUMT_LIBRARIES} ${SUPERLUMT_BLAS_LIBRARIES})
    mark_as_advanced(SUPERLUMT_BLAS_LIBRARIES)
endif(NOT LAPACK_ENABLE)