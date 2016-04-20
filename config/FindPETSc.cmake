# ---------------------------------------------------------------
# $Revision:  $
# $Date:  $
# ---------------------------------------------------------------
# Programmer:  Steven Smith @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2013, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# Find PETSC library.
# 

set(PRE "lib")
IF(WIN32)
  set(POST ".lib" ".dll")
else(WIN32)
  set(POST ".so")
endif(WIN32)

if (PETSC_LIBRARY)
    get_filename_component(PETSC_LIBRARY_DIR ${PETSC_LIBRARY} PATH)
else (PETSC_LIBRARY)
    set(PETSC_LIBRARY_NAME petsc)
    
    # find library path using potential names for static and/or shared libs
    set(temp_PETSC_LIBRARY_DIR ${PETSC_LIBRARY_DIR})
    unset(PETSC_LIBRARY_DIR CACHE)  
    find_path(PETSC_LIBRARY_DIR
        NAMES ${PRE}${PETSC_LIBRARY_NAME}${POST}
        PATHS ${temp_PETSC_LIBRARY_DIR}
        )

    mark_as_advanced(PETSC_LIBRARY)

    FIND_LIBRARY( PETSC_LIBRARY ${PRE}petsc${POST} ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH)
endif (PETSC_LIBRARY)


set(PETSC_LIBRARY ${PETSC_LIBRARY})
