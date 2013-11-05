# ---------------------------------------------------------------
# $Revision: 1.0 $
# $Date: 2013-02-14 $
# ---------------------------------------------------------------
# Programmer:  Steven Smith @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2013, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# Find KLU library.
# 

IF(MSVC)
  SET(PRE lib)
ENDIF(MSVC)

if (KLU_LIBRARIES)
    #print_warning("FindKLU.cmake KLU_LIBRARIES" "${KLU_LIBRARIES}")
    get_filename_component(KLU_LIBRARY_DIR ${KLU_LIBRARIES} PATH)
    #print_warning("FindKLU.cmake KLU_LIBRARY_DIR" "${KLU_LIBRARY_DIR}")
    
else (KLU_LIBRARIES)
    # SGS TODO Assumption here that all of SparseSuite is in the same dir
    # SGS TODO Not sure why this is convoluted.
    set(KLU_LIBRARY_NAME klu)
    
    # find library path using potential names for static and/or shared libs
    set(temp_KLU_LIBRARY_DIR ${KLU_LIBRARY_DIR})
    unset(KLU_LIBRARY_DIR CACHE)  
    find_path(KLU_LIBRARY_DIR
        NAMES lib${KLU_LIBRARY_NAME}.so lib${KLU_LIBRARY_NAME}.a
        PATHS ${temp_KLU_LIBRARY_DIR}
        )
    
    FIND_LIBRARY( KLU_LIBRARIES ${PRE}klu${POST} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)

    # AMD
    set(AMD_LIBRARY_NAME amd)

    # find library path using potential names for static and/or shared libs
    set(temp_AMD_LIBRARY_DIR ${KLU_LIBRARY_DIR})
    unset(AMD_LIBRARY_DIR CACHE)  
    find_path(AMD_LIBRARY_DIR
        NAMES lib${AMD_LIBRARY_NAME}.so lib${AMD_LIBRARY_NAME}.a
        PATHS ${temp_AMD_LIBRARY_DIR}
        )
    
    FIND_LIBRARY( AMD_LIBRARIES ${PRE}amd${POST} ${AMD_LIBRARY_DIR} NO_DEFAULT_PATH)

    # COLAMD
    set(COLAMD_LIBRARY_NAME colamd)
    
    # find library path using potential names for static and/or shared libs
    set(temp_COLAMD_LIBRARY_DIR ${KLU_LIBRARY_DIR})
    unset(COLAMD_LIBRARY_DIR CACHE)  
    find_path(COLAMD_LIBRARY_DIR
        NAMES lib${COLAMD_LIBRARY_NAME}.so lib${COLAMD_LIBRARY_NAME}.a
        PATHS ${temp_COLAMD_LIBRARY_DIR}
        )
    
    FIND_LIBRARY( COLAMD_LIBRARIES ${PRE}colamd${POST} ${COLAMD_LIBRARY_DIR} NO_DEFAULT_PATH)


    # BTF
    set(BTF_LIBRARY_NAME btf)
    
    # find library path using potential names for static and/or shared libs
    set(temp_BTF_LIBRARY_DIR ${KLU_LIBRARY_DIR})
    unset(BTF_LIBRARY_DIR CACHE)  
    find_path(BTF_LIBRARY_DIR
        NAMES lib${BTF_LIBRARY_NAME}.so lib${BTF_LIBRARY_NAME}.a
        PATHS ${temp_BTF_LIBRARY_DIR}
        )
    
    FIND_LIBRARY( BTF_LIBRARIES ${PRE}btf${POST} ${BTF_LIBRARY_DIR} NO_DEFAULT_PATH)

    LIST(APPEND KLU_LIBRARIES ${AMD_LIBRARIES})
    LIST(APPEND KLU_LIBRARIES ${COLAMD_LIBRARIES})
    LIST(APPEND KLU_LIBRARIES ${BTF_LIBRARIES})
    
    mark_as_advanced(KLU_LIBRARIES)
    
endif (KLU_LIBRARIES)

# add to KLU_LIBRARIES
set(KLU_LIBRARIES ${KLU_LIBRARIES})

