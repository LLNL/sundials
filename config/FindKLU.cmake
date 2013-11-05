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
    set(KLU_LIBRARY_NAME klu)
    
    # find library path using potential names for static and/or shared libs
    set(temp_KLU_LIBRARY_DIR ${KLU_LIBRARY_DIR})
    unset(KLU_LIBRARY_DIR CACHE)  
    find_path(KLU_LIBRARY_DIR
        NAMES lib${KLU_LIBRARY_NAME}.so lib${KLU_LIBRARY_NAME}.a
        PATHS ${temp_KLU_LIBRARY_DIR}
        )
    
    FIND_LIBRARY( KLU_LIBRARIES ${PRE}klu${POST} ${KLU_LIBRARY_DIR} NO_DEFAULT_PATH)
    mark_as_advanced(KLU_LIBRARIES)
    
endif (KLU_LIBRARIES)

# add to KLU_LIBRARIES
set(KLU_LIBRARIES ${KLU_LIBRARIES})

