# ---------------------------------------------------------------
# $Revision: 4713 $
# $Date: 2016-03-28 07:20:43 -0700 (Mon, 28 Mar 2016) $
# ---------------------------------------------------------------
# Programmer:  Slaven Peles @ LLNL, Jean Sexton @ SMU
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
# - Find hypre

#  HYPRE_INCLUDE_DIR = cached location of HYPRE.h
#  HYPRE_LIBRARY    = cached list of HYPRE library to link in
if (HYPRE_LIBRARY)
    get_filename_component(HYPRE_LIBRARY_DIR ${HYPRE_LIBRARY} PATH)
else (HYPRE_LIBRARY)

FIND_LIBRARY(HYPRE_LIBRARY
  NAMES HYPRE 
  PATHS /usr/lib /usr/local/lib /usr/local/hypre/lib
  "${CMAKE_SOURCE_DIR}/../hypre/src/hypre/lib"
  ${HYPRE_LIBRARY_DIR}
  )
get_filename_component(HYPRE_LIBRARY_DIR ${HYPRE_LIBRARY} PATH)

endif (HYPRE_LIBRARY)

FIND_PATH(HYPRE_INCLUDE_DIR HYPRE.h
  PATHS /usr/local/include 
  /usr/include 
  /usr/include/hypre
  /usr/local/hypre/include
  "${CMAKE_SOURCE_DIR}/../hypre/src/hypre/include"
  ${HYPRE_LIBRARY_DIR}/../include
  )
#MESSAGE(STATUS "HYPRE_LIBRARY_DIR:" ${HYPRE_LIBRARY_DIR})
#SET(HYPRE_INCLUDE_DIR "${HYPRE_LIBRARY_DIR}/../include")
#MESSAGE(STATUS "HYPRE_INCLUDE_DIR:" ${HYPRE_INCLUDE_DIR})
