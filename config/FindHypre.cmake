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
