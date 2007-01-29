# - Find LAPACK library
# This module finds an installed fortran library that implements the LAPACK
# linear-algebra interface (see http://www.netlib.org/lapack/).
#
# The approach follows that taken for the autoconf macro file, acx_lapack.m4
# (distributed at http://ac-archive.sourceforge.net/ac-archive/acx_lapack.html).
#
# This module sets the following variables:
#  LAPACK_FOUND - set to true if a library implementing the LAPACK interface
#    is found
#  LAPACK_LINKER_FLAGS - uncached list of required linker flags (excluding -l
#    and -L).
#  LAPACK_LIBRARIES - uncached list of libraries (using full path name) to 
#    link against to use LAPACK
#

include(CheckFortranFunctionExists)
set(LAPACK_FOUND FALSE)

if(LAPACK_FIND_QUIETLY OR NOT LAPACK_FIND_REQUIRED)
  find_package(BLAS)
else(LAPACK_FIND_QUIETLY OR NOT LAPACK_FIND_REQUIRED)
  find_package(BLAS REQUIRED)
endif(LAPACK_FIND_QUIETLY OR NOT LAPACK_FIND_REQUIRED)

if(BLAS_FOUND)
  set(LAPACK_LINKER_FLAGS ${BLAS_LINKER_FLAGS})
  # LAPACK linked to by default?  (is sometimes included in BLAS lib)
  set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES})
  check_fortran_function_exists(cheev LAPACK_BLAS_WORKS)
  mark_as_advanced(LAPACK_BLAS_WORKS)
  if(LAPACK_BLAS_WORKS)
    set(LAPACK_FOUND TRUE)
    set(LAPACK_LIBRARIES ${BLAS_LIBRARIES})
  endif(LAPACK_BLAS_WORKS)
  # Generic LAPACK library?
  if(NOT LAPACK_FOUND)
    find_library(LAPACK_LAPACK_LIBRARY
    NAMES lapack
    PATHS /usr/local/lib /usr/lib
    )
    mark_as_advanced(LAPACK_LAPACK_LIBRARY)
    if(LAPACK_LAPACK_LIBRARY)
      set(LAPACK_LIBRARIES ${LAPACK_LAPACK_LIBRARY} ${BLAS_LIBRARIES})
      # Test this combination of libraries.
      set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
      check_fortran_function_exists(cheev LAPACK_LAPACK_WORKS)
      mark_as_advanced(LAPACK_LAPACK_WORKS)
      set(CMAKE_REQUIRED_LIBRARIES)
      if(LAPACK_LAPACK_WORKS)
        set(LAPACK_FOUND TRUE)
      else(LAPACK_LAPACK_WORKS)
        set(LAPACK_LIBRARIES)
      endif(LAPACK_LAPACK_WORKS)
    endif(LAPACK_LAPACK_LIBRARY)
  endif(NOT LAPACK_FOUND)
  # Generic LAPACK rs6k library?
  if(NOT LAPACK_FOUND)
    find_library(LAPACK_RS6K_LIBRARY
    NAMES lapack_rs6k
    PATHS /usr/local/lib /usr/lib
    )
    mark_as_advanced(LAPACK_RS6K_LIBRARY)
    if(LAPACK_RS6K_LIBRARY)
      set(LAPACK_LIBRARIES ${LAPACK_RS6K_LIBRARY} ${BLAS_LIBRARIES})
      # Test this combination of libraries.
      set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
      check_fortran_function_exists(cheev LAPACK_RS6K_WORKS)
      mark_as_advanced(LAPACK_RS6K_WORKS)
      set(CMAKE_REQUIRED_LIBRARIES)
      if(LAPACK_RS6K_WORKS)
        set(LAPACK_FOUND TRUE)
      else(LAPACK_RS6K_WORKS)
        set(LAPACK_LIBRARIES)
      endif(LAPACK_RS6K_WORKS)
    endif(LAPACK_RS6K_LIBRARY)
  endif(NOT LAPACK_FOUND)
else(BLAS_FOUND)
  message(STATUS "LAPACK requires BLAS")
endif(BLAS_FOUND)

if(NOT LAPACK_FIND_QUIETLY)
  if(LAPACK_FOUND)
    message(STATUS "A library with LAPACK API found.")
  else(LAPACK_FOUND)
    if(LAPACK_FIND_REQUIRED)
      message(FATAL_ERROR 
      "A library with LAPACK API not found. Please specify library location."
      )
    else(LAPACK_FIND_REQUIRED)
      message(STATUS
      "A library with LAPACK API not found. Please specify library location."
      )
    endif(LAPACK_FIND_REQUIRED)
  endif(LAPACK_FOUND)
endif(NOT LAPACK_FIND_QUIETLY)
