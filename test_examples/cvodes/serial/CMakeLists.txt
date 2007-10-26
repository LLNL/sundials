# ---------------------------------------------------------------
# $Revision: 1.1 $
# $Date: 2007-10-26 16:36:34 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for CVODES serial examples


# Add variable cvodes_examples with the names of the serial CVODES examples

SET(cvodes_examples
  cvsminhx
  cvs2odex
  cvsnegsolx
  cvsdiscx
  cvsfwdswitchx
  cvsicx
  cvtestroot
  )

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE

IF(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODES_LIB sundials_cvodes_static)
  SET(CVODES_BL_LIB sundials_cvodes_bl_static)
  SET(NVECS_LIB sundials_nvecserial_static)
ELSE(LINK_LIBRARY_TYPE MATCHES "static")
  SET(CVODES_LIB sundials_cvodes_shared)
  SET(CVODES_BL_LIB sundials_cvodes_bl_shared)
  SET(NVECS_LIB sundials_nvecserial_shared)
ENDIF(LINK_LIBRARY_TYPE MATCHES "static")

# Add the build target for each CVODES example

FOREACH(example ${cvodes_examples})
  # add build target for example
  ADD_EXECUTABLE(${example} ${example}.c)
  # link to the sundials libraries
  TARGET_LINK_LIBRARIES(${example} ${CVODES_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})
ENDFOREACH(example ${cvodes_examples})

