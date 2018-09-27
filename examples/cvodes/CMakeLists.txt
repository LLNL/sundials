# ---------------------------------------------------------------
# Programmer:  David Gardner and Slaven Peles @ LLNL
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
# examples/cvodes level CMakeLists.txt for SUNDIALS (for CMake)
# ---------------------------------------------------------------

# C examples
IF(EXAMPLES_ENABLE_C)
  ADD_SUBDIRECTORY(serial)
  IF(MPI_C_FOUND)
    ADD_SUBDIRECTORY(parallel)
  ENDIF()
  IF(OPENMP_FOUND)
    ADD_SUBDIRECTORY(C_openmp)
  ENDIF()
ENDIF()
