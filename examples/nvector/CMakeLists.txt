# ---------------------------------------------------------------
# Programmer:  David Gardner, Slaven Peles, and Cody Balos @ LLNL
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
# examples/nvector level CMakeLists.txt for SUNDIALS (for CMake)
# ---------------------------------------------------------------

# Always add the nvector serial examples
ADD_SUBDIRECTORY(serial)

IF(MPI_C_FOUND)
  ADD_SUBDIRECTORY(parallel)
ENDIF()

IF(HYPRE_FOUND)
  ADD_SUBDIRECTORY(parhyp)
ENDIF()

IF(PTHREADS_FOUND)
  ADD_SUBDIRECTORY(pthreads)
ENDIF()

IF(OPENMP_FOUND)
  ADD_SUBDIRECTORY(C_openmp)
ENDIF()

IF(PETSC_FOUND)
  ADD_SUBDIRECTORY(petsc)
ENDIF()

IF(CUDA_FOUND)
  ADD_SUBDIRECTORY(cuda)
  IF(MPI_C_FOUND)
    ADD_SUBDIRECTORY(mpicuda)
  ENDIF()
ENDIF()

IF(RAJA_FOUND)
  ADD_SUBDIRECTORY(raja)
  IF(MPI_C_FOUND)
    ADD_SUBDIRECTORY(mpiraja)
  ENDIF()
ENDIF()
