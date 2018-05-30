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
# examples/arkode level CMakeLists.txt for SUNDIALS (for CMake)
# ---------------------------------------------------------------

# C examples
IF(EXAMPLES_ENABLE_C)
  ADD_SUBDIRECTORY(C_serial)
  IF(OPENMP_FOUND)
    ADD_SUBDIRECTORY(C_openmp)
  ENDIF()
  IF(MPI_C_FOUND)
    ADD_SUBDIRECTORY(C_parallel)
  ENDIF()
  IF(HYPRE_ENABLE AND HYPRE_FOUND)
    ADD_SUBDIRECTORY(C_parhyp)
  ENDIF()
ENDIF()

# C++ examples
IF(EXAMPLES_ENABLE_CXX AND CXX_FOUND)
  ADD_SUBDIRECTORY(CXX_serial)
  IF(MPI_CXX_FOUND)
    ADD_SUBDIRECTORY(CXX_parallel)
  ENDIF()
ENDIF()

# F77 examples
IF(EXAMPLES_ENABLE_F77 AND F77_FOUND)
  ADD_SUBDIRECTORY(F77_serial)
  IF(MPI_Fortran_FOUND)
    ADD_SUBDIRECTORY(F77_parallel)
  ENDIF()
ENDIF()

# F90 examples
IF(EXAMPLES_ENABLE_F90 AND F90_FOUND)
  ADD_SUBDIRECTORY(F90_serial)
  IF(MPI_Fortran_FOUND)
    ADD_SUBDIRECTORY(F90_parallel)
  ENDIF()
ENDIF()
