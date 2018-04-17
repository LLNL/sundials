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
# examples/sunlinsol level CMakeLists.txt for SUNDIALS (for CMake)
# ---------------------------------------------------------------

# Always add the serial sunlinearsolver dense/band/spils examples
ADD_SUBDIRECTORY(band)
ADD_SUBDIRECTORY(dense)

IF(KLU_FOUND)
  ADD_SUBDIRECTORY(klu)
ENDIF(KLU_FOUND)

IF(SUPERLUMT_FOUND)
  ADD_SUBDIRECTORY(superlumt)
ENDIF(SUPERLUMT_FOUND)

IF(LAPACK_FOUND)
  ADD_SUBDIRECTORY(lapackband)
  ADD_SUBDIRECTORY(lapackdense)
ENDIF(LAPACK_FOUND)

# Always add serial sunlinearsolver iterative examples
ADD_SUBDIRECTORY(spgmr/serial)
ADD_SUBDIRECTORY(spfgmr/serial)
ADD_SUBDIRECTORY(spbcgs/serial)
ADD_SUBDIRECTORY(sptfqmr/serial)
ADD_SUBDIRECTORY(pcg/serial)

IF(MPI_C_FOUND)
  ADD_SUBDIRECTORY(spgmr/parallel)
  ADD_SUBDIRECTORY(spfgmr/parallel)
  ADD_SUBDIRECTORY(spbcgs/parallel)
  ADD_SUBDIRECTORY(sptfqmr/parallel)
  #ADD_SUBDIRECTORY(pcg/parallel)
ENDIF()
