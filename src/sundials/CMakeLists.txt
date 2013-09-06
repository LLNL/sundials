# ---------------------------------------------------------------
# $Revision: 1.4 $
# $Date: 2009-02-17 02:52:53 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the generic SUNDIALS modules

# From here we only install the generic SUNDIALS headers. 
# The implementations themselves are incorporated in the individual SUNDIALS solver libraries.

INSTALL(CODE "MESSAGE(\"\nInstall shared components\n\")")

# Add variable sundials_HEADERS with the exported SUNDIALS header files
SET(sundials_HEADERS
  sundials_band.h
  sundials_dense.h
  sundials_direct.h
  sundials_iterative.h
  sundials_math.h
  sundials_nvector.h
  sundials_fnvector.h
  sundials_pcg.h
  sundials_sparse.h
  sundials_spbcgs.h
  sundials_spfgmr.h
  sundials_spgmr.h
  sundials_sptfqmr.h
  sundials_types.h
  )

# Add prefix with complete path to the SUNDIALS header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/sundials/ sundials_HEADERS)

# Install the SUNDIALS header files
INSTALL(FILES ${sundials_HEADERS} DESTINATION include/sundials)

# If Blas/Lapack support was enabled, install the Lapack interface headers
IF(LAPACK_FOUND)
  SET(sundials_BL_HEADERS sundials_lapack.h)
  ADD_PREFIX(${sundials_SOURCE_DIR}/include/sundials/ sundials_BL_HEADERS)
  INSTALL(FILES ${sundials_BL_HEADERS} DESTINATION include/sundials)
ENDIF(LAPACK_FOUND)
