/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the header file for a dense SUNMarix implementation using Kokkos.
 * ---------------------------------------------------------------------------*/

#ifndef _SUNMATRIX_KOKKOSDENSE_H
#define _SUNMATRIX_KOKKOSDENSE_H

#include <sundials/sundials_matrix.h>
#include <Kokkos_Core.hpp>

struct _SUNMatrixContent_KokkosDense
{
  sunindextype M;       /* number of rows in block    */
  sunindextype N;       /* number of columns in block */
  sunindextype nblocks; /* number of blocks in matrix */

  /* matrix data view */
  Kokkos::View<sunrealtype***> data_view;
};

typedef struct _SUNMatrixContent_KokkosDense *SUNMatrixContent_KokkosDense;

/* ----------------------------------
 * Implementation specific functions
 * ----------------------------------*/

SUNDIALS_EXPORT SUNMatrix SUNMatrix_KokkosDense(sunindextype M, sunindextype N,
                                                SUNContext sunctx);
SUNDIALS_EXPORT SUNMatrix SUNMatrix_KokkosDenseBlock(sunindextype nblocks,
                                                     sunindextype M,
                                                     sunindextype N,
                                                     SUNContext sunctx);
SUNDIALS_EXPORT sunindextype SUNMatrix_KokkosDense_Rows(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_KokkosDense_Columns(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_KokkosDense_BlockRows(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_KokkosDense_BlockColumns(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_KokkosDense_NumBlocks(SUNMatrix A);

/* ---------------------------------------
 * Implementation of SUNMatrix operations
 * ---------------------------------------*/

SUNDIALS_STATIC_INLINE
SUNMatrix_ID SUNMatGetID_KokkosDense(SUNMatrix A) { return SUNMATRIX_KOKKOSDENSE; }

SUNDIALS_EXPORT SUNMatrix SUNMatClone_KokkosDense(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatDestroy_KokkosDense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatZero_KokkosDense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatCopy_KokkosDense(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAddI_KokkosDense(sunrealtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatScaleAdd_KokkosDense(sunrealtype c, SUNMatrix A,
                                               SUNMatrix B);
SUNDIALS_EXPORT int SUNMatMatvec_KokkosDense(SUNMatrix A, N_Vector x,
                                             N_Vector y);

#endif
