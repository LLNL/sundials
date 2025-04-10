/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
 * Based on code sundials_direct.h by: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the band implementation of the
 * SUNMATRIX module, SUNMATRIX_BAND.
 *
 * Notes:
 *   - The definition of the generic SUNMatrix structure can be found
 *     in the header file sundials_matrix.h.
 *   - The definition of the type 'sunrealtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the
 *     configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type 'sunbooleantype' and 'indextype'.
 * -----------------------------------------------------------------
 */

#ifndef _SUNMATRIX_BAND_H
#define _SUNMATRIX_BAND_H

#include <stdio.h>
#include <sundials/sundials_matrix.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ---------------------------------
 * Band implementation of SUNMatrix
 * --------------------------------- */

struct _SUNMatrixContent_Band
{
  sunindextype M;
  sunindextype N;
  sunindextype ldim;
  sunindextype mu;
  sunindextype ml;
  sunindextype s_mu;
  sunrealtype* data;
  sunindextype ldata;
  sunrealtype** cols;
};

typedef struct _SUNMatrixContent_Band* SUNMatrixContent_Band;

/* ------------------------------------
 * Macros for access to SUNMATRIX_BAND
 * ------------------------------------ */

#define SM_CONTENT_B(A) ((SUNMatrixContent_Band)(A->content))

#define SM_ROWS_B(A) (SM_CONTENT_B(A)->M)

#define SM_COLUMNS_B(A) (SM_CONTENT_B(A)->N)

#define SM_LDATA_B(A) (SM_CONTENT_B(A)->ldata)

#define SM_UBAND_B(A) (SM_CONTENT_B(A)->mu)

#define SM_LBAND_B(A) (SM_CONTENT_B(A)->ml)

#define SM_SUBAND_B(A) (SM_CONTENT_B(A)->s_mu)

#define SM_LDIM_B(A) (SM_CONTENT_B(A)->ldim)

#define SM_DATA_B(A) (SM_CONTENT_B(A)->data)

#define SM_COLS_B(A) (SM_CONTENT_B(A)->cols)

#define SM_COLUMN_B(A, j) (((SM_CONTENT_B(A)->cols)[j]) + SM_SUBAND_B(A))

#define SM_COLUMN_ELEMENT_B(col_j, i, j) (col_j[(i) - (j)])

#define SM_ELEMENT_B(A, i, j) \
  ((SM_CONTENT_B(A)->cols)[j][(i) - (j) + SM_SUBAND_B(A)])

/* ----------------------------------------
 * Exported  Functions for SUNMATRIX_BAND
 * ---------------------------------------- */

SUNDIALS_EXPORT SUNMatrix SUNBandMatrix(sunindextype N, sunindextype mu,
                                        sunindextype ml, SUNContext sunctx);

SUNDIALS_EXPORT SUNMatrix SUNBandMatrixStorage(sunindextype N, sunindextype mu,
                                               sunindextype ml, sunindextype smu,
                                               SUNContext sunctx);

SUNDIALS_EXPORT void SUNBandMatrix_Print(SUNMatrix A, FILE* outfile);

SUNDIALS_EXPORT sunindextype SUNBandMatrix_Rows(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNBandMatrix_Columns(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNBandMatrix_LowerBandwidth(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNBandMatrix_UpperBandwidth(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNBandMatrix_StoredUpperBandwidth(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNBandMatrix_LDim(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNBandMatrix_LData(SUNMatrix A);
SUNDIALS_EXPORT sunrealtype* SUNBandMatrix_Data(SUNMatrix A);
SUNDIALS_EXPORT sunrealtype** SUNBandMatrix_Cols(SUNMatrix A);
SUNDIALS_EXPORT sunrealtype* SUNBandMatrix_Column(SUNMatrix A, sunindextype j);

SUNDIALS_EXPORT SUNMatrix_ID SUNMatGetID_Band(SUNMatrix A);
SUNDIALS_EXPORT SUNMatrix SUNMatClone_Band(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatDestroy_Band(SUNMatrix A);
SUNDIALS_EXPORT SUNErrCode SUNMatZero_Band(SUNMatrix A);
SUNDIALS_EXPORT SUNErrCode SUNMatCopy_Band(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT SUNErrCode SUNMatScaleAdd_Band(sunrealtype c, SUNMatrix A,
                                               SUNMatrix B);
SUNDIALS_EXPORT SUNErrCode SUNMatScaleAddI_Band(sunrealtype c, SUNMatrix A);
SUNDIALS_EXPORT SUNErrCode SUNMatMatvec_Band(SUNMatrix A, N_Vector x, N_Vector y);
SUNDIALS_EXPORT SUNErrCode SUNMatHermitianTransposeVec_Band(SUNMatrix A,
                                                            N_Vector x,
                                                            N_Vector y);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
SUNErrCode SUNMatSpace_Band(SUNMatrix A, long int* lenrw, long int* leniw);

#ifdef __cplusplus
}
#endif

#endif
