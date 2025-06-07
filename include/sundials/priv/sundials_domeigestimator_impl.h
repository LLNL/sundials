/* -----------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
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
 * This is the implementation header file for the eigenvalue
 * estimation of the DOMEIG module.
 * -----------------------------------------------------------------*/

#ifndef _DOMEIG_IMPL_H
#define _DOMEIG_IMPL_H

#include <sundials/sundials_domeigestimator.h>
#include <sundials/sundials_math.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  DOMEIG module private function prototypes
  ===============================================================*/

 SUNErrCode domeig_CheckNVector(N_Vector tmpl);
 sunrealtype domeig_Magnitude(const suncomplextype *c);
 int domeig_Compare(const void *a, const void *b);

 /*
  * -----------------------------------------------------------------
  * LAPACK function
  * -----------------------------------------------------------------
  */

 extern void dgeev_(char* jobvl, char* jobvr, int* n, sunrealtype* a, int* lda,
                    sunrealtype* wr, sunrealtype* wi, sunrealtype* vl, int* ldvl, sunrealtype* vr,
                    int* ldvr, sunrealtype* work, int* lwork, int* info);

/*===============================================================
  Reusable DOMEIG Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_DOMEIG_LAPACK_FAIL        "Error: LAPACK dgeev failed with info = %d\n"

#ifdef __cplusplus
}
#endif

#endif
