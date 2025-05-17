/* -----------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation header file for the eigenvalue
 * estimation of the SUNARNOLDI package.
 * -----------------------------------------------------------------*/

#ifndef _ARNOLDI_IMPL_H
#define _ARNOLDI_IMPL_H

#include <sundials/sundials_arnoldi.h>
#include <sundials/sundials_math.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ARNOLDI module private function prototypes
  ===============================================================*/

 sunbooleantype arnoldi_CheckNVector(N_Vector tmpl);
 sunrealtype arnoldi_Magnitude(const suncomplextype *c);
 int arnoldi_Compare(const void *a, const void *b);
 void* arnoldi_Create_Common(void* Adata, N_Vector q, int maxl, SUNContext sunctx);

 /*
  * -----------------------------------------------------------------
  * LAPACK function
  * -----------------------------------------------------------------
  */

 extern void dgeev_(char* jobvl, char* jobvr, int* n, sunrealtype* a, int* lda,
                    sunrealtype* wr, sunrealtype* wi, sunrealtype* vl, int* ldvl, sunrealtype* vr,
                    int* ldvr, sunrealtype* work, int* lwork, int* info);

/*===============================================================
  Reusable ARNOLDI Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_ARNOLDI_NULL_q             "q is null."
#define MSG_ARNOLDI_BAD_NVECTOR        "Bad NVector."
#define MSG_ARNOLDI_NULL_ATIMES        "ATimes is null."
#define MSG_ARNOLDI_NULL_SUNRHSFN      "SUNRhsFn is null."
#define MSG_ARNOLDI_ATIMES_FAIL_REC    "Atimes recoverable failure"
#define MSG_ARNOLDI_ATIMES_FAIL_UNREC  "Atimes unrecoverable failure"
#define MSG_ARNOLDI_NOT_ENOUGH_ITER    "Number of Krylov subspace is not enough (< 3)"
#define MSG_ARNOLDI_NULL_SUNCTX        "sunctx is null."
#define MSG_ARNOLDI_MEM_FAIL           "ARNOLDI memory fail."
#define MSG_ARNOLDI_GS_FAIL            "ARNOLDI Modified GS fail."
#define MSG_ARNOLDI_LAPACK_FAIL        "Error: LAPACK dgeev failed with info = %d\n"

#ifdef __cplusplus
}
#endif

#endif
