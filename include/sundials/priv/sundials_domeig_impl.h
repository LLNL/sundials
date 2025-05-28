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

#include <sundials/sundials_domeig.h>
#include <sundials/sundials_math.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  DOMEIG module private function prototypes
  ===============================================================*/

 sunbooleantype domeig_CheckNVector(N_Vector tmpl);
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
#define MSG_DOMEIG_NULL_q             "q is null."
#define MSG_DOMEIG_BAD_NVECTOR        "Bad NVector."
#define MSG_DOMEIG_NULL_ATIMES        "ATimes is null."
#define MSG_DOMEIG_ATIMES_FAIL_REC    "Atimes recoverable failure"
#define MSG_DOMEIG_ATIMES_FAIL_UNREC  "Atimes unrecoverable failure"
#define MSG_DOMEIG_NOT_ENOUGH_ITER    "Number of Krylov subspace is not enough (< 2)"
#define MSG_DOMEIG_NULL_SUNCTX        "sunctx is null."
#define MSG_DOMEIG_MEM_FAIL           "DOMEIG memory fail."
#define MSG_DOMEIG_GS_FAIL            "DOMEIG Modified GS fail."
#define MSG_DOMEIG_LAPACK_FAIL        "Error: LAPACK dgeev failed with info = %d\n"

#ifdef __cplusplus
}
#endif

#endif
