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

#include <sundials/sundials_math.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  DOMEIG module private function prototypes
  ===============================================================*/

int sundomeigest_Compare(const void* a, const void* b);

/*
  * -----------------------------------------------------------------
  * LAPACK function
  * -----------------------------------------------------------------
  */
#if defined(SUNDIALS_DOUBLE_PRECISION) || defined(SUNDIALS_EXTENDED_PRECISION)
extern void dgeev_(char* jobvl, char* jobvr, int* n, sunrealtype* a, int* lda,
                   sunrealtype* wr, sunrealtype* wi, sunrealtype* vl, int* ldvl,
                   sunrealtype* vr, int* ldvr, sunrealtype* work, int* lwork,
                   int* info);
#elif defined(SUNDIALS_SINGLE_PRECISION)
extern void sgeev_(char* jobvl, char* jobvr, int* n, sunrealtype* a, int* lda,
                   sunrealtype* wr, sunrealtype* wi, sunrealtype* vl, int* ldvl,
                   sunrealtype* vr, int* ldvr, sunrealtype* work, int* lwork,
                   int* info);
#endif

#ifdef __cplusplus
}
#endif

#endif
