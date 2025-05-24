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
 * This is the header file for the eigenvalue estimation of
 * the SUNARNOLDI package.
 * -----------------------------------------------------------------*/

#ifndef _ARNOLDI_H
#define _ARNOLDI_H

#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_math.h>
#include <sundials/sundials_iterative.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default ARNOLDI parameters */
#define ARNOLDI_MAXL_DEFAULT 3
#define DEFAULT_POWER_OF_A   0

/* SUNRhsFn type definition */
typedef int (*SUNRhsFn)(sunrealtype t, N_Vector y, N_Vector ydot,
                        void* user_data);

/*===============================================================
  ARNOLDI module data structure
  ===============================================================*/

/*---------------------------------------------------------------
  Types : struct ARNOLDIMemRec, ARNOLDIMem
  ---------------------------------------------------------------
  The type ARNOLDIMem is type pointer to struct
  ARNOLDIMemRec.  This structure contains fields to
  perform an ARNOLDI iteration.
  ---------------------------------------------------------------*/
typedef struct ARNOLDIMemRec
{
  /* ARNOLDI MEMORY specification */
  SUNATimesFn ATimes;   /* User provided ATimes function */
  void* Adata;          /* ATimes function data*/

  N_Vector *V, q;       /* Krylov subspace vectors */

  int maxl;             /* Krylov subspace dimension */
  int power_of_A;       /* Power of A in the preprocessing; initial q = A^{power_of_A}q/||A^{power_of_A}q|| */

  sunrealtype **Hes;    /* Hessenberg matrix Hes */
}* ARNOLDIMem;

// Struct to hold the real and imaginary parts
typedef struct {
    sunrealtype real;
    sunrealtype imag;
} suncomplextype;

/* -------------------------------------
 * Exported Functions for ARNOLDI
 * ------------------------------------- */

/* Creation and Estimation functions */

SUNDIALS_EXPORT void* ArnoldiCreate(SUNATimesFn ATimes, void* AData,
                N_Vector q, int maxl, SUNContext sunctx);

SUNDIALS_EXPORT int ArnoldiComputeHess(ARNOLDIMem arnoldi_mem);

SUNDIALS_EXPORT int ArnoldiPreProcess(ARNOLDIMem arnoldi_mem);

SUNDIALS_EXPORT suncomplextype ArnoldiEstimate(ARNOLDIMem arnoldi_mem);

#ifdef __cplusplus
}
#endif

#endif
