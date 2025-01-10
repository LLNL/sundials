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
 * This is the implementation file for the eigenvalue implementation of
 * the SUNLINSOL package.
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
#define MSG_ARK_NULL_ATIMES   "ATimes is null."
#define MSG_ARK_NULL_q        "q is null."
#define MSG_ARK_NULL_SUNCTX   "sunctx is null."
#define MSG_ARK_BAD_NVECTOR   "Bad NVector."
#define MSG_ARNOLDI_MEM_FAIL  "ARNOLDI memory fail."

#define DEFAULT_POWER_OF_A    0

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
  sunrealtype **Hes;    /* Hessenberg matrix Hes */

  int maxl;             /* Krylov subspace dimension */
  int power_of_A;       /* Power of A in the preprocessing; initial q = A^{power_of_A}q/||A^{power_of_A}q|| */

}* ARNOLDIMem;

// Struct to hold the real and imaginary parts
typedef struct {
    sunrealtype real;
    sunrealtype imag;
} suncomplextype;

/* -------------------------------------
 * Exported Functions for ARNOLDI
 * ------------------------------------- */

/* Creation and Reinitialization functions */

SUNDIALS_EXPORT void* ArnoldiCreate(SUNATimesFn ATimes, void* Adata,
                N_Vector q, int maxl, SUNContext sunctx);

SUNDIALS_EXPORT int ArnoldiComputeHess(ARNOLDIMem arnoldi_mem);

SUNDIALS_EXPORT int ArnoldiPreProcess(ARNOLDIMem arnoldi_mem);

SUNDIALS_EXPORT suncomplextype ArnoldiEstimate(ARNOLDIMem arnoldi_mem);

#ifdef __cplusplus
}
#endif

#endif
