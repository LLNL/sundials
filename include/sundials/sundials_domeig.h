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
 * This is the header file for the eigenvalue estimation of
 * the SUNDOMEIG package.
 * -----------------------------------------------------------------*/

#ifndef _DOMEIG_H
#define _DOMEIG_H

#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., macros */
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_math.h>
#include "sundials/sundials_errors.h"
#include <sundials/sundials_iterative.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default DOMEIG parameters */
#define DEFAULT_POWER_OF_A       0
#define DEFAULT_POWER_ITER_TOL   SUN_RCONST(0.01)
#define DEFAULT_MAX_POWER_ITER   100

/*===============================================================
  DOMEIG module data structure
  ===============================================================*/

// Struct to hold the real and imaginary parts
typedef struct {
    sunrealtype real;
    sunrealtype imag;
} suncomplextype;

/*---------------------------------------------------------------
  Types : struct DOMEIGMemRec, DOMEIGMem
  ---------------------------------------------------------------
  The type DOMEIGMem is type pointer to struct
  DOMEIGMemRec. This structure contains fields to
  perform an DOMEIG iteration.
  ---------------------------------------------------------------*/
typedef struct DOMEIGMemRec
{
  SUNContext sunctx;

  /* DOMEIG MEMORY specification */
  SUNATimesFn ATimes;   /* User provided ATimes function */
  void* Adata;          /* ATimes function data*/

  N_Vector *V, q;       /* Krylov subspace vectors */

  int maxl;             /* Krylov subspace dimension */
  int length;           /* Problem dimension */
  int power_of_A;       /* Power of A in the preprocessing; initial q = A^{power_of_A}q/||A^{power_of_A}q|| */

  sunrealtype powiter_tol; /* Convergence criteria for the power iteration */
  int max_powiter;         /* Maximum number of power iterations */

  sunrealtype* LAPACK_A;      /* The vector which holds rows of the Hessenberg matrix in the given order */
  sunrealtype* LAPACK_wr;     /* Real parts of eigenvalues */
  sunrealtype* LAPACK_wi;     /* Imaginary parts of eigenvalues */
  sunrealtype* LAPACK_work;   /* Workspace array */
  suncomplextype* LAPACK_arr; /* an array to sort eigenvalues*/

  sunrealtype **Hes;    /* Hessenberg matrix Hes */
}* DOMEIGMem;

/* -------------------------------------
 * Exported Functions for DOMEIG
 * ------------------------------------- */

/* Creation and Estimation functions */

SUNDIALS_EXPORT SUNErrCode DomEigCreate(SUNATimesFn ATimes, void* Adata,
                N_Vector q, int maxl, SUNContext sunctx, void** domeig_mem_out);

SUNDIALS_EXPORT SUNErrCode DomEigPreProcess(DOMEIGMem domeig_mem);

SUNDIALS_EXPORT SUNErrCode DomEigComputeHess(DOMEIGMem domeig_mem);

SUNDIALS_EXPORT SUNErrCode DomEigPowerIteration(DOMEIGMem domeig_mem, suncomplextype* dom_eig);

SUNDIALS_EXPORT SUNErrCode DomEigEstimate(DOMEIGMem domeig_mem, suncomplextype* dom_eig);

SUNDIALS_EXPORT void DomEigDestroy(void** domeig_mem);

#ifdef __cplusplus
}
#endif

#endif
