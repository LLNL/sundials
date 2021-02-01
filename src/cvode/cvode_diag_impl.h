/*
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Implementation header file for the diagonal linear solver, CVDIAG.
 * -----------------------------------------------------------------
 */

#ifndef _CVDIAG_IMPL_H
#define _CVDIAG_IMPL_H

#include <cvode/cvode_diag.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


typedef int (*CVDiagFormY)(const realtype h,
                           const realtype r,
                           const N_Vector fpred,
                           const N_Vector zn1,
                           const N_Vector ypred,
                           N_Vector ftemp,
                           N_Vector y);

typedef int (*CVDiagBuildM)(const realtype fract,
                            const realtype uround,
                            const realtype h,
                            const N_Vector ftemp,
                            const N_Vector fpred,
                            const N_Vector ewt,
                            N_Vector bit,
                            N_Vector bitcomp,
                            N_Vector y,
                            N_Vector M);

typedef int(*CVDiagUpdateM)(const realtype r, N_Vector M);

/*
 * -----------------------------------------------------------------
 * Types: CVDiagMemRec, CVDiagMem
 * -----------------------------------------------------------------
 * The type CVDiagMem is pointer to a CVDiagMemRec.
 * This structure contains CVDiag solver-specific data.
 * -----------------------------------------------------------------
 */

typedef struct {

  realtype di_gammasv; /* gammasv = gamma at the last call to setup */
                       /* or solve                                  */

  N_Vector di_M;       /* M = (I - gamma J)^{-1} , gamma = h / l1   */

  N_Vector di_bit;     /* temporary storage vector                  */

  N_Vector di_bitcomp; /* temporary storage vector                  */

  long int di_nfeDI;   /* no. of calls to f due to difference
                          quotient diagonal Jacobian approximation  */

  long int di_last_flag;    /* last error return flag               */

  booleantype   usefused;
  CVDiagFormY   formy_fused;
  CVDiagBuildM  buildM_fused;
  CVDiagUpdateM updateM_fused;

} CVDiagMemRec, *CVDiagMem;

/* Error Messages */

#define MSGDG_CVMEM_NULL "Integrator memory is NULL."
#define MSGDG_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGDG_MEM_FAIL "A memory request failed."
#define MSGDG_LMEM_NULL "CVDIAG memory is NULL."
#define MSGDG_RHSFUNC_FAILED "The right-hand side routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
