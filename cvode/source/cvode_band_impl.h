/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2006-06-15 15:38:55 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Implementation header file for the band linear solver, CVBAND.
 * -----------------------------------------------------------------
 */

#ifndef _CVBAND_IMPL_H
#define _CVBAND_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvode_band.h"

/*
 * -----------------------------------------------------------------
 * Types: CVBandMemRec, CVBandMem                                
 * -----------------------------------------------------------------
 * The type CVBandMem is pointer to a CVBandMemRec.
 * This structure contains CVBand solver-specific data.                
 * -----------------------------------------------------------------
 */                                                                

typedef struct {

  long int b_n;           /* N = problem dimension                    */

  CVBandJacFn b_jac;      /* jac = Jacobian routine to be called      */

  long int b_ml;          /* b_ml = lower bandwidth of savedJ         */
  
  long int b_mu;          /* b_mu = upper bandwidth of savedJ         */ 
  
  long int b_storage_mu;  /* upper bandwith of M = MIN(N-1,b_mu+b_ml) */
  
  BandMat b_M;            /* M = I - gamma J, gamma = h / l1          */
  
  long int *b_pivots;     /* pivots = pivot array for PM = LU         */
  
  BandMat b_savedJ;       /* savedJ = old Jacobian                    */
  
  long int b_nstlj;       /* nstlj = nst at last Jacobian eval.       */
  
  long int b_nje;         /* nje = no. of calls to jac                */
  
  long int b_nfeB;        /* nfeB = no. of calls to f due to difference
                             quotient band Jacobian approximation     */

  void *b_J_data;         /* J_data is passed to jac                  */

  int b_last_flag;        /* last error return flag                   */
  
} CVBandMemRec, *CVBandMem;

/* Error Messages */

#define MSGB_CVMEM_NULL "Integrator memory is NULL."

#define MSGB_MEM_FAIL "A memory request failed."
#define MSGB_BAD_SIZES "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGB_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGB_LMEM_NULL "CVBAND memory is NULL."
#define MSGB_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
