/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-02-02 00:36:31 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * KINBAND linear solver module header file (private version)
 * -----------------------------------------------------------------
 */

#ifndef _KINBAND_IMPL_H
#define _KINBAND_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "kinsol_band.h"

/*
 * -----------------------------------------------------------------
 * Types: KINBandMemRec, KINBandMem                                
 * -----------------------------------------------------------------
 * The type KINBandMem is pointer to a KINBandMemRec.
 * This structure contains KINBAND solver-specific data.                
 * -----------------------------------------------------------------
 */                                                                

typedef struct {

  long int b_n;           /* problem dimension                        */

  KINBandJacFn b_jac;      /* jac = Jacobian routine to be called     */

  long int b_ml;          /* b_ml = lower bandwidth of savedJ         */
  
  long int b_mu;          /* b_mu = upper bandwidth of savedJ         */ 
  
  long int b_storage_mu;  /* upper bandwith of M = MIN(N-1,b_mu+b_ml) */
  
  BandMat b_J;            /* problem Jacobian                         */
  
  long int *b_pivots;     /* pivots = pivot array for PM = LU         */
  
  long int b_nje;         /* nje = no. of calls to jac                */
  
  long int b_nfeB;        /* nfeB = no. of calls to f due to difference
                             quotient band Jacobian approximation     */

  void *b_J_data;         /* J_data is passed to jac                  */

  int b_last_flag;        /* last error return flag                   */
  
} KINBandMemRec, *KINBandMem;

/* Error Messages */

#define MSGB_MEM_FAIL    "A memory request failed."
#define MSGB_BAD_SIZES   "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGB_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGB_KINMEM_NULL "KINSOL memory is NULL."
#define MSGB_LMEM_NULL   "KINBAND memory is NULL."

#ifdef __cplusplus
}
#endif

#endif
