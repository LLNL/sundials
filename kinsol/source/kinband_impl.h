/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2005-08-12 19:30:11 $
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

#include <stdio.h>

#include "kinband.h"
#include "nvector.h"
#include "band.h"
#include "sundialstypes.h"

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

#define _KINBAND_         "KINBand-- "
#define MSGB_MEM_FAIL    _KINBAND_ "A memory request failed.\n\n"
#define MSGB_BAD_SIZES_1 _KINBAND_ "Illegal bandwidth parameter(s)."
#define MSGB_BAD_SIZES_2 "Must have 0 <=  ml, mu <= N-1.\n\n"
#define MSGB_BAD_SIZES   MSGB_BAD_SIZES_1 MSGB_BAD_SIZES_2
#define MSGB_BAD_NVECTOR _KINBAND_ "A required vector operation is not implemented.\n\n"
#define MSGB_KINMEM_NULL  _KINBAND_ "Solver memory is NULL.\n\n"

#define MSGB_SETGET_KINMEM_NULL "KINBandSet*/KINBandGet*-- Solver memory is NULL.\n\n"

#define MSGB_SETGET_LMEM_NULL "KINBandSet*/KINBandGet*-- kinband memory is NULL.\n\n"

#ifdef __cplusplus
}
#endif

#endif
