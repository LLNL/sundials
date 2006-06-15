/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2006-06-15 15:39:15 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh, and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see the LICENSE file
 * -----------------------------------------------------------------
 * This is the header file (private version) for the IDA band
 * linear solver module, IDABAND. It interfaces between the band
 * module and the integrator when a banded linear solver is
 * appropriate.
 * -----------------------------------------------------------------
 */

#ifndef _IDABAND_IMPL_H
#define _IDABAND_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "ida_band.h"

/*
 * -----------------------------------------------------------------
 * Types : IDABandMemRec, IDABandMem                             
 * -----------------------------------------------------------------
 */

typedef struct {

  long int b_neq;           /* Neq = problem size                           */

  IDABandJacFn b_jac;       /* jac = banded Jacobian routine to be called   */
  
  BandMat b_J;              /* J = dF/dy + cj*dF/dy', banded approximation. */
  
  long int b_mupper;        /* mupper = upper bandwidth of Jacobian matrix. */
  
  long int b_mlower;        /* mlower = lower bandwidth of Jacobian matrix. */
  
  long int b_storage_mu;    /* storage_mu = upper bandwidth with storage for
                               factoring = min(Neq-1, mupper+mlower).       */
  
  long int *b_pivots;       /* pivots = pivot array for PJ = LU             */
  
  long int b_nje;           /* nje = no. of calls to jac                    */
  
  long int b_nreB;          /* nreB = no. of calls to res due to 
                               difference quotient Jacobian evaluation      */

  void *b_jdata;            /* jdata = data structure required by jac.      */
  
  int b_last_flag;          /* last error return flag                       */

} IDABandMemRec, *IDABandMem;

/*
 * -----------------------------------------------------------------
 * Error Messages 
 * -----------------------------------------------------------------
 */

#define MSGB_IDAMEM_NULL "Integrator memory is NULL."
#define MSGB_MEM_FAIL    "A memory request failed."
#define MSGB_BAD_SIZES   "Illegal bandwidth parameter(s). Must have 0 <=  mlower, mupper <= N-1."
#define MSGB_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGB_LMEM_NULL   "IDABAND memory is NULL."
#define MSGB_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
