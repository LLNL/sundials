/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2004-10-08 15:27:24 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh, and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/idas/LICENSE
 * -----------------------------------------------------------------
 * This is the header file (private version) for the IDA/IDAS band
 * linear solver module, IDABAND. It interfaces between the band
 * module and the integrator when a banded linear solver is
 * appropriate.
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _idaband_impl_h
#define _idaband_impl_h

#include <stdio.h>

#include "idaband.h"

#include "band.h"
#include "nvector.h"
#include "sundialstypes.h"

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

#endif

#ifdef __cplusplus
}
#endif
