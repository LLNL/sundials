/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2004-11-05 23:35:41 $
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

/*
 * -----------------------------------------------------------------
 * Error Messages 
 * -----------------------------------------------------------------
 */

#define _IDABAND_              "IDABand-- "

#define MSGB_IDAMEM_NULL        _IDABAND_ "Integrator memory is NULL.\n\n"

#define MSGB_BAD_SIZES1         _IDABAND_ "illegal bandwidth parameter(s) "
#define MSGB_BAD_SIZES2         "Must have 0 <=  mlower, mupper <= N-1.\n\n"
#define MSGB_BAD_SIZES          MSGB_BAD_SIZES1 MSGB_BAD_SIZES2

#define MSGB_MEM_FAIL           _IDABAND_ "a memory request failed.\n\n"

#define MSGB_BAD_NVECTOR        _IDABAND_ "a required vector operation is not implemented.\n\n"

#define MSGB_WRONG_NVEC         _IDABAND_ "incompatible NVECTOR implementation.\n\n"

#define MSGB_SETGET_IDAMEM_NULL "IDABandSet*/IDABandGet*-- integrator memory is NULL. \n\n"

#define MSGB_SETGET_LMEM_NULL   "IDABandSet*/IDABandGet*-- IDABAND memory is NULL. \n\n"


#endif

#ifdef __cplusplus
}
#endif
