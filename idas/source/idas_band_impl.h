/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-03-24 15:57:25 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh, and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/idas/LICENSE
 * -----------------------------------------------------------------
 * This is the header file (private version) for the IDAS band
 * linear solver module, IDABAND. It interfaces between the band
 * module and the integrator when a banded linear solver is
 * appropriate.
 * -----------------------------------------------------------------
 */

#ifndef _IDASBAND_IMPL_H
#define _IDASBAND_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "idas_band.h"

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
   * Types : IDABandMemRecB, IDABandMemB       
   * -----------------------------------------------------------------
   * IDABandB attaches such a structure to the lmemB filed of IDAadjMem
   * -----------------------------------------------------------------
   */

  typedef struct {

    IDABandJacFnB b_bjacB;
    void *b_jac_dataB;

  } IDABandMemRecB, *IDABandMemB;


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

#define MSGB_CAMEM_NULL "idaadj_mem = NULL illegal."
#define MSGB_LMEMB_NULL "IDABAND memory is NULL for the backward integration."
#define MSGB_BAD_T "Bad t for interpolation."

#ifdef __cplusplus
}
#endif

#endif
