/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-01-12 22:53:38 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Implementation header file for the band linear solver, CVBAND.
 * -----------------------------------------------------------------
 */

#ifndef _CVSBAND_IMPL_H
#define _CVSBAND_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvodes_band.h"

  /*
   * -----------------------------------------------------------------
   * Types: CVBandMemRec, CVBandMem                                
   * -----------------------------------------------------------------
   * The type CVBandMem is pointer to a CVBandMemRec.
   * This structure contains CVBand solver-specific data.                
   *
   * CVBand attaches such a structure to the lmem field of CVodeMem
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


  /*
   * -----------------------------------------------------------------
   * Types : CVBandMemRecB, CVBandMemB       
   * -----------------------------------------------------------------
   * CVBandB attaches such a structure to the lmemB filed of CVadjMem
   * -----------------------------------------------------------------
   */

  typedef struct {

    CVBandJacFnB b_bjacB;
    void *b_jac_dataB;

  } CVBandMemRecB, *CVBandMemB;


  /*
   * -----------------------------------------------------------------
   * Error Messages 
   * -----------------------------------------------------------------
   */

#define _CVBAND_         "CVBand-- "
#define MSGB_MEM_FAIL    _CVBAND_ "A memory request failed.\n\n"
#define MSGB_BAD_SIZES_1 _CVBAND_ "Illegal bandwidth parameter(s)."
#define MSGB_BAD_SIZES_2 "Must have 0 <=  ml, mu <= N-1.\n\n"
#define MSGB_BAD_SIZES   MSGB_BAD_SIZES_1 MSGB_BAD_SIZES_2
#define MSGB_BAD_NVECTOR _CVBAND_ "A required vector operation is not implemented.\n\n"
#define MSGB_CVMEM_NULL  _CVBAND_ "Integrator memory is NULL.\n\n"

#define MSGB_SETGET_CVMEM_NULL "CVBandSet*/CVBandGet*-- Integrator memory is NULL.\n\n"

#define MSGB_SETGET_LMEM_NULL "CVBandSet*/CVBandGet*-- CVBAND memory is NULL.\n\n"

#ifdef __cplusplus
}
#endif

#endif
