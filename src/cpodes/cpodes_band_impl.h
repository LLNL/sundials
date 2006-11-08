/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:07:05 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Implementation header file for the band linear solver, CPBAND.
 * -----------------------------------------------------------------
 */

#ifndef _CPBAND_IMPL_H
#define _CPBAND_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cpodes/cpodes_band.h>

  /*
   * -----------------------------------------------------------------
   * CPBAND solver constants
   * -----------------------------------------------------------------
   * CPB_MSBJ : maximum number of steps between band Jacobian
   *            evaluations
   *
   * CPB_DGMAX : maximum change in gamma between band Jacobian
   *             evaluations
   * -----------------------------------------------------------------
   */
  
#define CPB_MSBJ  50
#define CPB_DGMAX RCONST(0.2)
  
  /*
   * -----------------------------------------------------------------
   * Types: CPBandMemRec, CPBandMem                                
   * -----------------------------------------------------------------
   * The type CPBandMem is pointer to a CPBandMemRec.
   * This structure contains CPBand solver-specific data.                
   * -----------------------------------------------------------------
   */                                                                

  typedef struct {
    
    long int b_n;           /* problem dimension                           */

    CPBandJacExplFn b_jacE; /* Jacobian routine to be called (CP_EXPL)     */
    CPBandJacImplFn b_jacI; /* Jacobian routine to be called (CP_IMPL)     */
    void *b_J_data;         /* J_data is passed to jacE or jacI            */

    long int b_ml;          /* b_ml = lower bandwidth of savedJ            */
    long int b_mu;          /* b_mu = upper bandwidth of savedJ            */ 
    long int b_storage_mu;  /* upper bandwith of M = MIN(N-1,b_mu+b_ml)    */
  
    BandMat b_M;            /* M = I-gamma*df/dy or M = dF/dy'+gamma*dF/dy */
    BandMat b_savedJ;       /* savedJ = old Jacobian (for ode=CP_EXPL)     */
    long int *b_pivots;     /* pivots = pivot array for PM = LU            */
  
    long int b_nstlj;       /* nstlj = nst at last Jacobian eval.          */
    long int b_nje;         /* nje = no. of calls to jac                   */  
    long int b_nfeB;        /* no. of calls to f due to DQ Jacobian approx.*/

    int b_last_flag;        /* last error return flag                      */
  
  } CPBandMemRec, *CPBandMem;

  /* Error Messages */

#define MSGB_CPMEM_NULL "Integrator memory is NULL."

#define MSGB_MEM_FAIL "A memory request failed."
#define MSGB_BAD_SIZES "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGB_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGB_LMEM_NULL "CPBAND memory is NULL."
#define MSGB_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
