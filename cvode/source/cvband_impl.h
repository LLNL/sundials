/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2004-05-26 19:54:25 $
 * ----------------------------------------------------------------- 
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, and         
 *              Radu Serban @ LLNL                              
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California 
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvode/LICENSE
 * -----------------------------------------------------------------
 * Implementation header file for the band linear solver, CVBAND.
 * -----------------------------------------------------------------
 */
 
#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvband_impl_h
#define _cvband_impl_h

#include <stdio.h>

#include "cvband.h"

#include "band.h"
#include "nvector.h"
#include "sundialstypes.h"

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
  
} CVBandMemRec, *CVBandMem;


#endif

#ifdef __cplusplus
}
#endif
