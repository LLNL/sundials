/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2004-08-25 16:19:23 $
 * ----------------------------------------------------------------- 
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, and         
 *              Radu Serban @ LLNL                              
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California 
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvodes/LICENSE
 * -----------------------------------------------------------------
 * Implementation header file for the dense linear solver, CVDENSE. 
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvdense_impl_h
#define _cvdense_impl_h

#include <stdio.h>

#include "cvdense.h"

#include "dense.h"
#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Types : CVDenseMemRec, CVDenseMem                             
 * -----------------------------------------------------------------
 * The type CVDenseMem is pointer to a CVDenseMemRec.
 * This structure contains CVDense solver-specific data. 
 * -----------------------------------------------------------------
 */

typedef struct {

  long int d_n;       /* problem dimension                      */

  CVDenseJacFn d_jac; /* jac = Jacobian routine to be called    */

  DenseMat d_M;       /* M = I - gamma J, gamma = h / l1        */
  
  long int *d_pivots; /* pivots = pivot array for PM = LU   */
  
  DenseMat d_savedJ;  /* savedJ = old Jacobian                  */
  
  long int  d_nstlj;  /* nstlj = nst at last Jacobian eval.     */
  
  long int d_nje;     /* nje = no. of calls to jac              */

  long int d_nfeD;    /* nfeD = no. of calls to f due to
                         difference quotient approximation of J */
  
  void *d_J_data;     /* J_data is passed to jac                */

  int d_last_flag;    /* last error return flag */
  
} CVDenseMemRec, *CVDenseMem;

#endif

#ifdef __cplusplus
}
#endif
