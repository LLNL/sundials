/*
 * -----------------------------------------------------------------
 * $Revision: 1.4.2.1 $
 * $Date: 2005-01-26 22:05:10 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Implementation header file for the dense linear solver, CVDENSE.
 * -----------------------------------------------------------------
 */

#ifndef _CVDENSE_IMPL_H
#define _CVDENSE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

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

/* Error Messages */

#define _CVDENSE_         "CVDense-- "
#define MSGDS_CVMEM_NULL  _CVDENSE_ "Integrator memory is NULL.\n\n"
#define MSGDS_BAD_NVECTOR _CVDENSE_ "A required vector operation is not implemented.\n\n"
#define MSGDS_MEM_FAIL    _CVDENSE_ "A memory request failed.\n\n"

#define MSGDS_SETGET_CVMEM_NULL "CVDenseSet*/CVDenseGet*-- Integrator memory is NULL.\n\n"

#define MSGDS_SETGET_LMEM_NULL "CVDenseSet*/CVDenseGet*-- cvdense memory is NULL.\n\n"

#ifdef __cplusplus
}
#endif

#endif
