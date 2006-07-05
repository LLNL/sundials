/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-07-05 15:32:35 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see the LICENSE file
 * -----------------------------------------------------------------
 * This is the header file (private version) for the IDA dense
 * linear solver module, IDADENSE.
 * -----------------------------------------------------------------
 */

#ifndef _IDADENSE_IMPL_H
#define _IDADENSE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <stdio.h>

#include <ida/ida_dense.h>

/*
 * -----------------------------------------------------------------
 * Types : IDADenseMemRec, IDADenseMem                             
 * -----------------------------------------------------------------
 */

typedef struct {

  long int d_neq;        /* Neq = problem dimension              */

  IDADenseJacFn d_jac;   /* jac = Jacobian routine to be called  */
  
  DenseMat d_J;          /* J = dF/dy + cj*dF/dy'                */
  
  long int *d_pivots;    /* pivots = pivot array for PJ = LU     */
  
  long int d_nje;        /* nje = no. of calls to jac            */
  
  long int d_nreD;       /* nreD = no. of calls to res due to 
                            diff. quotient Jacobian evaluation   */

  void *d_jdata;         /* jdata is passed to jac               */

  int d_last_flag;       /* last error return flag               */

} IDADenseMemRec, *IDADenseMem;

/*
 * -----------------------------------------------------------------
 * Error Messages
 * -----------------------------------------------------------------
 */

#define MSGD_IDAMEM_NULL "Integrator memory is NULL."
#define MSGD_MEM_FAIL    "A memory request failed."
#define MSGD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGD_LMEM_NULL   "IDADENSE memory is NULL."
#define MSGD_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
