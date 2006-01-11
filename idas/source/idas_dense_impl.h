/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-01-11 21:13:57 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/idas/LICENSE
 * -----------------------------------------------------------------
 * This is the header file (private version) for the IDA/IDAS dense
 * linear solver module, IDADENSE.
 * -----------------------------------------------------------------
 */

#ifndef _IDADENSE_IMPL_H
#define _IDADENSE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "idas_dense.h"

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

#define MSGD_IDAMEM_NULL        "IDADense-- integrator memory is NULL.\n\n"

#define MSGD_MEM_FAIL           "IDADense-- a memory request failed.\n\n"

#define MSGD_BAD_NVECTOR        "IDADense-- a required vector operation is not implemented.\n\n"

#define MSGD_WRONG_NVEC         "IDADense-- incompatible NVECTOR implementation.\n\n"

#define MSGD_SETGET_IDAMEM_NULL "IDADenseSet*/IDADenseGet*-- integrator memory is NULL. \n\n"

#define MSGD_SETGET_LMEM_NULL   "IDADenseSet*/IDADenseGet*-- IDADENSE memory is NULL. \n\n"

#ifdef __cplusplus
}
#endif

#endif
