/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2004-10-08 15:27:24 $
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

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _idadense_impl_h
#define _idadense_impl_h

#include <stdio.h>

#include "idadense.h"

#include "dense.h"
#include "nvector.h"
#include "sundialstypes.h"

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

#endif

#ifdef __cplusplus
}
#endif
