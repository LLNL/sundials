/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2004-06-02 23:22:22 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/ida/LICENSE
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
#include "sundialstypes.h"
#include "dense.h"
#include "nvector.h"
#include "idadense.h"

/******************************************************************
 *                                                                *
 * Types : IDADenseMemRec, IDADenseMem                            *
 *----------------------------------------------------------------*
 * The type IDADenseMem is pointer to an IDADenseMemRec. This     *
 * structure contains IDADense solver-specific data.              *
 *                                                                *
 ******************************************************************/

typedef struct {

  long int d_neq;        /* Neq = problem dimension              */

  IDADenseJacFn d_jac;   /* jac = Jacobian routine to be called  */
  
  DenseMat d_J;          /* J = dF/dy + cj*dF/dy'                */
  
  long int *d_pivots;    /* pivots = pivot array for PJ = LU     */
  
  long int d_nje;        /* nje = no. of calls to jac            */
  
  long int d_nreD;       /* nreD = no. of calls to res due to 
                            difference quotient Jacobian evaluation */

  void *d_jdata;         /* jdata is passed to jac               */

} IDADenseMemRec, *IDADenseMem;

#endif

#ifdef __cplusplus
}
#endif
