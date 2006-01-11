/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-01-11 21:14:00 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * KINDENSE linear solver module header file (private version)
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _KINDENSE_IMPL_H
#define _KINDENSE_IMPL_H

#include "kinsol_dense.h"

/*
 * -----------------------------------------------------------------
 * Types : KINDenseMemRec, KINDenseMem                             
 * -----------------------------------------------------------------
 * The type KINDenseMem is pointer to a KINDenseMemRec.
 * This structure contains KINDENSE solver-specific data. 
 * -----------------------------------------------------------------
 */

typedef struct {

  long int d_n;       /* problem dimension                      */

  KINDenseJacFn d_jac; /* jac = Jacobian routine to be called   */

  DenseMat d_J;       /* problem Jacobian                       */
  
  long int *d_pivots; /* pivots = pivot array for PM = LU       */
  
  long int d_nje;     /* nje = no. of calls to jac              */

  long int d_nfeD;    /* nfeD = no. of calls to f due to
                         difference quotient approximation of J */
  
  void *d_J_data;     /* J_data is passed to jac                */

  int d_last_flag;    /* last error return flag                 */
  
} KINDenseMemRec, *KINDenseMem;

/* Error Messages */

#define _KINDENSE_        "KINDense-- "
#define MSGDS_KINMEM_NULL  _KINDENSE_ "Solver memory is NULL.\n\n"
#define MSGDS_BAD_NVECTOR _KINDENSE_ "A required vector operation is not implemented.\n\n"
#define MSGDS_MEM_FAIL    _KINDENSE_ "A memory request failed.\n\n"

#define MSGDS_SETGET_KINMEM_NULL "KINDenseSet*/KINDenseGet*-- Solver memory is NULL.\n\n"

#define MSGDS_SETGET_LMEM_NULL "KINDenseSet*/KINDenseGet*-- kindense memory is NULL.\n\n"

#endif

#ifdef __cplusplus
}
#endif
