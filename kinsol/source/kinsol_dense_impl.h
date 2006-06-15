/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-06-15 15:39:45 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
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

#define MSGDS_KINMEM_NULL "KINSOL memory is NULL."
#define MSGDS_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGDS_MEM_FAIL    "A memory request failed."
#define MSGDS_LMEM_NULL   "KINDENSE memory is NULL."

#endif

#ifdef __cplusplus
}
#endif
