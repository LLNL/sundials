/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:01:30 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Implementation header file for the KINLAPACK linear solver.
 * -----------------------------------------------------------------
 */

#ifndef _KINLAPACK_IMPL_H
#define _KINLAPACK_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <kinsol/kinsol_lapack.h>

  /*
   * -----------------------------------------------------------------
   * Types : KINLapackMemRec, KINLapackMem                             
   * -----------------------------------------------------------------
   * KINLapackMem is pointer to a KINLapackMemRec structure.
   * This structure contains KINLapackDense solver-specific data. 
   * -----------------------------------------------------------------
   */

  typedef struct {

    int l_mtype;              /* Type of Jacobians (DENSE or BAND)            */

    int l_n;                  /* problem dimension                            */

    int b_ml;                 /* b_ml = lower bandwidth of savedJ             */
    int b_mu;                 /* b_mu = upper bandwidth of savedJ             */ 
    int b_smu;                /* upper bandwith of M = MIN(N-1,b_mu+b_ml)     */

    KINLapackDenseJacFn d_jac; /* dense Jacobian routine to be called         */

    KINLapackBandJacFn b_jac;  /* band Jacobian routine to be called          */

    void *l_J_data;           /* J_data is passed to jacE or jacI             */

    LapackMat l_J;            /* problem Jacobian                             */

    int *l_pivots;            /* pivots = pivot array for PM = LU             */
  
    long int l_nje;           /* nje = no. of calls to jac                    */

    long int l_nfeDQ;         /* no. of calls to f due to DQ Jacobian approx. */

    int l_last_flag;          /* last error return flag                       */
  
  } KINLapackMemRec, *KINLapackMem;

  /* Error Messages */

#define MSGLS_KINMEM_NULL "Solver memory is NULL."
#define MSGLS_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGLS_BAD_SIZES "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGLS_MEM_FAIL "A memory request failed."
#define MSGLS_LMEM_NULL "KINLAPACK memory is NULL."
#define MSGLS_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
