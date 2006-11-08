/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:01:18 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Implementation header file for the CVLAPACK linear solver.
 * -----------------------------------------------------------------
 */

#ifndef _CVLAPACK_IMPL_H
#define _CVLAPACK_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cvode/cvode_lapack.h>

  /*
   * -----------------------------------------------------------------
   * CVLAPACK solver constants
   * -----------------------------------------------------------------
   * CVL_MSBJ : maximum number of steps between dense Jacobian
   *            evaluations
   *
   * CVL_DGMAX : maximum change in gamma between dense Jacobian
   *             evaluations
   * -----------------------------------------------------------------
   */
  
#define CVL_MSBJ  50
#define CVL_DGMAX RCONST(0.2)

  /*
   * -----------------------------------------------------------------
   * Types : CVLapackMemRec, CVLapackMem                             
   * -----------------------------------------------------------------
   * CVLapackMem is pointer to a CVLapackMemRec structure.
   * This structure contains CVLapackDense solver-specific data. 
   * -----------------------------------------------------------------
   */

  typedef struct {

    int l_mtype;              /* Type of Jacobians (DENSE or BAND)            */

    int l_n;                  /* problem dimension                            */

    int b_ml;                 /* b_ml = lower bandwidth of savedJ             */
    int b_mu;                 /* b_mu = upper bandwidth of savedJ             */ 
    int b_smu;                /* upper bandwith of M = MIN(N-1,b_mu+b_ml)     */

    CVLapackDenseJacFn d_jac; /* Jacobian routine to be called                */

    CVLapackBandJacFn b_jac;  /* Jacobian routine to be called                */

    void *l_J_data;           /* J_data is passed to jac                      */

    LapackMat l_M;            /* M = I - gamma * df/dy                        */
    LapackMat l_savedJ;       /* savedJ = old Jacobian                        */
    int *l_pivots;            /* pivots = pivot array for PM = LU             */
  
    long int  l_nstlj;        /* nstlj = nst at last Jacobian eval.           */

    long int l_nje;           /* nje = no. of calls to jac                    */

    long int l_nfeDQ;         /* no. of calls to f due to DQ Jacobian approx. */

    int l_last_flag;          /* last error return flag                       */
  
  } CVLapackMemRec, *CVLapackMem;

  /* Error Messages */

#define MSGLS_CVMEM_NULL "Integrator memory is NULL."
#define MSGLS_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGLS_BAD_SIZES "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGLS_MEM_FAIL "A memory request failed."
#define MSGLS_LMEM_NULL "CVLAPACK memory is NULL."
#define MSGLS_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
