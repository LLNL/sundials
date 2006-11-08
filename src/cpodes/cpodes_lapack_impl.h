/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-08 01:07:06 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Implementation header file for the CPLAPACK linear solver.
 * -----------------------------------------------------------------
 */

#ifndef _CPLAPACK_IMPL_H
#define _CPLAPACK_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cpodes/cpodes_lapack.h>

  /*
   * -----------------------------------------------------------------
   * CPLAPACK solver constants
   * -----------------------------------------------------------------
   * CPL_MSBJ : maximum number of steps between dense Jacobian
   *            evaluations
   *
   * CPL_DGMAX : maximum change in gamma between dense Jacobian
   *             evaluations
   * -----------------------------------------------------------------
   */
  
#define CPL_MSBJ  50
#define CPL_DGMAX RCONST(0.2)

  /*
   * -----------------------------------------------------------------
   * Types : CPLapackMemRec, CPLapackMem                             
   * -----------------------------------------------------------------
   * CPLapackMem is pointer to a CPLapackMemRec structure.
   * This structure contains CPLapackDense solver-specific data. 
   * -----------------------------------------------------------------
   */

  typedef struct {

    int l_mtype;              /* Type of Jacobians (DENSE or BAND)            */

    int l_n;                  /* problem dimension                            */

    int b_ml;                 /* b_ml = lower bandwidth of savedJ             */
    int b_mu;                 /* b_mu = upper bandwidth of savedJ             */ 
    int b_smu;                /* upper bandwith of M = MIN(N-1,b_mu+b_ml)     */

    CPLapackDenseJacExplFn d_jacE; /* Jacobian routine to be called (CP_EXPL) */
    CPLapackDenseJacImplFn d_jacI; /* Jacobian routine to be called (CP_IMPL) */

    CPLapackBandJacExplFn b_jacE;  /* Jacobian routine to be called (CP_EXPL) */
    CPLapackBandJacImplFn b_jacI;  /* Jacobian routine to be called (CP_IMPL) */

    void *l_J_data;           /* J_data is passed to jacE or jacI             */

    LapackMat l_M;            /* M = I-gamma*df/dy or M = dF/dy'+gamma*dF/dy  */
    LapackMat l_savedJ;       /* savedJ = old Jacobian (for ode=CP_EXPL)      */
    int *l_pivots;            /* pivots = pivot array for PM = LU             */
  
    long int  l_nstlj;        /* nstlj = nst at last Jacobian eval.           */

    long int l_nje;           /* nje = no. of calls to jac                    */

    long int l_nfeDQ;         /* no. of calls to f due to DQ Jacobian approx. */

    int l_last_flag;          /* last error return flag                       */
  
  } CPLapackMemRec, *CPLapackMem;

  /*
   * -----------------------------------------------------------------
   * Types : CPLapackProjMemRec, CPLapackProjMem                             
   * -----------------------------------------------------------------
   * CPLapackProjMem is pointer to a CPLapackProjMemRec.
   * This structure contains CPLapackDenseProj solver-specific data. 
   * -----------------------------------------------------------------
   */

  typedef struct {

    int l_mtypeP;             /* type of Jacobians                            */

    int l_ny;                 /* number of states                             */

    int l_nc;                 /* number of constraints                        */

    CPLapackDenseProjJacFn l_jacP; /* Jacobian routine to be called           */

    void *l_JP_data;          /* J_data is passed to jacP                     */

    int l_ftype;              /* factorization type (LU, QR, or SC)           */

    int l_pnorm;              /* projection norm (L2 or WRMS)                 */

    LapackMat l_G;            /* G = (dc/dy)^T, transpose of cnstr. Jacobian  */
    LapackMat l_savedG;       /* saved Jacobian (before factorization)        */

    int l_nr;                 /* no. of independent constraints (ftype QRP)   */

    LapackMat l_K;            /* K matrix (s.p.d., form depends on ftype)     */

    int *l_pivotsP;           /* pivotsP = pivot array (for ftype LU or QRP)  */

    realtype *l_beta;         /* beta array (for ftype QR or QRP)             */

    realtype *l_wrk;          /* work array (for ftype QR or QRP)             */

    int l_len_wrk;            /* length of work array                         */

    long int  l_nstljP;       /* nstljP = nst at last Jacobian eval.          */

    long int l_njeP;          /* njeP = no. of calls to jacP                  */

    long int l_nceDQ;         /* no. of calls to c due to DQ Jacobian approx. */

  } CPLapackProjMemRec, *CPLapackProjMem;

  /* Error Messages */

#define MSGLS_CPMEM_NULL "Integrator memory is NULL."
#define MSGLS_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGLS_BAD_SIZES "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGLS_MEM_FAIL "A memory request failed."
#define MSGLS_LMEM_NULL "CPLAPACK memory is NULL."
#define MSGLS_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."

#define MSGLS_BAD_FACT "fact_type has an illegal value."

#ifdef __cplusplus
}
#endif

#endif
