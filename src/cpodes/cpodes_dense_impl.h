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
 * Implementation header file for the dense linear solver, CPDENSE.
 * -----------------------------------------------------------------
 */

#ifndef _CPDENSE_IMPL_H
#define _CPDENSE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cpodes/cpodes_dense.h>

  /*
   * -----------------------------------------------------------------
   * CPDENSE solver constants
   * -----------------------------------------------------------------
   * CPD_MSBJ : maximum number of steps between dense Jacobian
   *            evaluations
   *
   * CPD_DGMAX : maximum change in gamma between dense Jacobian
   *             evaluations
   * -----------------------------------------------------------------
   */
  
#define CPD_MSBJ  50
#define CPD_DGMAX RCONST(0.2)

  /*
   * -----------------------------------------------------------------
   * Types : CPDenseMemRec, CPDenseMem                             
   * -----------------------------------------------------------------
   * The type CPDenseMem is pointer to a CPDenseMemRec.
   * This structure contains CPDense solver-specific data. 
   * -----------------------------------------------------------------
   */
  
  typedef struct {

    long int d_n;             /* problem dimension                           */

    CPDenseJacExplFn d_jacE;  /* Jacobian routine to be called (CP_EXPL)     */
    CPDenseJacImplFn d_jacI;  /* Jacobian routine to be called (CP_IMPL)     */
    void *d_J_data;           /* J_data is passed to jacE or jacI            */

    DenseMat d_M;             /* M = I-gamma*df/dy or M = dF/dy'+gamma*dF/dy */

    DenseMat d_savedJ;        /* savedJ = old Jacobian (for ode=CP_EXPL)     */

    long int *d_pivots;       /* pivots = pivot array for PM = LU            */
  
    long int  d_nstlj;        /* nstlj = nst at last Jacobian eval.          */

    long int d_nje;           /* nje = no. of calls to jac                   */

    long int d_nfeD;          /* no. of calls to f due to DQ Jacobian approx.*/

    int d_last_flag;          /* last error return flag                      */
  
  } CPDenseMemRec, *CPDenseMem;

  /*
   * -----------------------------------------------------------------
   * Types : CPDenseProjMemRec, CPDenseProjMem                             
   * -----------------------------------------------------------------
   * The type CPDenseProjMem is pointer to a CPDenseProjMemRec.
   * This structure contains CPDenseProj solver-specific data. 
   * -----------------------------------------------------------------
   */
  
  typedef struct {

    long int d_nc;            /* number of constraints                       */
    long int d_ny;            /* number of states                            */

    CPDenseProjJacFn d_jacP;  /* Jacobian routine to be called               */
    void *d_JP_data;          /* J_data is passed to jacP                    */

    int d_ftype;              /* factorization type (LU, QR, or SC)          */
    int d_pnorm;              /* projection norm (L2 or WRMS)                */

    DenseMat d_G;             /* G = (dc/dy)^T, transpose of cnstr. Jacobian */
    DenseMat d_savedG;        /* saved Jacobian (before factorization)       */

    DenseMat d_K;             /* K matrix (s.p.d., form depends on ftype)    */
    long int *d_pivotsP;      /* pivotsP = pivot array (for ftype LU)        */

    realtype *d_beta;         /* beta array (for ftype QR)                   */
  
    long int  d_nstljP;       /* nstljP = nst at last Jacobian eval.         */

    long int d_njeP;          /* njeP = no. of calls to jacP                 */

    long int d_nceD;          /* no. of calls to c due to DQ Jacobian approx.*/

  } CPDenseProjMemRec, *CPDenseProjMem;



  /* Error Messages */

#define MSGDS_CPMEM_NULL "Integrator memory is NULL."
#define MSGDS_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGDS_MEM_FAIL "A memory request failed."
#define MSGDS_LMEM_NULL "CPDENSE memory is NULL."
#define MSGDS_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."

#define MSGDS_BAD_FACT "fact_type has an illegal value."

#ifdef __cplusplus
}
#endif

#endif
