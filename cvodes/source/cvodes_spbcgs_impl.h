/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-01-28 00:47:17 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Implementation header file for the scaled preconditioned Bi-CGSTAB
 * linear solver, CVSPBCG.
 * -----------------------------------------------------------------
 */

#ifndef _CVSSPBCG_IMPL_H
#define _CVSSPBCG_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvodes_spbcgs.h"

  /*
   * -----------------------------------------------------------------
   * Types : CVSpbcgMemRec, CVSpbcgMem
   * -----------------------------------------------------------------
   * The type CVSpbcgMem is pointer to a CVSpbcgMemRec.
   * This structure contains CVSpbcg solver-specific data.
   * -----------------------------------------------------------------
   */

  typedef struct {

    int b_pretype;        /* type of preconditioning                      */
    realtype b_sqrtN;     /* sqrt(N)                                      */
    realtype b_delt;      /* delt = user specified or DELT_DEFAULT        */
    realtype b_deltar;    /* deltar = delt * tq4                          */
    realtype b_delta;     /* delta = deltar * sqrtN                       */
    int b_maxl;           /* maxl = maximum dimension of the Krylov space */

    long int b_nstlpre;   /* value of nst at the last pset call           */
    long int b_npe;       /* npe = total number of pset calls             */
    long int b_nli;       /* nli = total number of linear iterations      */
    long int b_nps;       /* nps = total number of psolve calls           */
    long int b_ncfl;      /* ncfl = total number of convergence failures  */
    long int b_njtimes;   /* njtimes = total number of calls to jtimes    */
    long int b_nfeSB;     /* nfeSB = total number of calls to f for     
                             difference quotient Jacobian-vector products */

    N_Vector b_ytemp;     /* temp vector passed to jtimes and psolve      */
    N_Vector b_x;         /* temp vector used by CVSpbcgSolve             */
    N_Vector b_ycur;      /* CVODE current y vector in Newton Iteration   */
    N_Vector b_fcur;      /* fcur = f(tn, ycur)                           */

    CVSpilsPrecSetupFn b_pset; 
    /* pset = user-supplied routine to compute      */
    /* a preconditioner                             */

    CVSpilsPrecSolveFn b_psolve;   
    /* psolve = user-supplied routine to solve      */
    /* preconditioner linear system                 */

    void *b_P_data;       /* P_data passed to psolve and pset             */
    SpbcgMem b_spbcg_mem; /* spbcg_mem is memory used by the              */
    /* generic Spbcg solver                         */

    CVSpilsJacTimesVecFn b_jtimes;  
    /* jtimes = Jacobian * vector routine           */
    void *b_j_data;       /* j_data is passed to jtimes                   */

    int b_last_flag;      /* last error flag returned by any function     */

  } CVSpbcgMemRec, *CVSpbcgMem;


  /*
   * -----------------------------------------------------------------
   * CVSPBCG error messages
   * -----------------------------------------------------------------
   */

#define MSGBCG_CVMEM_NULL "Integrator memory is NULL."
#define MSGBCG_MEM_FAIL "A memory request failed."
#define MSGBCG_BAD_PRETYPE "Illegal value for pretype. Legal values are PREC_NONE, PREC_LEFT, PREC_RIGHT, and PREC_BOTH."
#define MSGBCG_PSOLVE_REQ "pretype != PREC_NONE, but PSOLVE = NULL is illegal."
#define MSGBCG_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGBCG_LMEM_NULL "CVSPBCG memory is NULL."
#define MSGBCG_BAD_DELT "delt < 0 illegal."

#define MSGBCG_CAMEM_NULL "cvadj_mem = NULL illegal."
#define MSGBCG_LMEMB_NULL "CVSPBCG memory is NULL for the backward integration."

#ifdef __cplusplus
}
#endif

#endif
