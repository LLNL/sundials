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
 * CVSPTFQMR linear solver module header file (private version)
 * -----------------------------------------------------------------
 */

#ifndef _CVSSPTFQMR_IMPL_H
#define _CVSSPTFQMR_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvodes_sptfqmr.h"

  /*
   * -----------------------------------------------------------------
   * Types : struct CVSptfqmrMemRec and struct *CVSptfqmrMem
   * -----------------------------------------------------------------
   * A variable declaration of type struct *CVSptfqmrMem denotes a
   * pointer to a data structure of type struct CVSptfqmrMemRec. The
   * CVSptfqmrMemRec structure contains fields that must be accessed
   * by CVSPTFQMR/SPTFQMR solver module routines.
   * -----------------------------------------------------------------
   */

  typedef struct {

    int q_pretype;        /* type of preconditioning                      */
    realtype q_sqrtN;     /* sqrt(N)                                      */
    realtype q_delt;      /* delt = user specified or DELT_DEFAULT        */
    realtype q_deltar;    /* deltar = delt * tq4                          */
    realtype q_delta;     /* delta = deltar * sqrtN                       */
    int q_maxl;           /* maxl = maximum dimension of the Krylov space */

    long int q_nstlpre;   /* value of nst at the last pset call           */
    long int q_npe;       /* npe = total number of pset calls             */
    long int q_nli;       /* nli = total number of linear iterations      */
    long int q_nps;       /* nps = total number of psolve calls           */
    long int q_ncfl;      /* ncfl = total number of convergence failures  */
    long int q_njtimes;   /* njtimes = total number of calls to jtimes    */
    long int q_nfeSQ;     /* nfeSQ = total number of calls to f for     
                             difference quotient Jacobian-vector products */

    N_Vector q_ytemp;     /* temp vector passed to jtimes and psolve      */
    N_Vector q_x;         /* temp vector used by CVSptfqmrSolve           */
    N_Vector q_ycur;      /* CVODE current y vector in Newton Iteration   */
    N_Vector q_fcur;      /* fcur = f(tn, ycur)                           */

    CVSpilsPrecSetupFn q_pset; 
    /* pset = user-supplied routine to compute      */
    /* a preconditioner                             */

    CVSpilsPrecSolveFn q_psolve;   
    /* psolve = user-supplied routine to solve      */
    /* preconditioner linear system                 */

    void *q_P_data;           /* P_data passed to psolve and pset         */
    SptfqmrMem q_sptfqmr_mem; /* sptfqmr_mem is memory used by the        */
    /* generic Sptfqmr solver                   */

    CVSpilsJacTimesVecFn q_jtimes;  
    /* jtimes = Jacobian * vector routine           */
    void *q_j_data;       /* j_data is passed to jtimes                   */

    int q_last_flag;      /* last error flag returned by any function     */

  } CVSptfqmrMemRec, *CVSptfqmrMem;


  /*
   * -----------------------------------------------------------------
   * CVSPTFQMR error messages
   * -----------------------------------------------------------------
   */

#define MSGTFQMR_CVMEM_NULL "Integrator memory is NULL."
#define MSGTFQMR_MEM_FAIL "A memory request failed."
#define MSGTFQMR_BAD_PRETYPE "Illegal value for pretype. Legal values are PREC_NONE, PREC_LEFT, PREC_RIGHT, and PREC_BOTH."
#define MSGTFQMR_PSOLVE_REQ "pretype != PREC_NONE, but PSOLVE = NULL is illegal."
#define MSGTFQMR_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGTFQMR_LMEM_NULL "CVSPTFQMR memory is NULL."
#define MSGTFQMR_BAD_DELT "delt < 0 illegal."

#define MSGTFQMR_CAMEM_NULL "cvadj_mem = NULL illegal."
#define MSGTFQMR_LMEMB_NULL "CVSPTFQMR memory is NULL for the backward integration."

#ifdef __cplusplus
}
#endif

#endif
