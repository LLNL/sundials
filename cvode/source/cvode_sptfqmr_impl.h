/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-01-11 21:13:48 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * CVSPTFQMR linear solver module header file (private version)
 * -----------------------------------------------------------------
 */

#ifndef _CVSPTFQMR_IMPL_H
#define _CVSPTFQMR_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvode_sptfqmr.h"

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

/* CVSptfqmr error messages */

#define _CVSPTFQMR_           "CVSptfqmr-- "
#define MSGTFQMR_CVMEM_NULL   _CVSPTFQMR_ "Integrator memory is NULL.\n\n"
#define MSGTFQMR_MEM_FAIL     _CVSPTFQMR_ "A memory request failed.\n\n"
#define MSGTFQMR_BAD_PRETYPE1 _CVSPTFQMR_ "Illegal value for pretype.\n"
#define MSGTFQMR_BAD_PRETYPE2 "The legal values are PREC_NONE, PREC_LEFT, "
#define MSGTFQMR_BAD_PRETYPE3 "PREC_RIGHT, and PREC_BOTH.\n\n"
#define MSGTFQMR_BAD_PRETYPE  MSGTFQMR_BAD_PRETYPE1 MSGTFQMR_BAD_PRETYPE2 MSGTFQMR_BAD_PRETYPE3
#define MSGTFQMR_PSOLVE_REQ   _CVSPTFQMR_ "pretype != PREC_NONE, but PSOLVE = NULL is illegal.\n\n"
#define MSGTFQMR_BAD_NVECTOR  _CVSPTFQMR_ "A required vector operation is not implemented.\n\n"

/* CVSptfqmrSet* and CVSptfqmrGet* error messages */

#define MSGTFQMR_SETGET_CVMEM_NULL "CVSptfqmrSet*/CVSptfqmrGet*-- Integrator memory is NULL.\n\n"
#define MSGTFQMR_SETGET_LMEM_NULL  "CVSptfqmrSet*/CVSptfqmrGet*-- cvsptfqmr memory is NULL.\n\n"

/* CVSptfqmrSetPrecType error messages */

#define MSGTFQMR_SET_BAD_PRETYPE1 "CVSptfqmrSetPrecType-- Illegal value for pretype.\n"
#define MSGTFQMR_SET_BAD_PRETYPE2 "The legal values are PREC_NONE, PREC_LEFT, "
#define MSGTFQMR_SET_BAD_PRETYPE3 "PREC_RIGHT, and PREC_BOTH.\n\n"
#define MSGTFQMR_SET_BAD_PRETYPE  MSGTFQMR_SET_BAD_PRETYPE1 MSGTFQMR_SET_BAD_PRETYPE2 MSGTFQMR_SET_BAD_PRETYPE3

/* CVSptfqmrSetDelt error message */

#define MSGTFQMR_SET_BAD_DELT "CVSptfqmrSetDelt-- delt < 0 illegal.\n\n"

#ifdef __cplusplus
}
#endif

#endif
