/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2005-10-04 21:59:34 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2004, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * Implementation header file for the scaled preconditioned Bi-CGSTAB
 * linear solver, CVSPBCG.
 * -----------------------------------------------------------------
 */

#ifndef _CVSPBCG_IMPL_H
#define _CVSPBCG_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <stdio.h>

#include "cvspbcg.h"
#include "iterative.h"

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

/* CVSpbcg error messages */

#define _CVSPBCG_         "CVSpbcg-- "
#define MSGBCG_CVMEM_NULL   _CVSPBCG_ "Integrator memory is NULL.\n\n"
#define MSGBCG_MEM_FAIL     _CVSPBCG_ "A memory request failed.\n\n"
#define MSGBCG_BAD_PRETYPE1 _CVSPBCG_ "Illegal value for pretype.\n"
#define MSGBCG_BAD_PRETYPE2 "The legal values are PREC_NONE, PREC_LEFT, "
#define MSGBCG_BAD_PRETYPE3 "PREC_RIGHT, and PREC_BOTH.\n\n"
#define MSGBCG_BAD_PRETYPE  MSGBCG_BAD_PRETYPE1 MSGBCG_BAD_PRETYPE2 MSGBCG_BAD_PRETYPE3
#define MSGBCG_PSOLVE_REQ   _CVSPBCG_ "pretype != PREC_NONE, but PSOLVE = NULL is illegal.\n\n"
#define MSGBCG_BAD_NVECTOR  _CVSPBCG_ "A required vector operation is not implemented.\n\n"

/* CVSpbcgSet* and CVSpbcgGet* error messages */

#define MSGBCG_SETGET_CVMEM_NULL "CVSpbcgSet*/CVSpbcgGet*-- Integrator memory is NULL.\n\n"
#define MSGBCG_SETGET_LMEM_NULL "CVSpbcgSet*/CVSpbcgGet*-- cvspbcg memory is NULL.\n\n"

/* CVSpbcgSetPrecType error messages */

#define MSGBCG_SET_BAD_PRETYPE1 "CVSpbcgSetPrecType-- Illegal value for pretype.\n"
#define MSGBCG_SET_BAD_PRETYPE2 "The legal values are PREC_NONE, PREC_LEFT, "
#define MSGBCG_SET_BAD_PRETYPE3 "PREC_RIGHT, and PREC_BOTH.\n\n"
#define MSGBCG_SET_BAD_PRETYPE  MSGBCG_SET_BAD_PRETYPE1 MSGBCG_SET_BAD_PRETYPE2 MSGBCG_SET_BAD_PRETYPE3

/* CVSpbcgSetDelt error message */

#define MSGBCG_SET_BAD_DELT "CVSpbcgSetDelt-- delt < 0 illegal.\n\n"

#ifdef __cplusplus
}
#endif

#endif
