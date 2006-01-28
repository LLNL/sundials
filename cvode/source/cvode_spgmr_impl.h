/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-01-28 00:47:27 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * Implementation header file for the scaled, preconditioned GMRES
 * linear solver, CVSPGMR.
 * -----------------------------------------------------------------
 */

#ifndef _CVSPGMR_IMPL_H
#define _CVSPGMR_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvode_spgmr.h"

/*
 * -----------------------------------------------------------------
 * Types : CVSpgmrMemRec, CVSpgmrMem
 * -----------------------------------------------------------------
 * The type CVSpgmrMem is pointer to a CVSpgmrMemRec.
 * This structure contains CVSpgmr solver-specific data.
 * -----------------------------------------------------------------
 */

typedef struct {

  int  g_pretype;       /* type of preconditioning                      */
  int  g_gstype;        /* type of Gram-Schmidt orthogonalization       */
  realtype g_sqrtN;     /* sqrt(N)                                      */
  realtype g_delt;      /* delt = user specified or DELT_DEFAULT        */
  realtype g_deltar;    /* deltar = delt * tq4                          */
  realtype g_delta;     /* delta = deltar * sqrtN                       */
  int  g_maxl;          /* maxl = maximum dimension of the Krylov space */

  long int g_nstlpre;   /* value of nst at the last pset call           */
  long int g_npe;       /* npe = total number of pset calls             */
  long int g_nli;       /* nli = total number of linear iterations      */
  long int g_nps;       /* nps = total number of psolve calls           */
  long int g_ncfl;      /* ncfl = total number of convergence failures  */
  long int g_njtimes;   /* njtimes = total number of calls to jtimes    */
  long int g_nfeSG;     /* nfeSG = total number of calls to f for     
                           difference quotient Jacobian-vector products */

  N_Vector g_ytemp;     /* temp vector passed to jtimes and psolve      */
  N_Vector g_x;         /* temp vector used by CVSpgmrSolve             */
  N_Vector g_ycur;      /* CVODE current y vector in Newton Iteration   */
  N_Vector g_fcur;      /* fcur = f(tn, ycur)                           */

  CVSpilsPrecSetupFn g_pset; 
                        /* pset = user-supplied routine to compute      */
                        /* a preconditioner                             */

  CVSpilsPrecSolveFn g_psolve;   
                        /* psolve = user-supplied routine to solve      */
                        /* preconditioner linear system                 */

  void *g_P_data;       /* P_data passed to psolve and pset             */
  SpgmrMem g_spgmr_mem; /* spgmr_mem is memory used by the              */
                        /* generic Spgmr solver                         */

  CVSpilsJacTimesVecFn g_jtimes;  
                        /* jtimes = Jacobian * vector routine           */
  void *g_j_data;       /* j_data is passed to jtimes                   */

  int g_last_flag;      /* last error flag returned by any function     */

} CVSpgmrMemRec, *CVSpgmrMem;

/* Error Messages */

#define MSGS_CVMEM_NULL "Integrator memory is NULL."

#define MSGS_MEM_FAIL "A memory request failed."

#define MSGS_BAD_PRETYPE "Illegal value for pretype. Legal values are PREC_NONE, PREC_LEFT, PREC_RIGHT, and PREC_BOTH."

#define MSGS_PSOLVE_REQ "pretype != PREC_NONE, but PSOLVE = NULL is illegal."

#define MSGS_BAD_NVECTOR "A required vector operation is not implemented."

#define MSGS_LMEM_NULL "CVSPGMR memory is NULL."

#define MSGS_BAD_GSTYPE "Illegal value for gstype. Legal values are MODIFIED_GS and CLASSICAL_GS."

#define MSGS_BAD_DELT "delt < 0 illegal."

#ifdef __cplusplus
}
#endif

#endif
