/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2005-01-26 22:23:31 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
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

#include <stdio.h>

#include "cvspgmr.h"

#include "spgmr.h"
#include "iterative.h"
#include "nvector.h"
#include "sundialstypes.h"

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

  CVSpgmrPrecSetupFn g_pset; 
                        /* pset = user-supplied routine to compute      */
                        /* a preconditioner                             */

  CVSpgmrPrecSolveFn g_psolve;   
                        /* psolve = user-supplied routine to solve      */
                        /* preconditioner linear system                 */

  void *g_P_data;       /* P_data passed to psolve and pset             */
  SpgmrMem g_spgmr_mem; /* spgmr_mem is memory used by the              */
                        /* generic Spgmr solver                         */

  CVSpgmrJacTimesVecFn g_jtimes;  
                        /* jtimes = Jacobian * vector routine           */
  void *g_j_data;       /* j_data is passed to jtimes                   */

  int g_last_flag;      /* last error flag returned by any function     */

} CVSpgmrMemRec, *CVSpgmrMem;

/* Error Messages */

#define _CVSPGMR_         "CVSpgmr-- "
#define MSGS_CVMEM_NULL   _CVSPGMR_ "Integrator memory is NULL.\n\n"
#define MSGS_MEM_FAIL     _CVSPGMR_ "A memory request failed.\n\n"
#define MSGS_BAD_PRETYPE1 _CVSPGMR_ "Illegal value for pretype.\n"
#define MSGS_BAD_PRETYPE2 "The legal values are PREC_NONE, PREC_LEFT, "
#define MSGS_BAD_PRETYPE3 "PREC_RIGHT, and PREC_BOTH.\n\n"
#define MSGS_BAD_PRETYPE  MSGS_BAD_PRETYPE1 MSGS_BAD_PRETYPE2 MSGS_BAD_PRETYPE3
#define MSGS_PSOLVE_REQ   _CVSPGMR_ "pretype != PREC_NONE, but PSOLVE = NULL is illegal.\n\n"
#define MSGS_BAD_NVECTOR  _CVSPGMR_ "A required vector operation is not implemented.\n\n"

#define MSGS_SETGET_CVMEM_NULL "CVSpgmrSet*/CVSpgmrGet*-- Integrator memory is NULL.\n\n"

#define MSGS_SETGET_LMEM_NULL "CVSpgmrSet*/CVSpgmrGet*-- cvspgmr memory is NULL.\n\n"

#define MSGS_SET_BAD_PRETYPE1 "CVSpgmrSetPrecType-- Illegal value for pretype.\n"
#define MSGS_SET_BAD_PRETYPE2 "The legal values are PREC_NONE, PREC_LEFT, "
#define MSGS_SET_BAD_PRETYPE3 "PREC_RIGHT, and PREC_BOTH.\n\n"
#define MSGS_SET_BAD_PRETYPE  MSGS_SET_BAD_PRETYPE1 MSGS_SET_BAD_PRETYPE2 MSGS_SET_BAD_PRETYPE3

#define MSGS_SET_BAD_GSTYPE1 "CVSpgmrSetGSType-- Illegal value for gstype.\n"
#define MSGS_SET_BAD_GSTYPE2 "The legal values are MODIFIED_GS and "
#define MSGS_SET_BAD_GSTYPE3 "CLASSICAL_GS.\n\n"
#define MSGS_SET_BAD_GSTYPE  MSGS_SET_BAD_GSTYPE1 MSGS_SET_BAD_GSTYPE2 MSGS_SET_BAD_GSTYPE3

#define MSGS_SET_BAD_DELT "CVSpgmrSetDelt-- delt < 0 illegal.\n\n"

#ifdef __cplusplus
}
#endif

#endif
