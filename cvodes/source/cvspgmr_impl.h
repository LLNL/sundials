/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2004-05-26 18:37:07 $
 * ----------------------------------------------------------------- 
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, and         
 *              Radu Serban @ LLNL                              
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California 
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvodes/LICENSE
 * -----------------------------------------------------------------
 * Implementation header file for the scaled, preconditioned GMRES 
 * linear solver, CVSPGMR.
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvspgmr_impl_h
#define _cvspgmr_impl_h

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

} CVSpgmrMemRec, *CVSpgmrMem;


#endif

#ifdef __cplusplus
}
#endif
