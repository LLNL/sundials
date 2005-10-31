/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2005-10-31 22:02:53 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2004, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/ida/LICENSE.
 * -----------------------------------------------------------------
 * This is the header file (private version) for the scaled
 * preconditioned Bi-CGSTAB linear solver module, IDASPBCG.
 * -----------------------------------------------------------------
 */

#ifndef _IDASPBCG_IMPL_H
#define _IDASPBCG_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "idaspbcg.h"
#include "iterative.h"
#include "nvector.h"
#include "spbcg.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Types : IDASpbcgMemRec, IDASpbcgMem
 * -----------------------------------------------------------------
 */

typedef struct {

  realtype b_sqrtN;    /* sqrt(N)                                       */
  int b_maxl;          /* maxl = maximum dimension of the Krylov space  */
  realtype b_eplifac;  /* eplifac = linear convergence factor           */
  realtype b_dqincfac; /* dqincfac = optional increment factor in Jv    */
  realtype b_epslin;   /* SpgrmSolve tolerance parameter                */

  int b_resflag;       /* flag from last res call                       */
  long int b_npe;      /* npe = total number of precond calls           */   
  long int b_nli;      /* nli = total number of linear iterations       */
  long int b_nps;      /* nps = total number of psolve calls            */
  long int b_ncfl;     /* ncfl = total number of convergence failures   */
  long int b_nreSB;    /* nreSB = total number of calls to res          */    
  long int b_njtimes;  /* njtimes = total number of calls to jtimes     */

  long int b_nst0;     /* nst0 = saved nst (for performance monitor)    */   
  long int b_nni0;     /* nni0 = saved nni (for performance monitor)    */   
  long int b_nli0;     /* nli0 = saved nli (for performance monitor)    */   
  long int b_ncfn0;    /* ncfn0 = saved ncfn (for performance monitor)  */   
  long int b_ncfl0;    /* ncfl0 = saved ncfl (for performance monitor)  */   
  long int b_nwarn;    /* nwarn = no. of warnings (for perf. monitor)   */   

  N_Vector b_ytemp;    /* temp vector used by IDAAtimesDQ               */ 
  N_Vector b_yptemp;   /* temp vector used by IDAAtimesDQ               */ 
  N_Vector b_xx;       /* temp vector used by IDASpbcgSolve             */
  N_Vector b_ycur;     /* current y vector in Newton iteration          */
  N_Vector b_ypcur;    /* current yp vector in Newton iteration         */
  N_Vector b_rcur;     /* rcur = F(tn, ycur, ypcur)                     */

  IDASpilsPrecSetupFn b_pset;     /* pset = user-supplied routine       */
                                  /* to compute a preconditioner        */

  IDASpilsPrecSolveFn b_psolve;   /* psolve = user-supplied routine to  */
                                  /* solve preconditioner linear system */

  void *b_pdata;                  /* pdata passed to psolve and precond */
  SpbcgMem b_spbcg_mem;           /* spbcg_mem is memory used by the    */
                                  /* generic Spbcg solver               */

  IDASpilsJacTimesVecFn b_jtimes; /* Jacobian*vector routine            */ 
  void *b_jdata;                  /* data passed to Jtimes              */

  int b_last_flag;                /* last error return flag             */

} IDASpbcgMemRec, *IDASpbcgMem;

/*
 * -----------------------------------------------------------------
 * Error and Warning Messages
 * -----------------------------------------------------------------
 */

#if defined(SUNDIALS_EXTENDED_PRECISION)

#define MSGBCG_TIME "at t = %Lg, "

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define MSGBCG_TIME "at t = %lg, "

#else

#define MSGBCG_TIME "at t = %g, "

#endif


/* Error Messages */

#define MSGBCG_IDAMEM_NULL        "IDASpbcg-- integrator memory is NULL.\n\n"

#define MSGBCG_MEM_FAIL           "IDASpbcg-- a memory request failed.\n\n"

#define MSGBCG_BAD_NVECTOR        "IDASpbcg-- a required vector operation is not implemented.\n\n"

#define MSGBCG_SETGET_IDAMEM_NULL "IDASpbcgSet*/IDASpbcgGet*-- integrator memory is NULL. \n\n"

#define MSGBCG_SETGET_LMEM_NULL   "IDASpbcgSet*/IDASpbcgGet*-- IDASPBCG memory is NULL. \n\n"

#define MSGBCG_IDAS_NEG_EPLIFAC   "IDASpbcgSetEpsLin-- eplifac < 0.0 illegal. \n\n"

#define MSGBCG_IDAS_NEG_DQINCFAC  "IDASpbcgSetIncrementFactor-- dqincfac < 0.0 illegal. \n\n"

/* Warning Messages */

#define MSGBCG_WARN1      "Warning. Poor iterative algorithm performance\n"
#define MSGBCG_WARN       "IDASpbcgPerf-- " MSGBCG_TIME MSGBCG_WARN1 

#if defined(SUNDIALS_EXTENDED_PRECISION)

#define MSGBCG_AVD_WARN1  "Average number of linear iterations is %Le.\n\n"
#define MSGBCG_AVD_WARN   MSGBCG_WARN MSGBCG_AVD_WARN1

#define MSGBCG_CFN_WARN1  "Nonlinear convergence failure rate is %Le.\n\n"
#define MSGBCG_CFN_WARN   MSGBCG_WARN MSGBCG_CFN_WARN1

#define MSGBCG_CFL_WARN1  "Linear convergence failure rate is %Le.\n\n"
#define MSGBCG_CFL_WARN   MSGBCG_WARN MSGBCG_CFL_WARN1

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define MSGBCG_AVD_WARN1  "Average number of linear iterations is %le.\n\n"
#define MSGBCG_AVD_WARN   MSGBCG_WARN MSGBCG_AVD_WARN1

#define MSGBCG_CFN_WARN1  "Nonlinear convergence failure rate is %le.\n\n"
#define MSGBCG_CFN_WARN   MSGBCG_WARN MSGBCG_CFN_WARN1

#define MSGBCG_CFL_WARN1  "Linear convergence failure rate is %le.\n\n"
#define MSGBCG_CFL_WARN   MSGBCG_WARN MSGBCG_CFL_WARN1

#else

#define MSGBCG_AVD_WARN1  "Average number of linear iterations is %e.\n\n"
#define MSGBCG_AVD_WARN   MSGBCG_WARN MSGBCG_AVD_WARN1

#define MSGBCG_CFN_WARN1  "Nonlinear convergence failure rate is %e.\n\n"
#define MSGBCG_CFN_WARN   MSGBCG_WARN MSGBCG_CFN_WARN1

#define MSGBCG_CFL_WARN1  "Linear convergence failure rate is %e.\n\n"
#define MSGBCG_CFL_WARN   MSGBCG_WARN MSGBCG_CFL_WARN1

#endif

#ifdef __cplusplus
}
#endif

#endif
