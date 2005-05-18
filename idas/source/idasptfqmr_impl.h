/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2005-05-18 18:17:33 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/idas/LICENSE.
 * -----------------------------------------------------------------
 * This is the header file (private version) for the scaled
 * preconditioned TFQMR linear solver module, IDASPTFQMR.
 * -----------------------------------------------------------------
 */

#ifndef _IDASPTFQMR_IMPL_H
#define _IDASPTFQMR_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "idasptfqmr.h"
#include "iterative.h"
#include "nvector.h"
#include "sptfqmr.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Types : IDASptfqmrMemRec, IDASptfqmrMem
 * -----------------------------------------------------------------
 */

typedef struct {

  realtype q_sqrtN;    /* sqrt(N)                                       */
  int q_maxl;          /* maxl = maximum dimension of the Krylov space  */
  realtype q_eplifac;  /* eplifac = linear convergence factor           */
  realtype q_dqincfac; /* dqincfac = optional increment factor in Jv    */
  realtype q_epslin;   /* SpgrmSolve tolerance parameter                */

  int q_resflag;       /* flag from last res call                       */
  long int q_npe;      /* npe = total number of precond calls           */   
  long int q_nli;      /* nli = total number of linear iterations       */
  long int q_nps;      /* nps = total number of psolve calls            */
  long int q_ncfl;     /* ncfl = total number of convergence failures   */
  long int q_nreSG;    /* nreSG = total number of calls to res          */    
  long int q_njtimes;  /* njtimes = total number of calls to jtimes     */

  long int q_nst0;     /* nst0 = saved nst (for performance monitor)    */   
  long int q_nni0;     /* nni0 = saved nni (for performance monitor)    */   
  long int q_nli0;     /* nli0 = saved nli (for performance monitor)    */   
  long int q_ncfn0;    /* ncfn0 = saved ncfn (for performance monitor)  */   
  long int q_ncfl0;    /* ncfl0 = saved ncfl (for performance monitor)  */   
  long int q_nwarn;    /* nwarn = no. of warnings (for perf. monitor)   */   

  N_Vector q_ytemp;    /* temp vector used by IDAAtimesDQ               */ 
  N_Vector q_yptemp;   /* temp vector used by IDAAtimesDQ               */ 
  N_Vector q_xx;       /* temp vector used by IDASptfqmrSolve           */
  N_Vector q_ycur;     /* current y vector in Newton iteration          */
  N_Vector q_ypcur;    /* current yp vector in Newton iteration         */
  N_Vector q_rcur;     /* rcur = F(tn, ycur, ypcur)                     */

  IDASpilsPrecSetupFn q_pset;     /* pset = user-supplied routine       */
                                  /* to compute a preconditioner        */

  IDASpilsPrecSolveFn q_psolve;   /* psolve = user-supplied routine to  */
                                  /* solve preconditioner linear system */

  void *q_pdata;                  /* pdata passed to psolve and precond */
  SptfqmrMem q_sptfqmr_mem;       /* sptfqmr_mem is memory used by the  */
                                  /* generic Sptfqmr solver             */

  IDASpilsJacTimesVecFn q_jtimes; /* Jacobian*vector routine            */ 
  void *q_jdata;                  /* data passed to Jtimes              */

  int q_last_flag;                /* last error return flag             */

} IDASptfqmrMemRec, *IDASptfqmrMem;

/*
 * -----------------------------------------------------------------
 * Error and Warning Messages
 * -----------------------------------------------------------------
 */

#if defined(SUNDIALS_EXTENDED_PRECISION)

#define MSGTFQMR_TIME "at t = %Lg, "

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define MSGTFQMR_TIME "at t = %lg, "

#else

#define MSGTFQMR_TIME "at t = %g, "

#endif


/* Error Messages */

#define MSGTFQMR_IDAMEM_NULL        "IDASptfqmr-- integrator memory is NULL.\n\n"

#define MSGTFQMR_MEM_FAIL           "IDASptfqmr-- a memory request failed.\n\n"

#define MSGTFQMR_BAD_NVECTOR        "IDASptfqmr-- a required vector operation is not implemented.\n\n"

#define MSGTFQMR_SETGET_IDAMEM_NULL "IDASptfqmrSet*/IDASptfqmrGet*-- integrator memory is NULL. \n\n"

#define MSGTFQMR_SETGET_LMEM_NULL   "IDASptfqmrSet*/IDASptfqmrGet*-- IDASPTFQMR memory is NULL. \n\n"

#define MSGTFQMR_IDAS_NEG_EPLIFAC   "IDASptfqmrSetEpsLin-- eplifac < 0.0 illegal. \n\n"

#define MSGTFQMR_IDAS_NEG_DQINCFAC  "IDASptfqmrSetIncrementFactor-- dqincfac < 0.0 illegal. \n\n"

/* Warning Messages */

#define MSGTFQMR_WARN1      "Warning. Poor iterative algorithm performance\n"
#define MSGTFQMR_WARN       "IDASptfqmrPerf-- " MSGTFQMR_TIME MSGTFQMR_WARN1 

#if defined(SUNDIALS_EXTENDED_PRECISION)

#define MSGTFQMR_AVD_WARN1  "Average number of linear iterations is %Le.\n\n"
#define MSGTFQMR_AVD_WARN   MSGTFQMR_WARN MSGTFQMR_AVD_WARN1

#define MSGTFQMR_CFN_WARN1  "Nonlinear convergence failure rate is %Le.\n\n"
#define MSGTFQMR_CFN_WARN   MSGTFQMR_WARN MSGTFQMR_CFN_WARN1

#define MSGTFQMR_CFL_WARN1  "Linear convergence failure rate is %Le.\n\n"
#define MSGTFQMR_CFL_WARN   MSGTFQMR_WARN MSGTFQMR_CFL_WARN1

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define MSGTFQMR_AVD_WARN1  "Average number of linear iterations is %le.\n\n"
#define MSGTFQMR_AVD_WARN   MSGTFQMR_WARN MSGTFQMR_AVD_WARN1

#define MSGTFQMR_CFN_WARN1  "Nonlinear convergence failure rate is %le.\n\n"
#define MSGTFQMR_CFN_WARN   MSGTFQMR_WARN MSGTFQMR_CFN_WARN1

#define MSGTFQMR_CFL_WARN1  "Linear convergence failure rate is %le.\n\n"
#define MSGTFQMR_CFL_WARN   MSGTFQMR_WARN MSGTFQMR_CFL_WARN1

#else

#define MSGTFQMR_AVD_WARN1  "Average number of linear iterations is %e.\n\n"
#define MSGTFQMR_AVD_WARN   MSGTFQMR_WARN MSGTFQMR_AVD_WARN1

#define MSGTFQMR_CFN_WARN1  "Nonlinear convergence failure rate is %e.\n\n"
#define MSGTFQMR_CFN_WARN   MSGTFQMR_WARN MSGTFQMR_CFN_WARN1

#define MSGTFQMR_CFL_WARN1  "Linear convergence failure rate is %e.\n\n"
#define MSGTFQMR_CFL_WARN   MSGTFQMR_WARN MSGTFQMR_CFL_WARN1

#endif

#ifdef __cplusplus
}
#endif

#endif
