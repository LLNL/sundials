/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2004-11-05 23:55:11 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/idas/LICENSE
 * -----------------------------------------------------------------
 * This is the header file (private version) for the Scaled
 * Preconditioned GMRES linear solver module, IDASPGMR.
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _idaspgmr_impl_h
#define _idaspgmr_impl_h

#include <stdio.h>

#include "idaspgmr.h"

#include "spgmr.h"
#include "iterative.h"
#include "sundialstypes.h"
#include "nvector.h"

/*
 * -----------------------------------------------------------------
 * Types : IDASpgmrMemRec, IDASpgmrMem                             
 * -----------------------------------------------------------------
 */

typedef struct {

  int  g_gstype;       /* type of Gram-Schmidt orthogonalization       */
  realtype g_sqrtN;    /* sqrt(N)                                      */
  int  g_maxl;         /* maxl = maximum dimension of the Krylov space */
  int  g_maxrs;        /* maxrs = max. number of GMRES restarts        */
  realtype g_eplifac;  /* eplifac = linear convergence factor          */
  realtype g_dqincfac; /* dqincfac = optional increment factor in Jv   */
  realtype g_epslin;   /* SpgrmSolve tolerance parameter               */

  int g_resflag;       /* flag from last res call                      */
  long int g_npe;      /* npe = total number of precond calls          */   
  long int g_nli;      /* nli = total number of linear iterations      */
  long int g_nps;      /* nps = total number of psolve calls           */
  long int g_ncfl;     /* ncfl = total number of convergence failures  */
  long int g_nreSG;    /* nreSG = total number of calls to res         */    
  long int g_njtimes;  /* njtimes = total number of calls to jtimes    */

  long int g_nst0;     /* nst0 = saved nst (for performance monitor)   */   
  long int g_nni0;     /* nni0 = saved nni (for performance monitor)   */   
  long int g_nli0;     /* nli0 = saved nli (for performance monitor)   */   
  long int g_ncfn0;    /* ncfn0 = saved ncfn (for performance monitor) */   
  long int g_ncfl0;    /* ncfl0 = saved ncfl (for performance monitor) */   
  long int g_nwarn;    /* nwarn = no. of warnings (for perf. monitor)  */   

  N_Vector g_ytemp;    /* temp vector used by IDAAtimesDQ              */ 
  N_Vector g_yptemp;   /* temp vector used by IDAAtimesDQ              */ 
  N_Vector g_xx;       /* temp vector used by IDASpgmrSolve            */
  N_Vector g_ycur;     /* current y vector in Newton iteration         */
  N_Vector g_ypcur;    /* current yp vector in Newton iteration        */
  N_Vector g_rcur;     /* rcur = F(tn, ycur, ypcur)                    */

  IDASpgmrPrecSetupFn g_pset;     /* pset = user-supplied routine      */
                                  /* to compute a preconditioner       */

  IDASpgmrPrecSolveFn g_psolve;   /* psolve = user-supplied routine to */
                                  /* solve preconditioner linear system*/

  void *g_pdata;                  /* pdata passed to psolve and precond*/
  SpgmrMem g_spgmr_mem;           /* spgmr_mem is memory used by the   */
                                  /* generic Spgmr solver              */

  IDASpgmrJacTimesVecFn g_jtimes; /* Jacobian*vector routine           */ 
  void *g_jdata;                  /* data passed to Jtimes             */

  int g_last_flag;                /* last error return flag            */

} IDASpgmrMemRec, *IDASpgmrMem;


/*
 * -----------------------------------------------------------------
 * Error and Warning Messages
 * -----------------------------------------------------------------
 */

#ifdef SUNDIALS_EXTENDED_PRECISION
#define MSGS_TIME "at t = %Lg, "
#else
#define MSGS_TIME "at t = %g, "
#endif


/* Error Messages */

#define MSGS_IDAMEM_NULL        "IDASpgmr-- integrator memory is NULL.\n\n"

#define MSGS_MEM_FAIL           "IDASpgmr-- a memory request failed.\n\n"

#define MSGS_BAD_NVECTOR        "IDASpgmr-- a required vector operation is not implemented.\n\n"

#define MSGS_SETGET_IDAMEM_NULL "IDASpgmrSet*/IDASpgmrGet*-- integrator memory is NULL. \n\n"

#define MSGS_SETGET_LMEM_NULL   "IDASpgmrSet*/IDASpgmrGet*-- IDASPGMR memory is NULL. \n\n"

#define MSGS_BAD_GSTYPE         "IDASpgmrSetGSType-- gstype has an illegal value.\n"

#define MSGS_IDAS_NEG_MAXRS     "IDASpgmrSetMaxRestarts-- maxrs < 0 illegal. \n\n"

#define MSGS_IDAS_NEG_EPLIFAC   "IDASpgmrSetEpsLin-- eplifac < 0.0 illegal. \n\n"

#define MSGS_IDAS_NEG_DQINCFAC  "IDASpgmrSetIncrementFactor-- dqincfac < 0.0 illegal. \n\n"

/* Warning Messages */

#define MSGS_WARN1      "Warning. Poor iterative algorithm performance\n"
#define MSGS_WARN       "IDASpgmrPerf-- " MSGS_TIME MSGS_WARN1 

#define MSGS_AVD_WARN1  "Average number of linear iterations is %e.\n\n"
#define MSGS_AVD_WARN   MSGS_WARN MSGS_AVD_WARN1

#define MSGS_CFN_WARN1  "Nonlinear convergence failure rate is %e.\n\n"
#define MSGS_CFN_WARN   MSGS_WARN MSGS_CFN_WARN1

#define MSGS_CFL_WARN1  "Linear convergence failure rate is %e.\n\n"
#define MSGS_CFL_WARN   MSGS_WARN MSGS_CFL_WARN1

#endif

#ifdef __cplusplus
}
#endif
