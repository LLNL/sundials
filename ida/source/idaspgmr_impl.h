/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2004-06-02 23:22:22 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/ida/LICENSE
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
#include "sundialstypes.h"
#include "spgmr.h"
#include "nvector.h"
#include "idaspgmr.h"

/******************************************************************
 *                                                                *           
 * Types : IDASpgmrMemRec, IDASpgmrMem                            *
 *----------------------------------------------------------------*
 * The type IDASpgmrMem is pointer to an IDASpgmrMemRec. This     *
 * structure contains IDASpgmr solver-specific data.              *
 *                                                                *
 ******************************************************************/

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

} IDASpgmrMemRec, *IDASpgmrMem;

#endif

#ifdef __cplusplus
}
#endif
