/*
 * -----------------------------------------------------------------
 * $Revision: 1.32 $
 * $Date: 2004-11-04 01:56:11 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the main KINSol solver.
 * It is independent of the KINSol linear solver in use.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "kinsol_impl.h"
#include "nvector.h"
#include "sundialsmath.h"
#include "sundialstypes.h"
#ifndef _SUNDIALS_CONFIG_H
#define _SUNDIALS_CONFIG_H
#include <sundials_config.h>
#endif

/*
 * -----------------------------------------------------------------
 * private constants
 * -----------------------------------------------------------------
 */

#define HALF           RCONST(0.5)
#define ZERO           RCONST(0.0)
#define ONE            RCONST(1.0)
#define ONEPT5         RCONST(1.5)
#define TWO            RCONST(2.0)
#define THREE          RCONST(3.0)
#define FIVE           RCONST(5.0)
#define TEN            RCONST(10.0)
#define POINT1         RCONST(0.1)
#define POINTOH1       RCONST(0.01)
#define POINT99        RCONST(0.99)
#define THOUSAND       RCONST(1000.0)
#define ONETHIRD       RCONST(.3333333333333333)
#define TWOTHIRDS      RCONST(.6666666666666667)
#define POINT6         RCONST(0.6)
#define POINT9         RCONST(0.9)
#define POINTOHOHONE   RCONST(0.001)
#define POINTOHOHOHONE RCONST(0.0001)

/*
 * -----------------------------------------------------------------
 * default constants
 * -----------------------------------------------------------------
 */
 
#define PRINTFL_DEFAULT 0
#define MXITER_DEFAULT  200
#define MXNBCF          10
#define MSBPRE          10

/* KINStop return value requesting more iterations */

#define CONTINUE_ITERATIONS -999

/*
 * -----------------------------------------------------------------
 * keys for KINPrintInfo
 * -----------------------------------------------------------------
 */

#define PRNT_RETVAL      1
#define PRNT_NNI         2
#define PRNT_TOL         3
#define PRNT_FMAX        4
#define PRNT_PNORM       5
#define PRNT_PNORM1      6
#define PRNT_FNORM       7
#define PRNT_LAM         8
#define PRNT_ALPHA       9
#define PRNT_BETA       10
#define PRNT_ALPHABETA  11
#define PRNT_ADJ        12

/*
 * -----------------------------------------------------------------
 * private helper function prototypes
 * -----------------------------------------------------------------
 */

static booleantype KINCheckNvector(N_Vector tmpl);
static booleantype KINAllocVectors(KINMem kin_mem, N_Vector tmpl);
static int KINSolInit(KINMem kin_mem);
static int KINConstraint(KINMem kin_mem );
static void KINForcingTerm(KINMem kin_mem, realtype fnormp);
static void KINFreeVectors(KINMem kin_mem);
static int  KINInexactNewton(KINMem kin_mem, realtype *fnormp, 
                             realtype *f1normp, booleantype *maxStepTaken);
static int  KINLineSearch(KINMem kin_mem, realtype *fnormp, 
                          realtype *f1normp, booleantype *maxSteptaken);
static booleantype KINInitialStop(KINMem kin_mem);
static int  KINLinSolDrv(KINMem kinmem , N_Vector bb , N_Vector xx );
static realtype KINScFNorm(N_Vector vv , N_Vector scale , N_Vector wrkv);
static realtype KINScSteplength(KINMem kin_mem, N_Vector ucur,
				N_Vector ss, N_Vector usc);
static int KINStop(KINMem kin_mem, booleantype maxStepTaken, int globalstratret);

static void KINPrintInfo(KINMem kin_mem, char *funcname, int key,...);

/* loop macro */
#define loop for(;;)

/*
 * -----------------------------------------------------------------
 * user-callable functions
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : KINCreate
 * -----------------------------------------------------------------
 * KINCreate creates an internal memory block for a problem to 
 * be solved by KINSOL. If successful, KINCreate returns a pointer
 * to the problem memory. This pointer should be passed to
 * KINMalloc. If an initialization error occurs, KINCreate prints
 * an error message to standard error and returns NULL. 
 * -----------------------------------------------------------------
 */

void *KINCreate(void)
{
  KINMem kin_mem;
  realtype uround;

  kin_mem = (KINMem) malloc(sizeof(struct KINMemRec));
  if (kin_mem == NULL) {
    fprintf(stderr, MSG_KINMEM_FAIL);
    return(NULL);
  }

  /* set uround (unit roundoff) */

  kin_mem->kin_uround = uround = UNIT_ROUNDOFF;
  
  /* set default values for solver optional inputs */

  kin_mem->kin_func         = NULL;
  kin_mem->kin_f_data       = NULL;
  kin_mem->kin_constraints  = NULL;
  kin_mem->kin_errfp        = stderr;
  kin_mem->kin_infofp       = stdout;
  kin_mem->kin_printfl      = PRINTFL_DEFAULT;
  kin_mem->kin_mxiter       = MXITER_DEFAULT;
  kin_mem->kin_noPrecInit   = FALSE;
  kin_mem->kin_msbpre       = MSBPRE;
  kin_mem->kin_pthrsh       = TWO;
  kin_mem->kin_noMinEps     = FALSE;
  kin_mem->kin_mxnewtstep   = ZERO;
  kin_mem->kin_sqrt_relfunc = RSqrt(uround);
  kin_mem->kin_scsteptol    = RPowerR(uround,TWOTHIRDS);
  kin_mem->kin_fnormtol     = RPowerR(uround,ONETHIRD);
  kin_mem->kin_etaflag      = KIN_ETACHOICE1;
  kin_mem->kin_eta          = POINT1;  /* default for KIN_ETACONSTANT */
  kin_mem->kin_eta_alpha    = TWO;  /* default for KIN_ETACHOICE2 */
  kin_mem->kin_eta_gamma    = POINT9;  /* default for KIN_ETACHOICE2 */
  kin_mem->kin_MallocDone   = FALSE;
  kin_mem->kin_setupNonNull = FALSE;

  return((void *) kin_mem);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetFdata
 * -----------------------------------------------------------------
 */

int KINSetFdata(void *kinmem, void *f_data)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_f_data = f_data;

  return(KIN_SUCCESS);
}

#define f_data (kin_mem->kin_f_data)

/*
 * -----------------------------------------------------------------
 * Function : KINSetErrFile
 * -----------------------------------------------------------------
 */

int KINSetErrFile(void *kinmem, FILE *errfp)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_errfp = errfp;

  return(KIN_SUCCESS);
}

#define errfp (kin_mem->kin_errfp)

/*
 * -----------------------------------------------------------------
 * Function : KINSetInfoFile
 * -----------------------------------------------------------------
 */

int KINSetInfoFile(void *kinmem, FILE *infofp)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_infofp = infofp;

  return(KIN_SUCCESS);
}

#define infofp (kin_mem->kin_infofp)

/*
 * -----------------------------------------------------------------
 * Function : KINSetPrintLevel
 * -----------------------------------------------------------------
 */

int KINSetPrintLevel(void *kinmem, int printfl)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if ((printfl < 0) || (printfl > 3)) {
    fprintf(errfp, MSG_BAD_PRINTFL);
    return(KIN_ILL_INPUT);
  }

  kin_mem->kin_printfl = printfl;

  return(KIN_SUCCESS);
}

#define printfl (kin_mem->kin_printfl)

/*
 * -----------------------------------------------------------------
 * Function : KINSetNumMaxIters
 * -----------------------------------------------------------------
 */

int KINSetNumMaxIters(void *kinmem, long int mxiter)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (mxiter <= 0) {
    fprintf(errfp, MSG_BAD_MXITER);
    return(KIN_ILL_INPUT);
  }

  kin_mem->kin_mxiter = mxiter;

  return(KIN_SUCCESS);
}

#define mxiter (kin_mem->kin_mxiter)

/*
 * -----------------------------------------------------------------
 * Function : KINSetNoPrecInit
 * -----------------------------------------------------------------
 */

int KINSetNoPrecInit(void *kinmem, booleantype noPrecInit)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_noPrecInit = noPrecInit;

  return(KIN_SUCCESS);
}

#define noPrecInit (kin_mem->kin_noPrecInit)

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxPrecCalls
 * -----------------------------------------------------------------
 */

int KINSetMaxPrecCalls(void *kinmem, long int msbpre)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (msbpre < 0) {
    fprintf(errfp, MSG_BAD_MSBPRE);
    return(KIN_ILL_INPUT);
  }

  kin_mem->kin_msbpre = msbpre;

  return(KIN_SUCCESS);
}

#define msbpre (kin_mem->kin_msbpre)

/*
 * -----------------------------------------------------------------
 * Function : KINSetEtaForm
 * -----------------------------------------------------------------
 */

int KINSetEtaForm(void *kinmem, int etachoice)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if ((etachoice != KIN_ETACONSTANT) && 
      (etachoice != KIN_ETACHOICE1)  && 
      (etachoice != KIN_ETACHOICE2)) {
    fprintf(errfp, MSG_BAD_ETACHOICE);
    return(KIN_ILL_INPUT);
  }

  kin_mem->kin_etaflag = etachoice;

  return(KIN_SUCCESS);
}

#define etaflag (kin_mem->kin_etaflag)

/*
 * -----------------------------------------------------------------
 * Function : KINSetEtaConstValue
 * -----------------------------------------------------------------
 */

int KINSetEtaConstValue(void *kinmem, realtype eta)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if ((eta <= ZERO) || (eta > ONE)) {
    fprintf(errfp, MSG_BAD_ETACONST);
    return(KIN_ILL_INPUT);
  }

  kin_mem->kin_eta = eta;

  return(KIN_SUCCESS);
}

#define eta (kin_mem->kin_eta)

/*
 * -----------------------------------------------------------------
 * Function : KINSetEtaParams
 * -----------------------------------------------------------------
 */

int KINSetEtaParams(void *kinmem, realtype egamma, realtype ealpha)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if ((ealpha <= ONE) || (ealpha > TWO))
    if (ealpha != ZERO) {
      fprintf(errfp, MSG_BAD_ALPHA);
      return(KIN_ILL_INPUT);
    }
  
  if (ealpha == ZERO) 
    kin_mem->kin_eta_alpha = TWO;
  else
    kin_mem->kin_eta_alpha = ealpha;

  if ((egamma <= ZERO) || (egamma > ONE))
    if (egamma != ZERO) {
      fprintf(errfp, MSG_BAD_GAMMA);
      return(KIN_ILL_INPUT);
    }

  if (egamma == ZERO)
    kin_mem->kin_eta_gamma = POINT9;
  else
    kin_mem->kin_eta_gamma = egamma;

  return(KIN_SUCCESS);
}

#define ealpha (kin_mem->kin_eta_alpha)
#define egamma (kin_mem->kin_eta_gamma)

/*
 * -----------------------------------------------------------------
 * Function : KINSetNoMinEps
 * -----------------------------------------------------------------
 */

int KINSetNoMinEps(void *kinmem, booleantype noMinEps)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_noMinEps = noMinEps;

  return(KIN_SUCCESS);
}

#define noMinEps (kin_mem->kin_noMinEps)

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxNewtonStep
 * -----------------------------------------------------------------
 */

int KINSetMaxNewtonStep(void *kinmem, realtype mxnewtstep)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (mxnewtstep <= ZERO) {
    fprintf(errfp, MSG_BAD_MXNEWTSTEP);
    return(KIN_ILL_INPUT);
  }

  kin_mem->kin_mxnewtstep = mxnewtstep;

  return(KIN_SUCCESS);
}

#define mxnewtstep (kin_mem->kin_mxnewtstep)

/*
 * -----------------------------------------------------------------
 * Function : KINSetRelErrFunc
 * -----------------------------------------------------------------
 */

int KINSetRelErrFunc(void *kinmem, realtype relfunc)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (relfunc <= ZERO) {
    fprintf(errfp, MSG_BAD_RELFUNC, relfunc);
    return(KIN_ILL_INPUT);
  }

  kin_mem->kin_sqrt_relfunc = RSqrt(relfunc);

  return(KIN_SUCCESS);
}

#define relfunc (kin_mem->kin_sqrt_relfunc)

/*
 * -----------------------------------------------------------------
 * Function : KINSetFuncNormTol
 * -----------------------------------------------------------------
 */

int KINSetFuncNormTol(void *kinmem, realtype fnormtol)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (fnormtol <= ZERO) {
    fprintf(errfp, MSG_BAD_FNORMTOL, fnormtol);
    return(KIN_ILL_INPUT);
  }

  kin_mem->kin_fnormtol = fnormtol;

  return(KIN_SUCCESS);
}

#define fnormtol (kin_mem->kin_fnormtol)

/*
 * -----------------------------------------------------------------
 * Function : KINSetScaledStepTol
 * -----------------------------------------------------------------
 */

int KINSetScaledStepTol(void *kinmem, realtype scsteptol)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (scsteptol <= ZERO) {
    fprintf(errfp, MSG_BAD_SCSTEPTOL, scsteptol);
    return(KIN_ILL_INPUT);
  }

  kin_mem->kin_scsteptol = scsteptol;

  return(KIN_SUCCESS);
}

#define scsteptol (kin_mem->kin_scsteptol)

/*
 * -----------------------------------------------------------------
 * Function : KINSetConstraints
 * -----------------------------------------------------------------
 */

int KINSetConstraints(void *kinmem, N_Vector constraints)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_constraints = constraints;

  return(KIN_SUCCESS);
}

#define constraints (kin_mem->kin_constraints)

/*
 * -----------------------------------------------------------------
 * Function : KINMalloc
 * -----------------------------------------------------------------
 * KINMalloc allocates memory for a problem or execution of KINSol. 
 * If memory is successfully allocated, KIN_SUCCESS is returned.
 * Otherwise, an error message is printed and an error flag
 * returned.
 * -----------------------------------------------------------------
 */

int KINMalloc(void *kinmem, KINSysFn func, N_Vector tmpl)
{
  long int liw1, lrw1;
  KINMem kin_mem;
  booleantype allocOK, nvectorOK;
  
  /* check kinmem */

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINM_NO_MEM);
    return(KIN_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (func == NULL) {
    fprintf(errfp, MSG_FUNC_NULL);
    return(KIN_ILL_INPUT);
  }

  /* check if all required vector operations are implemented */

  nvectorOK = KINCheckNvector(tmpl);
  if (!nvectorOK) {
    if (errfp != NULL) fprintf(errfp, MSG_BAD_NVECTOR);
    return(KIN_ILL_INPUT);
  }

  /* set space requirements for one N_Vector */

  if (tmpl->ops->nvspace != NULL) {
    N_VSpace(tmpl, &lrw1, &liw1);
    kin_mem->kin_lrw1 = lrw1;
    kin_mem->kin_liw1 = liw1;
  }
  else {
    kin_mem->kin_lrw1 = 0;
    kin_mem->kin_liw1 = 0;
  }

  /* allocate necessary vectors */

  allocOK = KINAllocVectors(kin_mem, tmpl);
  if (!allocOK) {
    fprintf(errfp, MSG_MEM_FAIL);
    free(kin_mem);
    return(KIN_MEM_FAIL);
  }

  /* copy the input parameter into KINSol state */

  kin_mem->kin_func = func;

  /* set the linear solver addresses to NULL */

  kin_mem->kin_linit  = NULL;
  kin_mem->kin_lsetup = NULL;
  kin_mem->kin_lsolve = NULL;
  kin_mem->kin_lfree  = NULL;
  kin_mem->kin_lmem   = NULL;
  
  /* problem memory has been successfully allocated */

  kin_mem->kin_MallocDone = TRUE;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetSysFunc
 * -----------------------------------------------------------------
 */

int KINSetSysFunc(void *kinmem, KINSysFn func)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (func == NULL) {
    fprintf(errfp, MSG_KINS_FUNC_NULL);
    return(KIN_ILL_INPUT);
  }

  kin_mem->kin_func = func;

  return(KIN_SUCCESS);
}

#define func (kin_mem->kin_func)

/*
 * -----------------------------------------------------------------
 * readability constants
 * -----------------------------------------------------------------
 */

#define uround (kin_mem->kin_uround)
#define nni (kin_mem->kin_nni)
#define nfe (kin_mem->kin_nfe)
#define nbcf (kin_mem->kin_nbcf)  
#define nbktrk (kin_mem->kin_nbktrk)
#define ncscmx (kin_mem->kin_ncscmx)
#define stepl (kin_mem->kin_stepl)
#define stepmul (kin_mem->kin_stepmul)
#define pthrsh (kin_mem->kin_pthrsh)
#define linit (kin_mem->kin_linit)
#define lsetup (kin_mem->kin_lsetup)
#define lsolve (kin_mem->kin_lsolve) 
#define lfree (kin_mem->kin_lfree)
#define constraintsSet (kin_mem->kin_constraintsSet) 
#define precondcurrent (kin_mem->kin_precondcurrent)          
#define nnilpre (kin_mem->kin_nnilpre)
#define lmem (kin_mem->kin_lmem)        
#define setupNonNull (kin_mem->kin_setupNonNull)
#define fval (kin_mem->kin_fval)      
#define fnorm (kin_mem->kin_fnorm)
#define f1norm (kin_mem->kin_f1norm)
#define etaflag (kin_mem->kin_etaflag)
#define callForcingTerm (kin_mem->kin_callForcingTerm)
#define uu (kin_mem->kin_uu)
#define uscale (kin_mem->kin_uscale)
#define fscale (kin_mem->kin_fscale)
#define globalstrategy (kin_mem->kin_globalstrategy)     
#define sJpnorm (kin_mem->kin_sJpnorm)
#define sfdotJp (kin_mem->kin_sfdotJp)
#define unew (kin_mem->kin_unew)
#define pp (kin_mem->kin_pp)
#define vtemp1 (kin_mem->kin_vtemp1)
#define vtemp2 (kin_mem->kin_vtemp2)
#define eps (kin_mem->kin_eps)
#define res_norm (kin_mem->kin_res_norm)
#define liw1 (kin_mem->kin_liw1)
#define lrw1 (kin_mem->kin_lrw1)
#define liw (kin_mem->kin_liw)
#define lrw (kin_mem->kin_lrw)

/*
 * -----------------------------------------------------------------
 * Function : KINSol
 * -----------------------------------------------------------------
 * KINSol (main KINSOL driver routine) manages the computational
 * process of computing an approximate solution of the nonlinear
 * system F(uu) = 0. The KINSol routine calls the following
 * subroutines:
 *
 *  KINSolInit  checks if initial guess satisfies user-supplied
 *              constraints and initializes linear solver
 *
 *  KINLinSolDrv  interfaces with linear solver to find a
 *                solution of the system J(uu)*x = b (calculate
 *                Newton step)
 *
 *  KINInexactNewton/KINLineSearch  implement the global strategy
 *
 *  KINForcingTerm  computes the forcing term (eta)
 *
 *  KINStop  determines if an approximate solution has been found
 * -----------------------------------------------------------------
 */

int KINSol(void *kinmem, N_Vector u, int strategy,  
           N_Vector u_scale, N_Vector f_scale)
{
  realtype fnormp, f1normp, epsmin;
  N_Vector bb, xx;
  KINMem kin_mem;
  int ret, globalstratret;
  booleantype maxStepTaken;

  globalstratret = 0;

  /* check for kinmem non-NULL */

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINSOL_NO_MEM);
    return(KIN_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if(kin_mem->kin_MallocDone == FALSE) {
    if (errfp != NULL) fprintf(errfp, MSG_KINSOL_NO_MALLOC);
    else fprintf(stderr, MSG_KINSOL_NO_MALLOC);
    return(KIN_NO_MALLOC);
  }

  /* load input arguments */

  globalstrategy = strategy;
  uu = u;
  uscale = u_scale;
  fscale = f_scale;

  /* initialize solver */

  ret = KINSolInit(kin_mem);
  if (ret != KIN_SUCCESS) return(ret);
  
  ncscmx = 0;

  /* Note: The following logic allows the choice of whether or not
     to force a call to the preconditioner upon a given call to
     KINSol. */

  if (noPrecInit == FALSE) pthrsh = TWO;
  if (noMinEps == FALSE) epsmin = POINTOH1 * fnormtol;

  loop{

    nni++;

    /* calculate the epsilon for this iteration based on eta from the
       routine KINForcingTerm */

    eps = (eta + uround) * fnorm;
    if(!noMinEps) eps = MAX(epsmin, eps);

    /* call KINLinSolDrv to calculate the approximate Newton step using
       the appropriate Krylov algorithm */

    /* use unew as work space for bb, the rhs vector, in the Newton
       steps and pp as xx, the x in Jx = b */

    bb = unew;
    xx = pp;

    /* call KINLinSolDrv and process the return flag */

    ret = KINLinSolDrv(kin_mem, bb, xx );
    if (ret != 0) break;

    /* call the appropriate routine to calculate an acceptable step pp */

    if (globalstrategy == KIN_INEXACT_NEWTON)
      globalstratret = KINInexactNewton(kin_mem, &fnormp, &f1normp, &maxStepTaken);
    else if (globalstrategy == KIN_LINESEARCH)
      globalstratret = KINLineSearch(kin_mem, &fnormp, &f1normp, &maxStepTaken);

    /* if too many beta condition failures, then stop iteration */

    if (nbcf > MXNBCF) {
      ret = KIN_LINESEARCH_BCFAIL;
      break;
    }

    /* evaluate eta by calling the forcing term routine */

    if (callForcingTerm) KINForcingTerm(kin_mem, fnormp);

    fnorm = fnormp;

    /* call KINStop to check if tolerances where met by this iteration */

    ret = KINStop(kin_mem, maxStepTaken, globalstratret); 

    /* update uu after the iteration */

    N_VScale(ONE, unew, uu);

    f1norm = f1normp;

    /* print the current nni, fnorm, and nfe values if printfl > 0 */

    if (printfl>0) 
      KINPrintInfo(kin_mem, "KINSol",PRNT_NNI,nni,nfe,fnorm);

    if (ret != CONTINUE_ITERATIONS) break; 

    fflush(errfp);

 }  /* end of loop so load optional outputs and return */

  if (printfl > 0)
    KINPrintInfo(kin_mem, "KINSol",PRNT_RETVAL,ret);

  fflush(infofp);

  return(ret);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetWorkSpace
 * -----------------------------------------------------------------
 */

int KINGetWorkSpace(void *kinmem, long int *lenrw, long int *leniw)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KING_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  *lenrw = lrw;
  *leniw = liw;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumNonlinSolvIters
 * -----------------------------------------------------------------
 */

int KINGetNumNonlinSolvIters(void *kinmem, long int *nniters)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KING_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  *nniters = nni;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumFuncEvals
 * -----------------------------------------------------------------
 */

int KINGetNumFuncEvals(void *kinmem, long int *nfevals)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KING_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  *nfevals = nfe;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumBetaCondFails
 * -----------------------------------------------------------------
 */

int KINGetNumBetaCondFails(void *kinmem, long int *nbcfails)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KING_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  *nbcfails = nbcf;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumBacktrackOps
 * -----------------------------------------------------------------
 */

int KINGetNumBacktrackOps(void *kinmem, long int *nbacktr)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KING_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  *nbacktr = nbktrk;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetFuncNorm
 * -----------------------------------------------------------------
 */

int KINGetFuncNorm(void *kinmem, realtype *funcnorm)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KING_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  *funcnorm = fnorm;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetStepLength
 * -----------------------------------------------------------------
 */

int KINGetStepLength(void *kinmem, realtype *steplength)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KING_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  *steplength = stepl;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINFree
 * -----------------------------------------------------------------
 * This routine frees the problem memory allocated by KINMalloc.
 * Such memory includes all the vectors allocated by
 * KINAllocVectors, and the memory lmem for the linear solver
 * (deallocated by a call to lfree).
 * -----------------------------------------------------------------
 */

void KINFree(void *kinmem)
{
  KINMem kin_mem;

  if (kinmem == NULL) return;

  kin_mem = (KINMem) kinmem;
  KINFreeVectors(kin_mem);

  /* call lfree if non-NULL */

  if (lfree != NULL) lfree(kin_mem);

  free(kin_mem);
}

/*
 * -----------------------------------------------------------------
 * private helper functions
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : KINCheckNvector
 * -----------------------------------------------------------------
 * This routine checks if all required vector operations are
 * implemented (excluding those required by KINConstraint). If all
 * necessary operations are present, then KINCheckNvector returns
 * TRUE. Otherwise, FALSE is returned.
 * -----------------------------------------------------------------
 */

static booleantype KINCheckNvector(N_Vector tmpl)
{
  if ((tmpl->ops->nvclone     == NULL) ||
      (tmpl->ops->nvdestroy   == NULL) ||
      (tmpl->ops->nvlinearsum == NULL) ||
      (tmpl->ops->nvprod      == NULL) ||
      (tmpl->ops->nvdiv       == NULL) ||
      (tmpl->ops->nvscale     == NULL) ||
      (tmpl->ops->nvabs       == NULL) ||
      (tmpl->ops->nvinv       == NULL) ||
      (tmpl->ops->nvmaxnorm   == NULL) ||
      (tmpl->ops->nvmin       == NULL) ||
      (tmpl->ops->nvwl2norm   == NULL)) return(FALSE);
  else return(TRUE);
}

/*
 * -----------------------------------------------------------------
 * Function : KINAllocVectors
 * -----------------------------------------------------------------
 * This routine allocates the KINSol vectors. If all memory
 * allocations are successful, KINAllocVectors returns TRUE.
 * Otherwise all allocated memory is freed and KINAllocVectors
 * returns FALSE.
 * -----------------------------------------------------------------
 */

static booleantype KINAllocVectors(KINMem kin_mem, N_Vector tmpl)
{
  /* allocate unew, fval, pp, vtemp1 and vtemp2 */
  
  unew = N_VClone(tmpl);
  if (unew == NULL) return(FALSE);

  fval = N_VClone(tmpl);
  if (fval == NULL) {
    N_VDestroy(unew);
    return(FALSE);
  }

  pp = N_VClone(tmpl);
  if (pp == NULL) {
    N_VDestroy(unew);
    N_VDestroy(fval);
    return(FALSE);
  }

  vtemp1 = N_VClone(tmpl);
  if (vtemp1 == NULL) {
    N_VDestroy(unew);
    N_VDestroy(fval);
    N_VDestroy(pp);
    return(FALSE);
  }

  vtemp2 = N_VClone(tmpl);
  if (vtemp2 == NULL) {
    N_VDestroy(unew);
    N_VDestroy(fval);
    N_VDestroy(pp);
    N_VDestroy(vtemp1);
    return(FALSE);
  }

  liw = 5*liw1;
  lrw = 5*lrw1;

  return(TRUE);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSolInit
 * -----------------------------------------------------------------
 * KINSolInit initializes the problem for the specific input
 * received in this call to KINSol (which calls KINSolInit). All
 * problem specification inputs are checked for errors. If any error
 * occurs during initialization, it is reported to the file whose
 * file pointer is errfp.
 *
 * The possible return values for KINSolInit are:
 *   KIN_SUCCESS : indicates a normal initialization
 *
 *   KINS_ILL_INPUT : indicates that an input error has been found
 *
 *   KIN_INITIAL_GUESS_OK : indicates that the guess uu
 *                          satisfied the system func(uu) = 0
 *                          within the tolerances specified
 * -----------------------------------------------------------------
 */

static int KINSolInit(KINMem kin_mem)
{
  int ret;
  
  /* check for illegal input parameters */

  if (uu == NULL) {
    fprintf(errfp, MSG_UU_NULL);   
    return(KIN_ILL_INPUT);
  }

  if ((globalstrategy != KIN_INEXACT_NEWTON) && (globalstrategy != KIN_LINESEARCH)) {
    fprintf(errfp, MSG_BAD_GLSTRAT, globalstrategy, KIN_INEXACT_NEWTON, KIN_LINESEARCH);
    return(KIN_ILL_INPUT);
  }

  if (uscale == NULL)  {
    fprintf(errfp, MSG_BAD_USCALE);
    return(KIN_ILL_INPUT);
  }

  if (N_VMin(uscale) <= ZERO){
    fprintf(errfp, MSG_USCALE_NONPOSITIVE); 
    return(KIN_ILL_INPUT);
  }

  if (fscale == NULL)  {
    fprintf(errfp, MSG_BAD_FSCALE);
    return(KIN_ILL_INPUT);
  }

  if (N_VMin(fscale) <= ZERO){
    fprintf(errfp, MSG_FSCALE_NONPOSITIVE);
    return(KIN_ILL_INPUT);
  }

  /* set the constraints flag */

  if (constraints == NULL) 
    constraintsSet = FALSE;
  else {
    constraintsSet = TRUE;
    if ((constraints->ops->nvconstrmask  == NULL) ||
	(constraints->ops->nvminquotient == NULL)) {
      if (errfp != NULL) fprintf(errfp, MSG_BAD_NVECTOR);
      return(KIN_ILL_INPUT);
    }
  }

  /* check the initial guess uu against the constraints */

  if (constraintsSet) {
    if (!N_VConstrMask(constraints, uu, vtemp1)) {
      fprintf(errfp, MSG_INITIAL_CNSTRNT);
      KINFreeVectors(kin_mem);
      free(kin_mem);
      return(KIN_ILL_INPUT);
    }
  }
  
  /* all error checking is complete at this point */

  if (printfl > 0)
    KINPrintInfo(kin_mem, "KINSolInit", PRNT_TOL, scsteptol, fnormtol);

  /* calculate the default value for mxnewtstep (maximum Newton step) */

  if (mxnewtstep == ZERO)
    mxnewtstep = THOUSAND * N_VWL2Norm(uu, uscale);
  if (mxnewtstep < ONE) mxnewtstep = ONE;

  /* set up the coefficients for the eta calculation */

  callForcingTerm = (etaflag != KIN_ETACONSTANT);

  /* this value is always used for choice #1 */

  if (etaflag == KIN_ETACHOICE1) ealpha = (ONE + RSqrt(FIVE)) * HALF;

  /* initial value for eta set to 0.5 for other than the KIN_ETACONSTANT option */

  if (etaflag != KIN_ETACONSTANT) eta = HALF;

  /* initialize counters */

  nfe = nnilpre = nni = nbcf = nbktrk = 0;

  /* see if the system func(uu) = 0 is satisfied by the initial guess uu */

  if (KINInitialStop(kin_mem)) return(KIN_INITIAL_GUESS_OK);

  /* initialize the linear solver if linit != NULL */

  if (linit != NULL) {
    ret = linit(kin_mem);
    if (ret != 0) {
      fprintf(errfp, MSG_LINIT_FAIL);
      return(KIN_LINIT_FAIL);
    }
  }

  /* initialize the L2 (Euclidean) norms of f for the linear iteration steps */
  fnorm = N_VWL2Norm(fval, fscale);
  f1norm = HALF * fnorm * fnorm;

  if (printfl > 0)
    KINPrintInfo(kin_mem, "KINSolInit", PRNT_NNI, nni, nfe, fnorm);

  /* problem has now been successfully initialized */

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINConstraint
 * -----------------------------------------------------------------
 * This routine checks if the proposed solution vector uu + pp
 * violates any constraints. If a constraint is violated, then the
 * scalar stepmul is determined such that uu + stepmul * pp does
 * not violate any constraints.
 *
 * Note: This routine is called by the global strategy routines
 * KINLineSearch and KINInexactNewton.
 * -----------------------------------------------------------------
 */

static int KINConstraint(KINMem kin_mem)
{
  N_VLinearSum(ONE, uu, ONE, pp, vtemp1);

  /* if vtemp1[i] violates constraint[i] then vtemp2[i] = 1
     else vtemp2[i] = 0 (vtemp2 is the mask vector) */

  if(N_VConstrMask(constraints, vtemp1, vtemp2)) return(0);

  /* vtemp1[i] = ABS(pp[i]) */

  N_VAbs(pp, vtemp1);

  /* consider vtemp1[i] only if vtemp2[i] = 1 (constraint violated) */

  N_VProd(vtemp2, vtemp1, vtemp1);

  N_VAbs(uu, vtemp2);
  stepmul = POINT9 * N_VMinQuotient(vtemp2, vtemp1);

  return(1);
}

/*
 * -----------------------------------------------------------------
 * Function : KINForcingTerm
 * -----------------------------------------------------------------
 * This routine computes eta, the scaling factor in the linear
 * convergence stopping tolerance eps when choice #1 or choice #2
 * forcing terms are used. Eta is computed here for all but the
 * first iterative step, which is set to the default in routine
 * KINSolInit.
 *
 * This routine was written by Homer Walker of Utah State
 * University with subsequent modifications by Allan Taylor @ LLNL.
 *
 * It is based on the concepts of the paper 'Choosing the forcing
 * terms in an inexact Newton method', SIAM J Sci Comput, 17
 * (1996), pp 16 - 32, or Utah State University Research Report
 * 6/94/75 of the same title.
 * -----------------------------------------------------------------
 */

static void KINForcingTerm(KINMem kin_mem, realtype fnormp)
{
  realtype eta_max, eta_min, eta_safe, linmodel_norm;

  eta_max  = POINT9;
  eta_min  = POINTOHOHOHONE;
  eta_safe = HALF;

  /* choice #1 forcing term */

  if (etaflag == KIN_ETACHOICE1) {

    /* compute the norm of f + Jp , scaled L2 norm */

    linmodel_norm = RSqrt((fnorm * fnorm) + (TWO * sfdotJp) + (sJpnorm * sJpnorm));

    /* form the safeguarded for choice #1 */ 

    eta_safe = RPowerR(eta, ealpha); 
    eta = ABS(fnormp - linmodel_norm) / fnorm; 
  }

  /* choice #2 forcing term */

  if (etaflag == KIN_ETACHOICE2) {
    eta_safe = egamma * RPowerR(eta, ealpha); 
    eta = egamma * RPowerR((fnormp / fnorm), ealpha); 
  }

  /* apply safeguards */
 
  if(eta_safe < POINT1) eta_safe = ZERO;
  eta = MAX(eta, eta_safe); 
  eta = MAX(eta, eta_min); 
  eta = MIN(eta, eta_max); 

  return; 
}

/*
 * -----------------------------------------------------------------
 * Function : KINFreeVectors
 * -----------------------------------------------------------------
 * This routine frees the KINSol vectors allocated by
 * KINAllocVectors.
 * -----------------------------------------------------------------
 */

static void KINFreeVectors(KINMem kin_mem)
{
  N_VDestroy(unew);
  N_VDestroy(fval);
  N_VDestroy(pp);
  N_VDestroy(vtemp1);
  N_VDestroy(vtemp2);
}

/*
 * -----------------------------------------------------------------
 * Function : KINInitialStop
 * -----------------------------------------------------------------
 * This routine checks the initial guess (uu) to see if the
 * system func(uu) = 0 is satisfied to the level of 0.01 * fnormtol
 * -----------------------------------------------------------------
 */

static booleantype KINInitialStop(KINMem kin_mem)
{
  realtype fmax;

  func(uu, fval, f_data); nfe++;
  fmax = KINScFNorm(fval, fscale, vtemp1);
  if (printfl > 1)
    KINPrintInfo(kin_mem, "KINInitialStop", PRNT_FMAX, fmax);

  return(fmax <= (POINTOH1 * fnormtol));
}

/*
 * -----------------------------------------------------------------
 * Function : KINInexactNewton
 * -----------------------------------------------------------------
 * This routine is the main driver for the inexact Newton
 * algorithm. Its purpose is to compute unew = uu + pp in the
 * direction pp from uu, taking the full inexact Newton step. The
 * step may be constrained if the constraint conditions are
 * violated, or if the norm of pp is greater than mxnewtstep. 
 * -----------------------------------------------------------------
 */

static int KINInexactNewton(KINMem kin_mem, realtype *fnormp, realtype *f1normp,
                            booleantype *maxStepTaken)
{
  realtype pnorm, ratio;
  int ret;

  *maxStepTaken = FALSE;
  pnorm = N_VWL2Norm(pp, uscale);
  ratio = ONE;
  if (pnorm > mxnewtstep) {
    ratio = mxnewtstep / pnorm;
    N_VScale(ratio, pp, pp);
    pnorm = mxnewtstep;
  }

  if (printfl > 0)
    KINPrintInfo(kin_mem,"KINInexactNewton",PRNT_PNORM,pnorm);

  /* if constraints are active, then constrain the step accordingly */

  stepl = pnorm;
  stepmul = ONE;
  if (constraintsSet) {
    ret = KINConstraint(kin_mem);  /* Note: This routine changes the step pp. */
    if (ret == 1) {
      ratio *= stepmul;
      N_VScale(stepmul, pp, pp);
      pnorm *= stepmul;
      stepl = pnorm;
      if (printfl > 0)
        KINPrintInfo(kin_mem,"KINInexactNewton",PRNT_PNORM,pnorm);
      if (pnorm <= scsteptol) return(1);
    }
  }
 
  /* scale these two expressions by ratio for later use in KINForcingTerm */

  sfdotJp *= ratio;
  sJpnorm *= ratio;
 
  /* compute the iterate unew = uu + pp */

  N_VLinearSum(ONE, uu, ONE, pp, unew);

  /* evaluate func(unew) and its norm, and return */
 
  func(unew, fval, f_data); nfe++;
  *fnormp = N_VWL2Norm(fval,fscale);
  *f1normp = HALF * (*fnormp) * (*fnormp);
  if (printfl > 1) 
    KINPrintInfo(kin_mem,"KINInexactNewton",PRNT_FNORM,*fnormp);
  if (pnorm > (POINT99 * mxnewtstep)) *maxStepTaken = TRUE; 
  return(0);

}

/*
 * -----------------------------------------------------------------
 * Function : KINLineSearch
 * -----------------------------------------------------------------
 * The routine KINLineSearch implements the LineSearch algorithm.
 * Its purpose is to find unew = uu + rl * pp in the direction pp
 * from uu so that:
 *                                    t
 *  func(unew) <= func(uu) + alpha * g  (unew - uu) (alpha = 1.e-4)
 *
 *    and
 *                                   t
 *  func(unew) >= func(uu) + beta * g  (unew - uu) (beta = 0.9)
 *
 * where 0 < rlmin <= rl <= rlmax.
 *
 * Note:
 *             mxnewtstep
 *  rlmax = ----------------   if uu+pp does not violateany constraints
 *          ||uscale*pp||_L2
 *
 *  rlmax = 1   otherwise
 *
 *    and
 *
 *                 scsteptol
 *  rlmin = ------------------------
 *          ||         pp         ||
 *          || ------------------ ||_L-infinity
 *          || (uscale + ABS(uu)) ||
 *
 * -----------------------------------------------------------------
 */

static int KINLineSearch(KINMem kin_mem, realtype *fnormp, realtype *f1normp,
                         booleantype *maxStepTaken)
{
  realtype pnorm, ratio, slpi, rlmin, rlength, rl, rlmax, rldiff;
  realtype rltmp, rlprev, pt1trl, f1nprv, rllo, rlinc, alpha, beta;
  realtype alpha_cond, beta_cond, rl_a, tmp1, rl_b, tmp2, disc;
  long int nfesave, rladjust;
  int ret, ivio;

  rladjust = 0;
  *maxStepTaken = FALSE;
  ratio = ONE;
  alpha = POINTOHOHOHONE;
  beta = POINT9;

  pnorm = N_VWL2Norm(pp, uscale);
  if(pnorm > mxnewtstep ) {
    ratio = mxnewtstep / pnorm;
    N_VScale(ratio, pp, pp);
    pnorm = mxnewtstep;
  }

  ivio = 0;

  /* check if constraints are active and if so constrain the step by the constraints */

  stepl = pnorm;
  stepmul = ONE;

  if(constraintsSet){
    ret = KINConstraint(kin_mem);
    if(ret == 1){
      ivio = 1;
      ratio *= stepmul;
      N_VScale(stepmul, pp, pp);
      pnorm *= stepmul;
      stepl = pnorm;
      if (printfl > 0)
        KINPrintInfo(kin_mem,"KINLineSearch",PRNT_PNORM1,pnorm);
    }
  }

  slpi = sfdotJp * ratio;
  rlength = KINScSteplength(kin_mem, uu, pp, uscale);
  rlmin = scsteptol / rlength;
  rl = ONE;

  if (printfl > 2)
    KINPrintInfo(kin_mem,"KINLineSearch",PRNT_LAM,rlmin,f1norm,pnorm);

  /* now begin the iteration to find an rl value which satisfies both the
     alpha and beta conditions */

  /* if rl < rlmin, then terminate and return 1 */

  nfesave = nfe;

  loop {

    N_VLinearSum(ONE, uu, rl, pp, unew);
    func(unew, fval, f_data); nfe++;
    *fnormp = N_VWL2Norm(fval, fscale);
    *f1normp = HALF * (*fnormp) * (*fnormp) ;
    alpha_cond = f1norm + (alpha * slpi * rl);

    if (printfl > 2 && rladjust > 0)
      KINPrintInfo(kin_mem,"KINLinesearch",PRNT_ALPHA,*fnormp,*f1normp,alpha_cond,rl);

    if ((*f1normp) <= alpha_cond) break;

    /* alpha condition not satisfied so perform backtracking to compute a new rl value */

    if (rl < rlmin) {

      /* no satisfactory unew can be found sufficiently distinct from uu
         so copy uu into unew and return */

      /* step remains unchanged */

      N_VScale(ONE, uu, unew);
      func(unew, fval, f_data); nfe++;
      *fnormp = N_VWL2Norm(fval, fscale);
      *f1normp = HALF * (*fnormp) * (*fnormp);
      nbktrk = nfe - nfesave - 1;
      return(1);
    }

    /* use cubic fit for all remaining backtracks */

    if (rl != ONE) {
      tmp1 = (*f1normp) - f1norm - (rl * slpi);
      tmp2 = f1nprv - f1norm - (rlprev * slpi);
      rl_a = ((ONE / (rl * rl)) * tmp1) - ((ONE / (rlprev * rlprev)) * tmp2);
      rl_b = ((-rlprev / (rl * rl)) * tmp1) + ((rl / (rlprev * rlprev)) * tmp2);
      tmp1 = ONE / (rl - rlprev);
      rl_a *= tmp1;
      rl_b *= tmp1;
      disc = (rl_b * rl_b) - (THREE * rl_a * slpi);

      /* cubic is actually just a quadratic (rl_a ~ 0) */

      if (ABS(rl_a) < uround) rltmp = -slpi / (TWO * rl_b);

      /* real cubic */

      else rltmp = (-rl_b + RSqrt(disc)) / (THREE * rl_a);

      if (rltmp > (HALF * rl)) rltmp = HALF * rl;
    }

    /* use quadratic fit only for initial backtrack */

    else if (rl == ONE) rltmp = -slpi / (TWO * ((*f1normp) - f1norm - slpi));

    rlprev = rl;
    f1nprv = (*f1normp);
    pt1trl = POINT1 * rl;
    rl = MAX(pt1trl, rltmp);
    rladjust++;
  }

  /* alpha condition is satisfied so now check the beta condition */

  beta_cond = f1norm + (beta * slpi * rl);
  if ((*f1normp) < beta_cond) {

    if ((rl == ONE) && (pnorm < mxnewtstep)) {
      rlmax = mxnewtstep / pnorm;
      if (ivio == 1) rlmax = ONE;
      do {
        rlprev = rl;
        f1nprv = *f1normp;
        rl = MIN((TWO * rl), rlmax);
        rladjust++;
        N_VLinearSum(ONE, uu, rl, pp, unew);
        func(unew, fval, f_data); nfe++;
        *fnormp = N_VWL2Norm(fval, fscale);
        *f1normp = HALF * (*fnormp) * (*fnormp);
        alpha_cond = f1norm + (alpha * slpi * rl);
        beta_cond = f1norm + (beta * slpi * rl);
        if (printfl > 2)
          KINPrintInfo(kin_mem,"KINLineSearch",PRNT_BETA,*f1normp, beta_cond, rl);
      } while (((*f1normp) <= alpha_cond) && 
	       ((*f1normp) < beta_cond) && (rl < rlmax));
    }
    if ((rl < ONE) || ((rl > ONE) && (*f1normp > alpha_cond))) {
      rllo = MIN(rl, rlprev);
      rldiff = ABS(rlprev - rl);

      do {
        rlinc = HALF * rldiff;
        rl = rllo + rlinc;
        rladjust++;
        N_VLinearSum(ONE, uu, rl, pp, unew);
        func(unew, fval, f_data); nfe++;
        *fnormp = N_VWL2Norm(fval, fscale);
        *f1normp = HALF * (*fnormp) * (*fnormp);
        alpha_cond = f1norm + (alpha * slpi * rl);
        beta_cond = f1norm + (beta * slpi * rl);
        if ((printfl > 2) && (rladjust > 0))
          KINPrintInfo(kin_mem, "KINLineSearch", PRNT_ALPHABETA, 
                       *f1normp, alpha_cond, beta_cond, rl);
        if ((*f1normp) > alpha_cond) rldiff = rlinc;
        else if (*f1normp < beta_cond) {
          rllo = rl;
          rldiff = rldiff - rlinc;
        }

      } while ((*f1normp > alpha_cond) ||
	       ((*f1normp < beta_cond) && (rldiff >= rlmin)));

      if ((*f1normp) < beta_cond) {

	/* beta condition could not be satisfied so set unew to last u value
	   that satisfied the alpha condition and continue */

	/* increment counter on number of steps beta condition not satisfied */

        N_VLinearSum(ONE, uu, rllo, pp, unew);
        func(unew, fval, f_data); nfe++;
        *fnormp = N_VWL2Norm(fval, fscale);
        *f1normp = HALF * (*fnormp) * (*fnormp);   
        nbcf++;
      }
    }  /* end of if (rl < ONE) block */
  }  /* end of (f1normp < beta_cond) loop */

  nbktrk= nfe - nfesave;
  if ((printfl > 1) && (rladjust > 0))
    KINPrintInfo(kin_mem, "KINLineSearch", PRNT_ADJ, rladjust);

  /* scale these two expressions by rl * ratio for later use in KINForcingTerm */

  sfdotJp = sfdotJp * rl * ratio;
  sJpnorm = sJpnorm * rl * ratio;

  if ((rl * pnorm) > (POINT99 * mxnewtstep)) *maxStepTaken = TRUE;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : KINLinSolvDrv
 * -----------------------------------------------------------------
 * This routine handles the process of solving for the approximate
 * solution xx of the Newton equations in the Newton iteration.
 * Subsequent routines handle the nonlinear aspects of its
 * application. It takes as arguments a KINSol memory pointer,
 * an N_Vector xx, which is used as scratch and when returned
 * contains the solution x!. It also takes an N_Vector bb which
 * is used as scratch but represents the rhs of the system Jx = b
 * whose solution xx is being approximated.
 * -----------------------------------------------------------------
 */

static int KINLinSolDrv(KINMem kin_mem , N_Vector bb , N_Vector xx )
{
  int ret;

  if ((nni - nnilpre) >= msbpre) pthrsh = TWO;

  loop{

    precondcurrent = FALSE;

    if ((pthrsh > ONEPT5) && setupNonNull) {
      ret = lsetup(kin_mem);  /* call pset routine if lsetup != NULL */
      precondcurrent = TRUE;
      nnilpre = nni;
      if (ret != 0) return(KIN_LSETUP_FAIL);
    }

  /* load bb with the current value of -fval */

  N_VScale(-ONE, fval, bb);

  /* call the generic 'lsolve' routine to solve the system Jx = b */

  ret = lsolve(kin_mem, xx, bb, &res_norm);

  if (ret == 0) return(0);
  else if (ret < 0) return(KIN_LSOLVE_FAIL);
  else if ((!setupNonNull) || (precondcurrent)) return(KIN_LINSOLV_NO_RECOVERY);

  /* loop back only if the preconditioner is in use and is not current */

  pthrsh = TWO;

  }
}

/*
 * -----------------------------------------------------------------
 * Function : KINScFNorm
 * -----------------------------------------------------------------
 * This routine computes the max norm for scaled vectors. The
 * scaling vector is scale, and the vector of which the norm is to
 * be determined is vv. The returned value, fnormval, is the
 * resulting scaled vector norm. This routine uses N_Vector
 * functions from the vector module.
 * -----------------------------------------------------------------
 */

static realtype KINScFNorm(N_Vector vv, N_Vector scale, N_Vector wrkv)
{
  N_VAbs(vv, wrkv);
  N_VProd(scale, wrkv, wrkv);
  return(N_VMaxNorm(wrkv));
}

/*
 * -----------------------------------------------------------------
 * Function : KINScStepLength
 * -----------------------------------------------------------------
 * This routine computes the max norm of the scaled steplength, ss.
 * Here ucur is the current step and usc is the u scale factor.
 * -----------------------------------------------------------------
 */

static realtype KINScSteplength(KINMem kin_mem, N_Vector ucur,
                                N_Vector ss, N_Vector usc)
{
  N_VInv(usc, vtemp1);
  N_VAbs(ucur, vtemp2);
  N_VLinearSum(ONE, vtemp1, ONE, vtemp2, vtemp1);
  N_VDiv(ss, vtemp1, vtemp1);
  return(N_VMaxNorm(vtemp1));
}

/*
 * -----------------------------------------------------------------
 * Function : KINStop
 * -----------------------------------------------------------------
 * This routine checks the current iterate unew to see if the
 * system func(unew) = 0 is satisfied by a variety of tests.
 * -----------------------------------------------------------------
 */

static int KINStop(KINMem kin_mem, booleantype maxStepTaken, int globalstratret)
{
  realtype fmax, rlength;

  if (globalstratret == 1){

    if (setupNonNull && !precondcurrent) {

      /* if the globalstratret was caused (potentially) by the 
	 preconditioner being out of date, then update the
	 preconditioner */

      pthrsh = TWO;

      return(CONTINUE_ITERATIONS);
    }

    else return( (globalstrategy == KIN_INEXACT_NEWTON) ?
		 KIN_STEP_LT_STPTOL : KIN_LINESEARCH_NONCONV );
  }

  /* check tolerance on scaled norm of func at the current iterate */

  fmax = KINScFNorm(fval, fscale, vtemp1);

  if (printfl > 1) 
    KINPrintInfo(kin_mem, "KINStop", PRNT_FMAX, fmax);

  if (fmax <= fnormtol) return(KIN_SUCCESS);

  /* check if the scaled distance between the last two steps is too small */

  N_VLinearSum(ONE, unew, -ONE, uu, vtemp1);
  rlength = KINScSteplength(kin_mem, unew, vtemp1, uscale);
  if (rlength <= scsteptol) {

    if (!precondcurrent) {

      /* for rlength too small and the preconditioner not current,
         try again with the preconditioner current */

      pthrsh = TWO;

      return(CONTINUE_ITERATIONS);
    }

    else

      /* otherwise return failure flag */

      return(KIN_STEP_LT_STPTOL);
  }

  /* if the maximum number of iterations is reached, return failure flag */

   if (nni >= mxiter) return(KIN_MAXITER_REACHED);

   /* check for consecutive number of steps taken of size mxnewtstep
      and if not maxStepTaken, then set ncscmx to 0 */
 
  if (maxStepTaken) ncscmx++;
  else ncscmx = 0;
 
  if (ncscmx == 5) return(KIN_MXNEWT_5X_EXCEEDED);

  /* load threshold for re-evaluating the preconditioner */

  pthrsh = rlength;

  /* if made it to here, then the iteration process is not finished
     so return CONTINUE_ITERATIONS flag */

  return(CONTINUE_ITERATIONS);
}


/*
 * -----------------------------------------------------------------
 * KINPrintInfo
 * -----------------------------------------------------------------
 */

static void KINPrintInfo(KINMem kin_mem, char *funcname, int key,...)
{
  va_list ap;
  int ret;

  fprintf(infofp, "---%s\n   ", funcname);

  /* initialize argument processing */
  va_start(ap, key); 

  switch(key) {

  case PRNT_RETVAL:
    ret = va_arg(ap,int);
    fprintf(infofp,"return value: %d ",ret);
    switch(ret) {
    case KIN_SUCCESS:
      fprintf(infofp, "(KIN_SUCCESS)\n");
      break;
    case KIN_STEP_LT_STPTOL:
      fprintf(infofp, "(KIN_STEP_LT_STPTOL)\n");
      break;
    case KIN_LINESEARCH_NONCONV:
      fprintf(infofp, "(KIN_LINESEARCH_NONCONV)\n");
      break;
    case KIN_LINESEARCH_BCFAIL:
      fprintf(infofp, "(KIN_LINESEARCH_BCFAIL)\n");
      break;
    case KIN_MAXITER_REACHED:
      fprintf(infofp, "(KIN_MAXITER_REACHED)\n");
      break;
    case KIN_MXNEWT_5X_EXCEEDED:
      fprintf(infofp, "(KIN_MXNEWT_5X_EXCEEDED)\n");
      break;
    case KIN_LINSOLV_NO_RECOVERY:
      fprintf(infofp, "(KIN_LINSOLV_NO_RECOVERY)\n");
      break;
    case KIN_LSETUP_FAIL:
      fprintf(infofp, "(KIN_PRECONDSET_FAILURE)\n");
      break;
    case KIN_LSOLVE_FAIL:
      fprintf(infofp, "(KIN_PRECONDSOLVE_FAILURE)\n");
      break;
    }
    break;

  case PRNT_NNI:
    fprintf(infofp, "nni = %4ld  ", va_arg(ap,long int));
    fprintf(infofp, "nfe = %6ld  ", va_arg(ap,long int));
          
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "fnorm = %26.16Lg\n", va_arg(ap,realtype));
#else
    fprintf(infofp, "fnorm = %26.16g\n", va_arg(ap,realtype));
#endif
    break;

  case PRNT_TOL:
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "scsteptol = %12.3Lg  fnormtol = %12.3Lg\n", 
            va_arg(ap,realtype), va_arg(ap,realtype));
#else
    fprintf(infofp, "scsteptol = %12.3g  fnormtol = %12.3g\n", 
            va_arg(ap,realtype), va_arg(ap,realtype));
#endif
    break;
    
  case PRNT_FMAX:
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "scaled f norm (for stopping) = %12.3Lg\n",
            va_arg(ap,realtype));
#else
    fprintf(infofp, "scaled f norm (for stopping) = %12.3g\n",
            va_arg(ap,realtype));
#endif
    break;
    
  case PRNT_PNORM:
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "pnorm = %12.4Le\n", va_arg(ap,realtype));
#else
    fprintf(infofp, "pnorm = %12.4e\n", va_arg(ap,realtype));
#endif
    break;

  case PRNT_PNORM1:
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "(ivio=1) pnorm = %12.4Le\n", va_arg(ap,realtype));
#else
    fprintf(infofp, "(ivio=1) pnorm = %12.4e\n", va_arg(ap,realtype));
#endif
    break;
     
  case PRNT_FNORM:
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "fnorm(L2) = %20.8Le\n", va_arg(ap,realtype));
#else
    fprintf(infofp, "fnorm(L2) = %20.8e\n", va_arg(ap,realtype));
#endif
    break;

  case PRNT_LAM:
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "min_lam = %11.4Le  ", va_arg(ap,realtype));
    fprintf(infofp, "f1norm = %11.4Le  ", va_arg(ap,realtype));
    fprintf(infofp, "pnorm = %11.4Le\n", va_arg(ap,realtype));
#else
    fprintf(infofp, "min_lam = %11.4e  ", va_arg(ap,realtype));
    fprintf(infofp, "f1norm = %11.4e  ", va_arg(ap,realtype));
    fprintf(infofp, "pnorm = %11.4e\n", va_arg(ap,realtype));
#endif
    break;

  case PRNT_ALPHA:
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "fnorm = %15.8Le  ", va_arg(ap,realtype));
    fprintf(infofp, "f1norm = %15.8Le  ", va_arg(ap,realtype));
    fprintf(infofp, "alpha_cond = %15.8Le  ", va_arg(ap,realtype));
    fprintf(infofp, "pnorm = %15.8Le\n", va_arg(ap,realtype));
#else
    fprintf(infofp, "fnorm = %15.8e  ", va_arg(ap,realtype));
    fprintf(infofp, "f1norm = %15.8e  ", va_arg(ap,realtype));
    fprintf(infofp, "alpha_cond = %15.8e  ", va_arg(ap,realtype));
    fprintf(infofp, "pnorm = %15.8e\n", va_arg(ap,realtype));
#endif
    break;

  case PRNT_BETA:
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "f1norm = %15.8Le  ", va_arg(ap,realtype));
    fprintf(infofp, "beta_cond = %15.8Le  ", va_arg(ap,realtype));
    fprintf(infofp, "lam = %15.8Le\n", va_arg(ap,realtype));
#else
    fprintf(infofp, "f1norm = %15.8e  ", va_arg(ap,realtype));
    fprintf(infofp, "beta_cond = %15.8e  ", va_arg(ap,realtype));
    fprintf(infofp, "lam = %15.8e\n", va_arg(ap,realtype));
#endif
    break;

  case PRNT_ALPHABETA:
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(infofp, "f1norm = %15.8Le  ", va_arg(ap,realtype));
    fprintf(infofp, "alpha_cond = %15.8Le  ", va_arg(ap,realtype));
    fprintf(infofp, "beta_cond = %15.8Le  ", va_arg(ap,realtype));
    fprintf(infofp, "lam = %15.8Le\n", va_arg(ap,realtype));
#else
    fprintf(infofp, "f1norm = %15.8e  ", va_arg(ap,realtype));
    fprintf(infofp, "alpha_cond = %15.8e  ", va_arg(ap,realtype));
    fprintf(infofp, "beta_cond = %15.8e  ", va_arg(ap,realtype));
    fprintf(infofp, "lam = %15.8e\n", va_arg(ap,realtype));
#endif
    break;

  case PRNT_ADJ:
    fprintf(infofp, "no. of lambda adjustments = %ld\n", va_arg(ap,long int));
    break;

  }

  /* finalize argument processing */
  va_end(ap);

}
