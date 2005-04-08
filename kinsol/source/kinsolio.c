/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2005-04-08 15:08:05 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the optional input and output
 * functions for the KINSOL solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "kinsol_impl.h"
#include "sundialstypes.h"
#include "sundialsmath.h"

#define ZERO      RCONST(0.0)
#define POINT1    RCONST(0.1)
#define ONETHIRD  RCONST(.3333333333333333)
#define HALF      RCONST(0.5)
#define TWOTHIRDS RCONST(.6666666666666667)
#define POINT9    RCONST(0.9)
#define ONE       RCONST(1.0)
#define TWO       RCONST(2.0)
#define TWOPT5    RCONST(2.5)

/* 
 * =================================================================
 * KINSOL optional input functions
 * =================================================================
 */

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

  if (mxiter < 0) {
    fprintf(errfp, MSG_BAD_MXITER);
    return(KIN_ILL_INPUT);
  }

  if (mxiter == 0)
    kin_mem->kin_mxiter = MXITER_DEFAULT;
  else
    kin_mem->kin_mxiter = mxiter;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetNoInitSetup
 * -----------------------------------------------------------------
 */

int KINSetNoInitSetup(void *kinmem, booleantype noInitSetup)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;
  kin_mem->kin_noInitSetup = noInitSetup;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxSetupCalls
 * -----------------------------------------------------------------
 */

int KINSetMaxSetupCalls(void *kinmem, long int msbset)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (msbset < 0) {
    fprintf(errfp, MSG_BAD_MSBSET);
    return(KIN_ILL_INPUT);
  }
  
  if (msbset == 0)
    kin_mem->kin_msbset = MSBSET_DEFAULT;
  else
    kin_mem->kin_msbset = msbset;

  return(KIN_SUCCESS);
}

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

  if ((eta < ZERO) || (eta > ONE)) {
    fprintf(errfp, MSG_BAD_ETACONST);
    return(KIN_ILL_INPUT);
  }

  if (eta == ZERO)
    kin_mem->kin_eta = POINT1;
  else
    kin_mem->kin_eta = eta;


  return(KIN_SUCCESS);
}

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

  if (mxnewtstep < ZERO) {
    fprintf(errfp, MSG_BAD_MXNEWTSTEP);
    return(KIN_ILL_INPUT);
  }

  /* Note: passing a value of 0.0 will use the default
     value (computed in KINSolinit) */

  kin_mem->kin_mxnewtstep = mxnewtstep;

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxBetaFails
 * -----------------------------------------------------------------
 */

int KINSetMaxBetaFails(void *kinmem, long int mxnbcf)
{
  KINMem kin_mem;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (mxnbcf < 0) {
    fprintf(errfp, MSG_BAD_MXNBCF);
    return(KIN_ILL_INPUT);
  }

  if (mxnbcf == 0)
    kin_mem->kin_mxnewtstep = MXNBCF_DEFAULT;
  else
    kin_mem->kin_mxnewtstep = mxnbcf;

  return(KIN_SUCCESS);

}

/*
 * -----------------------------------------------------------------
 * Function : KINSetRelErrFunc
 * -----------------------------------------------------------------
 */

int KINSetRelErrFunc(void *kinmem, realtype relfunc)
{
  KINMem kin_mem;
  realtype uround;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (relfunc < ZERO) {
    fprintf(errfp, MSG_BAD_RELFUNC);
    return(KIN_ILL_INPUT);
  }

  if (relfunc == ZERO) {
    uround = kin_mem->kin_uround;
    kin_mem->kin_sqrt_relfunc = RSqrt(uround);
  } else {
    kin_mem->kin_sqrt_relfunc = RSqrt(relfunc);
  }

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetFuncNormTol
 * -----------------------------------------------------------------
 */

int KINSetFuncNormTol(void *kinmem, realtype fnormtol)
{
  KINMem kin_mem;
  realtype uround;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (fnormtol < ZERO) {
    fprintf(errfp, MSG_BAD_FNORMTOL);
    return(KIN_ILL_INPUT);
  }

  if (fnormtol == ZERO) {
    uround = kin_mem->kin_uround;
    kin_mem->kin_fnormtol = RPowerR(uround,ONETHIRD);
  } else {
    kin_mem->kin_fnormtol = fnormtol;
  }

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetScaledStepTol
 * -----------------------------------------------------------------
 */

int KINSetScaledStepTol(void *kinmem, realtype scsteptol)
{
  KINMem kin_mem;
  realtype uround;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (scsteptol < ZERO) {
    fprintf(errfp, MSG_BAD_SCSTEPTOL);
    return(KIN_ILL_INPUT);
  }

  if (scsteptol == ZERO) {
    uround = kin_mem->kin_uround;
    kin_mem->kin_scsteptol = RPowerR(uround,TWOTHIRDS);
  } else {
    kin_mem->kin_scsteptol = scsteptol;
  }

  return(KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetConstraints
 * -----------------------------------------------------------------
 */

int KINSetConstraints(void *kinmem, N_Vector constraints)
{
  KINMem kin_mem;
  realtype temptest;

  if (kinmem == NULL) {
    fprintf(stderr, MSG_KINS_NO_MEM);
    return(KIN_MEM_NULL);
  }

  kin_mem = (KINMem) kinmem;

  if (constraints == NULL) {
    if (kin_mem->kin_constraintsSet) {
      N_VDestroy(kin_mem->kin_constraints);
    }
    kin_mem->kin_constraintsSet = FALSE;
    return(KIN_SUCCESS);
  }

  /*  Check the constraints vector */

  temptest = N_VMaxNorm(constraints);
  if((temptest > TWOPT5) || (temptest < HALF)){ 
    if(errfp!=NULL) fprintf(errfp, MSG_BAD_CONSTRAINTS); 
    return(KIN_ILL_INPUT); 
  }

  if (!kin_mem->kin_constraintsSet) {
    kin_mem->kin_constraints = N_VClone(constraints);
    kin_mem->kin_constraintsSet = TRUE;
  }

  /* Load the constraint vector */

  N_VScale(ONE, constraints, kin_mem->kin_constraints);

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

/* 
 * =================================================================
 * Readability constants
 * =================================================================
 */

#define nni (kin_mem->kin_nni)
#define nfe (kin_mem->kin_nfe)
#define nbcf (kin_mem->kin_nbcf)  
#define nbktrk (kin_mem->kin_nbktrk)
#define stepl (kin_mem->kin_stepl)
#define fnorm (kin_mem->kin_fnorm)
#define liw (kin_mem->kin_liw)
#define lrw (kin_mem->kin_lrw)

/* 
 * =================================================================
 * KINSOL optional input functions
 * =================================================================
 */

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

