/*******************************************************************
 *                                                                 *
 * File          : cvsspgmr.c                                      *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh and           *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 11 July 2003                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/cvodes/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the implementation file for the CVODES scaled,          *
 * preconditioned GMRES linear solver, CVSSPGMR.                   *
 *                                                                 *
 *******************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "cvsspgmr.h"
#include "cvodes.h"
#include "sundialstypes.h"
#include "nvector.h"
#include "sundialsmath.h"
#include "iterativ.h"
#include "spgmr.h"


/* Error Messages */

#define CVSPGMR               "CVSpgmr-- "

#define MSG_CVMEM_NULL        CVSPGMR "CVode Memory is NULL.\n\n"

#define MSG_MEM_FAIL          CVSPGMR "A memory request failed.\n\n"

#define MSG_BAD_PRETYPE1      CVSPGMR "pretype=%d illegal.\n"
#define MSG_BAD_PRETYPE2      "The legal values are NONE=%d, LEFT=%d, "
#define MSG_BAD_PRETYPE3      "RIGHT=%d, and BOTH=%d.\n\n"
#define MSG_BAD_PRETYPE       MSG_BAD_PRETYPE1 MSG_BAD_PRETYPE2 MSG_BAD_PRETYPE3

#define MSG_PSOLVE_REQ        CVSPGMR "pretype!=NONE, but PSOLVE=NULL is illegal. \n\n"

#define MSG_SETGET_CVMEM_NULL "CVSpgmrSet*/CVSpgmrGet*-- cvodes memory is NULL. \n\n"

#define MSG_SETGET_LMEM_NULL  "CVSpgmrSet*/CVSpgmrGet*-- cvsspgmr memory is NULL. \n\n"

#define MSG_CVS_BAD_PRETYPE1  "CVSpgmrSetPrecType-- pretype=%d illegal.\n"
#define MSG_CVS_BAD_PRETYPE2  "The legal values are NONE=%d, LEFT=%d, "
#define MSG_CVS_BAD_PRETYPE3  "RIGHT=%d, and BOTH=%d.\n\n"
#define MSG_CVS_BAD_PRETYPE   MSG_CVS_BAD_PRETYPE1 MSG_CVS_BAD_PRETYPE2 MSG_CVS_BAD_PRETYPE3

#define MSG_CVS_BAD_GSTYPE1   "CVSpgmrSetGSType-- gstype=%d illegal.\n"
#define MSG_CVS_BAD_GSTYPE2   "The legal values are MODIFIED_GS=%d and "
#define MSG_CVS_BAD_GSTYPE3   "CLASSICAL_GS=%d.\n\n"
#define MSG_CVS_BAD_GSTYPE    MSG_CVS_BAD_GSTYPE1 MSG_CVS_BAD_GSTYPE2 MSG_CVS_BAD_GSTYPE3

#define MSG_CVS_BAD_DELT      "CVSpgmrSetDelt-- delt<0.0 illegal. \n\n"

/* Other Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* CVSPGMR linit, lsetup, lsolve, lsolveS, and lfree routines */

static int CVSpgmrInit(CVodeMem cv_mem);

static int CVSpgmrSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                        N_Vector vtemp2, N_Vector vtemp3);

static int CVSpgmrSolve(CVodeMem cv_mem, N_Vector b, N_Vector ycur,
                        N_Vector fcur);

static int CVSpgmrSolveS(CVodeMem cv_mem, N_Vector b, N_Vector ycur,
                         N_Vector fcur, integertype is);

static void CVSpgmrFree(CVodeMem cv_mem);

/* CVSPGMR Atimes and PSolve routines called by generic SPGMR solver */

static int CVSpgmrAtimes(void *cv_mem, N_Vector v, N_Vector z);

static int CVSpgmrPSolve(void *cv_mem, N_Vector r, N_Vector z, int lr);

/* CVSPGMR difference quotient routine for J*v */

static int CVSpgmrDQJtimes(N_Vector v, N_Vector Jv, realtype t,
                           N_Vector y, N_Vector fy, void *jac_data,
                           N_Vector work);

/**********************************************************************/

/* Readability Replacements */

#define lrw1    (cv_mem->cv_lrw1)
#define liw1    (cv_mem->cv_liw1)
#define uround  (cv_mem->cv_uround)
#define tq      (cv_mem->cv_tq)
#define nst     (cv_mem->cv_nst)
#define tn      (cv_mem->cv_tn)
#define h       (cv_mem->cv_h)
#define gamma   (cv_mem->cv_gamma)
#define gammap  (cv_mem->cv_gammap)   
#define nfe     (cv_mem->cv_nfe)
#define f       (cv_mem->cv_f)
#define f_data  (cv_mem->cv_f_data)
#define ewt     (cv_mem->cv_ewt)
#define errfp   (cv_mem->cv_errfp)
#define mnewt   (cv_mem->cv_mnewt)
#define ropt    (cv_mem->cv_ropt)
#define linit   (cv_mem->cv_linit)
#define lsetup  (cv_mem->cv_lsetup)
#define lsolve  (cv_mem->cv_lsolve)
#define lsolveS (cv_mem->cv_lsolveS)
#define lfree   (cv_mem->cv_lfree)
#define lmem    (cv_mem->cv_lmem)
#define nvspec  (cv_mem->cv_nvspec)
#define setupNonNull (cv_mem->cv_setupNonNull)

#define ewtS    (cv_mem->cv_ewtS)

#define sqrtN   (cvspgmr_mem->g_sqrtN)   
#define ytemp   (cvspgmr_mem->g_ytemp)
#define x       (cvspgmr_mem->g_x)
#define ycur    (cvspgmr_mem->g_ycur)
#define fcur    (cvspgmr_mem->g_fcur)
#define delta   (cvspgmr_mem->g_delta)
#define deltar  (cvspgmr_mem->g_deltar)
#define npe     (cvspgmr_mem->g_npe)
#define nli     (cvspgmr_mem->g_nli)
#define nps     (cvspgmr_mem->g_nps)
#define ncfl    (cvspgmr_mem->g_ncfl)
#define nstlpre (cvspgmr_mem->g_nstlpre)
#define njtimes (cvspgmr_mem->g_njtimes)
#define nfeSG   (cvspgmr_mem->g_nfeSG)
#define spgmr_mem (cvspgmr_mem->g_spgmr_mem)


/*************** CVSpgmr *********************************************

 This routine initializes the memory record and sets various function
 fields specific to the Spgmr linear solver module. CVSpgmr first
 calls the existing lfree routine if this is not NULL.  It then sets
 the cv_linit, cv_lsetup, cv_lsolve, cv_lsolveS, cv_lfree fields in (*cvode_mem)
 to be CVSpgmrInit, CVSpgmrSetup, CVSpgmrSolve, CVSpgmrSolveS, and CVSpgmrFree,
 respectively.  It allocates memory for a structure of type
 CVSpgmrMemRec and sets the cv_lmem field in (*cvode_mem) to the
 address of this structure.  It sets setupNonNull in (*cvode_mem),
 and sets the following fields in the CVSpgmrMemRec structure:
   g_pretype = pretype                                       
   g_gstype  = gstype                                       
   g_maxl    = MIN(N,CVSPGMR_MAXL)  if maxl <= 0             
             = maxl                 if maxl > 0              
   g_delt    = CVSPGMR_DELT if delt == 0.0                     
             = delt         if delt != 0.0                     
   g_P_data  = P_data                                        
   g_pset    = pset                                       
   g_psolve  = psolve                                        
   g_jtimes  = input parameter jtimes  if jtimes != NULL
             = CVSpgmrDQJtimes         otherwise
   g_j_data  = input parameter jac_data
 Finally, CVSpgmr allocates memory for ytemp and x, and calls
 SpgmrMalloc to allocate memory for the Spgmr solver.  The CVSpgmr
 return value is SUCCESS = 0, LMEM_FAIL = -1, or LIN_ILL_INPUT = -2.

**********************************************************************/

int CVSpgmr(void *cvode_mem, int pretype, int maxl)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;
  int mxl;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lfree != NULL) lfree(cv_mem);

  /* Set five main function fields in cv_mem */
  linit  = CVSpgmrInit;
  lsetup = CVSpgmrSetup;
  lsolve = CVSpgmrSolve;
  lsolveS= CVSpgmrSolveS;
  lfree  = CVSpgmrFree;

  /* Get memory for CVSpgmrMemRec */
  cvspgmr_mem = (CVSpgmrMem) malloc(sizeof(CVSpgmrMemRec));
  if (cvspgmr_mem == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LMEM_FAIL);
  }

  /* Set Spgmr parameters that have been passed in call sequence */
  cvspgmr_mem->g_pretype    = pretype;
  mxl = cvspgmr_mem->g_maxl = (maxl <= 0) ? CVSPGMR_MAXL : maxl;

  /* Set default values for the rest of the Spgmr parameters */
  cvspgmr_mem->g_gstype     = MODIFIED_GS;
  cvspgmr_mem->g_delt       = CVSPGMR_DELT;
  cvspgmr_mem->g_P_data     = NULL;
  cvspgmr_mem->g_pset       = NULL;
  cvspgmr_mem->g_psolve     = NULL;
  cvspgmr_mem->g_jtimes     = CVSpgmrDQJtimes;
  cvspgmr_mem->g_j_data     = cvode_mem;

  setupNonNull = FALSE;

  /* Check for legal pretype */ 
  if ((pretype != NONE) && (pretype != LEFT) &&
      (pretype != RIGHT) && (pretype != BOTH)) {
    fprintf(errfp, MSG_BAD_PRETYPE, pretype, NONE, LEFT, RIGHT, BOTH);
    return(LIN_ILL_INPUT);
  }

  /* Allocate memory for ytemp and x */
  ytemp = N_VNew(nvspec);
  if (ytemp == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    return(LMEM_FAIL);
  }
  x = N_VNew(nvspec);
  if (x == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    N_VFree(ytemp);
    return(LMEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, ytemp);
  sqrtN = RSqrt( N_VDotProd(ytemp, ytemp) );

  /* Call SpgmrMalloc to allocate workspace for Spgmr */
  spgmr_mem = SpgmrMalloc(mxl, nvspec);
  if (spgmr_mem == NULL) {
    fprintf(errfp, MSG_MEM_FAIL);
    N_VFree(ytemp);
    N_VFree(x);
    return(LMEM_FAIL);
  }
  
  /* Attach linear solver memory to CVODES memory */
  lmem = cvspgmr_mem;

  return(SUCCESS);
}

/*************** CVSpgmrSetPrecType **********************************/

int CVSpgmrSetPrecType(void *cvode_mem, int pretype)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  /* Check for legal pretype */ 
  if ((pretype != NONE) && (pretype != LEFT) &&
      (pretype != RIGHT) && (pretype != BOTH)) {
    fprintf(errfp, MSG_CVS_BAD_PRETYPE, pretype, NONE, LEFT, RIGHT, BOTH);
    return(LIN_ILL_INPUT);
  }

  cvspgmr_mem->g_pretype = pretype;

  return(SUCCESS);
}

/*************** CVSpgmrSetGSType ************************************/

int CVSpgmrSetGSType(void *cvode_mem, int gstype)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  /* Check for legal gstype */
  if ((gstype != MODIFIED_GS) && (gstype != CLASSICAL_GS)) {
    fprintf(errfp, MSG_CVS_BAD_GSTYPE, gstype, MODIFIED_GS, CLASSICAL_GS);
    return(LIN_ILL_INPUT);
  }

  cvspgmr_mem->g_gstype = gstype;

  return(SUCCESS);
}

/*************** CVSpgmrSetDelt **************************************/

int CVSpgmrSetDelt(void *cvode_mem, realtype delt)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  /* Check for legal delt */
  if(delt < ZERO) {
    fprintf(errfp, MSG_CVS_BAD_DELT);
    return(LIN_ILL_INPUT);
  }

  cvspgmr_mem->g_delt = (delt == ZERO) ? CVSPGMR_DELT : delt;

  return(SUCCESS);
}

/*************** CVSpgmrSetPrecSetupFn *******************************/

int CVSpgmrSetPrecSetupFn(void *cvode_mem, CVSpgmrPrecSetupFn pset)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  cvspgmr_mem->g_pset = pset;

  return(SUCCESS);
}

/*************** CVSpgmrSetPrecSolveFn *******************************/

int CVSpgmrSetPrecSolveFn(void *cvode_mem, CVSpgmrPrecSolveFn psolve)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  cvspgmr_mem->g_psolve = psolve;

  return(SUCCESS);
}

/*************** CVSpgmrSetPrecData **********************************/

int CVSpgmrSetPrecData(void *cvode_mem, void *P_data)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  cvspgmr_mem->g_P_data = P_data;

  return(SUCCESS);
}

/*************** CVSpgmrSetJacTimesVecFn *****************************/

int CVSpgmrSetJacTimesVecFn(void *cvode_mem, CVSpgmrJacTimesVecFn jtimes)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  cvspgmr_mem->g_jtimes = jtimes;

  return(SUCCESS);
}

/*************** CVSpgmrSetJacData ***********************************/

int CVSpgmrSetJacData(void *cvode_mem, void *jac_data)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  cvspgmr_mem->g_j_data = jac_data;

  return(SUCCESS);
}


/*************** CVSpgmrGetIntWorkSpace ******************************/

int CVSpgmrGetIntWorkSpace(void *cvode_mem, long int *leniwSG)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;
  int maxl;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  maxl = cvspgmr_mem->g_maxl;
  *leniwSG = liw1*(maxl + 5);

  return(OKAY);
}

/*************** CVSpgmrGetRealWorkSpace *****************************/

int CVSpgmrGetRealWorkSpace(void *cvode_mem, long int *lenrwSG)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;
  int maxl;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  maxl = cvspgmr_mem->g_maxl;
  *lenrwSG = lrw1*(maxl + 5) + maxl*(maxl + 4) + 1;

  return(OKAY);
}

/*************** CVSpgmrGetNumPrecEvals ******************************/

int CVSpgmrGetNumPrecEvals(void *cvode_mem, int *npevals)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  *npevals = npe;

  return(OKAY);
}

/*************** CVSpgmrGetNumPrecSolves *****************************/

int CVSpgmrGetNumPrecSolves(void *cvode_mem, int *npsolves)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  *npsolves = nps;

  return(OKAY);
}

/*************** CVSpgmrGetNumLinIters *******************************/

int CVSpgmrGetNumLinIters(void *cvode_mem, int *nliters)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  *nliters = nli;

  return(OKAY);
}

/*************** CVSpgmrGetNumConvFails ******************************/

int CVSpgmrGetNumConvFails(void *cvode_mem, int *nlcfails)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  *nlcfails = ncfl;

  return(OKAY);
}

/*************** CVSpgmrGetNumJtimesEvals ****************************/

int CVSpgmrGetNumJtimesEvals(void *cvode_mem, int *njvevals)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  *njvevals = njtimes;

  return(OKAY);
}

/*************** CVSpgmrGetNumRhsEvals *******************************/

int CVSpgmrGetNumRhsEvals(void *cvode_mem, int *nfevalsSG)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_SETGET_CVMEM_NULL);
    return(LIN_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    fprintf(errfp, MSG_SETGET_LMEM_NULL);
    return(LIN_NO_LMEM);
  }
  cvspgmr_mem = (CVSpgmrMem) lmem;

  *nfevalsSG = nfeSG;

  return(OKAY);
}


/* Additional readability Replacements */

#define pretype (cvspgmr_mem->g_pretype)
#define gstype  (cvspgmr_mem->g_gstype)
#define delt    (cvspgmr_mem->g_delt)
#define maxl    (cvspgmr_mem->g_maxl)
#define psolve  (cvspgmr_mem->g_psolve)
#define pset    (cvspgmr_mem->g_pset)
#define P_data  (cvspgmr_mem->g_P_data)
#define jtimes  (cvspgmr_mem->g_jtimes)
#define j_data  (cvspgmr_mem->g_j_data)


/*************** CVSpgmrInit *****************************************

 This routine does remaining initializations specific to the Spgmr 
 linear solver.

**********************************************************************/

static int CVSpgmrInit(CVodeMem cv_mem)
{
  CVSpgmrMem cvspgmr_mem;
  cvspgmr_mem = (CVSpgmrMem) lmem;

  /* Initialize counters */
  npe = nli = nps = ncfl = nstlpre = 0;
  njtimes = nfeSG = 0;

  /* Check for legal combination pretype - psolve */
  if ((pretype != NONE) && (psolve == NULL)) {
    fprintf(errfp, MSG_PSOLVE_REQ);
    return(LINIT_ERR);
  }

  /* Set setupNonNull = TRUE iff there is preconditioning (pretype != NONE)
     and there is a preconditioning setup phase (pset != NULL)             */
  setupNonNull = (pretype != NONE) && (pset != NULL);

  /* If jtimes is NULL at this time, set it to DQ */
  if (jtimes == NULL) {
    jtimes = CVSpgmrDQJtimes;
    j_data = cv_mem;
  }

  return(LINIT_OK);
}

/*************** CVSpgmrSetup ****************************************

 This routine does the setup operations for the Spgmr linear solver.
 It makes a decision as to whether or not to signal for re-evaluation
 of Jacobian data in the pset routine, based on various state
 variables, then it calls pset.  If we signal for re-evaluation,
 then we reset jcur = *jcurPtr to TRUE, regardless of the pset output.
 In any case, if jcur == TRUE, we increment npe and save nst in nstlpre.

**********************************************************************/

static int CVSpgmrSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                        N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  int  ier;
  CVSpgmrMem cvspgmr_mem;

  cvspgmr_mem = (CVSpgmrMem) lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = ABS((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlpre + CVSPGMR_MSBPRE) ||
         ((convfail == FAIL_BAD_J) && (dgamma < CVSPGMR_DGMAX)) ||
    (convfail == FAIL_OTHER);
  *jcurPtr = jbad;
  jok = !jbad;

  /* Call pset routine and possibly reset jcur */
  ier = pset(tn, ypred, fpred, jok, jcurPtr, gamma, P_data, 
             vtemp1, vtemp2, vtemp3);
  if (jbad) *jcurPtr = TRUE;

  /* If jcur = TRUE, increment npe and save nst value */
  if (*jcurPtr) {
    npe++;
    nstlpre = nst;
  }

  /* Return the same value ier that pset returned */
  return(ier);
}

/*************** CVSpgmrSolve ****************************************

 This routine handles the call to the generic solver SpgmrSolve
 for the solution of the linear system Ax = b with the SPGMR method,
 without restarts.  The solution x is returned in the vector b.

 If the WRMS norm of b is small, we return x = b (if this is the first
 Newton iteration) or x = 0 (if a later Newton iteration).

 Otherwise, we set the tolerance parameter and initial guess (x = 0),
 call SpgmrSolve, and copy the solution x into b.  The x-scaling and
 b-scaling arrays are both equal to ewt, and no restarts are allowed.

 The counters nli, nps, and ncfl are incremented, and the return value
 is set according to the success of SpgmrSolve.  The success flag is
 returned if SpgmrSolve converged, or if this is the first Newton
 iteration and the residual norm was reduced below its initial value.

**********************************************************************/

static int CVSpgmrSolve(CVodeMem cv_mem, N_Vector b, N_Vector ynow,
                        N_Vector fnow)
{
  realtype bnorm, res_norm;
  CVSpgmrMem cvspgmr_mem;
  int nli_inc, nps_inc, ier;
  
  cvspgmr_mem = (CVSpgmrMem) lmem;

  /* Test norm(b); if small, return x = 0 or x = b */
  deltar = delt*tq[4]; 

  bnorm = N_VWrmsNorm(b, ewt);
  if (bnorm <= deltar) {
    if (mnewt > 0) N_VConst(ZERO, b); 
    return(0);
  }

  /* Set vectors ycur and fcur for use by the Atimes and Psolve routines */
  ycur = ynow;
  fcur = fnow;

  /* Set inputs delta and initial guess x = 0 to SpgmrSolve */  
  delta = deltar * sqrtN;
  N_VConst(ZERO, x);
  
  /* Call SpgmrSolve and copy x to b */
  ier = SpgmrSolve(spgmr_mem, cv_mem, x, b, pretype, gstype, delta, 0,
                   cv_mem, ewt, ewt, CVSpgmrAtimes, CVSpgmrPSolve,
                   &res_norm, &nli_inc, &nps_inc);

  N_VScale(ONE, x, b);
  
  /* Increment counters nli, nps, and ncfl */
  nli += nli_inc;
  nps += nps_inc;
  if (ier != 0) ncfl++;

  /* Set return value to -1, 0, or 1 */
  if (ier < 0) return(-1);  
  if ((ier == SPGMR_SUCCESS) || 
      ((ier == SPGMR_RES_REDUCED) && (mnewt == 0)))
    return(0);
  return(1);  
}

/*************** CVSpgmrSolveS ***************************************

 This routine handles the call to the generic solver SpgmrSolve
 for the solution of the linear system Ax = b for sensitivity 'is' 
 with the SPGMR method, without restarts.  
 The solution x is returned in the vector b.

**********************************************************************/

static int CVSpgmrSolveS(CVodeMem cv_mem, N_Vector b, N_Vector ynow,
                         N_Vector fnow, integertype is)
{
  realtype bnorm, res_norm;
  CVSpgmrMem cvspgmr_mem;
  int nli_inc, nps_inc, ier;
  
  cvspgmr_mem = (CVSpgmrMem) lmem;

  /* Test norm(b); if small, return x = 0 or x = b */
  deltar = delt*tq[4]; 

  bnorm = N_VWrmsNorm(b, ewtS[is]);
  if (bnorm <= deltar) {
    if (mnewt > 0) N_VConst(ZERO, b); 
    return(0);
  }
  
  /* Set vectors ycur and fcur for use by the Atimes and Psolve routines */
  ycur = ynow;
  fcur = fnow;

  /* Set inputs delta and initial guess x = 0 to SpgmrSolve */  
  delta = deltar * sqrtN;
  N_VConst(ZERO, x);
  
  /* Call SpgmrSolve and copy x to b */
  ier = SpgmrSolve(spgmr_mem, cv_mem, x, b, pretype, gstype, delta, 0,
                   cv_mem, ewtS[is], ewtS[is], CVSpgmrAtimes, CVSpgmrPSolve,
                   &res_norm, &nli_inc, &nps_inc);

  N_VScale(ONE, x, b);
  
  /* Increment counters nli, nps, and ncfl */
  nli += nli_inc;
  nps += nps_inc;
  if (ier != 0) ncfl++;

  /* Set return value to -1, 0, or 1 */
  if (ier < 0) return(-1);  
  if ((ier == SPGMR_SUCCESS) || 
      ((ier == SPGMR_RES_REDUCED) && (mnewt == 0)))
    return(0);
  return(1);  
}

/*************** CVSpgmrFree *****************************************

 This routine frees memory specific to the Spgmr linear solver.

**********************************************************************/

static void CVSpgmrFree(CVodeMem cv_mem)
{
  CVSpgmrMem cvspgmr_mem;

  cvspgmr_mem = (CVSpgmrMem) lmem;
  
  N_VFree(ytemp);
  N_VFree(x);
  SpgmrFree(spgmr_mem);
  free(cvspgmr_mem);
}

/***************** CVSpgmrAtimes *************************************

 This routine generates the matrix-vector product z = Mv, where
 M = I - gamma*J. The product J*v is obtained by calling the jtimes 
 routine. It is then scaled by -gamma and added to v to obtain M*v.
 The return value is the same as the value returned by jtimes --
 0 if successful, nonzero otherwise.

**********************************************************************/

static int CVSpgmrAtimes(void *cvode_mem, N_Vector v, N_Vector z)
{
  CVodeMem   cv_mem;
  CVSpgmrMem cvspgmr_mem;
  int jtflag;

  cv_mem = (CVodeMem) cvode_mem;
  cvspgmr_mem = (CVSpgmrMem) lmem;

  jtflag = jtimes(v, z, tn, ycur, fcur, j_data, ytemp);
  njtimes++;
  if (jtflag != 0) return(jtflag);

  N_VLinearSum(ONE, v, -gamma, z, z);

  return(0);
}

/*************** CVSpgmrPSolve ***************************************

 This routine interfaces between the generic SpgmrSolve routine and
 the user's psolve routine.  It passes to psolve all required state 
 information from cvode_mem.  Its return value is the same as that
 returned by psolve. Note that the generic SPGMR solver guarantees
 that CVSpgmrPSolve will not be called in the case in which
 preconditioning is not done. This is the only case in which the
 user's psolve routine is allowed to be NULL.

**********************************************************************/

static int CVSpgmrPSolve(void *cvode_mem, N_Vector r, N_Vector z, int lr)
{
  CVodeMem   cv_mem;
  CVSpgmrMem cvspgmr_mem;
  int ier;

  cv_mem = (CVodeMem) cvode_mem;
  cvspgmr_mem = (CVSpgmrMem)lmem;

  ier = psolve(tn, ycur, fcur, r, z, gamma, delta, lr, P_data, ytemp);
  /* This call is counted in nps within the CVSpgmrSolve routine */

  return(ier);     
}


/*************** CVSpgmrDQJtimes *************************************

 This routine generates a difference quotient approximation to
 the Jacobian times vector f_y(t,y) * v. The approximation is 
 Jv = vnrm[f(y + v/vnrm) - f(y)], where vnrm = (WRMS norm of v) is
 input, i.e. the WRMS norm of v/vnrm is 1.

**********************************************************************/

static int CVSpgmrDQJtimes(N_Vector v, N_Vector Jv, realtype t, 
                           N_Vector y, N_Vector fy,
                           void *jac_data, N_Vector work)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;
  realtype vnrm;

  /* jac_data is cvode_mem */
  cv_mem = (CVodeMem) jac_data;
  cvspgmr_mem = (CVSpgmrMem) lmem;

  /* Evaluate norm of v */
  vnrm = N_VWrmsNorm(v, ewt);

  /* Set work = y + (1/vnrm) v */
  N_VLinearSum(ONE/vnrm, v, ONE, y, work);

  /* Set Jv = f(tn, work) */
  f(t, work, Jv, f_data); 
  nfeSG++;

  /* Replace Jv by vnrm*(Jv - fy) */
  N_VLinearSum(ONE, Jv, -ONE, fy, Jv);
  N_VScale(vnrm, Jv, Jv);

  return(0);
}
