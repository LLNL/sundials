/*
 * -----------------------------------------------------------------
 * $Revision: 1.48 $
 * $Date: 2005-08-12 23:34:28 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh, Radu Serban and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the Fortran interface to
 * the CVODE package.  See fcvode.h for usage.
 * NOTE: some routines are necessarily stored elsewhere to avoid
 * linking problems.  Therefore, see also fcvpreco.c, fcvpsol.c,
 * and fcvjtimes.c for all the options available.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cvband.h"         /* prototypes for CVBAND interface routines        */
#include "cvdense.h"        /* prototypes for CVDENSE interface routines       */
#include "cvdiag.h"         /* prototypes for CVDIAG interface routines        */
#include "cvode.h"          /* CVODE constants and prototypes                  */
#include "cvsptfqmr.h"      /* prototypes for CVSPTFQMR interface routines     */
#include "cvspbcg.h"        /* prototypes for CVSPBCG interface routines       */
#include "cvspgmr.h"        /* prototypes for CVSPGMR interface routines       */
#include "fcvode.h"         /* actual function names, prototypes, global vars. */
#include "nvector.h"        /* definitions of type N_Vector and vector macros  */
#include "sundialstypes.h"  /* definition of type realtype                     */

/***************************************************************************/

/* Definitions for global variables shared amongst various routines */

N_Vector CV_ewt;

void *CV_cvodemem;
long int *CV_iout;
realtype *CV_rout;
int CV_nrtfn;
int CV_ls;

/***************************************************************************/

/* private constant(s) */
#define ZERO RCONST(0.0)

/***************************************************************************/

/* Prototypes of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_FUN(realtype*, realtype*, realtype*);
#ifdef __cplusplus
}
#endif

/**************************************************************************/

void FCV_MALLOC(realtype *t0, realtype *y0, 
                int *meth, int *itmeth, int *iatol, 
                realtype *rtol, realtype *atol,
                long int *iout, realtype *rout, 
                int *ier)
{
  int lmm, iter, itol;
  N_Vector Vatol;
  void *atolptr;

  Vatol = NULL;
  atolptr = NULL;

  /* Initialize itol to avoid compiler warning message */
  itol = -1;

  if(F2C_CVODE_vec->ops->nvgetarraypointer == NULL ||
     F2C_CVODE_vec->ops->nvsetarraypointer == NULL) {
    *ier = -1;
    printf("A required vector operation is not implemented.\n\n");
    return;
  }

  /* Set data in F2C_CVODE_vec to y0 */
  N_VSetArrayPointer(y0, F2C_CVODE_vec);

  lmm = (*meth == 1) ? CV_ADAMS : CV_BDF;
  iter = (*itmeth == 1) ? CV_FUNCTIONAL : CV_NEWTON;
  switch (*iatol) {
  case 1:
    itol = CV_SS; 
    atolptr = (void *) atol; 
    break;
  case 2:
    itol = CV_SV; 
    Vatol = N_VCloneEmpty(F2C_CVODE_vec);
    N_VSetArrayPointer(atol, Vatol);
    atolptr = (void *) Vatol; 
    break;
  case 3:
    itol = CV_WF;
    break;
  }

  /* 
     Call CVodeCreate and CVodeMalloc to initialize CVODE: 
     lmm           is the method specifier
     iter          is the iteration method specifier
     CVf           is the user's right-hand side function in y'=f(t,y)
     *t0           is the initial time
     F2C_CVODE_vec is the initial dependent variable vector
     itol          is the tolerance type
     rtol          is the scalar relative tolerance
     atolptr       is the absolute tolerance pointer (to scalar or vector)

     A pointer to CVODE problem memory is createded and stored in CV_cvodemem. 
  */

  *ier = 0;

  CV_cvodemem = CVodeCreate(lmm, iter);

  if (CV_cvodemem == NULL) {
    *ier = -1;
    return;
  }

  *ier = CVodeMalloc(CV_cvodemem, FCVf, *t0, F2C_CVODE_vec, itol, *rtol, atolptr);

  /* destroy Vatol if allocated */
  if (itol == CV_SV)
    N_VDestroy(Vatol);

  if(*ier != CV_SUCCESS) {
    *ier = -1;
    return;
  }

  /* Grab optional output arrays and store them in global vriables */

  CV_iout = iout;
  CV_rout = rout;

  /* Store the unit roundoff in rout for user access */
  CV_rout[5] = UNIT_ROUNDOFF;

  CV_ewt = NULL;

  return;
}

/***************************************************************************/

void FCV_REINIT(realtype *t0, realtype *y0, 
                int *iatol, realtype *rtol, realtype *atol, 
                int *ier)
{
  int itol;
  N_Vector Vatol;
  void *atolptr;

  Vatol = NULL;
  atolptr = NULL;

  /* Initialize itol to avoid compiler warning message */
  itol = -1;

  /* Set data in F2C_CVODE_vec to y0 */
  N_VSetArrayPointer(y0, F2C_CVODE_vec);

  switch (*iatol) {
  case 1:
    itol = CV_SS; 
    atolptr = (void *) atol; 
    break;
  case 2:
    itol = CV_SV; 
    Vatol = N_VCloneEmpty(F2C_CVODE_vec);
    N_VSetArrayPointer(atol, Vatol);
    atolptr = (void *) Vatol; 
    break;
  case 3:
    itol = CV_WF;
    break;
  }

  /* 
     Call CVReInit to re-initialize CVODE: 
     CVf           is the user's right-hand side function in y'=f(t,y)
     t0            is the initial time
     F2C_CVODE_vec is the initial dependent variable vector
     itol          is the tolerance type
     rtol          is the scalar relative tolerance
     atolptr       is the absolute tolerance pointer (to scalar or vector)
  */

  *ier = CVodeReInit(CV_cvodemem, FCVf, *t0, F2C_CVODE_vec, itol, *rtol, atolptr);

  /* destroy Vatol if allocated */
  if (itol == CV_SV)
    N_VDestroy(Vatol);

  if (*ier != CV_SUCCESS) {
    *ier = -1;
    return;
  }

  return;
}

/***************************************************************************/

void FCV_SETIIN(char key_name[], long int *ival, int *ier, int key_len)
{

  char key[20];

  strncpy(key, key_name, key_len);

  printf("KEY: >>%s<<    LEN:%d     NEW_KEY: >>%s<<\n",key_name,key_len,key);
  printf("VAL: %ld\n",*ival);

  if (!strcmp(key,"MAX_ORD")) 
    *ier = CVodeSetMaxOrd(CV_cvodemem, (int) *ival);
  else if (!strcmp(key,"MAX_NSTEPS")) 
    *ier = CVodeSetMaxNumSteps(CV_cvodemem, (int) *ival);
  else if (!strcmp(key,"MAX_ERRFAIL")) 
    *ier = CVodeSetMaxErrTestFails(CV_cvodemem, (int) *ival);
  else if (!strcmp(key,"MAX_NITERS")) 
    *ier = CVodeSetMaxNonlinIters(CV_cvodemem, (int) *ival);
  else if (!strcmp(key,"MAX_CONVFAIL")) 
    *ier = CVodeSetMaxConvFails(CV_cvodemem, (int) *ival);
  else if (!strcmp(key,"HNIL_WARNS")) 
    *ier = CVodeSetMaxHnilWarns(CV_cvodemem, (int) *ival);
  else if (!strcmp(key,"STAB_LIM")) 
    *ier = CVodeSetStabLimDet(CV_cvodemem, (int) *ival);
  else {
    *ier = -99;
    printf("FCVSETIIN: Unrecognized key.\n\n");
  }

}

/***************************************************************************/

void FCV_SETRIN(char key_name[], realtype *rval, int *ier, int key_len)
{
  char key[20];

  strncpy(key, key_name, key_len);

  if (!strcmp(key,"INIT_STEP")) 
    *ier = CVodeSetInitStep(CV_cvodemem, *rval);
  else if (!strcmp(key,"MAX_STEP")) 
    *ier = CVodeSetMaxStep(CV_cvodemem, *rval);
  else if (!strcmp(key,"MIN_STEP")) 
    *ier = CVodeSetMinStep(CV_cvodemem, *rval);
  else if (!strcmp(key,"STOP_TIME")) 
    *ier = CVodeSetStopTime(CV_cvodemem, *rval);
  else if (!strcmp(key,"NLCONV_COEF")) 
    *ier = CVodeSetNonlinConvCoef(CV_cvodemem, *rval);
  else {
    *ier = -99;
    printf("FCVSETRIN: Unrecognized key.\n\n");
  }

}

/***************************************************************************/

void FCV_DENSE(long int *neq, int *ier)
{
  /* neq  is the problem size */

  *ier = CVDense(CV_cvodemem, *neq);

  CV_ls = CV_LS_DENSE;
}

/***************************************************************************/

void FCV_BAND(long int *neq, long int *mupper, long int *mlower, int *ier)
{
  /* 
     neq        is the problem size
     mupper     is the upper bandwidth
     mlower     is the lower bandwidth 
  */

  *ier = CVBand(CV_cvodemem, *neq, *mupper, *mlower);

  CV_ls = CV_LS_BAND;
}

/***************************************************************************/

void FCV_DIAG(int *ier)
{
  *ier = CVDiag(CV_cvodemem);

  CV_ls = CV_LS_DIAG;
}

/***************************************************************************/

void FCV_SPTFQMR(int *pretype, int *maxl, realtype *delt, int *ier)
{
  /* 
     pretype    the preconditioner type
     maxl       the maximum Krylov dimension
     delt       the linear convergence tolerance factor 
  */

  *ier = CVSptfqmr(CV_cvodemem, *pretype, *maxl);
  if (*ier != CVSPTFQMR_SUCCESS) return;

  *ier = CVSptfqmrSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPTFQMR_SUCCESS) return;

  CV_ls = CV_LS_SPTFQMR;
}

/***************************************************************************/

void FCV_SPBCG(int *pretype, int *maxl, realtype *delt, int *ier)
{
  /* 
     pretype    the preconditioner type
     maxl       the maximum Krylov dimension
     delt       the linear convergence tolerance factor 
  */

  *ier = CVSpbcg(CV_cvodemem, *pretype, *maxl);
  if (*ier != CVSPBCG_SUCCESS) return;

  *ier = CVSpbcgSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPBCG_SUCCESS) return;

  CV_ls = CV_LS_SPBCG;
}

/***************************************************************************/

void FCV_SPGMR(int *pretype, int *gstype, int *maxl, realtype *delt, int *ier)
{
  /* 
     pretype    the preconditioner type
     maxl       the maximum Krylov dimension
     gstype     the Gram-Schmidt process type
     delt       the linear convergence tolerance factor 
  */

  *ier = CVSpgmr(CV_cvodemem, *pretype, *maxl);
  if (*ier != CVSPGMR_SUCCESS) return;

  *ier = CVSpgmrSetGSType(CV_cvodemem, *gstype);
  if (*ier != CVSPGMR_SUCCESS) return;

  *ier = CVSpgmrSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPGMR_SUCCESS) return;

  CV_ls = CV_LS_SPGMR;
}

/***************************************************************************/

void FCV_SPTFQMRREINIT(int *pretype, realtype *delt, int *ier)
{
  /* 
     pretype    the preconditioner type
     delt       the linear convergence tolerance factor 
  */

  *ier = CVSptfqmrSetPrecType(CV_cvodemem, *pretype);
  if (*ier != CVSPTFQMR_SUCCESS) return;

  *ier = CVSptfqmrSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPTFQMR_SUCCESS) return;

  CV_ls = CV_LS_SPTFQMR;
}

/***************************************************************************/

void FCV_SPBCGREINIT(int *pretype, realtype *delt, int *ier)
{
  /* 
     pretype    the preconditioner type
     delt       the linear convergence tolerance factor 
  */

  *ier = CVSpbcgSetPrecType(CV_cvodemem, *pretype);
  if (*ier != CVSPBCG_SUCCESS) return;

  *ier = CVSpbcgSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPBCG_SUCCESS) return;

  CV_ls = CV_LS_SPBCG;
}

/***************************************************************************/

void FCV_SPGMRREINIT(int *pretype, int *gstype, realtype *delt, int *ier)
{
  /* 
     pretype    the preconditioner type
     gstype     the Gram-Schmidt process type
     delt       the linear convergence tolerance factor 
  */

  *ier = CVSpgmrSetPrecType(CV_cvodemem, *pretype);
  if (*ier != CVSPGMR_SUCCESS) return;

  *ier = CVSpgmrSetGSType(CV_cvodemem, *gstype);
  if (*ier != CVSPGMR_SUCCESS) return;

  *ier = CVSpgmrSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPGMR_SUCCESS) return;

  CV_ls = CV_LS_SPGMR;
}

/***************************************************************************/

void FCV_CVODE(realtype *tout, realtype *t, realtype *y, int *itask, int *ier)
{
  realtype h0u;

  /* 
     tout          is the t value where output is desired
     F2C_CVODE_vec is the N_Vector containing the solution on return
     t             is the returned independent variable value
     itask         is the task indicator (1 = CV_NORMAL, 2 = CV_ONE_STEP, 
                                          3 = CV_NORMAL_TSTOP, 4 = CV_ONE_STEP_TSTOP) 
  */

  N_VSetArrayPointer(y, F2C_CVODE_vec);

  *ier = CVode(CV_cvodemem, *tout, F2C_CVODE_vec, t, *itask);

  /* Load optional outputs in iout & rout */
  CVodeGetWorkSpace(CV_cvodemem,
                    &CV_iout[0],                          /* LENRW   */
                    &CV_iout[1]);                         /* LENIW   */
  CVodeGetIntegratorStats(CV_cvodemem, 
                          &CV_iout[2],                    /* NST     */
                          &CV_iout[3],                    /* NFE     */ 
                          &CV_iout[7],                    /* NSETUPS */ 
                          &CV_iout[4],                    /* NETF    */ 
                          (int *) &CV_iout[8],            /* QU      */
                          (int *) &CV_iout[9],            /* QCUR    */
                          &CV_rout[0],                    /* H0U     */
                          &CV_rout[1],                    /* HU      */ 
                          &CV_rout[2],                    /* HCUR    */ 
                          &CV_rout[3]);                   /* TCUR    */ 
  CVodeGetTolScaleFactor(CV_cvodemem, &CV_rout[4]);       /* TOLSFAC */
  CVodeGetNonlinSolvStats(CV_cvodemem,
                          &CV_iout[6],                    /* NNI     */
                          &CV_iout[5]);                   /* NCFN    */
  CVodeGetNumStabLimOrderReds(CV_cvodemem, &CV_iout[10]); /* NOR     */
  
  /* Root finding is on */
  if (CV_nrtfn != 0)
    CVodeGetNumGEvals(CV_cvodemem, &CV_iout[11]);         /* NGE     */
  
  switch(CV_ls) {
  case CV_LS_DENSE:
    CVDenseGetWorkSpace(CV_cvodemem, &CV_iout[12], &CV_iout[13]);  /* LRW and LIW */
    CVDenseGetLastFlag(CV_cvodemem, (int *) &CV_iout[14]);         /* last linear solver flag */
    CVDenseGetNumRhsEvals(CV_cvodemem, &CV_iout[15]);              /* NFE */
    CVDenseGetNumJacEvals(CV_cvodemem, &CV_iout[16]);              /* NJE */
    break;
  case CV_LS_BAND:
    CVBandGetWorkSpace(CV_cvodemem, &CV_iout[12], &CV_iout[13]);   /* LRW and LIW */
    CVBandGetLastFlag(CV_cvodemem, (int *) &CV_iout[14]);          /* last linear solver flag */
    CVBandGetNumRhsEvals(CV_cvodemem, &CV_iout[15]);               /* NFE */
    CVBandGetNumJacEvals(CV_cvodemem, &CV_iout[16]);               /* NJE */
    break;
  case CV_LS_DIAG:
    CVDiagGetWorkSpace(CV_cvodemem, &CV_iout[12], &CV_iout[13]);   /* LRW and LIW */
    CVDiagGetLastFlag(CV_cvodemem, (int *) &CV_iout[14]);          /* last linear solver flag */
    CVDiagGetNumRhsEvals(CV_cvodemem, &CV_iout[15]);               /* NFE */
    break;
  case CV_LS_SPGMR:
    CVSpgmrGetWorkSpace(CV_cvodemem, &CV_iout[12], &CV_iout[13]);  /* LRW and LIW */
    CVSpgmrGetLastFlag(CV_cvodemem, (int *) &CV_iout[14]);         /* last linear solver flag */
    CVSpgmrGetNumRhsEvals(CV_cvodemem, &CV_iout[15]);              /* NFE */
    CVSpgmrGetNumJtimesEvals(CV_cvodemem, &CV_iout[16]);           /* NJTV */
    CVSpgmrGetNumPrecEvals(CV_cvodemem, &CV_iout[17]);             /* NPE */
    CVSpgmrGetNumPrecSolves(CV_cvodemem, &CV_iout[18]);            /* NPS */
    CVSpgmrGetNumLinIters(CV_cvodemem, &CV_iout[19]);              /* NLI */
    CVSpgmrGetNumConvFails(CV_cvodemem, &CV_iout[20]);             /* NCFL */
    break;
  case CV_LS_SPBCG:
    CVSpbcgGetWorkSpace(CV_cvodemem, &CV_iout[12], &CV_iout[13]);  /* LRW and LIW */
    CVSpbcgGetLastFlag(CV_cvodemem, (int *) &CV_iout[14]);         /* last linear solver flag */
    CVSpbcgGetNumRhsEvals(CV_cvodemem, &CV_iout[15]);              /* NFE */
    CVSpbcgGetNumJtimesEvals(CV_cvodemem, &CV_iout[16]);           /* NJTV */
    CVSpbcgGetNumPrecEvals(CV_cvodemem, &CV_iout[17]);             /* NPE */
    CVSpbcgGetNumPrecSolves(CV_cvodemem, &CV_iout[18]);            /* NPS */
    CVSpbcgGetNumLinIters(CV_cvodemem, &CV_iout[19]);              /* NLI */
    CVSpbcgGetNumConvFails(CV_cvodemem, &CV_iout[20]);             /* NCFL */
    break;
  case CV_LS_SPTFQMR:
    CVSptfqmrGetWorkSpace(CV_cvodemem, &CV_iout[12], &CV_iout[13]); /* LRW and LIW */
    CVSptfqmrGetLastFlag(CV_cvodemem, (int *) &CV_iout[14]);        /* last linear solver flag */
    CVSptfqmrGetNumRhsEvals(CV_cvodemem, &CV_iout[15]);             /* NFE */
    CVSptfqmrGetNumJtimesEvals(CV_cvodemem, &CV_iout[16]);          /* NJTV */
    CVSptfqmrGetNumPrecEvals(CV_cvodemem, &CV_iout[17]);            /* NPE */
    CVSptfqmrGetNumPrecSolves(CV_cvodemem, &CV_iout[18]);           /* NPS */
    CVSptfqmrGetNumLinIters(CV_cvodemem, &CV_iout[19]);             /* NLI */
    CVSptfqmrGetNumConvFails(CV_cvodemem, &CV_iout[20]);            /* NCFL */
    break;
  }
}

/***************************************************************************/

void FCV_DKY (realtype *t, int *k, realtype *dky, int *ier)
{
  /* 
     t             is the t value where output is desired
     k             is the derivative order
     F2C_CVODE_vec is the N_Vector containing the solution derivative on return 
  */

  N_VSetArrayPointer(dky, F2C_CVODE_vec);

  *ier = CVodeGetDky(CV_cvodemem, *t, *k, F2C_CVODE_vec);

}

/***************************************************************************/

void FCV_FREE ()
{
  CVodeFree(CV_cvodemem);
  N_VDestroy(F2C_CVODE_vec);
  if (CV_ewt != NULL) N_VDestroy(CV_ewt);
}

/***************************************************************************/

/* 
 * C function CVf to interface between CVODE and a Fortran subroutine FCVFUN.
 * Addresses of t, y, and ydot are passed to CVFUN, using the
 * routine N_VGetArrayPointer from the NVECTOR module.
 * Auxiliary data is assumed to be communicated by Common. 
 */

void FCVf(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  realtype *ydata, *dydata;

  ydata  = N_VGetArrayPointer(y);
  dydata = N_VGetArrayPointer(ydot);

  FCV_FUN(&t, ydata, dydata);
}
