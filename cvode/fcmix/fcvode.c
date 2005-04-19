/*
 * -----------------------------------------------------------------
 * $Revision: 1.46 $
 * $Date: 2005-04-19 21:13:49 $
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

#include "cvband.h"         /* prototypes for CVBAND interface routines        */
#include "cvdense.h"        /* prototypes for CVDENSE interface routines       */
#include "cvdiag.h"         /* prototypes for CVDIAG interface routines        */
#include "cvode.h"          /* CVODE constants and prototypes                  */
#include "cvspbcg.h"        /* prototypes for CVSPBCG interface routines       */
#include "cvspgmr.h"        /* prototypes for CVSPGMR interface routines       */
#include "fcvode.h"         /* actual function names, prototypes, global vars. */
#include "nvector.h"        /* definitions of type N_Vector and vector macros  */
#include "sundialstypes.h"  /* definition of type realtype                     */

/***************************************************************************/

/* Definitions for global variables shared amongst various routines */

N_Vector CV_ewt;

void *CV_cvodemem;
booleantype CV_optin;
long int *CV_iopt;
realtype *CV_ropt;
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
                int *optin, long int *iopt, realtype *ropt, 
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
     Call CVodeCreate, CVodeSet*, and CVodeMalloc to initialize CVODE: 
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

  if (*optin == 1) {
    CV_optin = TRUE;
    if (iopt[0] > 0)     CVodeSetMaxOrd(CV_cvodemem, (int)iopt[0]);
    if (iopt[1] > 0)     CVodeSetMaxNumSteps(CV_cvodemem, iopt[1]);
    if (iopt[2] > 0)     CVodeSetMaxHnilWarns(CV_cvodemem, (int)iopt[2]);
    if (iopt[13] > 0)    CVodeSetStabLimDet(CV_cvodemem, TRUE);
    if (iopt[21] > 0)    CVodeSetMaxErrTestFails(CV_cvodemem, (int)iopt[21]);
    if (iopt[22] > 0)    CVodeSetMaxNonlinIters(CV_cvodemem, (int)iopt[22]);
    if (iopt[23] > 0)    CVodeSetMaxConvFails(CV_cvodemem, (int)iopt[23]);
    if (ropt[0] != ZERO) CVodeSetInitStep(CV_cvodemem, ropt[0]);
    if (ropt[1] > ZERO)  CVodeSetMaxStep(CV_cvodemem, ropt[1]);
    if (ropt[2] > ZERO)  CVodeSetMinStep(CV_cvodemem, ropt[2]);
    if (ropt[7] != ZERO) CVodeSetStopTime(CV_cvodemem, ropt[7]);
    if (ropt[8] > ZERO)  CVodeSetNonlinConvCoef(CV_cvodemem, ropt[8]);
  } else {
    CV_optin = FALSE;
  }

  *ier = CVodeMalloc(CV_cvodemem, FCVf, *t0, F2C_CVODE_vec, itol, *rtol, atolptr);

  /* destroy Vatol if allocated */
  if (itol == CV_SV)
    N_VDestroy(Vatol);

  if(*ier != CV_SUCCESS) {
    *ier = -1;
    return;
  }

  /* Store the unit roundoff in ropt for user access */
  ropt[9] = UNIT_ROUNDOFF;

  CV_iopt = iopt;
  CV_ropt = ropt;

  CV_ewt = NULL;

  return;
}

/***************************************************************************/

void FCV_REINIT(realtype *t0, realtype *y0, int *iatol, realtype *rtol,
                realtype *atol, int *optin, long int *iopt,
                realtype *ropt, int *ier)
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
     Call CVodeSet* and CVReInit to re-initialize CVODE: 
     CVf           is the user's right-hand side function in y'=f(t,y)
     t0            is the initial time
     F2C_CVODE_vec is the initial dependent variable vector
     itol          is the tolerance type
     rtol          is the scalar relative tolerance
     atolptr       is the absolute tolerance pointer (to scalar or vector)
  */

  if (*optin == 1) {
    CV_optin = TRUE;
    if (iopt[0] > 0)     CVodeSetMaxOrd(CV_cvodemem, (int)iopt[0]);
    if (iopt[1] > 0)     CVodeSetMaxNumSteps(CV_cvodemem, iopt[1]);
    if (iopt[2] > 0)     CVodeSetMaxHnilWarns(CV_cvodemem, (int)iopt[2]);
    if (iopt[13] > 0)    CVodeSetStabLimDet(CV_cvodemem, TRUE);
    if (iopt[21] > 0)    CVodeSetMaxErrTestFails(CV_cvodemem, (int)iopt[21]);
    if (iopt[22] > 0)    CVodeSetMaxNonlinIters(CV_cvodemem, (int)iopt[22]);
    if (iopt[23] > 0)    CVodeSetMaxConvFails(CV_cvodemem, (int)iopt[23]);
    if (ropt[0] != ZERO) CVodeSetInitStep(CV_cvodemem, ropt[0]);
    if (ropt[1] > ZERO)  CVodeSetMaxStep(CV_cvodemem, ropt[1]);
    if (ropt[2] > ZERO)  CVodeSetMinStep(CV_cvodemem, ropt[2]);
    if (ropt[7] != ZERO) CVodeSetStopTime(CV_cvodemem, ropt[7]);
    if (ropt[8] > ZERO)  CVodeSetNonlinConvCoef(CV_cvodemem, ropt[8]);
  } else {
    CV_optin = FALSE;
  }

  *ier = CVodeReInit(CV_cvodemem, FCVf, *t0, F2C_CVODE_vec, itol, *rtol, atolptr);

  /* destroy Vatol if allocated */
  if (itol == CV_SV)
    N_VDestroy(Vatol);

  if (*ier != CV_SUCCESS) {
    *ier = -1;
    return;
  }

  CV_iopt = iopt;
  CV_ropt = ropt;

  return;
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
  int qu, qcur;

  /* 
     tout          is the t value where output is desired
     F2C_CVODE_vec is the N_Vector containing the solution on return
     t             is the returned independent variable value
     itask         is the task indicator (1 = CV_NORMAL, 2 = CV_ONE_STEP, 
                                          3 = CV_NORMAL_TSTOP, 4 = CV_ONE_STEP_TSTOP) 
  */

  N_VSetArrayPointer(y, F2C_CVODE_vec);

  *ier = CVode(CV_cvodemem, *tout, F2C_CVODE_vec, t, *itask);

  /* Load optional outputs in iopt & ropt */
  if ((CV_iopt != NULL) && (CV_ropt != NULL)) {

    CVodeGetIntegratorStats(CV_cvodemem, 
                            &CV_iopt[3],  /* NST */
                            &CV_iopt[4],  /* NFE */ 
                            &CV_iopt[5],  /* NSETUPS */ 
                            &CV_iopt[8],  /* NETF */ 
                            &qu,
                            &qcur,
                            &h0u,
                            &CV_ropt[3],  /* HU */ 
                            &CV_ropt[4],  /* HCUR */ 
                            &CV_ropt[5]);  /* TCUR */ 
    CV_iopt[9]  = (long int)qu;  /* QU */ 
    CV_iopt[10] = (long int)qcur;  /* QCUR */ 
    CVodeGetTolScaleFactor(CV_cvodemem, &CV_ropt[6]);
    CVodeGetNonlinSolvStats(CV_cvodemem,
                            &CV_iopt[6],  /* NNI */
                            &CV_iopt[7]);  /* NCFN */
    CVodeGetWorkSpace(CV_cvodemem,
                      &CV_iopt[11],  /* LENRW */
                      &CV_iopt[12]);  /* LENIW */
    if (CV_optin && (CV_iopt[13] > 0))
      CVodeGetNumStabLimOrderReds(CV_cvodemem, &CV_iopt[14]);  /* NOR */

    /* Root finding is on */
    if (CV_nrtfn != 0)
      CVodeGetNumGEvals(CV_cvodemem, &CV_iopt[24]);

    switch(CV_ls) {
    case 1:
      CVDenseGetWorkSpace(CV_cvodemem, &CV_iopt[15], &CV_iopt[16]);  /* LRW and LIW */
      CVDenseGetNumJacEvals(CV_cvodemem, &CV_iopt[17]);  /* NJE */
      CVDenseGetLastFlag(CV_cvodemem, (int *) &CV_iopt[25]);  /* last linear solver flag */
      break;
    case 2:
      CVBandGetWorkSpace(CV_cvodemem, &CV_iopt[15], &CV_iopt[16]);  /* LRW and LIW */
      CVBandGetNumJacEvals(CV_cvodemem, &CV_iopt[17]);  /* NJE */
      CVBandGetLastFlag(CV_cvodemem, (int *) &CV_iopt[25]);  /* last linear solver flag */
      break;
    case 3:
      CVDiagGetWorkSpace(CV_cvodemem, &CV_iopt[15], &CV_iopt[16]);  /* LRW and LIW */
      CVDiagGetLastFlag(CV_cvodemem, (int *) &CV_iopt[25]);  /* last linear solver flag */
      break;
    case 4:
      CVSpgmrGetWorkSpace(CV_cvodemem, &CV_iopt[15], &CV_iopt[16]);  /* LRW and LIW */
      CVSpgmrGetNumPrecEvals(CV_cvodemem, &CV_iopt[17]);  /* NPE */
      CVSpgmrGetNumLinIters(CV_cvodemem, &CV_iopt[18]);  /* NLI */
      CVSpgmrGetNumPrecSolves(CV_cvodemem, &CV_iopt[19]);  /* NPS */
      CVSpgmrGetNumConvFails(CV_cvodemem, &CV_iopt[20]);  /* NCFL */
      CVSpgmrGetLastFlag(CV_cvodemem, (int *) &CV_iopt[25]);  /* last linear solver flag */
      break;
    case 5:
      CVSpbcgGetWorkSpace(CV_cvodemem, &CV_iopt[15], &CV_iopt[16]);  /* LRW and LIW */
      CVSpbcgGetNumPrecEvals(CV_cvodemem, &CV_iopt[17]);  /* NPE */
      CVSpbcgGetNumLinIters(CV_cvodemem, &CV_iopt[18]);  /* NLI */
      CVSpbcgGetNumPrecSolves(CV_cvodemem, &CV_iopt[19]);  /* NPS */
      CVSpbcgGetNumConvFails(CV_cvodemem, &CV_iopt[20]);  /* NCFL */
      CVSpbcgGetLastFlag(CV_cvodemem, (int *) &CV_iopt[25]);  /* last linear solver flag */
      break;
    }
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
