/******************************************************************
 * File          : fcvode.c                                       *
 * Programmers   : Alan C. Hindmarsh and Radu Serban @ LLNL       *
 * Version of    : 31 July 2003                                   *
 *----------------------------------------------------------------*
 * This is the implementation file for the Fortran interface to   *
 * the CVODE package.  See fcvode.h for usage.                    *
 * NOTE: some routines are necessarily stored elsewhere to avoid  *
 * linking problems.  Therefore, see also fcvpreco.c, fcvpsol.c,  *
 * fcvjtimes.c, and the five files fcvspgmr**.c (where ** is one  *
 * 01, 10, 11, 20, 21) for all the options available.             *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of types realtype and integertype   */
#include "nvector.h"       /* definitions of type N_Vector and vector macros  */
#include "cvode.h"         /* CVODE constants and prototypes                  */
#include "cvdiag.h"        /* prototypes for CVDIAG interface routines        */
#include "cvdense.h"       /* prototypes for CVDENSE interface routines       */
#include "cvband.h"        /* prototypes for CVBAND interface routines        */
#include "cvspgmr.h"       /* prototypes for CVSPGMR interface routines       */
#include "fcmixpar.h"      /* global F2C_machEnv variable                     */
#include "fcvode.h"        /* actual function names, prototypes, global vars. */

/***************************************************************************/

/* Prototypes of the Fortran routines */
void FCV_FUN(realtype*, realtype*, realtype*);

/**************************************************************************/

void FCV_MALLOC(realtype *t0, realtype *y0, 
                int *meth, int *itmeth, int *iatol, 
                realtype *rtol, realtype *atol,
                int *optin, int *iopt, realtype *ropt, int *ier)
{
  int lmm, iter, itol;
  N_Vector atolvec;
  void *atolptr;

  CV_yvec = N_VMake(y0, F2C_nvspec);
  lmm = (*meth == 1) ? ADAMS : BDF;
  iter = (*itmeth == 1) ? FUNCTIONAL : NEWTON;
  if (*iatol == 1)
    { itol = SS; atolptr = atol; }
  else
    { atolvec = N_VMake(atol, F2C_nvspec);
      itol = SV; atolptr = atolvec; }

  /* Call CVodeCreate, CVodeSet*, and CVodeMalloc to initialize CVODE: 
     lmm     is the method specifier
     iter    is the iteration method specifier
     CVf     is the user's right-hand side function in y'=f(t,y)
     *t0     is the initial time
     CV_yvec is the initial dependent variable vector
     itol    specifies tolerance type
     rtol    is a pointer to the scalar relative tolerance
     atolptr is the absolute tolerance pointer (to scalar or vector)
     F2C_nvspec is the pointer to the vector specification

     A pointer to CVODE problem memory is createded and stored in CV_cvodemem. */


  *ier = 0;

  CV_cvodemem = CVodeCreate(lmm, iter);

  if (CV_cvodemem == NULL) {
    *ier = -1;
    return;
  }

  if (*optin == 1) {
    if (iopt[0]>0)      CVodeSetMaxOrd(CV_cvodemem, iopt[0]);
    if (iopt[1]>0)      CVodeSetMaxNumSteps(CV_cvodemem, iopt[1]);
    if (iopt[2]>0)      CVodeSetMaxHnilWarns(CV_cvodemem, iopt[2]);
    if (iopt[13]>0)     CVodeSetStabLimDet(CV_cvodemem, TRUE);
    if (ropt[0] != 0.0) CVodeSetInitStep(CV_cvodemem, ropt[0]);
    if (ropt[1] > 0.0)  CVodeSetMaxStep(CV_cvodemem, ropt[1]);
    if (ropt[2] > 0.0)  CVodeSetMinStep(CV_cvodemem, ropt[2]);
  }

  *ier = CVodeMalloc(CV_cvodemem, CVf, *t0, CV_yvec,
                     itol, rtol, atolptr, F2C_nvspec);

  CV_iopt = iopt;
  CV_ropt = ropt;

  return;

}

/***************************************************************************/

void FCV_REINIT(realtype *t0, realtype *y0, int *iatol, realtype *rtol,
                realtype *atol, int *optin, int *iopt,
                realtype *ropt, int *ier)
{
  int itol;
  N_Vector atolvec;
  void *atolptr;

  N_VSetData(y0, CV_yvec);

  if (*iatol == 1)
    { itol = SS; atolptr = atol; }
  else
    { atolvec = N_VMake(atol, F2C_nvspec);
      itol = SV; atolptr = atolvec; }

  /* Call CVodeSet* and CVReInit to re-initialize CVODE: 
     CV_cvodemem is the pointer to the CVODE memory block
     CVf     is the user's right-hand side function in y'=f(t,y)
     *t0     is the initial time
     CV_yvec is the initial dependent variable vector
     itol    specifies tolerance type
     rtol    is a pointer to the scalar relative tolerance
     atolptr is the absolute tolerance pointer (to scalar or vector)
     F2C_nvspec is the pointer to the vector specification */

  if (*optin == 1) {
    if (iopt[0]>0)      CVodeSetMaxOrd(CV_cvodemem, iopt[0]);
    if (iopt[1]>0)      CVodeSetMaxNumSteps(CV_cvodemem, iopt[1]);
    if (iopt[2]>0)      CVodeSetMaxHnilWarns(CV_cvodemem, iopt[2]);
    if (iopt[13]>0)     CVodeSetStabLimDet(CV_cvodemem, TRUE);
    if (ropt[0] != 0.0) CVodeSetInitStep(CV_cvodemem, ropt[0]);
    if (ropt[1] > 0.0)  CVodeSetMaxStep(CV_cvodemem, ropt[1]);
    if (ropt[2] > 0.0)  CVodeSetMinStep(CV_cvodemem, ropt[2]);
  }  


  *ier = CVodeReInit(CV_cvodemem, CVf, *t0, CV_yvec,
                     itol, rtol, atolptr);

  CV_iopt = iopt;
  CV_ropt = ropt;

  return;
}

/***************************************************************************/

void FCV_DIAG(int *ier)
{
  /* Call CVDiag:
     CV_cvodemem is the pointer to the CVODE memory block  */

  *ier = CVDiag(CV_cvodemem);

  CV_ls = 3;
}

/***************************************************************************/

void FCV_DENSE(integertype *neq, int *ier)
{
  /* Call CVDense:
     *neq        is the problem size
     CV_cvodemem is the pointer to the CVODE memory block */

  *ier = CVDense(CV_cvodemem, *neq);

  CV_ls = 1;
}

/***************************************************************************/

void FCV_BAND(integertype *neq, integertype *mupper, integertype *mlower, int *ier)
{
  /* Call CVBand:
     CV_cvodemem is the pointer to the CVODE memory block 
     *neq        is the problem size
     *mupper     is the upper bandwidth
     *mlower     is the lower bandwidth */

  *ier = CVBand(CV_cvodemem, *neq, *mupper, *mlower);

  CV_ls = 2;
}

/***************************************************************************/

void FCV_SPGMR(int *pretype, int *gstype, int *maxl, realtype *delt, int *ier)
{
  /* Call CVSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     0           is the preconditioner type (none)
     *maxl       is the maximum Krylov dimension
     *gstype     is the Gram-Schmidt process type
     *delt       is the linear convergence tolerance factor */

  *ier = CVSpgmr(CV_cvodemem, *pretype, *maxl);
  if (*ier != 0) return;

  *ier = CVSpgmrSetGSType(CV_cvodemem, *gstype);
  if (*ier != 0) return;

  *ier = CVSpgmrSetDelt(CV_cvodemem, *delt);
  if (*ier != 0) return;

  CV_ls = 4;

}

/***************************************************************************/

void FCV_REINSPGMR(int *pretype, int *gstype, realtype *delt, int *ier)
{
  /* Call CVSpgmrSet* to specify 
     *gstype     is the Gram-Schmidt process type
     *delt       is the linear convergence tolerance factor */

  *ier = CVSpgmrResetPrecType(CV_cvodemem, *pretype);
  if (*ier != 0) return;

  *ier = CVSpgmrSetGSType(CV_cvodemem, *gstype);
  if (*ier != 0) return;

  *ier = CVSpgmrSetDelt(CV_cvodemem, *delt);
  if (*ier != 0) return;

  CV_ls = 4;

}

/***************************************************************************/

void FCV_CVODE(realtype *tout, realtype *t, realtype *y, int *itask, int *ier)
{
  CVodeMem CV_cvmem;
  realtype h0u;
  integertype liw, lrw;

  /* Call CVode:
     CV_cvodemem is the pointer to the CVODE memory block
     *tout       is the t value where output is desired
     CV_yvec     is the N_Vector containing the solution on return
     t           is the returned independent variable value
     *itask      is the task indicator (NORMAL or ONE_STEP) */

  *ier = CVode(CV_cvodemem, *tout, CV_yvec, t, *itask);
  y = N_VGetData(CV_yvec);

  /* Load optional outputs in iopt & ropt */
  if ( (CV_iopt != NULL) && (CV_ropt != NULL) ) {

    CV_cvmem = (CVodeMem) CV_cvodemem;

    CVodeGetIntegratorStats(CV_cvodemem, 
                            &CV_iopt[3],  /* NST */
                            &CV_iopt[4],  /* NFE  */ 
                            &CV_iopt[5],  /* NSETUPS  */ 
                            &CV_iopt[8],  /* NETF  */ 
                            &CV_iopt[9],  /* QU  */ 
                            &CV_iopt[11], /* QCUR  */ 
                            &h0u,
                            &CV_ropt[3],  /* HU */ 
                            &CV_ropt[4],  /* HCUR  */ 
                            &CV_ropt[5]); /* TCUR  */ 
    CVodeGetTolScaleFactor(CV_cvodemem, &CV_ropt[6]);
    CVodeGetNonlinSolvStats(CV_cvodemem,
                            &CV_iopt[6],  /* NNI */
                            &CV_iopt[7]); /* NCFN */
    CVodeGetWorkSpace(CV_cvodemem, &liw, &lrw);
    CV_iopt[12] = (int) lrw;              /* LENRW */
    CV_iopt[13] = (int) liw;              /* LENIW */
    if (CV_cvmem->cv_sldeton)
      CVodeGetNumStabLimOrderReds(CV_cvodemem, &CV_iopt[14]); /* NOR */

    switch(CV_ls) {
    case 1:
      CVDenseGetNumJacEvals(CV_cvodemem, &CV_iopt[15]);   /* NJE */
      CVDenseGetRealWorkSpace(CV_cvodemem, &lrw); /* LRW */
      CV_iopt[16] = (int) lrw;
      CVDenseGetIntWorkSpace(CV_cvodemem, &liw);  /* LIW */
      CV_iopt[17] = (int) liw;
      break;
    case 2:
      CVBandGetNumJacEvals(CV_cvodemem, &CV_iopt[15]);    /* NJE */
      CVBandGetRealWorkSpace(CV_cvodemem, &lrw); /* LRW */
      CV_iopt[16] = (int) lrw;
      CVBandGetIntWorkSpace(CV_cvodemem, &liw);  /* LIW */
      CV_iopt[17] = (int) liw;
      break;
    case 3:
      CVDiagGetRealWorkSpace(CV_cvodemem, &lrw); /* LRW */
      CV_iopt[16] = (int) lrw;
      CVDiagGetIntWorkSpace(CV_cvodemem, &liw);  /* LIW */
      CV_iopt[17] = (int) liw;
      break;
    case 4:
      CVSpgmrGetNumPrecEvals(CV_cvodemem, &CV_iopt[15]);  /* NPE */
      CVSpgmrGetNumLinIters(CV_cvodemem, &CV_iopt[16]);   /* NLI */
      CVSpgmrGetNumPrecSolves(CV_cvodemem, &CV_iopt[17]); /* NPS */
      CVSpgmrGetNumConvFails(CV_cvodemem, &CV_iopt[18]);  /* NCFL */
      CVSpgmrGetRealWorkSpace(CV_cvodemem, &lrw); /* LRW */
      CV_iopt[19] = (int) lrw;
      CVSpgmrGetIntWorkSpace(CV_cvodemem, &liw);  /* LIW */
      CV_iopt[20] = (int) liw;
      break;
    }

  }

}

/***************************************************************************/

void FCV_DKY (realtype *t, int *k, realtype *dky, int *ier)
{
  /* Call CVodeDky:
     CV_cvodemem is the pointer to the CVODE memory block
     *t          is the t value where output is desired
     *k          is the derivative order
     CV_yvec     is the N_Vector containing the solution derivative on return */

  *ier = CVodeGetDky(CV_cvodemem, *t, *k, CV_yvec);

  dky = N_VGetData(CV_yvec);

}

/***************************************************************************/

void FCV_FREE ()
{
  /* Call CVodeFree:
     CV_cvodemem is the pointer to the CVODE memory block */

  CVodeFree (CV_cvodemem);
}

/***************************************************************************/

/* C function CVf to interface between CVODE and a Fortran subroutine CVFUN.
   Addresses of t, y, and ydot are passed to CVFUN, using the
   routine N_VGetData from the NVECTOR module.
   Auxiliary data is assumed to be communicated by Common. */

void CVf(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  realtype *ydata, *dydata;

  ydata = N_VGetData(y);
  dydata = N_VGetData(ydot);

  FCV_FUN (&t, ydata, dydata);

  N_VSetData(dydata, ydot);

}

