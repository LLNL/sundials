/*
 * -----------------------------------------------------------------
 * $Revision: 1.25 $
 * $Date: 2004-06-21 23:07:12 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvode/LICENSE
 * -----------------------------------------------------------------
 * This is the implementation file for the Fortran interface to   
 * the CVODE package.  See fcvode.h for usage.                    
 * NOTE: some routines are necessarily stored elsewhere to avoid  
 * linking problems.  Therefore, see also fcvpreco.c, fcvpsol.c,  
 * fcvjtimes.c, and fcvspgmr.c for all the options available.    
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of type realtype                    */
#include "nvector.h"       /* definitions of type N_Vector and vector macros  */
#include "cvode.h"         /* CVODE constants and prototypes                  */
#include "cvdiag.h"        /* prototypes for CVDIAG interface routines        */
#include "cvdense.h"       /* prototypes for CVDENSE interface routines       */
#include "cvband.h"        /* prototypes for CVBAND interface routines        */
#include "cvspgmr.h"       /* prototypes for CVSPGMR interface routines       */
#include "fcvode.h"        /* actual function names, prototypes, global vars. */

/***************************************************************************/

/* Prototypes of the Fortran routines */
void FCV_FUN(realtype*, realtype*, realtype*);

/**************************************************************************/

void FCV_MALLOC(realtype *t0, realtype *y0, 
                int *meth, int *itmeth, int *iatol, 
                realtype *rtol, realtype *atol,
                int *optin, long int *iopt, realtype *ropt, 
                int *ier)
{
  int lmm, iter, itol;
  N_Vector atolvec;
  void *atolptr;

  CV_yvec = N_VMake((void *)y0, F2C_nvspec);
  lmm = (*meth == 1) ? ADAMS : BDF;
  iter = (*itmeth == 1) ? FUNCTIONAL : NEWTON;
  if (*iatol == 1)
    { itol = SS; atolptr = atol; }
  else
    { atolvec = N_VMake((void *)atol, F2C_nvspec);
      itol = SV; atolptr = atolvec; }

  /* 
     Call CVodeCreate, CVodeSet*, and CVodeMalloc to initialize CVODE: 
     lmm     is the method specifier
     iter    is the iteration method specifier
     CVf     is the user's right-hand side function in y'=f(t,y)
     *t0     is the initial time
     CV_yvec is the initial dependent variable vector
     itol    specifies tolerance type
     rtol    is a pointer to the scalar relative tolerance
     atolptr is the absolute tolerance pointer (to scalar or vector)
     F2C_nvspec is the pointer to the vector specification

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
    if (iopt[0]>0)      CVodeSetMaxOrd(CV_cvodemem, (int)iopt[0]);
    if (iopt[1]>0)      CVodeSetMaxNumSteps(CV_cvodemem, iopt[1]);
    if (iopt[2]>0)      CVodeSetMaxHnilWarns(CV_cvodemem, (int)iopt[2]);
    if (iopt[13]>0)     CVodeSetStabLimDet(CV_cvodemem, TRUE);
    if (iopt[21]>0)     CVodeSetMaxErrTestFails(CV_cvodemem, (int)iopt[21]);
    if (iopt[22]>0)     CVodeSetMaxNonlinIters(CV_cvodemem, (int)iopt[22]);
    if (iopt[23]>0)     CVodeSetMaxConvFails(CV_cvodemem, (int)iopt[23]);
    if (ropt[0] != 0.0) CVodeSetInitStep(CV_cvodemem, ropt[0]);
    if (ropt[1] > 0.0)  CVodeSetMaxStep(CV_cvodemem, ropt[1]);
    if (ropt[2] > 0.0)  CVodeSetMinStep(CV_cvodemem, ropt[2]);
    if (ropt[7] != 0.0) CVodeSetStopTime(CV_cvodemem, ropt[7]);
    if (ropt[8] > 0.0)  CVodeSetNonlinConvCoef(CV_cvodemem, ropt[8]);
  } else {
    CV_optin = FALSE;
  }

  *ier = CVodeMalloc(CV_cvodemem, FCVf, *t0, CV_yvec,
                     itol, rtol, atolptr, F2C_nvspec);

  /* Store the unit roundoff in ropt for user access */
  ropt[9] = UNIT_ROUNDOFF;

  CV_iopt = iopt;
  CV_ropt = ropt;

  return;

}

/***************************************************************************/

void FCV_REINIT(realtype *t0, realtype *y0, int *iatol, realtype *rtol,
                realtype *atol, int *optin, long int *iopt,
                realtype *ropt, int *ier)
{
  int itol;
  N_Vector atolvec;
  void *atolptr;

  N_VSetData((void *)y0, CV_yvec);

  if (*iatol == 1)
    { itol = SS; atolptr = atol; }
  else
    { atolvec = N_VMake((void *)atol, F2C_nvspec);
      itol = SV; atolptr = atolvec; }

  /* 
     Call CVodeSet* and CVReInit to re-initialize CVODE: 
     CVf     is the user's right-hand side function in y'=f(t,y)
     t0      is the initial time
     CV_yvec is the initial dependent variable vector
     itol    specifies tolerance type
     rtol    is a pointer to the scalar relative tolerance
     atolptr is the absolute tolerance pointer (to scalar or vector)
     F2C_nvspec is the pointer to the vector specification 
  */

  if (*optin == 1) {
    CV_optin = TRUE;
    if (iopt[0]>0)      CVodeSetMaxOrd(CV_cvodemem, (int)iopt[0]);
    if (iopt[1]>0)      CVodeSetMaxNumSteps(CV_cvodemem, iopt[1]);
    if (iopt[2]>0)      CVodeSetMaxHnilWarns(CV_cvodemem, (int)iopt[2]);
    if (iopt[13]>0)     CVodeSetStabLimDet(CV_cvodemem, TRUE);
    if (iopt[21]>0)     CVodeSetMaxErrTestFails(CV_cvodemem, (int)iopt[21]);
    if (iopt[22]>0)     CVodeSetMaxNonlinIters(CV_cvodemem, (int)iopt[22]);
    if (iopt[23]>0)     CVodeSetMaxConvFails(CV_cvodemem, (int)iopt[23]);
    if (ropt[0] != 0.0) CVodeSetInitStep(CV_cvodemem, ropt[0]);
    if (ropt[1] > 0.0)  CVodeSetMaxStep(CV_cvodemem, ropt[1]);
    if (ropt[2] > 0.0)  CVodeSetMinStep(CV_cvodemem, ropt[2]);
    if (ropt[7] != 0.0) CVodeSetStopTime(CV_cvodemem, ropt[7]);
    if (ropt[8] > 0.0)  CVodeSetNonlinConvCoef(CV_cvodemem, ropt[8]);
  } else {
    CV_optin = FALSE;
  }


  *ier = CVodeReInit(CV_cvodemem, FCVf, *t0, CV_yvec,
                     itol, rtol, atolptr);

  CV_iopt = iopt;
  CV_ropt = ropt;

  return;
}

/***************************************************************************/

void FCV_DIAG(int *ier)
{

  *ier = CVDiag(CV_cvodemem);

  CV_ls = 3;
}

/***************************************************************************/

void FCV_DENSE(long int *neq, int *ier)
{
  /* 
     neq  is the problem size
  */

  *ier = CVDense(CV_cvodemem, *neq);

  CV_ls = 1;
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

  CV_ls = 2;
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
  if (*ier != 0) return;

  *ier = CVSpgmrSetGSType(CV_cvodemem, *gstype);
  if (*ier != 0) return;

  *ier = CVSpgmrSetDelt(CV_cvodemem, *delt);
  if (*ier != 0) return;

  CV_ls = 4;

}

/***************************************************************************/

void FCV_SPGMRREINIT(int *pretype, int *gstype, realtype *delt, int *ier)
{
  /* 
     pretype    the preconditioner type
     gstype     the Gram-Schmidt process type
     delt       the linear convergence tolerance factor 
  */

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
  realtype h0u;
  int i, qu, qcur, flag;
  int *rootsfound;

  /* 
     tout      is the t value where output is desired
     CV_yvec   is the N_Vector containing the solution on return
     t         is the returned independent variable value
     itask     is the task indicator (1=NORMAL, 2=ONE_STEP, 
                                      3=NORMAL_TSTOP, 4=ONE_STEP_TSTOP) 
  */

  *ier = CVode(CV_cvodemem, *tout, CV_yvec, t, *itask);
  if (*ier < 0) return;

  y = (realtype *) N_VGetData(CV_yvec);

  /* CVode() succeeded and found at least one root */
  if (*ier == ROOT_RETURN) {
    CVodeGetNumGEvals(CV_cvodemem, &CV_iopt[24]);
    flag = CVodeGetRootInfo(CV_cvodemem, &rootsfound);
    if (flag == SUCCESS) {
      printf("   rootsfound[] = ");
      for (i = 0; i < CV_nrtfn; ++i) printf("%d ", *(rootsfound + i));
      printf("\n");
    }
    else {
      *ier = flag;
      return;
    }
  }

  /* Load optional outputs in iopt & ropt */
  if ( (CV_iopt != NULL) && (CV_ropt != NULL) ) {

    CVodeGetIntegratorStats(CV_cvodemem, 
                            &CV_iopt[3],  /* NST */
                            &CV_iopt[4],  /* NFE  */ 
                            &CV_iopt[5],  /* NSETUPS  */ 
                            &CV_iopt[8],  /* NETF  */ 
                            &qu,
                            &qcur,
                            &h0u,
                            &CV_ropt[3],  /* HU */ 
                            &CV_ropt[4],  /* HCUR  */ 
                            &CV_ropt[5]); /* TCUR  */ 
    CV_iopt[9]  = (long int)qu;    /* QU  */ 
    CV_iopt[10] = (long int)qcur;  /* QCUR  */ 
    CVodeGetTolScaleFactor(CV_cvodemem, &CV_ropt[6]);
    CVodeGetNonlinSolvStats(CV_cvodemem,
                            &CV_iopt[6],  /* NNI */
                            &CV_iopt[7]); /* NCFN */
    CVodeGetWorkSpace(CV_cvodemem, 
                      &CV_iopt[12],       /* LENIW */
                      &CV_iopt[11]);      /* LENRW */
    if ( CV_optin && (CV_iopt[13]>0) )
      CVodeGetNumStabLimOrderReds(CV_cvodemem, &CV_iopt[14]); /* NOR */

    switch(CV_ls) {
    case 1:
      CVDenseGetRealWorkSpace(CV_cvodemem, &CV_iopt[15]); /* LRW */
      CVDenseGetIntWorkSpace(CV_cvodemem, &CV_iopt[16]);  /* LIW */
      CVDenseGetNumJacEvals(CV_cvodemem, &CV_iopt[17]);   /* NJE */
      break;
    case 2:
      CVBandGetRealWorkSpace(CV_cvodemem, &CV_iopt[15]);  /* LRW */
      CVBandGetIntWorkSpace(CV_cvodemem, &CV_iopt[16]);   /* LIW */
      CVBandGetNumJacEvals(CV_cvodemem, &CV_iopt[17]);    /* NJE */
      break;
    case 3:
      CVDiagGetRealWorkSpace(CV_cvodemem, &CV_iopt[15]);  /* LRW */
      CVDiagGetIntWorkSpace(CV_cvodemem, &CV_iopt[16]);   /* LIW */
      break;
    case 4:
      CVSpgmrGetRealWorkSpace(CV_cvodemem, &CV_iopt[15]); /* LRW */
      CVSpgmrGetIntWorkSpace(CV_cvodemem, &CV_iopt[16]);  /* LIW */
      CVSpgmrGetNumPrecEvals(CV_cvodemem, &CV_iopt[17]);  /* NPE */
      CVSpgmrGetNumLinIters(CV_cvodemem, &CV_iopt[18]);   /* NLI */
      CVSpgmrGetNumPrecSolves(CV_cvodemem, &CV_iopt[19]); /* NPS */
      CVSpgmrGetNumConvFails(CV_cvodemem, &CV_iopt[20]);  /* NCFL */
      break;
    }

  }

}

/***************************************************************************/

void FCV_DKY (realtype *t, int *k, realtype *dky, int *ier)
{
  /* 
     t        is the t value where output is desired
     k        is the derivative order
     CV_yvec  is the N_Vector containing the solution derivative on return 
  */

  *ier = CVodeGetDky(CV_cvodemem, *t, *k, CV_yvec);

  dky = (realtype *) N_VGetData(CV_yvec);

}

/***************************************************************************/

void FCV_FREE ()
{
  CVodeFree(CV_cvodemem);
}

/***************************************************************************/

/* 
   C function CVf to interface between CVODE and a Fortran subroutine CVFUN.
   Addresses of t, y, and ydot are passed to CVFUN, using the
   routine N_VGetData from the NVECTOR module.
   Auxiliary data is assumed to be communicated by Common. 
*/

void FCVf(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  realtype *ydata, *dydata;

  ydata  = (realtype *) N_VGetData(y);
  dydata = (realtype *) N_VGetData(ydot);

  FCV_FUN(&t, ydata, dydata);

  N_VSetData((void *)dydata, ydot);

}

