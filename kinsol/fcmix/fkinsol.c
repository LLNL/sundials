/******************************************************************
 * File          : fkinsol.c                                      *
 * Programmers   : Allan G Taylor, Alan C. Hindmarsh, and         * 
 *                 Radu Serban @ LLNL                             *
 * Version of    : 07 February 2004                               *
 *----------------------------------------------------------------*
 * This is the implementation file for the Fortran interface to   *
 * the KINSOL package. See fkinsol.h for usage.                   *
 * NOTE: some routines are necessarily stored elsewhere to avoid  *
 * linking problems. See also, therefore, fkinpreco.c, fkinpsol.c *
 * fkinuatimes.c, and fkinspgmr.c                                 *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of type realtype                 */
#include "nvector.h"       /* def's of type N_Vector and related routines  */
#include "kinsol.h"        /* KINSOL constants and prototypes              */
#include "kinspgmr.h"      /* prototypes of KINSPGMR interface routines    */
#include "fkinsol.h"       /* prototypes of interfaces, global variables   */

/**************************************************************************/

/* Prototype of the user-supplied Fortran routine */
void FK_FUN(realtype*, realtype*);

/**************************************************************************/

void FKIN_MALLOC(long int *msbpre, realtype *fnormtol, realtype *scsteptol,
                 realtype *constraints, 
                 int *optin, long int *iopt, realtype *ropt,
                 int *ier)
{
  N_Vector constr_vec;

  KIN_mem = KINCreate();

  if (KIN_mem == NULL) {
    *ier = -1;
    return;
  }

  constr_vec = N_VMake(constraints, F2C_nvspec);

  KINSetMaxPrecCalls(KIN_mem, *msbpre);
  KINSetFuncNormTol(KIN_mem, *fnormtol);
  KINSetScaledStepTol(KIN_mem, *scsteptol);
  KINSetConstraints(KIN_mem, constr_vec);

  if (*optin == 1) {

    if(iopt[0]>0) KINSetPrintLevel(KIN_mem, (int) iopt[0]);
    if(iopt[1]>0) KINSetNumMaxIters(KIN_mem, iopt[1]);
    if(iopt[2]>0) KINSetNoPrecInit(KIN_mem, TRUE);
    if(iopt[7]>=0) KINSetEtaForm(KIN_mem, (int) iopt[7]);
    if(iopt[8]>0) KINSetNoMinEps(KIN_mem, TRUE);

    if(ropt[0]>0.0) KINSetMaxNewtonStep(KIN_mem, ropt[0]);
    if(ropt[1]>0.0) KINSetRelErrFunc(KIN_mem, ropt[1]);
    if(ropt[4]>0.0) KINSetEtaConstValue(KIN_mem, ropt[4]);
    if(ropt[5]>0.0 || ropt[6]>0.0) 
      KINSetEtaParams(KIN_mem, ropt[5], ropt[6]);
  }

  *ier = KINMalloc(KIN_mem, FKINfunc, F2C_nvspec);

  KIN_iopt = iopt;
  KIN_ropt = ropt;
}


/***************************************************************************/

void FKIN_SPGMR(int *maxl, int *maxlrst, int *ier)
{
  *ier = KINSpgmr(KIN_mem, *maxl);
  KINSpgmrSetMaxRestarts(KIN_mem, *maxlrst);

  KIN_ls = 1;
}

/***************************************************************************/

void FKIN_SOL(realtype *uu, int *globalstrategy, 
              realtype *uscale , realtype *fscale, int *ier)

{ 
  N_Vector uuvec, uscalevec, fscalevec;
  long int nniters, nfevals, nbcfails, nbacktr;
  long int nliters, npevals, npsolves, nlcfails;

  uuvec     = N_VMake(uu, F2C_nvspec);
  uscalevec = N_VMake(uscale, F2C_nvspec);
  fscalevec = N_VMake(fscale, F2C_nvspec);

  *ier = KINSol(KIN_mem, uuvec, *globalstrategy, uscalevec, fscalevec);

  N_VDispose(uuvec);
  N_VDispose(uscalevec);
  N_VDispose(fscalevec);

  /* Load optional outputs in iopt & ropt */
  if ( (KIN_iopt != NULL) && (KIN_ropt != NULL) ) {

    KINGetNumNonlinSolvIters(KIN_mem, &nniters);
    KINGetNumFuncEvals(KIN_mem, &nfevals);
    KINGetNumBetaCondFails(KIN_mem, &nbcfails);
    KINGetNumBacktrackOps(KIN_mem, &nbacktr);

    KIN_iopt[3] = nniters;
    KIN_iopt[4] = nfevals;
    KIN_iopt[5] = nbcfails;
    KIN_iopt[6] = nbacktr;

    KINGetFuncNorm(KIN_mem, &KIN_ropt[3]);
    KINGetStepLength(KIN_mem, &KIN_ropt[4]);

    switch(KIN_ls) {
    case 1:
      KINSpgmrGetNumLinIters(KIN_mem, &nliters);
      KINSpgmrGetNumPrecEvals(KIN_mem, &npevals);
      KINSpgmrGetNumPrecSolves(KIN_mem, &npsolves);
      KINSpgmrGetNumConvFails(KIN_mem, &nlcfails);
      break;
    }

    KIN_iopt[10] = nliters;
    KIN_iopt[11] = npevals;
    KIN_iopt[12] = npsolves;
    KIN_iopt[13] = nlcfails;

  }
}

/***************************************************************************/

void FKIN_FREE()
{
  /* Call KINFree:
     KIN_mem is the pointer to the KINSOL memory block */

  KINFree(KIN_mem);
}

/***************************************************************************/


/* C function KINfunc acts as an interface between KINSol and the Fortran 
   user-supplied subroutine KFUN.
   Addresses of Neq and the data for uu and fval are passed to KFUN,
   using the routine N_VGetData from the NVECTOR module.
   The data in the returned N_Vector fval is set using N_VSetData. 
   Auxiliary data is assumed to be communicated by Common. */

void FKINfunc(N_Vector uu, N_Vector fval, void *f_data)
{
  realtype *udata, *fdata;

  udata = N_VGetData(uu);
  fdata = N_VGetData(fval);

  FK_FUN(udata, fdata);

  N_VSetData(fdata, fval);

}
