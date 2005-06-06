/*
 * -----------------------------------------------------------------
 * $Revision: 1.36 $
 * $Date: 2005-06-06 21:34:51 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the Fortran interface to
 * the KINSOL package. See fkinsol.h for usage.
 *
 * Note: Some routines are necessarily stored elsewhere to avoid
 * linking problems. See also, therefore, fkinpreco.c, fkinjtimes.c,
 * and fkinbbd.c.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fkinsol.h"        /* prototypes of interfaces and global variables */
#include "kinsol.h"         /* KINSOL constants and prototypes               */
#include "kindense.h"       /* prototypes of KINDENSE interface routines     */
#include "kinsptfqmr.h"     /* prototypes of KINSPTFQMR interface routines   */
#include "kinspbcg.h"       /* prototypes of KINSPBCG interface routines     */
#include "kinspgmr.h"       /* prototypes of KINSPGMR interface routines     */
#include "nvector.h"        /* definition of type N_Vector and prototypes
                               of related routines                           */
#include "sundialstypes.h"  /* definition of type realtype                   */

/*
 * ----------------------------------------------------------------
 * definitions of global variables shared amongst various routines
 * ----------------------------------------------------------------
 */

void *KIN_kinmem;
booleantype KIN_optin;
long int *KIN_iopt;
realtype *KIN_ropt;
int KIN_ls;

/*
 * ----------------------------------------------------------------
 * private constants
 * ----------------------------------------------------------------
 */

#define ZERO RCONST(0.0)

/*
 * ----------------------------------------------------------------
 * prototype of user-supplied fortran routine
 * ----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FK_FUN(realtype*, realtype*);

#ifdef __cplusplus
}
#endif

/*
 * ----------------------------------------------------------------
 * Function : FKIN_MALLOC
 * ----------------------------------------------------------------
 */

void FKIN_MALLOC(long int *msbpre, realtype *fnormtol, realtype *scsteptol,
                 realtype *constraints, int *optin, long int *iopt,
		 realtype *ropt, int *ier)
{
  /* check for required vector operations */

  if ((F2C_KINSOL_vec->ops->nvgetarraypointer == NULL) ||
      (F2C_KINSOL_vec->ops->nvsetarraypointer == NULL)) {
    *ier = -1;
    printf("A required vector operation is not implemented.\n\n");
    return;
  }

  *ier = 0;
  KIN_kinmem = NULL;
  KIN_iopt   = NULL;
  KIN_ropt   = NULL;

  KIN_kinmem = KINCreate();

  if (KIN_kinmem == NULL) {
    *ier = -1;
    return;
  }

  /* F2C_KINSOL_vec used as template vector */

  *ier = KINMalloc(KIN_kinmem, FKINfunc, F2C_KINSOL_vec);

  if (*ier != KIN_SUCCESS) {
    *ier = -1;
    return;
  }

  /* set data pointer in F2C_KINSOL_vec to constraints */

  N_VSetArrayPointer(constraints, F2C_KINSOL_vec);
  KINSetConstraints(KIN_kinmem, F2C_KINSOL_vec);

  /* set optional inputs controlled by arguments to FKIN_MALLOC routine */

  KINSetMaxSetupCalls(KIN_kinmem, *msbpre);
  KINSetFuncNormTol(KIN_kinmem, *fnormtol);
  KINSetScaledStepTol(KIN_kinmem, *scsteptol);

  if (*optin == 1) {
    KIN_optin = TRUE;

    if (iopt[0]  > 0) KINSetPrintLevel(KIN_kinmem, (int) iopt[0]);
    if (iopt[1]  > 0) KINSetNumMaxIters(KIN_kinmem, iopt[1]);
    if (iopt[2]  > 0) KINSetNoInitSetup(KIN_kinmem, TRUE);
    if (iopt[7]  > 0) KINSetEtaForm(KIN_kinmem, (int) iopt[7]);
    if (iopt[8]  > 0) KINSetNoMinEps(KIN_kinmem, TRUE);
    if (iopt[9]  > 0) KINSetNoResMon(KIN_kinmem, TRUE);
    if (iopt[15] > 0) KINSetMaxSubSetupCalls(KIN_kinmem, iopt[15]);

    if (ropt[0] > ZERO) KINSetMaxNewtonStep(KIN_kinmem, ropt[0]);
    if (ropt[1] > ZERO) KINSetRelErrFunc(KIN_kinmem, ropt[1]);
    if (ropt[4] > ZERO) KINSetEtaConstValue(KIN_kinmem, ropt[4]);
    if ((ropt[5] > ZERO) || (ropt[6] > ZERO)) 
      KINSetEtaParams(KIN_kinmem, ropt[5], ropt[6]);
    if (ropt[7] > ZERO) KINSetResMonConstValue(KIN_kinmem, ropt[7]);
    if ((ropt[8] > ZERO) || (ropt[9] > ZERO))
      KINSetResMonParams(KIN_kinmem, ropt[8], ropt[9]);
  }
  else KIN_optin = FALSE;

  /* reset data pointer */

  N_VSetArrayPointer(NULL, F2C_KINSOL_vec);

  KIN_iopt = iopt;
  KIN_ropt = ropt;

  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_DENSE
 * ----------------------------------------------------------------
 */

void FKIN_DENSE(long int *neq, int *ier)
{
  *ier = KINDense(KIN_kinmem, *neq);

  KIN_ls = KIN_DENSE;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SPTFQMR
 * ----------------------------------------------------------------
 */

void FKIN_SPTFQMR(int *maxl, int *ier)
{
  *ier = KINSptfqmr(KIN_kinmem, *maxl);

  KIN_ls = KIN_SPTFQMR;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SPBCG
 * ----------------------------------------------------------------
 */

void FKIN_SPBCG(int *maxl, int *ier)
{
  *ier = KINSpbcg(KIN_kinmem, *maxl);

  KIN_ls = KIN_SPBCG;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SPGMR
 * ----------------------------------------------------------------
 */

void FKIN_SPGMR(int *maxl, int *maxlrst, int *ier)
{
  *ier = KINSpgmr(KIN_kinmem, *maxl);
  KINSpgmrSetMaxRestarts(KIN_kinmem, *maxlrst);

  KIN_ls = KIN_SPGMR;

  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SOL
 * ----------------------------------------------------------------
 */

void FKIN_SOL(realtype *uu, int *globalstrategy, 
              realtype *uscale , realtype *fscale, int *ier)

{
  N_Vector uuvec, uscalevec, fscalevec;

  *ier = 0;
  uuvec = uscalevec = fscalevec = NULL;

  uuvec = F2C_KINSOL_vec;
  N_VSetArrayPointer(uu, uuvec);

  uscalevec = N_VCloneEmpty(F2C_KINSOL_vec);
  N_VSetArrayPointer(uscale, uscalevec);

  fscalevec = N_VCloneEmpty(F2C_KINSOL_vec);
  N_VSetArrayPointer(fscale, fscalevec);

  *ier = KINSol(KIN_kinmem, uuvec, *globalstrategy, uscalevec, fscalevec);

  N_VSetArrayPointer(NULL, uuvec);

  N_VSetArrayPointer(NULL, uscalevec);
  N_VDestroy(uscalevec);

  N_VSetArrayPointer(NULL, fscalevec);
  N_VDestroy(fscalevec);

  /* load optional outputs into iopt[] and ropt[] */

  if ((KIN_iopt != NULL) && (KIN_ropt != NULL)) {

    KINGetNumNonlinSolvIters(KIN_kinmem, &KIN_iopt[3]);
    KINGetNumFuncEvals(KIN_kinmem, &KIN_iopt[4]);
    KINGetNumBetaCondFails(KIN_kinmem, &KIN_iopt[5]);
    KINGetNumBacktrackOps(KIN_kinmem, &KIN_iopt[6]);

    KINGetFuncNorm(KIN_kinmem, &KIN_ropt[2]);
    KINGetStepLength(KIN_kinmem, &KIN_ropt[3]);

    switch(KIN_ls) {

    case KIN_DENSE:
      KINDenseGetNumJacEvals(KIN_kinmem, &KIN_iopt[10]);
      KINDenseGetNumFuncEvals(KIN_kinmem, &KIN_iopt[11]);
      KINDenseGetLastFlag(KIN_kinmem, (int *) &KIN_iopt[12]);

    case KIN_SPTFQMR:
      KINSptfqmrGetNumLinIters(KIN_kinmem, &KIN_iopt[10]);
      KINSptfqmrGetNumPrecEvals(KIN_kinmem, &KIN_iopt[11]);
      KINSptfqmrGetNumPrecSolves(KIN_kinmem, &KIN_iopt[12]);
      KINSptfqmrGetNumConvFails(KIN_kinmem, &KIN_iopt[13]);
      KINSptfqmrGetLastFlag(KIN_kinmem, (int *) &KIN_iopt[14]);
      break;

    case KIN_SPBCG:
      KINSpbcgGetNumLinIters(KIN_kinmem, &KIN_iopt[10]);
      KINSpbcgGetNumPrecEvals(KIN_kinmem, &KIN_iopt[11]);
      KINSpbcgGetNumPrecSolves(KIN_kinmem, &KIN_iopt[12]);
      KINSpbcgGetNumConvFails(KIN_kinmem, &KIN_iopt[13]);
      KINSpbcgGetLastFlag(KIN_kinmem, (int *) &KIN_iopt[14]);
      break;

    case KIN_SPGMR:
      KINSpgmrGetNumLinIters(KIN_kinmem, &KIN_iopt[10]);
      KINSpgmrGetNumPrecEvals(KIN_kinmem, &KIN_iopt[11]);
      KINSpgmrGetNumPrecSolves(KIN_kinmem, &KIN_iopt[12]);
      KINSpgmrGetNumConvFails(KIN_kinmem, &KIN_iopt[13]);
      KINSpgmrGetLastFlag(KIN_kinmem, (int *) &KIN_iopt[14]);
      break;

    }

  }

  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_FREE
 * ----------------------------------------------------------------
 */

void FKIN_FREE(void)
{

  /* call KINFree: KIN_kinmem is the pointer to the KINSOL memory block */

  KINFree(KIN_kinmem);

  N_VSetArrayPointer(NULL , F2C_KINSOL_vec);
  N_VDestroy(F2C_KINSOL_vec);

  return;
}


/*
 * ----------------------------------------------------------------
 * Function : FKINfunc
 * ----------------------------------------------------------------
 * The C function FKINfunc acts as an interface between KINSOL and
 * the Fortran user-supplied subroutine FKFUN. Addresses of the
 * data uu and fdata are passed to FKFUN, using the routine
 * N_VGetArrayPointer from the NVECTOR module. The data in the
 * returned N_Vector fval is set using N_VSetArrayPointer. Auxiliary
 * data is assumed to be communicated by 'Common'.
 * ----------------------------------------------------------------
 */

void FKINfunc(N_Vector uu, N_Vector fval, void *f_data)
{
  realtype *udata, *fdata;

  udata = N_VGetArrayPointer(uu);
  fdata = N_VGetArrayPointer(fval);

  FK_FUN(udata, fdata);

  return;
}
