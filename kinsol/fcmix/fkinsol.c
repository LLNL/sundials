/*
 * -----------------------------------------------------------------
 * $Revision: 1.28.2.3 $
 * $Date: 2005-04-05 20:23:33 $
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
 * linking problems. See also, therefore, fkinpreco.c, fkinpsol.c
 * fkinjtimes.c, and fkinbbd.c.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fkinsol.h"        /* prototypes of interfaces and global variables */
#include "kinsol.h"         /* KINSOL constants and prototypes               */
#include "kinspgmr.h"       /* prototypes of KINSPGMR interface routines     */
#include "nvector.h"        /* definition of type N_Vector and prototypes
                               of related routines                           */
#include "sundialstypes.h"  /* definition of type realtype                   */

/*
 * ----------------------------------------------------------------
 * definitions of global variables shared amongst various routines
 * ----------------------------------------------------------------
 */

realtype *data_F2C_vec;
void *KIN_mem;
long int *KIN_iopt;
realtype *KIN_ropt;
int KIN_ls;

enum { KIN_SPGMR = 1 };

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
  KIN_mem = KINCreate();

  if (KIN_mem == NULL) {
    *ier = -1;
    return;
  }

  /* check for required vector operations */

  if ((F2C_vec->ops->nvgetarraypointer == NULL) ||
      (F2C_vec->ops->nvsetarraypointer == NULL)) {
    *ier = -1;
    printf("A required vector operation is not implemented.\n\n");
    return;
  }

  /* set pointer data_F2C_vec to point to data array from global F2C_vec
     and then overwrite it with the user-supplied constraints */

  data_F2C_vec = N_VGetArrayPointer(F2C_vec);
  N_VSetArrayPointer(constraints, F2C_vec);

  KINSetMaxPrecCalls(KIN_mem, *msbpre);
  KINSetFuncNormTol(KIN_mem, *fnormtol);
  KINSetScaledStepTol(KIN_mem, *scsteptol);
  KINSetConstraints(KIN_mem, F2C_vec);

  if (*optin == 1) {

    if (iopt[0] > 0) KINSetPrintLevel(KIN_mem, (int) iopt[0]);
    if (iopt[1] > 0) KINSetNumMaxIters(KIN_mem, iopt[1]);
    if (iopt[2] > 0) KINSetNoPrecInit(KIN_mem, TRUE);
    if (iopt[7] > 0) KINSetEtaForm(KIN_mem, (int) iopt[7]);
    if (iopt[8] > 0) KINSetNoMinEps(KIN_mem, TRUE);

    if (ropt[0] > ZERO) KINSetMaxNewtonStep(KIN_mem, ropt[0]);
    if (ropt[1] > ZERO) KINSetRelErrFunc(KIN_mem, ropt[1]);
    if (ropt[4] > ZERO) KINSetEtaConstValue(KIN_mem, ropt[4]);
    if (ropt[5] > ZERO || ropt[6] > ZERO) 
      KINSetEtaParams(KIN_mem, ropt[5], ropt[6]);
  }

  /* F2C_vec used as template vector */

  *ier = KINMalloc(KIN_mem, FKINfunc, F2C_vec);

  /* reset data pointer */

  N_VSetArrayPointer(data_F2C_vec, F2C_vec);

  KIN_iopt = iopt;
  KIN_ropt = ropt;

  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SPGMR
 * ----------------------------------------------------------------
 */

void FKIN_SPGMR(int *maxl, int *maxlrst, int *ier)
{
  *ier = KINSpgmr(KIN_mem, *maxl);
  KINSpgmrSetMaxRestarts(KIN_mem, *maxlrst);

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
  long int nniters, nfevals, nbcfails, nbacktr;
  long int nliters, npevals, npsolves, nlcfails;
  int lsflag;
  N_Vector uuvec, uscalevec, fscalevec;
  realtype *data_uuvec, *data_uscalevec, *data_fscalevec;

  uuvec = F2C_vec;
  data_uuvec = N_VGetArrayPointer(uuvec);
  N_VSetArrayPointer(uu, uuvec);

  uscalevec = N_VClone(F2C_vec);
  data_uscalevec = N_VGetArrayPointer(uscalevec);
  N_VSetArrayPointer(uscale, uscalevec);

  fscalevec = N_VClone(F2C_vec);
  data_fscalevec = N_VGetArrayPointer(fscalevec);
  N_VSetArrayPointer(fscale, fscalevec);

  *ier = KINSol(KIN_mem, uuvec, *globalstrategy, uscalevec, fscalevec);

  N_VSetArrayPointer(data_uuvec, uuvec);

  N_VSetArrayPointer(data_uscalevec, uscalevec);
  N_VDestroy(uscalevec);

  N_VSetArrayPointer(data_fscalevec, fscalevec);
  N_VDestroy(fscalevec);

  /* load optional outputs into iopt[] and ropt[] */

  if ( (KIN_iopt != NULL) && (KIN_ropt != NULL) ) {

    KINGetNumNonlinSolvIters(KIN_mem, &nniters);
    KINGetNumFuncEvals(KIN_mem, &nfevals);
    KINGetNumBetaCondFails(KIN_mem, &nbcfails);
    KINGetNumBacktrackOps(KIN_mem, &nbacktr);

    KIN_iopt[3] = nniters;
    KIN_iopt[4] = nfevals;
    KIN_iopt[5] = nbcfails;
    KIN_iopt[6] = nbacktr;

    KINGetFuncNorm(KIN_mem, &KIN_ropt[2]);
    KINGetStepLength(KIN_mem, &KIN_ropt[3]);

    switch(KIN_ls) {

    case KIN_SPGMR:
      KINSpgmrGetNumLinIters(KIN_mem, &nliters);
      KINSpgmrGetNumPrecEvals(KIN_mem, &npevals);
      KINSpgmrGetNumPrecSolves(KIN_mem, &npsolves);
      KINSpgmrGetNumConvFails(KIN_mem, &nlcfails);
      KINSpgmrGetLastFlag(KIN_mem, &lsflag);
      break;

    }

    KIN_iopt[10] = nliters;
    KIN_iopt[11] = npevals;
    KIN_iopt[12] = npsolves;
    KIN_iopt[13] = nlcfails;
    KIN_iopt[14] = lsflag;
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

  /* call KINFree: KIN_mem is the pointer to the KINSOL memory block */

  KINFree(KIN_mem);

  /* restore data array in F2C_vec */

  N_VSetArrayPointer(data_F2C_vec, F2C_vec);

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
