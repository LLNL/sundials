/*
 * -----------------------------------------------------------------
 * $Revision: 1.43 $
 * $Date: 2006-02-02 00:36:20 $
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
#include <string.h>

#include "fkinsol.h"          /* prototypes of interfaces and global variables  */
#include "kinsol.h"           /* KINSOL constants and prototypes                */

#include "kinsol_band.h"      /* prototypes of KINBAND interface routines       */
#include "kinsol_dense.h"     /* prototypes of KINDENSE interface routines      */
#include "kinsol_sptfqmr.h"   /* prototypes of KINSPTFQMR interface routines    */
#include "kinsol_spbcgs.h"    /* prototypes of KINSPBCG interface routines      */
#include "kinsol_spgmr.h"     /* prototypes of KINSPGMR interface routines      */

#include "kinsol_impl.h"      /* definition of KINMem type                      */

#include "sundials_nvector.h" /* definitions of type N_Vector and vector macros */
#include "sundials_types.h"   /* definition of type realtype                    */

/*
 * ----------------------------------------------------------------
 * definitions of global variables shared amongst various routines
 * ----------------------------------------------------------------
 */

void *KIN_kinmem;
long int *KIN_iout;
realtype *KIN_rout;
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

extern void FK_FUN(realtype*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/*
 * ----------------------------------------------------------------
 * Function : FKIN_MALLOC
 * ----------------------------------------------------------------
 */

void FKIN_MALLOC(long int *iout, realtype *rout, int *ier)
{
  
  /* check for required vector operations */
  if ((F2C_KINSOL_vec->ops->nvgetarraypointer == NULL) ||
      (F2C_KINSOL_vec->ops->nvsetarraypointer == NULL)) {
    *ier = -1;
    printf("A required vector operation is not implemented.\n\n");
    return;
  }

  /* Initialize pointers to NULL */
  KIN_kinmem = NULL;

  /* Create KINSOL object */
  KIN_kinmem = KINCreate();
  if (KIN_kinmem == NULL) {
    *ier = -1;
    return;
  }

  /* Call KINMalloc */
  *ier = 0;
  *ier = KINMalloc(KIN_kinmem, FKINfunc, F2C_KINSOL_vec);

  /* On failure, exit */
  if (*ier != KIN_SUCCESS) {
    *ier = -1;
    return;
  }

  /* Grab optional output arrays and store them in global variables */
  KIN_iout = iout;
  KIN_rout = rout;

  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SETIIN
 * ----------------------------------------------------------------
 */

void FKIN_SETIIN(char key_name[], long int *ival, int *ier, int key_len)
{
  if (!strncmp(key_name,"PRNT_LEVEL", (size_t)key_len)) 
    *ier = KINSetPrintLevel(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NITERS", (size_t)key_len)) 
    *ier = KINSetNumMaxIters(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"ETA_FORM", (size_t)key_len)) 
    *ier = KINSetEtaForm(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"MAX_SETUPS", (size_t)key_len)) 
    *ier = KINSetMaxSetupCalls(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"MAX_SP_SETUPS", (size_t)key_len)) 
    *ier = KINSetMaxSubSetupCalls(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"NO_INIT_SETUP", (size_t)key_len)) 
    *ier = KINSetNoInitSetup(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"NO_MIN_EPS", (size_t)key_len)) 
    *ier = KINSetNoMinEps(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"NO_RES_MON", (size_t)key_len)) 
    *ier = KINSetNoResMon(KIN_kinmem, (int) *ival);
  else {
    *ier = -99;
    printf("FKINSETIIN: Unrecognized key.\n\n");
  }

}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SETRIN
 * ----------------------------------------------------------------
 */

void FKIN_SETRIN(char key_name[], realtype *rval, int *ier, int key_len)
{

  if (!strncmp(key_name,"FNORM_TOL", (size_t)key_len)) 
    *ier = KINSetFuncNormTol(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"SSTEP_TOL", (size_t)key_len)) 
    *ier = KINSetScaledStepTol(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"MAX_STEP", (size_t)key_len)) 
    *ier = KINSetMaxNewtonStep(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"RERR_FUNC", (size_t)key_len)) 
    *ier = KINSetRelErrFunc(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"ETA_CONST", (size_t)key_len)) 
    *ier = KINSetEtaConstValue(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"ETA_PARAMS", (size_t)key_len)) 
    *ier = KINSetEtaParams(KIN_kinmem, rval[0], rval[1]);
  else if (!strncmp(key_name,"RMON_CONST", (size_t)key_len)) 
    *ier = KINSetResMonConstValue(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"RMON_PARAMS", (size_t)key_len)) 
    *ier = KINSetResMonParams(KIN_kinmem, rval[0], rval[1]);
  else {
    *ier = -99;
    printf("FKINSETRIN: Unrecognized key.\n\n");
  }

}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SETVIN
 * ----------------------------------------------------------------
 */

void FKIN_SETVIN(char key_name[], realtype *vval, int *ier, int key_len)
{
  N_Vector Vec;

  if (!strncmp(key_name,"CONSTR_VEC", (size_t)key_len)) {
    Vec = NULL;
    Vec = N_VCloneEmpty(F2C_KINSOL_vec);
    if (Vec == NULL) {
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(vval, Vec);
    KINSetConstraints(KIN_kinmem, Vec);
    N_VDestroy(Vec);
  } else {
    *ier = -99;
    printf("FKINSETVIN: Unrecognized key.\n\n");
  }

}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_DENSE
 * ----------------------------------------------------------------
 */

void FKIN_DENSE(long int *neq, int *ier)
{
  *ier = KINDense(KIN_kinmem, *neq);

  KIN_ls = KIN_LS_DENSE;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BAND
 * ----------------------------------------------------------------
 */

void FKIN_BAND(long int *neq, long int *mupper, long int *mlower, int *ier)
{
  *ier = KINBand(KIN_kinmem, *neq, *mupper, *mlower);

  KIN_ls = KIN_LS_BAND;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SPTFQMR
 * ----------------------------------------------------------------
 */

void FKIN_SPTFQMR(int *maxl, int *ier)
{
  *ier = KINSptfqmr(KIN_kinmem, *maxl);

  KIN_ls = KIN_LS_SPTFQMR;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SPBCG
 * ----------------------------------------------------------------
 */

void FKIN_SPBCG(int *maxl, int *ier)
{
  *ier = KINSpbcg(KIN_kinmem, *maxl);

  KIN_ls = KIN_LS_SPBCG;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SPGMR
 * ----------------------------------------------------------------
 */

void FKIN_SPGMR(int *maxl, int *maxlrst, int *ier)
{
  *ier = KINSpgmr(KIN_kinmem, *maxl);
  KINSpilsSetMaxRestarts(KIN_kinmem, *maxlrst);

  KIN_ls = KIN_LS_SPGMR;

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

  uscalevec = NULL;
  uscalevec = N_VCloneEmpty(F2C_KINSOL_vec);
  if (uscalevec == NULL) {
    *ier = -4;  /* KIN_MEM_FAIL */
    return;
  }
  N_VSetArrayPointer(uscale, uscalevec);

  fscalevec = NULL;
  fscalevec = N_VCloneEmpty(F2C_KINSOL_vec);
  if (fscalevec == NULL) {
    N_VDestroy(uscalevec);
    *ier = -4;  /* KIN_MEM_FAIL */
    return;
  }
  N_VSetArrayPointer(fscale, fscalevec);

  /* Call main solver function */
  *ier = KINSol(KIN_kinmem, uuvec, *globalstrategy, uscalevec, fscalevec);

  N_VSetArrayPointer(NULL, uuvec);

  N_VSetArrayPointer(NULL, uscalevec);
  N_VDestroy(uscalevec);

  N_VSetArrayPointer(NULL, fscalevec);
  N_VDestroy(fscalevec);

  /* load optional outputs into iout[] and rout[] */
  KINGetWorkSpace(KIN_kinmem, &KIN_iout[0], &KIN_iout[1]);   /* LENRW & LENIW */
  KINGetNumNonlinSolvIters(KIN_kinmem, &KIN_iout[2]);        /* NNI */
  KINGetNumFuncEvals(KIN_kinmem, &KIN_iout[3]);              /* NFE */
  KINGetNumBetaCondFails(KIN_kinmem, &KIN_iout[4]);          /* NBCF */
  KINGetNumBacktrackOps(KIN_kinmem, &KIN_iout[5]);           /* NBCKTRK */

  KINGetFuncNorm(KIN_kinmem, &KIN_rout[0]);                  /* FNORM */
  KINGetStepLength(KIN_kinmem, &KIN_rout[1]);                /* SSTEP */

  switch(KIN_ls) {

  case KIN_LS_DENSE:
    KINDenseGetWorkSpace(KIN_kinmem, &KIN_iout[6], &KIN_iout[7]); /* LRW & LIW */
    KINDenseGetLastFlag(KIN_kinmem, (int *) &KIN_iout[8]);        /* LSTF */
    KINDenseGetNumFuncEvals(KIN_kinmem, &KIN_iout[9]);            /* NFE */
    KINDenseGetNumJacEvals(KIN_kinmem, &KIN_iout[10]);            /* NJE */
    
  case KIN_LS_BAND:
    KINBandGetWorkSpace(KIN_kinmem, &KIN_iout[6], &KIN_iout[7]);  /* LRW & LIW */
    KINBandGetLastFlag(KIN_kinmem, (int *) &KIN_iout[8]);         /* LSTF */
    KINBandGetNumFuncEvals(KIN_kinmem, &KIN_iout[9]);             /* NFE */
    KINBandGetNumJacEvals(KIN_kinmem, &KIN_iout[10]);             /* NJE */
    
  case KIN_LS_SPTFQMR:
  case KIN_LS_SPBCG:
  case KIN_LS_SPGMR:
    KINSpilsGetWorkSpace(KIN_kinmem, &KIN_iout[6], &KIN_iout[7]); /* LRW & LIW */
    KINSpilsGetLastFlag(KIN_kinmem, (int *) &KIN_iout[8]);        /* LSTF */
    KINSpilsGetNumFuncEvals(KIN_kinmem, &KIN_iout[9]);            /* NFE */
    KINSpilsGetNumJtimesEvals(KIN_kinmem, &KIN_iout[10]);         /* NJE */
    KINSpilsGetNumPrecEvals(KIN_kinmem, &KIN_iout[11]);           /* NPE */
    KINSpilsGetNumPrecSolves(KIN_kinmem, &KIN_iout[12]);          /* NPS */
    KINSpilsGetNumLinIters(KIN_kinmem, &KIN_iout[13]);            /* NLI */
    KINSpilsGetNumConvFails(KIN_kinmem, &KIN_iout[14]);           /* NCFL */
    break;

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

  KINFree(&KIN_kinmem);

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

int FKINfunc(N_Vector uu, N_Vector fval, void *f_data)
{
  realtype *udata, *fdata;
  int ier;

  udata = N_VGetArrayPointer(uu);
  fdata = N_VGetArrayPointer(fval);

  FK_FUN(udata, fdata, &ier);

  return(0);
}
