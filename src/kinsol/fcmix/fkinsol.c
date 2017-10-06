/*
 * -----------------------------------------------------------------
 * $Revision: 4906 $
 * $Date: 2016-09-14 16:05:17 -0700 (Wed, 14 Sep 2016) $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
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

#include "fkinsol.h"               /* prototypes of interfaces and global vars.   */
#include "kinsol_impl.h"           /* definition of KINMem type                   */

#include <kinsol/kinsol_band.h>    /* prototypes of KINBAND interface routines    */
#include <kinsol/kinsol_dense.h>   /* prototypes of KINDENSE interface routines   */
#include <kinsol/kinsol_klu.h>     /* prototypes of KINKLU interface routines     */
#include <kinsol/kinsol_superlumt.h> /* prototypes of KINSUPERLU interface routines */
#include <kinsol/kinsol_sptfqmr.h> /* prototypes of KINSPTFQMR interface routines */
#include <kinsol/kinsol_spbcgs.h>  /* prototypes of KINSPBCG interface routines   */
#include <kinsol/kinsol_spgmr.h>   /* prototypes of KINSPGMR interface routines   */
#include <kinsol/kinsol_spfgmr.h>  /* prototypes of KINSPFGMR interface routines  */

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
 * Function : FKIN_CREATE
 * ----------------------------------------------------------------
 */

void FKIN_CREATE(int *ier)
{
  
  *ier = 0;
  /* check for required vector operations */
  if ((F2C_KINSOL_vec->ops->nvgetarraypointer == NULL) ||
      (F2C_KINSOL_vec->ops->nvsetarraypointer == NULL)) {
    *ier = -1;
    printf("FKINCREATE: A required vector operation is not implemented.\n\n");
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
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_INIT
 * ----------------------------------------------------------------
 */

void FKIN_INIT(long int *iout, realtype *rout, int *ier)
{
  
  /* Call KINInit */
  *ier = 0;
  *ier = KINInit(KIN_kinmem, FKINfunc, F2C_KINSOL_vec);

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

  /* Call KINInit */
  *ier = 0;
  *ier = KINInit(KIN_kinmem, FKINfunc, F2C_KINSOL_vec);

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

void FKIN_SETIIN(char key_name[], long int *ival, int *ier)
{
  if (!strncmp(key_name,"PRNT_LEVEL",10))
    *ier = KINSetPrintLevel(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NITERS",10))
    *ier = KINSetNumMaxIters(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"ETA_FORM",8))
    *ier = KINSetEtaForm(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"MAA",3))
    *ier = KINSetMAA(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"MAX_SETUPS",10))
    *ier = KINSetMaxSetupCalls(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"MAX_SP_SETUPS",13))
    *ier = KINSetMaxSubSetupCalls(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"NO_INIT_SETUP",13))
    *ier = KINSetNoInitSetup(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"NO_MIN_EPS",10))
    *ier = KINSetNoMinEps(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"NO_RES_MON",10))
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

void FKIN_SETRIN(char key_name[], realtype *rval, int *ier)
{

  if (!strncmp(key_name,"FNORM_TOL",9))
    *ier = KINSetFuncNormTol(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"SSTEP_TOL",9))
    *ier = KINSetScaledStepTol(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"MAX_STEP",8))
    *ier = KINSetMaxNewtonStep(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"RERR_FUNC",9))
    *ier = KINSetRelErrFunc(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"ETA_CONST",9))
    *ier = KINSetEtaConstValue(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"ETA_PARAMS",10))
    *ier = KINSetEtaParams(KIN_kinmem, rval[0], rval[1]);
  else if (!strncmp(key_name,"RMON_CONST",10))
    *ier = KINSetResMonConstValue(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"RMON_PARAMS",11))
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

void FKIN_SETVIN(char key_name[], realtype *vval, int *ier)
{
  N_Vector Vec;

  if (!strncmp(key_name,"CONSTR_VEC",10)) {
    Vec = NULL;
    Vec = N_VCloneEmpty(F2C_KINSOL_vec);
    if (Vec == NULL) {
      *ier = -1;
      return;
    }
    *ier = 0;
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
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SPFGMR
 * ----------------------------------------------------------------
 */

void FKIN_SPFGMR(int *maxl, int *maxlrst, int *ier)
{
  *ier = KINSpfgmr(KIN_kinmem, *maxl);
  KINSpilsSetMaxRestarts(KIN_kinmem, *maxlrst);
  KIN_ls = KIN_LS_SPFGMR;
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
  case KIN_LS_BAND:
  case KIN_LS_LAPACKDENSE:
  case KIN_LS_LAPACKBAND:
    KINDlsGetWorkSpace(KIN_kinmem, &KIN_iout[6], &KIN_iout[7]); /* LRW & LIW */
    KINDlsGetLastFlag(KIN_kinmem, &KIN_iout[8]);                /* LSTF */
    KINDlsGetNumFuncEvals(KIN_kinmem, &KIN_iout[9]);            /* NFE */
    KINDlsGetNumJacEvals(KIN_kinmem, &KIN_iout[10]);            /* NJE */    
  case KIN_LS_KLU:
  case KIN_LS_SUPERLUMT:
    KINSlsGetLastFlag(KIN_kinmem, &KIN_iout[8]);                /* LSTF  */
    KINSlsGetNumJacEvals(KIN_kinmem, &KIN_iout[10]);            /* NJE   */
    break;
  case KIN_LS_SPTFQMR:
  case KIN_LS_SPBCG:
  case KIN_LS_SPFGMR:
  case KIN_LS_SPGMR:
    KINSpilsGetWorkSpace(KIN_kinmem, &KIN_iout[6], &KIN_iout[7]); /* LRW & LIW */
    KINSpilsGetLastFlag(KIN_kinmem, &KIN_iout[8]);                /* LSTF */
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

int FKINfunc(N_Vector uu, N_Vector fval, void *user_data)
{
  realtype *udata, *fdata;
  int ier;

  udata = N_VGetArrayPointer(uu);
  fdata = N_VGetArrayPointer(fval);

  FK_FUN(udata, fdata, &ier);

  return(ier);
}
