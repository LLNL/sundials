/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-01-06 19:00:23 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * MEX implementation for KINSOL Matlab interface.
 * -----------------------------------------------------------------
 */

#include <string.h>
#include <stdlib.h>
#include "kim.h"
#include "nvm.h"

/*
 * ---------------------------------------------------------------------------------
 * Definitions for global variables (declared in kim.h)
 * ---------------------------------------------------------------------------------
 */

/* KINSOL data */

void *kin_mem;     /* KINSOL solver memory */
void *bbd_data;    /* BBD preconditioner data */
N_Vector y;        /* solution vector */
int N;             /* problem dimension */
int ls;            /* linear solver type */
int pm;            /* preconditioner module */

/* Matlab data */

mxArray *mx_mtlb_SYSfct;
mxArray *mx_mtlb_JACfct;
mxArray *mx_mtlb_PSETfct;
mxArray *mx_mtlb_PSOLfct;
mxArray *mx_mtlb_GLOCfct;
mxArray *mx_mtlb_GCOMfct;

mxArray *mx_mtlb_data;

/*
 * ---------------------------------------------------------------------------------
 * static function prototypes
 * ---------------------------------------------------------------------------------
 */

static void KIM_init();
static void KIM_final();

static void KIM_Malloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static void KIM_Solve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static void KIM_Stats(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static void KIM_Get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static void KIM_Set(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static void KIM_Free(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/*
 * ---------------------------------------------------------------------------------
 * Main entry point
 * ---------------------------------------------------------------------------------
 */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] )
{
  int mode;
  /* 
     Modes:
     
     1 - initialize KINSOL solver
     2 - solve problem
     3 - get solver stats
     4 - extract data from kin_mem
     5 - set one optional input at a time
     6 - finalize
  */

  mode = (int)mxGetScalar(prhs[0]);
  switch(mode) {
  case 1:
    KIM_init();
    KIM_Malloc(nlhs, plhs,nrhs-1,&prhs[1]);
    mexLock();
    mexMakeArrayPersistent(mx_mtlb_data);
    mexMakeArrayPersistent(mx_mtlb_SYSfct);
    mexMakeArrayPersistent(mx_mtlb_JACfct);
    mexMakeArrayPersistent(mx_mtlb_PSETfct);
    mexMakeArrayPersistent(mx_mtlb_PSOLfct);
    mexMakeArrayPersistent(mx_mtlb_GLOCfct);
    mexMakeArrayPersistent(mx_mtlb_GCOMfct);
    break;
  case 2:
    KIM_Solve(nlhs, plhs,nrhs-1,&prhs[1]);
    break;
  case 3:
    KIM_Stats(nlhs, plhs,nrhs-1,&prhs[1]);
    break;
  case 4:
    KIM_Get(nlhs, plhs,nrhs-1,&prhs[1]);
    break;
  case 5:
    KIM_Set(nlhs, plhs,nrhs-1,&prhs[1]);
    break;
  case 6:
    KIM_Free(nlhs, plhs,nrhs-1,&prhs[1]);
    KIM_final();
    break;
  }

}

/*
 * ---------------------------------------------------------------------------------
 * Exported procedures
 * ---------------------------------------------------------------------------------
 */

static void KIM_Malloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int i, status;
  double *tmp;

  int vec_type;
  mxArray *mx_in[1], *mx_out[2];
  mxArray *mx_comm;

  int mxiter, msbset, etachoice, mxnbcf;
  double eta, egamma, ealpha, mxnewtstep, relfunc, fnormtol, scsteptol;
  booleantype noInitSetup, noMinEps;

  double *constraints;
  N_Vector NVconstraints;

  int ptype;
  int mudq, mldq, mukeep, mlkeep;
  int maxl, maxrs;
  double dq;

  mxArray *options;

  /* 
   * -----------------------------
   * Find out the vector type and
   * then pass it to the vector
   * library.
   * -----------------------------
   */

  /* Call the MEX file nvm with mode=3 */
  
  mx_in[0] = mxCreateScalarDouble(3);
  mexCallMATLAB(2,mx_out,1,mx_in,"nvm");
  vec_type = (int)*mxGetPr(mx_out[0]);
  mx_comm = mxDuplicateArray(mx_out[1]);

  /* Send vec_type and mx_comm */
  
  InitVectors(vec_type, mx_comm);

  /* 
   * -----------------------------
   * Extract stuff from arguments:
   * - SYS function
   * - problem dimension
   * - solver options
   * - user data
   * -----------------------------
   */

  /* Matlab user-provided function */

  mx_mtlb_SYSfct = mxDuplicateArray(prhs[0]);
  
  /* problem dimension */

  N = (int) mxGetScalar(prhs[1]);

  /* Solver Options -- optional argument */

  status = get_SolverOptions(prhs[2],
                             &mxiter, &msbset, &etachoice, &mxnbcf,
                             &eta, &egamma, &ealpha, &mxnewtstep, 
                             &relfunc, &fnormtol, &scsteptol,
                             &constraints,
                             &noInitSetup, &noMinEps);
                             

  /* User data -- optional argument */

  mx_mtlb_data = mxDuplicateArray(prhs[3]);

  /* 
   * -----------------------------------------------------
   * Set solution vector (used as a template to KINMAlloc)
   * -----------------------------------------------------
   */

  y = NewVector(N);

  /* 
   * ---------------------------------------
   * Create kinsol object and allocate memory
   * ---------------------------------------
   */

  kin_mem = KINCreate();

  /* Call KINMalloc */

  status = KINMalloc(kin_mem, mtlb_KINSys, y);

  /* Redirect output */
  status = KINSetErrFile(kin_mem, stdout);

  /* Optional inputs */

  status = KINSetNumMaxIters(kin_mem,mxiter);
  status = KINSetNoInitSetup(kin_mem,noInitSetup);
  status = KINSetNoMinEps(kin_mem,noMinEps);
  status = KINSetMaxSetupCalls(kin_mem,msbset);
  status = KINSetMaxBetaFails(kin_mem,mxnbcf);
  status = KINSetEtaForm(kin_mem,etachoice);
  status = KINSetEtaConstValue(kin_mem,eta);
  status = KINSetEtaParams(kin_mem,egamma,ealpha);
  status = KINSetMaxNewtonStep(kin_mem,mxnewtstep);
  status = KINSetRelErrFunc(kin_mem,relfunc);
  status = KINSetFuncNormTol(kin_mem,fnormtol);
  status = KINSetScaledStepTol(kin_mem,scsteptol);
  if (constraints != NULL) {
    NVconstraints = N_VCloneEmpty(y);
    N_VSetArrayPointer(constraints, NVconstraints);
    status = KINSetConstraints(kin_mem,NVconstraints);
    N_VDestroy(NVconstraints);
  }

  status = get_LinSolvOptions(prhs[2], &ls, 
                              &maxl, &maxrs, &ptype, &pm,
                              &mudq, &mldq, &mukeep, &mlkeep, &dq,
                              &mx_mtlb_JACfct, &mx_mtlb_PSETfct, &mx_mtlb_PSOLfct,
                              &mx_mtlb_GLOCfct, &mx_mtlb_GCOMfct);

  switch (ls) {
  case LS_NONE:
    mexErrMsgTxt("KINMalloc:: no linear solver specified.");
  case LS_DENSE:
    status = KINDense(kin_mem, N);
    if (!mxIsEmpty(mx_mtlb_JACfct))
      status = KINDenseSetJacFn(kin_mem, mtlb_KINDenseJac, NULL);
    break;
  case LS_SPGMR:

    if (pm == PM_NONE) {
      status = KINSpgmr(kin_mem, maxl);
      if (!mxIsEmpty(mx_mtlb_PSOLfct)) {
        if (!mxIsEmpty(mx_mtlb_PSETfct))
          status = KINSpgmrSetPreconditioner(kin_mem, mtlb_KINSpilsPset, mtlb_KINSpilsPsol, NULL);
        else
          status = KINSpgmrSetPreconditioner(kin_mem, NULL, mtlb_KINSpilsPsol, NULL);
      }
    } else if (pm == PM_BBDPRE) {
      bbd_data = KINBBDPrecAlloc(kin_mem, N, mudq, mldq, mukeep, mlkeep, dq, mtlb_KINGloc, mtlb_KINGcom);
      status = KINBBDSpgmr(kin_mem, maxl, bbd_data);
    }
    status = KINSpgmrSetMaxRestarts(kin_mem, maxrs);
    if (!mxIsEmpty(mx_mtlb_JACfct))
      status = KINSpgmrSetJacTimesVecFn(kin_mem, mtlb_KINSpilsJac, NULL);
    break;
  case LS_SPBCG:
    if (pm == PM_NONE) {
      status = KINSpbcg(kin_mem, maxl);
      if (!mxIsEmpty(mx_mtlb_PSOLfct)) {
        if (!mxIsEmpty(mx_mtlb_PSETfct))
          status = KINSpbcgSetPreconditioner(kin_mem, mtlb_KINSpilsPset, mtlb_KINSpilsPsol, NULL);
        else
          status = KINSpbcgSetPreconditioner(kin_mem, NULL, mtlb_KINSpilsPsol, NULL);
      }
    } else if (pm == PM_BBDPRE) {
      bbd_data = KINBBDPrecAlloc(kin_mem, N, mudq, mldq, mukeep, mlkeep, dq, mtlb_KINGloc, mtlb_KINGcom);
      status = KINBBDSpbcg(kin_mem, maxl, bbd_data);
    }
    if (!mxIsEmpty(mx_mtlb_JACfct))
      status = KINSpbcgSetJacTimesVecFn(kin_mem, mtlb_KINSpilsJac, NULL);
    break;
  }
  

}

static void KIM_Solve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  

  double *y0, *ys, *fs;
  N_Vector yscale, fscale;
  int buflen, status, strategy;
  char *bufval;
 
  int i;
  double *tmp;

  /* Exract y0 */

  y0 = mxGetPr(prhs[0]);

  /* Extract strategy */

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"None")) strategy = KIN_NONE;
  if(!strcmp(bufval,"LineSearch")) strategy = KIN_LINESEARCH;

  /* Extract yscale */

  ys = mxGetPr(prhs[2]);
  yscale = N_VCloneEmpty(y);
  N_VSetArrayPointer(ys, yscale);

  /* Extract fscale */

  fs = mxGetPr(prhs[3]);
  fscale = N_VCloneEmpty(y);
  N_VSetArrayPointer(fs, fscale);

  /* Call KINSol */

  /* Solution vector */
  plhs[1] = mxCreateDoubleMatrix(N,1,mxREAL);
  N_VSetArrayPointer(mxGetPr(plhs[1]), y);

  PutData(y, y0, N);

  status = KINSol(kin_mem, y, strategy, yscale, fscale);

  /* KINSOL return flag */
  plhs[0] = mxCreateScalarDouble((double)status);

  N_VDestroy(yscale);
  N_VDestroy(fscale);

}

static void KIM_Stats(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const char *fnames_solver[]={
    "nfe",
    "nni",
    "nbcf",
    "nbops",
    "fnorm",
    "step",
    "LSInfo",
  };
  const char *fnames_dense[]={
    "name",
    "njeD",
    "nfeD"
  };
  const char *fnames_spils[]={
    "name",
    "nli",
    "npe",
    "nps",
    "ncfl",
  };

  long int nfe, nni, nbcf, nbops;
  double fnorm, step;

  long int njeD, nfeD;
  long int nli, npe, nps, ncfl;

  int flag;
  mxArray *mx_ls;
  int nfields;

  flag = KINGetNumNonlinSolvIters(kin_mem, &nni);
  flag = KINGetNumFuncEvals(kin_mem, &nfe);
  flag = KINGetNumBetaCondFails(kin_mem, &nbcf);
  flag = KINGetNumBacktrackOps(kin_mem, &nbops);
  
  flag = KINGetFuncNorm(kin_mem, &fnorm);
  flag = KINGetStepLength(kin_mem, &step);

  nfields = sizeof(fnames_solver)/sizeof(*fnames_solver);
  plhs[0] = mxCreateStructMatrix(1, 1, nfields, fnames_solver);
 
  mxSetField(plhs[0], 0, "nfe",   mxCreateScalarDouble((double)nfe));
  mxSetField(plhs[0], 0, "nni",   mxCreateScalarDouble((double)nni));
  mxSetField(plhs[0], 0, "nbcf",  mxCreateScalarDouble((double)nbcf));
  mxSetField(plhs[0], 0, "nbops", mxCreateScalarDouble((double)nbops));
  mxSetField(plhs[0], 0, "fnorm", mxCreateScalarDouble(fnorm));
  mxSetField(plhs[0], 0, "step",  mxCreateScalarDouble(step));

  /* Linear Solver Statistics */

  switch(ls){

  case LS_DENSE:
    
    flag = KINDenseGetNumJacEvals(kin_mem, &njeD);
    flag = KINDenseGetNumFuncEvals(kin_mem, &nfeD);
    
    nfields = sizeof(fnames_dense)/sizeof(*fnames_dense);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_dense);

    mxSetField(mx_ls, 0, "name", mxCreateString("Dense"));
    mxSetField(mx_ls, 0, "njeD", mxCreateScalarDouble((double)njeD));
    mxSetField(mx_ls, 0, "nfeD", mxCreateScalarDouble((double)nfeD));
    
    break;

  case LS_SPGMR:

    flag = KINSpgmrGetNumLinIters(kin_mem, &nli);
    flag = KINSpgmrGetNumPrecEvals(kin_mem, &npe);
    flag = KINSpgmrGetNumPrecSolves(kin_mem, &nps);
    flag = KINSpgmrGetNumConvFails(kin_mem, &ncfl);
    
    nfields = sizeof(fnames_spils)/sizeof(*fnames_spils);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_spils);
 
    mxSetField(mx_ls, 0, "name",  mxCreateString("GMRES"));
    mxSetField(mx_ls, 0, "nli",   mxCreateScalarDouble((double)nli));
    mxSetField(mx_ls, 0, "npe",   mxCreateScalarDouble((double)npe));
    mxSetField(mx_ls, 0, "nps",   mxCreateScalarDouble((double)nps));
    mxSetField(mx_ls, 0, "ncfl",  mxCreateScalarDouble((double)ncfl));
    
    break;

  case LS_SPBCG:

    flag = KINSpbcgGetNumLinIters(kin_mem, &nli);
    flag = KINSpbcgGetNumPrecEvals(kin_mem, &npe);
    flag = KINSpbcgGetNumPrecSolves(kin_mem, &nps);
    flag = KINSpbcgGetNumConvFails(kin_mem, &ncfl);
    
    nfields = sizeof(fnames_spils)/sizeof(*fnames_spils);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_spils);
 
    mxSetField(mx_ls, 0, "name",  mxCreateString("BiCGStab"));
    mxSetField(mx_ls, 0, "nli",   mxCreateScalarDouble((double)nli));
    mxSetField(mx_ls, 0, "npe",   mxCreateScalarDouble((double)npe));
    mxSetField(mx_ls, 0, "nps",   mxCreateScalarDouble((double)nps));
    mxSetField(mx_ls, 0, "ncfl",  mxCreateScalarDouble((double)ncfl));
    
    break;
    
  }

  mxSetField(plhs[0], 0, "LSInfo", mx_ls);

}


static void KIM_Set(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
}

static void KIM_Get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
}

static void KIM_Free(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  N_VDestroy(y);
  if (pm == PM_BBDPRE) KINBBDPrecFree(bbd_data);

  KINFree(kin_mem);

}


static void KIM_init()
{
  mxArray *empty;

  /* Initialize global dimensions to zero */
  N = 0;

  empty = mxCreateDoubleMatrix(0,0,mxREAL);

  mx_mtlb_data     = mxDuplicateArray(empty);

  mx_mtlb_SYSfct   = mxDuplicateArray(empty);
  mx_mtlb_JACfct   = mxDuplicateArray(empty);
  mx_mtlb_PSETfct  = mxDuplicateArray(empty);
  mx_mtlb_PSOLfct  = mxDuplicateArray(empty);
  mx_mtlb_GLOCfct  = mxDuplicateArray(empty);
  mx_mtlb_GCOMfct  = mxDuplicateArray(empty);
  
  mxDestroyArray(empty);

}

static void KIM_final()
{
  mexUnlock();

  mxDestroyArray(mx_mtlb_data);
  
  mxDestroyArray(mx_mtlb_SYSfct);
  mxDestroyArray(mx_mtlb_JACfct);
  mxDestroyArray(mx_mtlb_PSETfct);
  mxDestroyArray(mx_mtlb_PSOLfct);
  mxDestroyArray(mx_mtlb_GLOCfct);
  mxDestroyArray(mx_mtlb_GCOMfct);

}
