/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-07-05 16:00:44 $
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

kim_KINSOLdata kim_Kdata;  /* KINSOL data */
kim_MATLABdata kim_Mdata;  /* MATLAB data */

/*
 * ---------------------------------------------------------------------------------
 * static function prototypes
 * ---------------------------------------------------------------------------------
 */

static void KIM_init();
static void KIM_makePersistent();
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

  mexUnlock();

  switch(mode) {
  case 1:
    KIM_init();
    KIM_Malloc(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 2:
    KIM_Solve(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 3:
    KIM_Stats(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 4:
    KIM_Get(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 5:
    KIM_Set(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 6:
    KIM_Free(nlhs, plhs, nrhs-1, &prhs[1]);
    KIM_final();
    break;
  }

  KIM_makePersistent();
  mexLock();

  return;

}

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */

#define kin_mem    (kim_Kdata->kin_mem)
#define bbd_data   (kim_Kdata->bbd_data)
#define y          (kim_Kdata->y)
#define N          (kim_Kdata->N)
#define ls         (kim_Kdata->ls)
#define pm         (kim_Kdata->pm)

#define mx_data    (kim_Mdata->mx_data)

#define mx_SYSfct  (kim_Mdata->mx_SYSfct)
#define mx_JACfct  (kim_Mdata->mx_JACfct)
#define mx_PSETfct (kim_Mdata->mx_PSETfct)
#define mx_PSOLfct (kim_Mdata->mx_PSOLfct)
#define mx_GLOCfct (kim_Mdata->mx_GLOCfct)
#define mx_GCOMfct (kim_Mdata->mx_GCOMfct)

#define fig_handle (kim_Mdata->fig_handle)

/*
 * ---------------------------------------------------------------------------------
 * Exported procedures
 * ---------------------------------------------------------------------------------
 */

static void KIM_Malloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int i, status;
  double *tmp;

  mxArray *mx_in[3], *mx_out[2];

  int mxiter, msbset, msbsetsub, etachoice, mxnbcf;
  double eta, egamma, ealpha, mxnewtstep, relfunc, fnormtol, scsteptol;
  booleantype verbose, noInitSetup, noMinEps;

  double *constraints;
  N_Vector NVconstraints;

  int ptype;
  int mudq, mldq, mupper, mlower;
  int maxl, maxrs;
  double dqrely;

  mxArray *options;

  /* 
   * -----------------------------
   * Find out the vector type and
   * then pass it to the vector
   * library.
   * -----------------------------
   */

  /* Send vec_type and mx_comm */
  
  InitVectors();

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

  mx_SYSfct = mxDuplicateArray(prhs[0]);
  
  /* problem dimension */

  N = (int) mxGetScalar(prhs[1]);

  /* Solver Options -- optional argument */

  status = get_SolverOptions(prhs[2],
                             &verbose,
                             &mxiter, &msbset, &msbsetsub, &etachoice, &mxnbcf,
                             &eta, &egamma, &ealpha, &mxnewtstep, 
                             &relfunc, &fnormtol, &scsteptol,
                             &constraints,
                             &noInitSetup, &noMinEps);
                             

  /* User data -- optional argument */

  mx_data = mxDuplicateArray(prhs[3]);

  /* 
   * -----------------------------------------------------
   * Set solution vector (used as a template to KINMAlloc)
   * -----------------------------------------------------
   */

  y = NewVector(N);

  /* 
   * ----------------------------------------
   * Create kinsol object and allocate memory
   * ----------------------------------------
   */

  kin_mem = KINCreate();

  /* attach error handler function */
  status = KINSetErrHandlerFn(kin_mem, mtlb_KINErrHandler, NULL);

  if (verbose) {
    status = KINSetPrintLevel(kin_mem,3);
    /* attach info handler function */
    status = KINSetInfoHandlerFn(kin_mem, mtlb_KINInfoHandler, NULL);
    /* initialize the output window */
    mx_in[0] = mxCreateScalarDouble(0);
    mx_in[1] = mxCreateScalarDouble(0); /* ignored */
    mx_in[2] = mxCreateScalarDouble(0); /* ignored */
    mexCallMATLAB(1,mx_out,3,mx_in,"kim_info");
    fig_handle = (int)*mxGetPr(mx_out[0]);
  }

  /* Call KINMalloc */

  status = KINMalloc(kin_mem, mtlb_KINSys, y);

  /* Redirect output */
  status = KINSetErrFile(kin_mem, stdout);

  /* Optional inputs */

  status = KINSetNumMaxIters(kin_mem,mxiter);
  status = KINSetNoInitSetup(kin_mem,noInitSetup);
  status = KINSetNoMinEps(kin_mem,noMinEps);
  status = KINSetMaxSetupCalls(kin_mem,msbset);
  status = KINSetMaxSubSetupCalls(kin_mem,msbsetsub);
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

  status = get_LinSolvOptions(prhs[2],
                              &mupper, &mlower,
                              &mudq, &mldq, &dqrely,
                              &ptype, &maxrs, &maxl);

  switch (ls) {

  case LS_NONE:

    mexErrMsgTxt("KINMalloc:: no linear solver specified.");

    break;

  case LS_DENSE:

    status = KINDense(kin_mem, N);
    if (!mxIsEmpty(mx_JACfct))
      status = KINDenseSetJacFn(kin_mem, mtlb_KINDenseJac, NULL);

    break;

  case LS_BAND:

    status = KINBand(kin_mem, N, mupper, mlower);
    if (!mxIsEmpty(mx_JACfct))
        status = KINBandSetJacFn(kin_mem, mtlb_KINBandJac, NULL);

    break;

  case LS_SPGMR:

    switch(pm) {
    case PM_NONE:
      status = KINSpgmr(kin_mem, maxl);
      if (!mxIsEmpty(mx_PSOLfct)) {
        if (!mxIsEmpty(mx_PSETfct))
          status = KINSpilsSetPreconditioner(kin_mem, mtlb_KINSpilsPset, mtlb_KINSpilsPsol, NULL);
        else
          status = KINSpilsSetPreconditioner(kin_mem, NULL, mtlb_KINSpilsPsol, NULL);
      }
      break;
    case PM_BBDPRE:
      if (!mxIsEmpty(mx_GCOMfct))
        bbd_data = KINBBDPrecAlloc(kin_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_KINGloc, mtlb_KINGcom);
      else
        bbd_data = KINBBDPrecAlloc(kin_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_KINGloc, NULL);
      status = KINBBDSpgmr(kin_mem, maxl, bbd_data);
      break;
    }

    status = KINSpilsSetMaxRestarts(kin_mem, maxrs);

    if (!mxIsEmpty(mx_JACfct))
      status = KINSpilsSetJacTimesVecFn(kin_mem, mtlb_KINSpilsJac, NULL);

    break;

  case LS_SPBCG:

    switch(pm) {
    case PM_NONE:
      status = KINSpbcg(kin_mem, maxl);
      if (!mxIsEmpty(mx_PSOLfct)) {
        if (!mxIsEmpty(mx_PSETfct))
          status = KINSpilsSetPreconditioner(kin_mem, mtlb_KINSpilsPset, mtlb_KINSpilsPsol, NULL);
        else
          status = KINSpilsSetPreconditioner(kin_mem, NULL, mtlb_KINSpilsPsol, NULL);
      }
      break;
    case PM_BBDPRE:
      if (!mxIsEmpty(mx_GCOMfct))
        bbd_data = KINBBDPrecAlloc(kin_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_KINGloc, mtlb_KINGcom);
      else
        bbd_data = KINBBDPrecAlloc(kin_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_KINGloc, NULL);
      status = KINBBDSpbcg(kin_mem, maxl, bbd_data);
      break;
    }

    if (!mxIsEmpty(mx_JACfct))
      status = KINSpilsSetJacTimesVecFn(kin_mem, mtlb_KINSpilsJac, NULL);

    break;

  case LS_SPTFQMR:

    switch(pm) {
    case PM_NONE:
      status = KINSptfqmr(kin_mem, maxl);
      if (!mxIsEmpty(mx_PSOLfct)) {
        if (!mxIsEmpty(mx_PSETfct))
          status = KINSpilsSetPreconditioner(kin_mem, mtlb_KINSpilsPset, mtlb_KINSpilsPsol, NULL);
        else
          status = KINSpilsSetPreconditioner(kin_mem, NULL, mtlb_KINSpilsPsol, NULL);
      }
      break;
    case PM_BBDPRE:
      if (!mxIsEmpty(mx_GCOMfct))
        bbd_data = KINBBDPrecAlloc(kin_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_KINGloc, mtlb_KINGcom);
      else
        bbd_data = KINBBDPrecAlloc(kin_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_KINGloc, NULL);
      status = KINBBDSptfqmr(kin_mem, maxl, bbd_data);
      break;
    }

    if (!mxIsEmpty(mx_JACfct))
      status = KINSpilsSetJacTimesVecFn(kin_mem, mtlb_KINSpilsJac, NULL);

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

  /* Exract y0 and load initial guess in y */
  y0 = mxGetPr(prhs[0]);
  PutData(y, y0, N);

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

  /* call KINSol() */
  status = KINSol(kin_mem, y, strategy, yscale, fscale);

  /* KINSOL return flag */
  plhs[0] = mxCreateScalarDouble((double)status);

  /* Solution vector */
  plhs[1] = mxCreateDoubleMatrix(N,1,mxREAL);
  GetData(y, mxGetPr(plhs[1]), N);

  /* Free temporary N_Vectors */
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
  case LS_SPBCG:
  case LS_SPTFQMR:

    flag = KINSpilsGetNumLinIters(kin_mem, &nli);
    flag = KINSpilsGetNumPrecEvals(kin_mem, &npe);
    flag = KINSpilsGetNumPrecSolves(kin_mem, &nps);
    flag = KINSpilsGetNumConvFails(kin_mem, &ncfl);
    
    nfields = sizeof(fnames_spils)/sizeof(*fnames_spils);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_spils);

    if (ls == LS_SPGMR)
      mxSetField(mx_ls, 0, "name",  mxCreateString("GMRES"));
    else if (ls == LS_SPBCG)
      mxSetField(mx_ls, 0, "name",  mxCreateString("BiCGStab"));
    else
      mxSetField(mx_ls, 0, "name",  mxCreateString("TFQMR"));

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

  /* Allocate space for global KINSOL and MATLAB data structures */
  
  kim_Kdata = (kim_KINSOLdata) mxMalloc(sizeof(struct kim_KINSOLdataStruct));
  kim_Mdata = (kim_MATLABdata) mxMalloc(sizeof(struct kim_MATLABdataStruct));

  /* Initialize global KINSOL data */

  kin_mem  = NULL;
  bbd_data = NULL;

  y = NULL;

  N = 0;

  ls  = LS_DENSE;
  pm  = PM_NONE;

  /* Initialize global MATLAB data */

  empty = mxCreateDoubleMatrix(0,0,mxREAL);

  mx_data     = mxDuplicateArray(empty);

  mx_SYSfct   = mxDuplicateArray(empty);
  mx_JACfct   = mxDuplicateArray(empty);
  mx_PSETfct  = mxDuplicateArray(empty);
  mx_PSOLfct  = mxDuplicateArray(empty);
  mx_GLOCfct  = mxDuplicateArray(empty);
  mx_GCOMfct  = mxDuplicateArray(empty);
  
  mxDestroyArray(empty);

  return;
}

static void KIM_makePersistent() {

  mexMakeArrayPersistent(mx_data);
  mexMakeArrayPersistent(mx_SYSfct);
  mexMakeArrayPersistent(mx_JACfct);
  mexMakeArrayPersistent(mx_PSETfct);
  mexMakeArrayPersistent(mx_PSOLfct);
  mexMakeArrayPersistent(mx_GLOCfct);
  mexMakeArrayPersistent(mx_GCOMfct);

  mexMakeMemoryPersistent(kim_Kdata);
  mexMakeMemoryPersistent(kim_Mdata);

}


static void KIM_final()
{
  mxDestroyArray(mx_data);
  
  mxDestroyArray(mx_SYSfct);
  mxDestroyArray(mx_JACfct);
  mxDestroyArray(mx_PSETfct);
  mxDestroyArray(mx_PSOLfct);
  mxDestroyArray(mx_GLOCfct);
  mxDestroyArray(mx_GCOMfct);

  mxFree(kim_Kdata);
  mxFree(kim_Mdata);

}
