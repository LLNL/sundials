/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
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
 * Global interface data variable
 * ---------------------------------------------------------------------------------
 */

kimInterfaceData kimData = NULL;

/*
 * ---------------------------------------------------------------------------------
 * Static function prototypes
 * ---------------------------------------------------------------------------------
 */

static void kimInitKINSOLdata();
static void kimPersistKINSOLdata();
static void kimFinalKINSOLdata();

static int KIM_Initialization(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int KIM_Solve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int KIM_Stats(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int KIM_Get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int KIM_Set(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int KIM_Free(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

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

  if ( (mode != 1) && (kimData == NULL) ) {
    mexErrMsgTxt("KINSOL - Illegal attempt to call before KINInit.");
  }


  switch(mode) {

  case 1:
    if (kimData != NULL) {
      KIM_Free(nlhs, plhs, nrhs-1, &prhs[1]);
      kimFinalKINSOLdata();
    }
    kimInitKINSOLdata();
    KIM_Initialization(nlhs, plhs, nrhs-1, &prhs[1]);
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
    kimFinalKINSOLdata();
    break;

  }

  /* Unless this was the KINFree call,
   * make data persistent and lock the MEX file */
  if (mode != 6) {
    kimPersistKINSOLdata();
    mexLock();
  }

  return;
}

/*
 * ---------------------------------------------------------------------------------
 * Private functions to deal with the global data
 * ---------------------------------------------------------------------------------
 */

static void kimInitKINSOLdata()
{
  mxArray *empty;

  /* Allocate space for global KINSOL data structure */

  kimData = (kimInterfaceData) mxMalloc(sizeof(struct kimInterfaceData_));

  /* Initialize global KINSOL data */

  kimData->kin_mem = NULL;

  kimData->n  = 0;

  kimData->Y  = NULL;

  kimData->LS = LS_DENSE;
  kimData->PM = PM_NONE;

  kimData->errMsg = SUNTRUE;

  /* Initialize Matlab mex arrays to empty */

  empty = mxCreateDoubleMatrix(0,0,mxREAL);

  kimData->SYSfct   = mxDuplicateArray(empty);
  kimData->JACfct   = mxDuplicateArray(empty);
  kimData->PSETfct  = mxDuplicateArray(empty);
  kimData->PSOLfct  = mxDuplicateArray(empty);
  kimData->GLOCfct  = mxDuplicateArray(empty);
  kimData->GCOMfct  = mxDuplicateArray(empty);

  kimData->mtlb_data = mxDuplicateArray(empty);

  mxDestroyArray(empty);

  return;
}

static void kimPersistKINSOLdata()
{
  /* Make global memory persistent */

  mexMakeArrayPersistent(kimData->mtlb_data);

  mexMakeArrayPersistent(kimData->SYSfct);
  mexMakeArrayPersistent(kimData->JACfct);
  mexMakeArrayPersistent(kimData->PSETfct);
  mexMakeArrayPersistent(kimData->PSOLfct);
  mexMakeArrayPersistent(kimData->GLOCfct);
  mexMakeArrayPersistent(kimData->GCOMfct);

  mexMakeMemoryPersistent(kimData);

  return;
}

static void kimFinalKINSOLdata()
{  
  if (kimData == NULL) return;

  if (kimData->Y != NULL) N_VDestroy(kimData->Y);

  mxDestroyArray(kimData->mtlb_data);

  mxDestroyArray(kimData->SYSfct);
  mxDestroyArray(kimData->JACfct);
  mxDestroyArray(kimData->PSETfct);
  mxDestroyArray(kimData->PSOLfct);
  mxDestroyArray(kimData->GLOCfct);
  mxDestroyArray(kimData->GCOMfct);

  mxFree(kimData);
  kimData = NULL;

  return;
}

/*
 * ---------------------------------------------------------------------------------
 * Error handler function.
 *
 * This function is both passed as the KINSOL error handler and used throughout
 * the Matlab interface.
 *
 * If called directly by one of the interface functions, error_code = -999 to
 * indicate an error and err_code = +999 to indicate a warning. Otherwise,
 * err_code is set by the calling KINSOL function.
 *
 * NOTE: mexErrMsgTxt will end the execution of the MEX file. Therefore we do
 *       not have to intercept any of the KINSOL error return flags.
 *       The only return flags we intercept are those from CVode() and CVodeB()
 *       which are passed back to the user (only positive values will make it). 
 * ---------------------------------------------------------------------------------
 */

void kimErrHandler(int error_code, 
                   const char *module, const char *function, 
                   char *msg, void *eh_data)
{
  char err_msg[256];

  if (!(kimData->errMsg)) return;

  if (error_code > 0) {
    sprintf(err_msg,"Warning in ==> %s\n%s",function,msg);
    mexWarnMsgTxt(err_msg);    
  } else if (error_code < 0) {
    /*mexUnlock();
      kimFinalKINSOLdata();*/
    sprintf(err_msg,"Error using ==> %s\n%s",function,msg);
    mexErrMsgTxt(err_msg);
  }

  return;
}

/*
 * ---------------------------------------------------------------------------------
 * Info handler function
 * 
 * This function is passed as the KINSOL info handler if verbose output was
 * requested. It is a wrapper around the Matlab m-file kim_info.m which posts
 * all info messages to a separate Matlab figure.
 * ---------------------------------------------------------------------------------
 */

void kimInfoHandler(const char *module, const char *function, 
                    char *msg, void *ih_data)
{
  char my_msg[400];
  mxArray *mx_in[3];

  sprintf(my_msg,"[%s] %s\n  %s\n",module,function,msg);

  /* action=1 -> append */
  mx_in[0] = mxCreateDoubleScalar(1);
  mx_in[1] = mxCreateDoubleScalar((double)kimData->fig_handle);
  mx_in[2] = mxCreateString(my_msg);

  mexCallMATLAB(0,NULL,3,mx_in,"kim_info");

}

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */

#define kin_mem      (kimData->kin_mem)

#define y            (kimData->Y)
#define N            (kimData->n)
#define ls           (kimData->LS)
#define pm           (kimData->PM)

#define mtlb_data    (kimData->mtlb_data)

#define mtlb_SYSfct  (kimData->SYSfct)
#define mtlb_JACfct  (kimData->JACfct)
#define mtlb_PSETfct (kimData->PSETfct)
#define mtlb_PSOLfct (kimData->PSOLfct)
#define mtlb_GLOCfct (kimData->GLOCfct)
#define mtlb_GCOMfct (kimData->GCOMfct)

/*
 * ---------------------------------------------------------------------------------
 * Interface procedures
 * ---------------------------------------------------------------------------------
 */

static int KIM_Initialization(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray *options;
  mxArray *mx_in[3], *mx_out[2];

  int status;

  int mxiter, msbset, msbsetsub, etachoice, mxnbcf;
  double eta, egamma, ealpha, mxnewtstep, relfunc, fnormtol, scsteptol;
  booleantype verbose, errmsg, noInitSetup, noMinEps;

  double *constraints;
  N_Vector NVconstraints;

  int ptype;
  sunindextype mudq, mldq, mupper, mlower;
  int maxl, maxrs;
  double dqrely;

  /*
   * ------------------------------------
   * Initialize appropriate vector module
   * ------------------------------------
   */

  InitVectors();

  /* 
   * ----------------------------
   * Extract stuff from arguments
   * ----------------------------
   */

  /* Matlab user-provided function */

  mxDestroyArray(mtlb_SYSfct);
  mtlb_SYSfct = mxDuplicateArray(prhs[0]);
  
  /* Problem dimension */

  N = (sunindextype) mxGetScalar(prhs[1]);

  /* Solver Options (may be empty) */

  options = prhs[2];

  /* Create the solution N_Vector */

  y = NewVector(N);

  /*
   * ----------------------------- 
   * Process the options structure 
   * -----------------------------
   */

  status = get_SolverOptions(options,
                             &verbose, &errmsg,
                             &mxiter, &msbset, &msbsetsub, &etachoice, &mxnbcf,
                             &eta, &egamma, &ealpha, &mxnewtstep, 
                             &relfunc, &fnormtol, &scsteptol,
                             &constraints,
                             &noInitSetup, &noMinEps);
  if (status != 0) goto error_return;     

  /* 
   * ------------------------------------------------------
   * Create KINSOL object, allocate memory, and set options
   * ------------------------------------------------------
   */

  kin_mem = KINCreate();
  if (kin_mem == NULL) goto error_return;

  /* Attach the global KINSOL data as 'user_data' */
  status = KINSetUserData(kin_mem, kimData);
  if (status != KIN_SUCCESS) goto error_return;

  /* Attach error handler function */
  status = KINSetErrHandlerFn(kin_mem, kimErrHandler, NULL);
  if (status != KIN_SUCCESS) goto error_return;

  /* If verbose was set to SUNTRUE */
  if (verbose) {
    /* Set print level to its highest value */
    status = KINSetPrintLevel(kin_mem,3);
    if (status != KIN_SUCCESS) goto error_return;
    /* Attach info handler function */
    status = KINSetInfoHandlerFn(kin_mem, kimInfoHandler, NULL);
    if (status != KIN_SUCCESS) goto error_return;
    /* Initialize the output window and store the figure handle */
    mx_in[0] = mxCreateDoubleScalar(0); /* action=0, initialize */
    mx_in[1] = mxCreateDoubleScalar(0); /* ignored */
    mx_in[2] = mxCreateDoubleScalar(0); /* ignored */
    mexCallMATLAB(1,mx_out,3,mx_in,"kim_info");
    kimData->fig_handle = (int)*mxGetPr(mx_out[0]);
  }

  /* Call KINInit */
  status = KINInit(kin_mem, mxW_KINSys, y);
  if (status != KIN_SUCCESS) goto error_return;

  /* Redirect output */
  status = KINSetErrFile(kin_mem, stdout);
  if (status != KIN_SUCCESS) goto error_return;

  /* Optional inputs */
  status = KINSetNumMaxIters(kin_mem, mxiter);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINSetNoInitSetup(kin_mem, noInitSetup);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINSetNoMinEps(kin_mem, noMinEps);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINSetMaxSetupCalls(kin_mem, msbset);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINSetMaxSubSetupCalls(kin_mem, msbsetsub);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINSetMaxBetaFails(kin_mem, mxnbcf);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINSetEtaForm(kin_mem, etachoice);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINSetEtaConstValue(kin_mem, eta);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINSetEtaParams(kin_mem, egamma, ealpha);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINSetMaxNewtonStep(kin_mem, mxnewtstep);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINSetRelErrFunc(kin_mem, relfunc);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINSetFuncNormTol(kin_mem, fnormtol);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINSetScaledStepTol(kin_mem, scsteptol);
  if (status != KIN_SUCCESS) goto error_return;

  /* Constraints */
  if (constraints != NULL) {
    NVconstraints = N_VCloneEmpty(y);
    N_VSetArrayPointer(constraints, NVconstraints);
    status = KINSetConstraints(kin_mem, NVconstraints);
    if (status != KIN_SUCCESS) goto error_return;
    N_VDestroy(NVconstraints);
  }

  /*
   * --------------------
   * Attach linear solver
   * --------------------
   */

  status = get_LinSolvOptions(options,
                              &mupper, &mlower,
                              &mudq, &mldq, &dqrely,
                              &ptype, &maxrs, &maxl);
  if (status != 0) goto error_return;

  switch (ls) {

  case LS_NONE:

    kimErrHandler(-999, "KINSOL", "KINInit",
                  "No linear solver was specified.", NULL);
    goto error_return;

  case LS_DENSE:

    status = KINDense(kin_mem,  N);
    if (status != KIN_SUCCESS) goto error_return;
    if (!mxIsEmpty(mtlb_JACfct)) {
      status = KINDlsSetDenseJacFn(kin_mem, mxW_KINDenseJac);
      if (status != KIN_SUCCESS) goto error_return;
    }

    break;

  case LS_BAND:

    status = KINBand(kin_mem, N, mupper, mlower);
    if (status != KIN_SUCCESS) goto error_return;
    if (!mxIsEmpty(mtlb_JACfct)) {
      status = KINDlsSetBandJacFn(kin_mem, mxW_KINBandJac);
      if (status != KIN_SUCCESS) goto error_return;
    }

    break;

  case LS_SPGMR:

    status = KINSpgmr(kin_mem, maxl);
    if (status != KIN_SUCCESS) goto error_return;
    status = KINSpilsSetMaxRestarts(kin_mem, maxrs);
    if (status != KIN_SUCCESS) goto error_return;

    break;

  case LS_SPBCG:

    status = KINSpbcg(kin_mem, maxl);
    if (status != KIN_SUCCESS) goto error_return;

    break;

  case LS_SPTFQMR:

    status = KINSptfqmr(kin_mem, maxl);
    if (status != KIN_SUCCESS) goto error_return;

    break;

  }

  /* Jacobian * vector and preconditioner for SPILS linear solvers */

  if ( (ls==LS_SPGMR) || (ls==LS_SPBCG) || (ls==LS_SPTFQMR) ) {

    if (!mxIsEmpty(mtlb_JACfct)) {
      status = KINSpilsSetJacTimesVecFn(kin_mem, mxW_KINSpilsJac);
      if (status != KIN_SUCCESS) goto error_return;
    }

    switch (pm) {

    case PM_NONE:

      if (!mxIsEmpty(mtlb_PSOLfct)) {
        
        if (!mxIsEmpty(mtlb_PSETfct)) status = KINSpilsSetPreconditioner(kin_mem, mxW_KINSpilsPset, mxW_KINSpilsPsol);
        else                          status = KINSpilsSetPreconditioner(kin_mem, NULL, mxW_KINSpilsPsol);
        if (status != KIN_SUCCESS) goto error_return;
        
      }
      
      break;

    case PM_BBDPRE:

      if (!mxIsEmpty(mtlb_GCOMfct)) status = KINBBDPrecInit(kin_mem, N, mudq, mldq, mupper, mlower, dqrely, mxW_KINGloc, mxW_KINGcom);
      else                          status = KINBBDPrecInit(kin_mem, N, mudq, mldq, mupper, mlower, dqrely, mxW_KINGloc, NULL);
      if (status != KIN_SUCCESS) goto error_return;
      
      break;

    }


  }

  /* Set errMsg field in global data 
   * (all error messages from here on will respect this) */
  
  kimData->errMsg = errmsg;

  /* Successfull return */

  status = 0;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(0);
  
  /* Error return */

 error_return:
  status = -1;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(-1);

}

static int KIM_Solve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *y0, *ys, *fs;
  N_Vector yscale = NULL, fscale = NULL;
  int buflen, status, strategy;
  char *bufval;

  /*
   * ----------------------------------------------------------------
   * Extract input arguments
   * ----------------------------------------------------------------
   */

  /* Exract y0 and load initial guess in y */

  y0 = mxGetPr(prhs[0]);
  PutData(y, y0, N);

  /* Extract strategy */

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"None")) {
    strategy = KIN_NONE;
  } else if(!strcmp(bufval,"LineSearch")) {
    strategy = KIN_LINESEARCH;
  } else {
    kimErrHandler(-999, "KINSOL", "KINSol",
                  "Illegal value for strategy.", NULL); 
    goto error_return;
  }

  /* Extract yscale */

  ys = mxGetPr(prhs[2]);
  yscale = N_VCloneEmpty(y);
  N_VSetArrayPointer(ys, yscale);

  /* Extract fscale */

  fs = mxGetPr(prhs[3]);
  fscale = N_VCloneEmpty(y);
  N_VSetArrayPointer(fs, fscale);

  /*
   * -------------------------
   * Call main solver function
   * -------------------------
   */

  status = KINSol(kin_mem, y, strategy, yscale, fscale);
  if (status < 0) goto error_return;

  /* Extract solution vector */
  plhs[1] = mxCreateDoubleMatrix(N,1,mxREAL);
  GetData(y, mxGetPr(plhs[1]), N);

  /* Free temporary vectors */
  N_VDestroy(yscale);
  N_VDestroy(fscale);

  /* KINSOL return flag (only non-negative values make it here) */
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(0);

 error_return:
  status = -1;
  plhs[0] = mxCreateDoubleScalar((double)status);
  plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL);
  if (yscale != NULL) N_VDestroy(yscale);
  if (fscale != NULL) N_VDestroy(fscale);
  return(-1);

}

static int KIM_Stats(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
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
  const char *fnames_band[]={
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

  sunindextype nfe, nni, nbcf, nbops;
  double fnorm, step;

  sunindextype njeD, nfeD;
  sunindextype nli, npe, nps, ncfl;

  mxArray *mx_ls;
  int nfields;

  int status;


  if (kimData == NULL) return;

  status = KINGetNumNonlinSolvIters(kin_mem, &nni);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINGetNumFuncEvals(kin_mem, &nfe);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINGetNumBetaCondFails(kin_mem, &nbcf);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINGetNumBacktrackOps(kin_mem, &nbops);
  if (status != KIN_SUCCESS) goto error_return;
  
  status = KINGetFuncNorm(kin_mem, &fnorm);
  if (status != KIN_SUCCESS) goto error_return;
  status = KINGetStepLength(kin_mem, &step);
  if (status != KIN_SUCCESS) goto error_return;

  nfields = sizeof(fnames_solver)/sizeof(*fnames_solver);
  plhs[0] = mxCreateStructMatrix(1, 1, nfields, fnames_solver);
 
  mxSetField(plhs[0], 0, "nfe",   mxCreateDoubleScalar((double)nfe));
  mxSetField(plhs[0], 0, "nni",   mxCreateDoubleScalar((double)nni));
  mxSetField(plhs[0], 0, "nbcf",  mxCreateDoubleScalar((double)nbcf));
  mxSetField(plhs[0], 0, "nbops", mxCreateDoubleScalar((double)nbops));
  mxSetField(plhs[0], 0, "fnorm", mxCreateDoubleScalar(fnorm));
  mxSetField(plhs[0], 0, "step",  mxCreateDoubleScalar(step));

  /* Linear Solver Statistics */

  switch(ls){

  case LS_DENSE:
    
    status = KINDlsGetNumJacEvals(kin_mem, &njeD);
    if (status != KIN_SUCCESS) goto error_return;
    status = KINDlsGetNumFuncEvals(kin_mem, &nfeD);
    if (status != KIN_SUCCESS) goto error_return;

    nfields = sizeof(fnames_dense)/sizeof(*fnames_dense);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_dense);

    mxSetField(mx_ls, 0, "name", mxCreateString("Dense"));
    mxSetField(mx_ls, 0, "njeD", mxCreateDoubleScalar((double)njeD));
    mxSetField(mx_ls, 0, "nfeD", mxCreateDoubleScalar((double)nfeD));
    
    break;

  case LS_BAND:
    
    status = KINDlsGetNumJacEvals(kin_mem, &njeD);
    if (status != KIN_SUCCESS) goto error_return;
    status = KINDlsGetNumFuncEvals(kin_mem, &nfeD);
    if (status != KIN_SUCCESS) goto error_return;
    
    nfields = sizeof(fnames_band)/sizeof(*fnames_band);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_band);

    mxSetField(mx_ls, 0, "name", mxCreateString("Band"));
    mxSetField(mx_ls, 0, "njeD", mxCreateDoubleScalar((double)njeD));
    mxSetField(mx_ls, 0, "nfeD", mxCreateDoubleScalar((double)nfeD));
    
    break;

  case LS_SPGMR:
  case LS_SPBCG:
  case LS_SPTFQMR:

    status = KINSpilsGetNumLinIters(kin_mem, &nli);
    if (status != KIN_SUCCESS) goto error_return;
    status = KINSpilsGetNumPrecEvals(kin_mem, &npe);
    if (status != KIN_SUCCESS) goto error_return;
    status = KINSpilsGetNumPrecSolves(kin_mem, &nps);
    if (status != KIN_SUCCESS) goto error_return;
    status = KINSpilsGetNumConvFails(kin_mem, &ncfl);
    if (status != KIN_SUCCESS) goto error_return;
    
    nfields = sizeof(fnames_spils)/sizeof(*fnames_spils);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_spils);

    if (ls == LS_SPGMR)
      mxSetField(mx_ls, 0, "name",  mxCreateString("GMRES"));
    else if (ls == LS_SPBCG)
      mxSetField(mx_ls, 0, "name",  mxCreateString("BiCGStab"));
    else
      mxSetField(mx_ls, 0, "name",  mxCreateString("TFQMR"));

    mxSetField(mx_ls, 0, "nli",   mxCreateDoubleScalar((double)nli));
    mxSetField(mx_ls, 0, "npe",   mxCreateDoubleScalar((double)npe));
    mxSetField(mx_ls, 0, "nps",   mxCreateDoubleScalar((double)nps));
    mxSetField(mx_ls, 0, "ncfl",  mxCreateDoubleScalar((double)ncfl));
    
    break;

  }

  mxSetField(plhs[0], 0, "LSInfo", mx_ls);

  /* Successfull return */

  status = 0;
  plhs[1] = mxCreateDoubleScalar((double)status);
  return(0);

  /* Error return */

 error_return:
  status = -1;
  plhs[1] = mxCreateDoubleScalar((double)status);
  return(-1);

}


static int KIM_Set(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}

static int KIM_Get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}

static int KIM_Free(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (kimData == NULL) return(0);

  KINFree(&kin_mem);

  return;
}
