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
 * MEX implementation for CVODES Matlab interface.
 * -----------------------------------------------------------------
 */

#include <string.h>
#include <stdlib.h>
#include "cvm.h"
#include "nvm.h"

/*
 * ---------------------------------------------------------------------------------
 * Global interface data variable
 * ---------------------------------------------------------------------------------
 */

cvmInterfaceData cvmData = NULL;

/*
 * ---------------------------------------------------------------------------------
 * Static function prototypes
 * ---------------------------------------------------------------------------------
 */

static void cvmInitCVODESdata();
static void cvmPersistCVODESdata();
static void cvmFinalCVODESdata();

static void cvmInitPbData(cvmPbData pb);
static void cvmPersistPbData(cvmPbData pb);
static void cvmFinalPbData(cvmPbData pb);


static int CVM_Initialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_QuadInitialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_SensInitialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_AdjInitialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

static int CVM_InitializationB(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_QuadInitializationB(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

static int CVM_SensToggleOff(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

static int CVM_Solve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_SolveB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

static int cvmSolveB_one(mxArray *plhs[], int NtoutB, double *toutB, int itaskB);
static int cvmSolveB_more(mxArray *plhs[], int NtoutB, double *toutB, int itaskB,
                          booleantype any_quadrB, booleantype any_monB);


static int CVM_Stats(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_StatsB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);


static int CVM_Set(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_SetB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

static int CVM_Get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);


static int CVM_Free(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

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
     
     1 - initialize CVODES solver
     2 - initialize quadratures
     3 - initialize forward sensitivity calculations
     4 - initialize adjoint sensitivity calculations
     5 - initialize backward solver
     6 - initialize backward quadratures

    11 - reinitialize CVODES solver
    12 - reinitialize quadratures
    13 - reinitialize forward sensitivity calculations
    14 - reinitialize adjoint sensitivity calculations
    15 - reinitialize backward solver
    16 - reinitialize backward quadratures

    18 - toggle FSA off

    20 - solve problem
    21 - solve backward problem

    30 - get integrator stats
    31 - get backward integrator stats
    32 - extract data from cvode_mem

    33 - set one optional input at a time
    34 - set one optional input at a time for backward problems

    40 - finalize
  */

  mode = (int)mxGetScalar(prhs[0]);

  mexUnlock();


  if ( (mode != 1) && (cvmData == NULL) ) {
    mexErrMsgTxt("CVODES - Illegal attempt to call before CVodeInit.");
  }


  switch(mode) {

    /* Initialization functions */

  case 1:
    if (cvmData != NULL) {
      CVM_Free(nlhs, plhs, nrhs-1, &prhs[1]);
      cvmFinalCVODESdata();
    }
    cvmInitCVODESdata();
    CVM_Initialization(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 2:
    CVM_QuadInitialization(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 3:
    CVM_SensInitialization(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 4:
    CVM_AdjInitialization(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 5:
    CVM_InitializationB(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 6:
    CVM_QuadInitializationB(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

    /* Re-initialization functions */

  case 11:
    CVM_Initialization(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 12:
    CVM_QuadInitialization(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 13:
    CVM_SensInitialization(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 14:
    CVM_AdjInitialization(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 15:
    CVM_InitializationB(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 16:
    CVM_QuadInitializationB(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

    /* Sensitivity toggle function */

  case 18:
    CVM_SensToggleOff(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
    
    /* Solve functions */

  case 20:
    CVM_Solve(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 21:
    CVM_SolveB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

    /* Optional output extraction functions */

  case 30:
    CVM_Stats(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 31:
    CVM_StatsB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 32:
    CVM_Get(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 33:
    CVM_Set(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 34:
    CVM_SetB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

    /* Memory deallocation function */

  case 40:
    CVM_Free(nlhs, plhs, nrhs-1, &prhs[1]);
    cvmFinalCVODESdata();
    return;

  }

  /* Unless this was the CVodeFree call,
   * make data persistent and lock the MEX file */
  if (mode != 40) {
    cvmPersistCVODESdata();
    mexLock();
  }

  return;
}

/*
 * ---------------------------------------------------------------------------------
 * Private functions
 * ---------------------------------------------------------------------------------
 */

static void cvmInitCVODESdata()
{
  /* Allocate space for global CVODES data structure */

  cvmData = (cvmInterfaceData) mxMalloc(sizeof(struct cvmInterfaceData_));

  /* Initialize global CVODES data */

  cvmData->cvode_mem = NULL;

  cvmData->fwdPb     = NULL;
  cvmData->bckPb     = NULL;

  cvmData->NbckPb    = 0;

  cvmData->Nd        = 0;
  cvmData->Nc        = 0;
  cvmData->asa       = SUNFALSE;

  cvmData->errMsg    = SUNTRUE;

  return;
}


static void cvmInitPbData(cvmPbData pb)
{
  mxArray *empty;

  pb->n  = 0;
  pb->nq = 0;
  pb->ng = 0;
  pb->ns = 0;

  pb->Y  = NULL;
  pb->YQ = NULL;
  pb->YS = NULL;

  pb->Quadr  = SUNFALSE;
  pb->Fsa    = SUNFALSE;
  pb->Mon    = SUNFALSE;

  pb->LS = LS_DENSE;
  pb->PM = PM_NONE;

  empty = mxCreateDoubleMatrix(0,0,mxREAL);

  pb->RHSfct   = mxDuplicateArray(empty);
  pb->Gfct     = mxDuplicateArray(empty);
  pb->QUADfct  = mxDuplicateArray(empty);
  pb->SRHSfct  = mxDuplicateArray(empty);
  pb->JACfct   = mxDuplicateArray(empty);
  pb->PSETfct  = mxDuplicateArray(empty);
  pb->PSOLfct  = mxDuplicateArray(empty);
  pb->GLOCfct  = mxDuplicateArray(empty);
  pb->GCOMfct  = mxDuplicateArray(empty);

  pb->MONfct   = mxDuplicateArray(empty);
  pb->MONdata  = mxDuplicateArray(empty);

  pb->mtlb_data = mxDuplicateArray(empty);

  pb->fwd = cvmData->fwdPb;

  pb->index = 0;
  pb->next  = NULL;

  mxDestroyArray(empty);
}


static void cvmPersistCVODESdata()
{
  cvmPbData tmpPb;

  /* Make global memory persistent */

  if (cvmData->fwdPb != NULL) {
    cvmPersistPbData(cvmData->fwdPb);
    mexMakeMemoryPersistent(cvmData->fwdPb);
  }

  tmpPb = cvmData->bckPb;
  while(tmpPb != NULL) {
    cvmPersistPbData(tmpPb);
    mexMakeMemoryPersistent(tmpPb);
    tmpPb = tmpPb->next;
  }
  
  mexMakeMemoryPersistent(cvmData);

  return;
}


static void cvmPersistPbData(cvmPbData pb)
{
  mexMakeArrayPersistent(pb->mtlb_data);

  mexMakeArrayPersistent(pb->RHSfct);
  mexMakeArrayPersistent(pb->Gfct);
  mexMakeArrayPersistent(pb->QUADfct);
  mexMakeArrayPersistent(pb->SRHSfct);
  mexMakeArrayPersistent(pb->JACfct);
  mexMakeArrayPersistent(pb->PSETfct);
  mexMakeArrayPersistent(pb->PSOLfct);
  mexMakeArrayPersistent(pb->GLOCfct);
  mexMakeArrayPersistent(pb->GCOMfct);

  mexMakeArrayPersistent(pb->MONfct);
  mexMakeArrayPersistent(pb->MONdata);
}

static void cvmFinalCVODESdata()
{  
  cvmPbData tmpPb;

  if (cvmData == NULL) return;

  if (cvmData->fwdPb != NULL) {
    cvmFinalPbData(cvmData->fwdPb);
    mxFree(cvmData->fwdPb);
    cvmData->fwdPb = NULL;
  }

  while(cvmData->bckPb != NULL) {
    tmpPb = cvmData->bckPb->next;
    mxFree(cvmData->bckPb);
    cvmData->bckPb = tmpPb;
  }

  mxFree(cvmData);
  cvmData = NULL;

  return;
}


static void cvmFinalPbData(cvmPbData pb)
{
  if (pb->Y != NULL) N_VDestroy(pb->Y);
  if (pb->YQ != NULL) N_VDestroy(pb->YQ);
  if (pb->YS != NULL) N_VDestroyVectorArray(pb->YS, pb->ns);

  mxDestroyArray(pb->mtlb_data);

  mxDestroyArray(pb->RHSfct);
  mxDestroyArray(pb->Gfct);
  mxDestroyArray(pb->QUADfct);
  mxDestroyArray(pb->SRHSfct);
  mxDestroyArray(pb->JACfct);
  mxDestroyArray(pb->PSETfct);
  mxDestroyArray(pb->PSOLfct);
  mxDestroyArray(pb->GLOCfct);
  mxDestroyArray(pb->GCOMfct);

  mxDestroyArray(pb->MONfct);
  mxDestroyArray(pb->MONdata);
}

/*
 * ---------------------------------------------------------------------------------
 * Error handler function.
 *
 * This function is both passed as the CVODES error handler and used throughout
 * the Matlab interface.
 *
 * If called directly by one of the interface functions, error_code = -999 to
 * indicate an error and err_code = +999 to indicate a warning. Otherwise,
 * err_code is set by the calling CVODES function.
 *
 * NOTE: mexErrMsgTxt will end the execution of the MEX file. Therefore we do
 *       not have to intercept any of the CVODES error return flags.
 *       The only return flags we intercept are those from CVode() and CVodeB()
 *       which are passed back to the user (only positive values will make it). 
 * ---------------------------------------------------------------------------------
 */

void cvmErrHandler(int error_code, 
                   const char *module, const char *function, 
                   char *msg, void *eh_data)
{
  char err_msg[256];

  if (!(cvmData->errMsg)) return;

  if (error_code > 0) {
    sprintf(err_msg,"Warning in ==> %s\n%s",function,msg);
    mexWarnMsgTxt(err_msg);    
  } else if (error_code < 0) {
    /*mexUnlock();
      cvmFinalCVODESdata();*/
    sprintf(err_msg,"Error using ==> %s\n%s",function,msg);
    mexErrMsgTxt(err_msg);
  }

  return;
}

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */

#define cvode_mem   (cvmData->cvode_mem)

#define asa         (cvmData->asa)
#define Nd          (cvmData->Nd) 
#define Nc          (cvmData->Nc) 
#define NbckPb      (cvmData->NbckPb)

#define fsa         (fwdPb->Fsa)
#define quadr       (fwdPb->Quadr)
#define mon         (fwdPb->Mon)
#define rootSet     (fwdPb->RootSet)
#define tstopSet    (fwdPb->TstopSet)

#define y           (fwdPb->Y) 
#define yQ          (fwdPb->YQ) 
#define yS          (fwdPb->YS) 
#define N           (fwdPb->n) 
#define Nq          (fwdPb->nq) 
#define Ng          (fwdPb->ng) 
#define Ns          (fwdPb->ns) 
#define ls          (fwdPb->LS) 
#define pm          (fwdPb->PM)

#define mtlb_data     (fwdPb->mtlb_data)

#define mtlb_RHSfct   (fwdPb->RHSfct)
#define mtlb_QUADfct  (fwdPb->QUADfct)
#define mtlb_JACfct   (fwdPb->JACfct)
#define mtlb_PSETfct  (fwdPb->PSETfct)
#define mtlb_PSOLfct  (fwdPb->PSOLfct)
#define mtlb_GLOCfct  (fwdPb->GLOCfct)
#define mtlb_GCOMfct  (fwdPb->GCOMfct)
#define mtlb_Gfct     (fwdPb->Gfct)
#define mtlb_SRHSfct  (fwdPb->SRHSfct)

#define mtlb_MONfct   (fwdPb->MONfct)
#define mtlb_MONdata  (fwdPb->MONdata)



#define indexB      (bckPb->index)

#define quadrB      (bckPb->Quadr)
#define monB        (bckPb->Mon)

#define yB          (bckPb->Y) 
#define yQB         (bckPb->YQ) 
#define NB          (bckPb->n) 
#define NqB         (bckPb->nq) 
#define lsB         (bckPb->LS) 
#define pmB         (bckPb->PM) 

#define mtlb_dataB    (bckPb->mtlb_data)

#define mtlb_RHSfctB  (bckPb->RHSfct)
#define mtlb_QUADfctB (bckPb->QUADfct)
#define mtlb_JACfctB  (bckPb->JACfct)
#define mtlb_PSETfctB (bckPb->PSETfct)
#define mtlb_PSOLfctB (bckPb->PSOLfct)
#define mtlb_GLOCfctB (bckPb->GLOCfct)
#define mtlb_GCOMfctB (bckPb->GCOMfct)

#define mtlb_MONfctB  (bckPb->MONfct)
#define mtlb_MONdataB (bckPb->MONdata)



/*
 * ---------------------------------------------------------------------------------
 * Exported procedures
 * ---------------------------------------------------------------------------------
 */

/* CVM_Initialization
 *
 * action = 0   -> CVodeCreate + CVodeInit
 * prhs contains:
 *   fct 
 *   lmm
 *   iter
 *   t0
 *   y0
 *   options
 *   data
 *
 * action = 1   -> CVodeReInit
 * prhs contains:
 *   t0
 *   y0
 *   options
 *
 */

static int CVM_Initialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData fwdPb;

  const mxArray *options;

  double t0, *y0;

  int lmm, iter, maxord;

  sunindextype mxsteps;

  int itol;
  realtype reltol, Sabstol, *Vabstol;
  N_Vector NV_abstol;

  double hin, hmax, hmin;

  double tstop;

  booleantype sld;

  booleantype errmsg;

  booleantype rhs_s; /* ignored */

  sunindextype mupper, mlower;
  int ptype, gstype, maxl;
  sunindextype mudq, mldq;
  double dqrely;

  char *bufval;
  int buflen, status;

  /* 
   * ------------------------------------
   * Process inputs based on action
   * ------------------------------------
   */

  switch (action) {

  case 0:     /* SOLVER INITIALIZATION */

    /* Create and initialize a new problem */

    fwdPb = (cvmPbData) mxMalloc(sizeof(struct cvmPbData_));
    cvmInitPbData(fwdPb);

    cvmData->fwdPb = fwdPb;

    /* Initialize appropriate vector module */

    InitVectors();

    /* Extract user-provided RHS function */

    mxDestroyArray(mtlb_RHSfct);
    mtlb_RHSfct = mxDuplicateArray(prhs[0]);

    /* Extract lmm */

    buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(prhs[1], bufval, buflen);
    if(status != 0) {
      cvmErrHandler(-999, "CVODES", "CVodeInit",
                    "Cannot parse LMM input argument.", NULL);
      goto error_return;
    }
    if(!strcmp(bufval,"Adams"))    {
      lmm = CV_ADAMS;
    } else if(!strcmp(bufval,"BDF")) {
      lmm = CV_BDF;
    } else {
      cvmErrHandler(-999, "CVODES", "CVodeInit",
                    "LMM has an illegal value.", NULL);
      goto error_return;
    }
    mxFree(bufval);

    /* Extract iter */

    buflen = mxGetM(prhs[2]) * mxGetN(prhs[2]) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(prhs[2], bufval, buflen);
    if(status != 0) {
      cvmErrHandler(-999, "CVODES", "CVodeInit",
                    "Cannot parse NLS input argument.", NULL);
      goto error_return;
    }
    if(!strcmp(bufval,"Functional"))  {
      iter = CV_FUNCTIONAL;
    } else if(!strcmp(bufval,"Newton")) {
      iter = CV_NEWTON;
    } else {
      cvmErrHandler(-999, "CVODES", "CVodeInit",
                    "NLS has an illegal value.", NULL);
      goto error_return;
    }
    mxFree(bufval);

    /* Extract initial time */

    t0 = (double)mxGetScalar(prhs[3]);

    /* Extract initial conditions */

    y0 = mxGetPr(prhs[4]);
    N = mxGetM(prhs[4]);

    /* Create the solution N_Vector */

    y = NewVector(N);

    /* Load initial conditions */

    PutData(y, y0, N);

    /* Extract options structure */
    
    options = prhs[5];

    break;

  case 1:    /* SOLVER RE-INITIALIZATION */

    fwdPb = cvmData->fwdPb;

    /* If monitoring was enabled, finalize it now. */

    if (mon) mxW_CVodeMonitor(2, 0.0, NULL, NULL, NULL, fwdPb);

    /* Extract initial time */

    t0 = (double)mxGetScalar(prhs[0]);

    /* Extract initial conditions */

    y0 = mxGetPr(prhs[1]);

    if (mxGetM(prhs[1]) != N) {
      cvmErrHandler(-999, "CVODES", "CVodeReInit",
                    "Size of y0 changed from CVodeInit call.", NULL);
      goto error_return;
    }

    /* Load initial conditions */

    PutData(y, y0, N);

    /* Extract options structure */
    
    options = prhs[2];

    break;

  }

  /* Process the options structure */

  status = get_IntgrOptions(options, fwdPb, SUNTRUE, lmm,
                            &maxord, &sld, &errmsg, &mxsteps,
                            &itol, &reltol, &Sabstol, &Vabstol,
                            &hin, &hmax, &hmin, &tstop, &rhs_s);
  if (status != 0) goto error_return;

  /* 
   * ----------------------------------------
   * Call appropriate CVODES functions
   *
   * If action = 0
   *    Create CVODES object and allocate memory
   *    Attach error handler function
   *    Redirect output
   * If action = 1
   *    Reinitialize solver
   * ----------------------------------------
   */

  switch (action) {

  case 0:

    /* Create CVODES object */
    cvode_mem = CVodeCreate(lmm, iter);
    if (cvode_mem == NULL) goto error_return;

    /* Attach the global CVODES data as 'user-data' */
    status = CVodeSetUserData(cvode_mem, fwdPb);
    if (status != CV_SUCCESS) goto error_return;

    /* Attach error handler function */
    status = CVodeSetErrHandlerFn(cvode_mem, cvmErrHandler, NULL);
    if (status != CV_SUCCESS) goto error_return;

    /* Call CVodeInit */
    status = CVodeInit(cvode_mem, mxW_CVodeRhs, t0, y);
    if (status != CV_SUCCESS) goto error_return;

    /* Redirect output */
    status = CVodeSetErrFile(cvode_mem, stdout);
    if (status != CV_SUCCESS) goto error_return;

    break;

  case 1:

    /* Reinitialize solver */
    status = CVodeReInit(cvode_mem, t0, y);
    if (status != CV_SUCCESS) goto error_return;

    break;

  }

  /*
   * ----------------------------------------
   * Set tolerances
   * ----------------------------------------
   */

  switch (itol) {
    case CV_SS:
      status = CVodeSStolerances(cvode_mem, reltol, Sabstol);
      if (status != CV_SUCCESS) goto error_return;
      break;
    case CV_SV:
      NV_abstol = N_VClone(y);
      PutData(NV_abstol, Vabstol, N);
      status = CVodeSVtolerances(cvode_mem, reltol, NV_abstol);
      if (status != CV_SUCCESS) goto error_return;
      N_VDestroy(NV_abstol);
      break;
    }
  
  /*
   * --------------------------------
   * Set various optional inputs
   * --------------------------------
   */

  /* set maxorder (default is consistent with LMM) */
  status = CVodeSetMaxOrd(cvode_mem, maxord);
  if (status != CV_SUCCESS) goto error_return;

  /* set initial step size (the default value of 0.0 is ignored by CVODES) */
  status = CVodeSetInitStep(cvode_mem, hin);
  if (status != CV_SUCCESS) goto error_return;

  /* set max step (default is infinity) */
  status = CVodeSetMaxStep(cvode_mem, hmax);
  if (status != CV_SUCCESS) goto error_return;

  /* set min step (default is 0) */
  status = CVodeSetMinStep(cvode_mem, hmin);
  if (status != CV_SUCCESS) goto error_return;
 
  /* set number of max steps */
  status = CVodeSetMaxNumSteps(cvode_mem, mxsteps);
  if (status != CV_SUCCESS) goto error_return;

  /* set tstop? */
  if (tstopSet) {
    status = CVodeSetStopTime(cvode_mem, tstop);
    if (status != CV_SUCCESS) goto error_return;
  }

  /* set stability limit detection (default is SUNFALSE) */
  status = CVodeSetStabLimDet(cvode_mem, sld);
  if (status != CV_SUCCESS) goto error_return;

  /* Rootfinding? */
  if ( !mxIsEmpty(mtlb_Gfct) && (Ng > 0) ) {
    status = CVodeRootInit(cvode_mem, Ng, mxW_CVodeGfct);
    if (status != CV_SUCCESS) goto error_return;
    rootSet = SUNTRUE;
  } else {
    rootSet = SUNFALSE;
  }

  /*
   * ----------------------------------------
   * Need a linear solver?
   * ----------------------------------------
   */

  if (iter == CV_NEWTON) {

    status = get_LinSolvOptions(options, fwdPb, SUNTRUE,
                                &mupper, &mlower,
                                &mudq, &mldq, &dqrely,
                                &ptype, &gstype, &maxl);
    if (status != 0) goto error_return;

    switch (ls) {

    case LS_DENSE:

      status = CVDense(cvode_mem, N);
      if (status != CV_SUCCESS) goto error_return;
      if (!mxIsEmpty(mtlb_JACfct)) {
        status = CVDlsSetDenseJacFn(cvode_mem, mxW_CVodeDenseJac);
        if (status != CV_SUCCESS) goto error_return;
      }

      break;

    case LS_DIAG:

      status = CVDiag(cvode_mem);
      if (status != CV_SUCCESS) goto error_return;

      break;

    case LS_BAND:

      status = CVBand(cvode_mem, N, mupper, mlower);
      if (status != CV_SUCCESS) goto error_return;
      if (!mxIsEmpty(mtlb_JACfct)) {
        status = CVDlsSetBandJacFn(cvode_mem, mxW_CVodeBandJac);
        if (status != CV_SUCCESS) goto error_return;
      }

      break;

    case LS_SPGMR:

      status = CVSpgmr(cvode_mem, ptype, maxl);
      if (status != CV_SUCCESS) goto error_return;
      status = CVSpilsSetGSType(cvode_mem, gstype);
      if (status != CV_SUCCESS) goto error_return;

      break;

    case LS_SPBCG:

      status = CVSpbcg(cvode_mem, ptype, maxl);
      if (status != CV_SUCCESS) goto error_return;

      break;

    case LS_SPTFQMR:

      status = CVSptfqmr(cvode_mem, ptype, maxl);
      if (status != CV_SUCCESS) goto error_return;

      break;

    }


    /* Jacobian * vector and preconditioner for SPILS linear solvers */

    if ( (ls==LS_SPGMR) || (ls==LS_SPBCG) || (ls==LS_SPTFQMR) ) {

      if (!mxIsEmpty(mtlb_JACfct)) {
        status = CVSpilsSetJacTimesVecFn(cvode_mem, mxW_CVodeSpilsJac);
        if (status != CV_SUCCESS) goto error_return;
      }

      switch (pm) {

      case PM_NONE:

        if (!mxIsEmpty(mtlb_PSOLfct)) {

          if (!mxIsEmpty(mtlb_PSETfct)) status = CVSpilsSetPreconditioner(cvode_mem, mxW_CVodeSpilsPset, mxW_CVodeSpilsPsol);
          else                          status = CVSpilsSetPreconditioner(cvode_mem, NULL, mxW_CVodeSpilsPsol);
          if (status != CV_SUCCESS) goto error_return;

        }

        break;

      case PM_BANDPRE:

        status = CVBandPrecInit(cvode_mem, N, mupper, mlower);
        if (status != CV_SUCCESS) goto error_return;

        break;

      case PM_BBDPRE:

        if (!mxIsEmpty(mtlb_GCOMfct)) status = CVBBDPrecInit(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mxW_CVodeBBDgloc, mxW_CVodeBBDgcom);
        else                          status = CVBBDPrecInit(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mxW_CVodeBBDgloc, NULL);
        if (status != CV_SUCCESS) goto error_return;

        break;
      }

    }

  } else {

    ls = LS_NONE;

  }

  /* Do we monitor? */
  
  if (mon) mxW_CVodeMonitor(0, t0, NULL, NULL, NULL, fwdPb);

  /* Set errMsg field in global data 
   * (all error messages from here on will respect this) */

  cvmData->errMsg = errmsg;

  /* Successful return */

  status = 0;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(0);

  /* Error return */

 error_return:
  status = -1;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(-1);

}

/* CVM_QuadInitialization
 *
 * action = 0   -> CVodeQuadInit
 * prhs contains:
 *   fQ
 *   y0
 *   options
 *
 * action = 1   -> CVodeQuadReInit
 * prhs contains:
 *   y0
 *   options
 *
 */

static int CVM_QuadInitialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData fwdPb;

  const mxArray *options;

  double *yQ0;

  booleantype rhs_s; /* ignored */

  booleantype errconQ;
  int itolQ;
  realtype reltolQ, SabstolQ, *VabstolQ;
  N_Vector NV_abstolQ;

  int status;

  fwdPb = cvmData->fwdPb;

  /* 
   * ------------------------------------
   * Process inputs based on action
   * ------------------------------------
   */

  switch (action) {

  case 0:     /* QUADRATURE INITIALIZATION */

    /* Extract user-provided quadrature RHS function */

    mxDestroyArray(mtlb_QUADfct);
    mtlb_QUADfct = mxDuplicateArray(prhs[0]);
  
    /* Extract quadrature initial conditions */

    yQ0 = mxGetPr(prhs[1]);
    Nq = mxGetM(prhs[1]);

    /* Create the quadrature N_Vector */

    yQ = NewVector(Nq);

    /* Load quadrature initial conditions */
    
    PutData(yQ, yQ0, Nq);

    /* Extract quadrature options structure */

    options = prhs[2];

    break;

  case 1:     /* QUADRATURE RE-INITIALIZATION */

    /* Extract quadrature initial conditions */

    yQ0 = mxGetPr(prhs[0]);

    if (mxGetM(prhs[0]) != Nq) {
      cvmErrHandler(-999, "CVODES", "CVodeQuadReInit",
                    "Size of yQ0 changed from CVodeQuadInit call.", NULL);
      goto error_return;
    }

    /* Load quadrature initial conditions */
    
    PutData(yQ, yQ0, Nq);

    /* Extract quadrature options structure */

    options = prhs[1];

    break;

  }

  /* Process the options structure */

  status = get_QuadOptions(options, fwdPb, SUNTRUE,
                           Nq, &rhs_s,
                           &errconQ, 
                           &itolQ, &reltolQ, &SabstolQ, &VabstolQ);
  if (status != 0) goto error_return;

  /* 
   * ----------------------------------------
   * Call appropriate CVODES functions
   *
   * If action = 0
   *    Initialize quadratures
   * If action = 1
   *    Reinitialize quadratures
   * ----------------------------------------
   */

  switch (action) {
  case 0:
    status = CVodeQuadInit(cvode_mem, mxW_CVodeQUADfct, yQ);
    if (status != CV_SUCCESS) goto error_return;
    break;
  case 1:
    status = CVodeQuadReInit(cvode_mem, yQ);
    if (status != CV_SUCCESS) goto error_return;
    break;
  }

  /*
   * ----------------------------------------
   * Set tolerances for quadrature variables
   * ----------------------------------------
   */

  status = CVodeSetQuadErrCon(cvode_mem, errconQ);
  if (status != CV_SUCCESS) goto error_return;

  if (errconQ) {
    
    switch (itolQ) {
    case CV_SS:
      status = CVodeQuadSStolerances(cvode_mem, reltolQ, SabstolQ);
      if (status != CV_SUCCESS) goto error_return;
      break;
    case CV_SV:
      NV_abstolQ = N_VClone(yQ);
      PutData(NV_abstolQ, VabstolQ, Nq);
      status = CVodeQuadSVtolerances(cvode_mem, reltolQ, NV_abstolQ);
      if (status != CV_SUCCESS) goto error_return;
      N_VDestroy(NV_abstolQ);
      break;
    }
    
  }

  /* Quadratures will be integrated */

  quadr = SUNTRUE;

  /* Successful return */

  status = 0;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(0);

  /* Error return */

 error_return:
  status = -1;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(-1);

}

/* CVM_SensInitialization
 *
 * action = 0 -> CVodeSensInit
 * prhs contains:
 *   Ns
 *   fS
 *   yS0
 *   options
 * action = 1 -> CVodeSensReInit
 *   yS0
 *   options
 *
 */

static int CVM_SensInitialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData fwdPb;

  const mxArray *options;

  booleantype fS_DQ;
  CVSensRhsFn rhsS;

  double *yS0;

  int ism;

  mxArray *pfield;
  char *pfield_name;

  booleantype errconS;
  int itolS;
  realtype reltolS;
  realtype *SabstolS, *VabstolS;
  N_Vector *NV_abstolS;

  int *plist, dqtype;
  double *p, *pbar, rho;

  int is, status;

  p = NULL;
  plist = NULL;
  pbar = NULL;

  fwdPb = cvmData->fwdPb;

  /* 
   * ------------------------------------
   * Process inputs based on action
   * ------------------------------------
   */

  switch (action) {

  case 0:     /* FSA INITIALIZATION */

    /* Extract number of sensitivities */

    Ns = (int)mxGetScalar(prhs[0]);

    /* Extract user-provided sensitivity RHS function */

    if ( mxIsEmpty(prhs[1]) ) {
      rhsS = NULL;
      fS_DQ = SUNTRUE;
    } else {
      mxDestroyArray(mtlb_SRHSfct);
      mtlb_SRHSfct = mxDuplicateArray(prhs[1]);
      rhsS = mxW_CVodeSensRhs;
      fS_DQ = SUNFALSE;
    }

    /* Extract sensitivity initial condition */

    yS0 = mxGetPr(prhs[2]);

    /* Create the sensitivity N_Vectors */

    yS = N_VCloneVectorArray(Ns, y);

    /* Load sensitivity initial conditions */

    for (is=0;is<Ns;is++)
      PutData(yS[is], &yS0[is*N], N);

    /* Extract FSA options structure */

    options = prhs[3];

    break;

  case 1:     /* FSA RE-INITIALIZATION */

    /* Extract sensitivity initial condition */

    yS0 = mxGetPr(prhs[0]);

    if ( (mxGetM(prhs[0]) != N) || (mxGetN(prhs[0]) != Ns) )  {
      cvmErrHandler(-999, "CVODES", "CVodeSensReInit",
                    "Size of yS0 changed from CVodeSensInit call.", NULL);
      goto error_return;
    }

    /* Load sensitivity initial conditions */

    for (is=0;is<Ns;is++)
      PutData(yS[is], &yS0[is*N], N);

    /* Extract qFSA options structure */

    options = prhs[1];

    break;

  }

  /* Process the options structure */

  status = get_FSAOptions(options, fwdPb, 
                          &ism,
                          &pfield_name, &plist, &pbar,
                          &dqtype, &rho,
                          &errconS, &itolS, &reltolS, &SabstolS, &VabstolS);
  if (status != 0) goto error_return;

  /* 
   * ----------------------------------------
   * Call appropriate CVODES functions
   *
   * If action = 0
   *    Check if required inputs are available
   *    Initialize FSA
   * If action = 1
   *    Reinitialize FSA
   * ----------------------------------------
   */

  switch (action) {

  case 0:

    /* Test required inputs */
    if ( fS_DQ && (pfield_name == NULL) ) {
      cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                    "pfield required but was not provided.", NULL);
      goto error_return;
    }

    status = CVodeSensInit(cvode_mem, Ns, ism, rhsS, yS);
    if (status != CV_SUCCESS) goto error_return;

    break;

  case 1:

    status = CVodeSensReInit(cvode_mem, ism, yS);
    if (status != CV_SUCCESS) goto error_return;

    break;

  }

  /*
   * ----------------------------------------
   * Set tolerances for sensitivity variables
   * ----------------------------------------
   */

  switch (itolS) {
  case CV_SS:
    status = CVodeSensSStolerances(cvode_mem, reltolS, SabstolS);
    if (status != CV_SUCCESS) goto error_return;
    break;
  case CV_SV:
    NV_abstolS = N_VCloneVectorArray(Ns, y);
    for (is=0;is<Ns;is++)
      PutData(NV_abstolS[is], &VabstolS[is*N], N);
    status = CVodeSensSVtolerances(cvode_mem, reltolS, NV_abstolS);
    if (status != CV_SUCCESS) goto error_return;
    N_VDestroyVectorArray(NV_abstolS, Ns);
    break;
  case CV_EE:
    status = CVodeSensEEtolerances(cvode_mem);
    if (status != CV_SUCCESS) goto error_return;
    break;
  }

  /*
   * --------------------------------
   * Set various optional inputs
   * --------------------------------
   */

  if (pfield_name != NULL) {
    pfield = mxGetField(mtlb_data,0,pfield_name);
    if (pfield == NULL) {
      cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                    "illegal pfield input.", NULL);
      goto error_return;
    }
    p = mxGetPr(pfield);
  }
  
  status = CVodeSetSensParams(cvode_mem, p, pbar, plist);
  if (status != CV_SUCCESS) goto error_return;

  status = CVodeSetSensDQMethod(cvode_mem, dqtype, rho);
  if (status != CV_SUCCESS) goto error_return;

  status = CVodeSetSensErrCon(cvode_mem, errconS);
  if (status != CV_SUCCESS) goto error_return;

  fsa = SUNTRUE;

  /* Successful return */

  status = 0;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(0);

  /* Error return */

 error_return:
  status = -1;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(-1);

}

/*
 * CVM_SensToggleOff
 *
 * deactivates FSA
 */

static int CVM_SensToggleOff(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData fwdPb;
  int status;

  fwdPb = cvmData->fwdPb;

  status = CVodeSensToggleOff(cvode_mem);
  if (status != CV_SUCCESS) {
    status = -1;
    plhs[0] = mxCreateDoubleScalar((double)status);
    return(-1);
  }

  fsa = SUNFALSE;

  status = 0;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(0);
}


/* CVM_AdjInitialization
 *
 * prhs contains:
 *   Nd - number of interpolatin data points (i.e. steps between check points)
 *   interp - type of interpolation
 *
 */

static int CVM_AdjInitialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int interp;

  int buflen, status;
  char *bufval;

  switch (action) {

  case 0:

    /* Number of steps */

    Nd = (int)mxGetScalar(prhs[0]);

    /* Interpolation method */

    buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(prhs[1], bufval, buflen);
    if(status != 0) {
      cvmErrHandler(-999, "CVODES", "CVodeAdjInit", 
                    "Could not parse InterpType.", NULL);
      goto error_return;
    }
    if(!strcmp(bufval,"Hermite")) {
      interp = CV_HERMITE;
    } else if(!strcmp(bufval,"Polynomial")) {
      interp = CV_POLYNOMIAL;
    } else {
      cvmErrHandler(-999, "CVODES", "CVodeAdjInit",
                    "Interp. type has an illegal value.", NULL);
      goto error_return;
    }

    status = CVodeAdjInit(cvode_mem, Nd, interp);
    if (status != CV_SUCCESS) goto error_return;

    break;

  case 1:

    status = CVodeAdjReInit(cvode_mem);
    if (status != CV_SUCCESS) goto error_return;

    break;

  }

  asa = SUNTRUE;
  
  /* Successful return */

  status = 0;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(0);

  /* Error return */

 error_return:
  status = -1;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(-1);

}


/* CVM_InitializationB
 *
 * action = 0   -> CVodeCreateB + CVodeInitB
 * prhs contains:
 *   fctB
 *   lmmB
 *   iterB
 *   tF
 *   yB0
 *   options
 * plhs contains:
 *   indexB
 *
 * action = 1   -> CVodeReInitB
 * prhs contains:
 *   indexB
 *   tF
 *   yB0
 *   options
 *
 */

static int CVM_InitializationB(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData bckPb;

  const mxArray *options;

  int idxB;

  double tB0, *yB0;

  int lmmB, iterB, maxordB;
  sunindextype mxstepsB;

  int itolB;
  realtype reltolB, SabstolB, *VabstolB;
  N_Vector NV_abstolB;

  double hinB, hmaxB, hminB;

  double tstopB;            /* ignored */
  booleantype sldB;         /* ignored */
  booleantype errmsgB;      /* ignored */

  booleantype rhs_s;

  sunindextype mupperB, mlowerB;
  int ptypeB, gstypeB, maxlB;
  sunindextype mudqB, mldqB;
  double dqrelyB;

  booleantype found_bck;

  char *bufval;
  int buflen;

  int status;
  int i_status;

  /* Set output containing status */

  i_status = (action == 0) ? 1 : 0;

  /* 
   * -----------------------------
   * Finalize Forward monitoring
   * -----------------------------
   */

  if (cvmData->fwdPb->Mon) {
    mxW_CVodeMonitor(2, 0.0, NULL, NULL, NULL, cvmData->fwdPb);
    cvmData->fwdPb->Mon = SUNFALSE;
  }

  /* 
   * ------------------------------------
   * Process inputs based on action
   * ------------------------------------
   */

  switch (action) {

  case 0:     /* BACKWARD SOLVER INITIALIZATION */

    /* Create and initialize a new problem */

    bckPb = (cvmPbData) mxMalloc(sizeof(struct cvmPbData_));
    cvmInitPbData(bckPb);

    bckPb->next = cvmData->bckPb;
    cvmData->bckPb = bckPb;

    /* Extract user-provided RHS function */

    mxDestroyArray(mtlb_RHSfctB);
    mtlb_RHSfctB = mxDuplicateArray(prhs[0]);

    /* Extract lmmB */

    buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(prhs[1], bufval, buflen);
    if(status != 0) {
      cvmErrHandler(-999, "CVODES", "CVodeInitB",
                    "Cannot parse LMM input argument.", NULL);
      goto error_return;
    }
    if(!strcmp(bufval,"Adams")) {
      lmmB = CV_ADAMS;
    } else if(!strcmp(bufval,"BDF")) {
      lmmB = CV_BDF;
    } else {
      cvmErrHandler(-999, "CVODES", "CVodeInitB",
                    "LMM has an illegal value.", NULL);
      goto error_return;
    }
    mxFree(bufval);

    /* Extract iterB */

    buflen = mxGetM(prhs[2]) * mxGetN(prhs[2]) + 1;
    bufval = mxCalloc(buflen, sizeof(char));
    status = mxGetString(prhs[2], bufval, buflen);
    if(status != 0) {
      cvmErrHandler(-999, "CVODES", "CVodeInitB",
                    "Cannot parse NLS input argument.", NULL);
      goto error_return;
    }
    if(!strcmp(bufval,"Functional")) {
      iterB = CV_FUNCTIONAL;
    } else if(!strcmp(bufval,"Newton")) {
      iterB = CV_NEWTON;
    } else {
      cvmErrHandler(-999, "CVODES", "CVodeInitB",
                    "NLS has an illegal value.", NULL);
      goto error_return;
    }
    mxFree(bufval);

    /* Extract final time */

    tB0 = (double)mxGetScalar(prhs[3]);

    /* Extract final conditions */

    yB0 = mxGetPr(prhs[4]);
    NB = mxGetM(prhs[4]);

    /* Create the solution N_Vector */

    yB = NewVector(NB);

    /* Load final conditions */

    PutData(yB, yB0, NB);

    /* Extract options structure */
    
    options = prhs[5];

    break;

  case 1:     /* BACKWARD SOLVER RE-INITIALIZATION */

    /* Extract index of current backward problem */
    
    idxB = (int)mxGetScalar(prhs[0]);

    /* Find current backward problem */

    found_bck = SUNFALSE;
    bckPb = cvmData->bckPb;
    while (bckPb != NULL) {
      if (indexB == idxB) {
        found_bck = SUNTRUE;
        break;
      }
      bckPb = bckPb->next;
    }

    if (!found_bck) {
      cvmErrHandler(-999, "CVODES", "CVodeReInitB",
                    "idxB has an illegal value.", NULL);
      goto error_return;
    }

    /* If backward monitoring was enabled, finalize it now. */

    if (monB) mxW_CVodeMonitorB(2, indexB, 0.0, NULL, NULL, bckPb);

    /* Extract final time */

    tB0 = (double)mxGetScalar(prhs[1]);

    /* Extract final conditions */

    yB0 = mxGetPr(prhs[2]);

    if (mxGetM(prhs[2]) != NB) {
      cvmErrHandler(-999, "CVODES", "CVodeReInitB",
                    "Size of yB0 changed from CVodeInitB call.", NULL);
      goto error_return;
    }
    
    /* Load final conditions */

    PutData(yB, yB0, NB);

    /* Extract options structure */
    
    options = prhs[3];

    break;

  }

  /* Process the options structure */

  status = get_IntgrOptions(options, bckPb, SUNFALSE, lmmB,
                            &maxordB, &sldB, &errmsgB, &mxstepsB,
                            &itolB, &reltolB, &SabstolB, &VabstolB,
                            &hinB, &hmaxB, &hminB, &tstopB, &rhs_s);
  if (status != 0) goto error_return;

  /* 
   * ----------------------------------------
   * Call appropriate CVODES functions
   *
   * If action = 0
   *    Create CVODES object and allocate memory
   *    Initialize and allocate memory
   * If action = 1
   *    Reinitialize solver
   * ----------------------------------------
   */

  switch (action) {

  case 0:

    status = CVodeCreateB(cvode_mem, lmmB, iterB, &idxB);
    if (status != CV_SUCCESS) goto error_return;

    status = CVodeSetUserDataB(cvode_mem, idxB, bckPb);
    if (status != CV_SUCCESS) goto error_return;

    if (rhs_s) {
      status = CVodeInitBS(cvode_mem, idxB, mxW_CVodeRhsBS, tB0, yB);
    } else {
      status = CVodeInitB(cvode_mem, idxB, mxW_CVodeRhsB, tB0, yB);
    }
    if (status != CV_SUCCESS) goto error_return;

    /* Load idxB in the 1st output (status is 2nd one for this action) */

    plhs[0] = mxCreateDoubleScalar((double)idxB);

    indexB = idxB;

    NbckPb++;

    break;

  case 1:

    status = CVodeReInitB(cvode_mem, idxB, tB0, yB);
    if (status != CV_SUCCESS) goto error_return;

    break;

  }

  /*
   * ----------------------------------------
   * Set tolerances
   * ----------------------------------------
   */

  switch (itolB) {
  case CV_SS:
    status = CVodeSStolerancesB(cvode_mem, idxB, reltolB, SabstolB);
    if (status != CV_SUCCESS) goto error_return;
    break;
  case CV_SV:
    NV_abstolB = N_VClone(yB);
    PutData(NV_abstolB, VabstolB, NB);
    status = CVodeSVtolerancesB(cvode_mem, idxB, reltolB, NV_abstolB);
    if (status != CV_SUCCESS) goto error_return;
    N_VDestroy(NV_abstolB);
    break;
  }

  /*
   * --------------------------------
   * Set various optional inputs
   * --------------------------------
   */

  /* set maxorder (default is consistent with LMM) */
  status = CVodeSetMaxOrdB(cvode_mem, idxB, maxordB);
  if (status != CV_SUCCESS) goto error_return;

  /* set initial step size (the default value of 0.0 is ignored by CVODES) */
  status = CVodeSetInitStepB(cvode_mem, idxB, hinB);
  if (status != CV_SUCCESS) goto error_return;

  /* set max step (default is infinity) */
  status = CVodeSetMaxStepB(cvode_mem, idxB, hmaxB);
  if (status != CV_SUCCESS) goto error_return;

  /* set min step (default is 0) */
  status = CVodeSetMinStepB(cvode_mem, idxB, hminB);
  if (status != CV_SUCCESS) goto error_return;
 
  /* set number of max steps */
  status = CVodeSetMaxNumStepsB(cvode_mem, idxB, mxstepsB);
  if (status != CV_SUCCESS) goto error_return;

  /*
   * ----------------------------------------
   * Need a linear solver?
   * ----------------------------------------
   */

  if (iterB == CV_NEWTON) {

    status = get_LinSolvOptions(options, bckPb, SUNFALSE,
                                &mupperB, &mlowerB,
                                &mudqB, &mldqB, &dqrelyB,
                                &ptypeB, &gstypeB, &maxlB);
    if (status != 0) goto error_return;

    switch(lsB) {

    case LS_DENSE:

      status = CVDenseB(cvode_mem, idxB, NB);
      if (status != CV_SUCCESS) goto error_return;
      if (!mxIsEmpty(mtlb_JACfctB)) {
        status = CVDlsSetDenseJacFnB(cvode_mem, idxB, mxW_CVodeDenseJacB);
        if (status != CV_SUCCESS) goto error_return;
      }

      break;

    case LS_DIAG:

      status = CVDiagB(cvode_mem, idxB);
      if (status != CV_SUCCESS) goto error_return;

      break;

    case LS_BAND:

      status = CVBandB(cvode_mem, idxB, NB, mupperB, mlowerB);
      if (status != CV_SUCCESS) goto error_return;
      if (!mxIsEmpty(mtlb_JACfctB)) {
        status = CVDlsSetBandJacFnB(cvode_mem, idxB, mxW_CVodeBandJacB);
        if (status != CV_SUCCESS) goto error_return;
      }

      break;

    case LS_SPGMR:
      
      status = CVSpgmrB(cvode_mem, idxB, ptypeB, maxlB);
      if (status != CV_SUCCESS) goto error_return;
      status = CVSpilsSetGSTypeB(cvode_mem, idxB, gstypeB);
      if (status != CV_SUCCESS) goto error_return;

      break;

    case LS_SPBCG:

      status = CVSpbcgB(cvode_mem, idxB, ptypeB, maxlB);
      if (status != CV_SUCCESS) goto error_return;

      break;

    case LS_SPTFQMR:

      status = CVSptfqmrB(cvode_mem, idxB, ptypeB, maxlB);
      if (status != CV_SUCCESS) goto error_return;

      break;

    }

    /* Jacobian * vector and preconditioner for SPILS linear solvers */

    if ( (lsB==LS_SPGMR) || (lsB==LS_SPBCG) || (lsB==LS_SPTFQMR) ) {

      if (!mxIsEmpty(mtlb_JACfctB)) {
        status = CVSpilsSetJacTimesVecFnB(cvode_mem, idxB, mxW_CVodeSpilsJacB);
        if (status != CV_SUCCESS) goto error_return;
      }

      switch (pmB) {

      case PM_NONE:

        if (!mxIsEmpty(mtlb_PSOLfctB)) {
          if (!mxIsEmpty(mtlb_PSETfctB)) status = CVSpilsSetPreconditionerB(cvode_mem, idxB, mxW_CVodeSpilsPsetB, mxW_CVodeSpilsPsolB);
          else                           status = CVSpilsSetPreconditionerB(cvode_mem, idxB, NULL, mxW_CVodeSpilsPsolB);
        }
        if (status != CV_SUCCESS) goto error_return;

        break;

      case PM_BANDPRE:

        status = CVBandPrecInitB(cvode_mem, idxB, NB, mupperB, mlowerB);
        if (status != CV_SUCCESS) goto error_return;

        break;

      case PM_BBDPRE:

        if (!mxIsEmpty(mtlb_GCOMfctB)) status = CVBBDPrecInitB(cvode_mem, idxB, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mxW_CVodeBBDglocB, mxW_CVodeBBDgcomB);
        else                           status = CVBBDPrecInitB(cvode_mem, idxB, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mxW_CVodeBBDglocB, NULL);
        if (status != CV_SUCCESS) goto error_return;

        break;

      }

    }

  } else {

    lsB = LS_NONE;

  }

  /* Do we monitor? */

  if (monB) mxW_CVodeMonitorB(0, idxB, tB0, NULL, NULL, bckPb);

  /* Successful return */

  status = 0;
  plhs[i_status] = mxCreateDoubleScalar((double)status);
  return(0);

  /* Error return */

 error_return:
  status = -1;
  plhs[i_status] = mxCreateDoubleScalar((double)status);
  return(-1);

}


/* CVM_QuadInitializationB
 *
 * action = 0   -> CVodeQuadInitB
 * prhs contains:
 *   idxB
 *   fQB
 *   yQB0
 *   options
 *
 * action = 1   -> CVodeQuadReInitB
 *   idxB
 *   yQB0
 *   options
 *
 */

static int CVM_QuadInitializationB(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData bckPb;

  const mxArray *options;

  int idxB;

  double *yQB0;

  booleantype rhs_s;
  
  booleantype errconQB;
  int itolQB;
  realtype reltolQB, SabstolQB, *VabstolQB;
  N_Vector NV_abstolQB;

  booleantype found_bck;

  int status;

  /* Extract index of current backward problem */
    
  idxB = (int)mxGetScalar(prhs[0]);

  /* Find current backward problem */

  found_bck = SUNFALSE;
  bckPb = cvmData->bckPb;
  while (bckPb != NULL) {
    if (indexB == idxB) {
      found_bck = SUNTRUE;
      break;
    }
    bckPb = bckPb->next;
  }

  if (!found_bck) {
    cvmErrHandler(-999, "CVODES", "CVodeQuadInitB/CVodeQuadReInitB",
                  "idxB has an illegal value.", NULL);
    goto error_return;
  }

  /* 
   * ------------------------------------
   * Process inputs based on action
   * ------------------------------------
   */

  switch (action) {

  case 0:     /* BACKWARD QUADRATURE INITIALIZATION */

    /* Extract user-provided quadrature RHS function */

    mxDestroyArray(mtlb_QUADfctB);
    mtlb_QUADfctB = mxDuplicateArray(prhs[1]);

    /* Extract quadrature final conditions */

    yQB0 = mxGetPr(prhs[2]);
    NqB = mxGetM(prhs[2]);

    /* Create the backward quadrature N_Vector */

    yQB = NewVector(NqB);

    /* Load quadrature final conditions */

    PutData(yQB, yQB0, NqB);

    /* Extract quadrature options structure */

    options = prhs[3];

    break;

  case 1:     /* BACKWARD QUADRATURE RE-INITIALIZATION */


    /* Extract quadrature final conditions */

    yQB0 = mxGetPr(prhs[1]);

    if (mxGetM(prhs[1]) != NqB) cvmErrHandler(-999, "CVODES", "CVodeQuadReInitB",
                                              "Size of yQB0 changed from CVodeQuadInitB call.", NULL);

    /* Load quadrature final conditions */

    PutData(yQB, yQB0, NqB);

    /* Extract quadrature options structure */

    options = prhs[2];

    break;

  }

  /* Process the options structure */

  status = get_QuadOptions(options, bckPb, SUNFALSE,
                           NqB, &rhs_s,
                           &errconQB, 
                           &itolQB, &reltolQB, &SabstolQB, &VabstolQB);
  if (status != 0) goto error_return;

  /* 
   * ----------------------------------------
   * Call appropriate CVODES functions
   *
   * If action = 0
   *    Initialize backward quadratures
   * If action = 1
   *    Reinitialize backward quadratures
   * ----------------------------------------
   */

  switch (action) {
  case 0:
    if (rhs_s) status = CVodeQuadInitBS(cvode_mem, idxB, mxW_CVodeQUADfctBS, yQB);
    else       status = CVodeQuadInitB(cvode_mem, idxB, mxW_CVodeQUADfctB, yQB);
    if (status != CV_SUCCESS) goto error_return;
    break;
  case 1:
    status = CVodeQuadReInitB(cvode_mem, idxB, yQB);
    if (status != CV_SUCCESS) goto error_return;
    break;
  }

  /*
   * ----------------------------------------
   * Set tolerances for quadrature variables
   * ----------------------------------------
   */
  
  status = CVodeSetQuadErrConB(cvode_mem, idxB, errconQB);
  if (status != CV_SUCCESS) goto error_return;

  if (errconQB) {

    switch (itolQB) {
    case CV_SS:
      status = CVodeQuadSStolerancesB(cvode_mem, idxB, reltolQB, SabstolQB);
      if (status != CV_SUCCESS) goto error_return;
      break;
    case CV_SV:
      NV_abstolQB = N_VClone(yQB);
      PutData(NV_abstolQB, VabstolQB, NqB);
      status = CVodeQuadSVtolerancesB(cvode_mem, idxB, reltolQB, NV_abstolQB);
      if (status != CV_SUCCESS) goto error_return;
      N_VDestroy(NV_abstolQB);
      break;
    }
    
  }

  quadrB = SUNTRUE;

  /* Successful return */

  status = 0;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(0);

  /* Error return */

 error_return:
  status = -1;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(-1);

}


/*
 * CVM_Solve  - Main solution function
 *
 */

static int CVM_Solve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData fwdPb;

  int buflen;
  char *bufval;

  int nlhs_bad, dims[3];

  int status, cv_status;
  int itask, is, Ntout, itout, s_idx;
  double *tout, tret, h;
  double *tdata, *ydata, *yQdata, *ySdata;
  sunindextype nst;

  fwdPb = cvmData->fwdPb;


  /* Set index of output corresponding to FSA */

  if (fsa) {
    s_idx = quadr ? 4 : 3;
  }


  /*
   * ----------------------------------------------------------------
   * Verify if number of output arguments agrees with current options
   * ----------------------------------------------------------------
   */

  nlhs_bad = 0;

  if (nlhs < 3) nlhs_bad = -1;
  if (nlhs > 5) nlhs_bad = 1;
  if ( (nlhs == 3) && (quadr || fsa) ) nlhs_bad = -1;
  if ( (nlhs == 4) && (quadr && fsa) ) nlhs_bad = -1;
  if ( (nlhs == 5) && (!quadr || !fsa) ) nlhs_bad = 1;

  if (nlhs_bad < 0) {
    cvmErrHandler(-999, "CVODES", "CVode",
                  "Too few output arguments.", NULL);
    goto error_return;
  }

  if (nlhs_bad > 0) {
    cvmErrHandler(-999, "CVODES", "CVode",
                  "Too many output arguments.", NULL);
    goto error_return;
  }

  /*
   * ----------------------------------------------------------------
   * Extract input arguments
   * ----------------------------------------------------------------
   */

  /* Extract tout */

  Ntout = mxGetM(prhs[0]) * mxGetN(prhs[0]);
  tout = mxGetPr(prhs[0]);

  /* If rootfinding or tstop are enabled, we do not allow multiple output times */

  if (rootSet && (Ntout>1)) {
    cvmErrHandler(-999, "CVODES", "CVode",
                  "More than one tout value prohibited with rootfinding enabled.", NULL);
    goto error_return;
  }

  if (tstopSet && (Ntout>1)) {
    cvmErrHandler(-999, "CVODES", "CVode",
                  "More than one tout value prohibited with tstop enabled.", NULL);
    goto error_return;
  }

  /* Extract itask */

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Normal")) {
    itask = CV_NORMAL;
  } else if(!strcmp(bufval,"OneStep")) {
    itask = CV_ONE_STEP;
  } else {
    cvmErrHandler(-999, "CVODES", "CVode",
                  "Illegal value for itask.", NULL); 
    goto error_return;
  }

  if (itask == CV_ONE_STEP) {

    /* If itask==CV_ONE_STEP, we do not allow multiple output times and we do not monitor */

    if (Ntout > 1) {
      cvmErrHandler(-999, "CVODES", "CVode",
                    "More than one tout value prohibited in ONE_STEP mode.", NULL); 
      goto error_return;
    }

    if (mon) {
      cvmErrHandler(+999, "CVODES", "CVode",
                    "Monitoring disabled in ONE_STEP mode.", NULL);
      mon = SUNFALSE;
    }

  } else {

    /* Check if tout values are legal */

    status = CVodeGetCurrentTime(cvode_mem, &tret);
    if (status != CV_SUCCESS) goto error_return;
    status = CVodeGetNumSteps(cvode_mem, &nst);
    if (status != CV_SUCCESS) goto error_return;

    /* h is used throughout this function as integration direction only */
    if (nst == 0) {
      h = tout[0] - tret;
    } else {
      status = CVodeGetLastStep(cvode_mem, &h);
      if (status != CV_SUCCESS) goto error_return;
      if ( (tout[0] - tret + h)*h < 0.0 ) {
        cvmErrHandler(-999, "CVODES", "CVode",
                      "Illegal value of tout.", NULL);
        goto error_return;
      }
    }
    
    for (itout=1; itout<Ntout; itout++) 
      if ( (tout[itout] - tout[itout-1])*h < 0.0 ) {
        cvmErrHandler(-999, "CVODES", "CVode",
                      "tout values are not monotonic.", NULL);
        goto error_return;
      }

  }

  /*
   * ----------------------------------------------------------------
   * Prepare the output arrays
   * ----------------------------------------------------------------
   */

  /* Return time(s) */

  plhs[1] = mxCreateDoubleMatrix(1,Ntout,mxREAL);
  tdata = mxGetPr(plhs[1]);

  /* Solution vector(s) */

  plhs[2] = mxCreateDoubleMatrix(N,Ntout,mxREAL);
  ydata = mxGetPr(plhs[2]);

  /* Quadrature vector(s) */

  if (quadr) {
    plhs[3] = mxCreateDoubleMatrix(Nq,Ntout,mxREAL);
    yQdata = mxGetPr(plhs[3]);
  }

  /* Sensitivity vectors */

  if (fsa) {
    dims[0] = N;
    dims[1] = Ns;
    dims[2] = Ntout;
    plhs[s_idx] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    ySdata = mxGetPr(plhs[s_idx]);
  }

  /*
   * ----------------------------------------------------------------
   * Call the CVODES main solver function
   * ----------------------------------------------------------------
   */
 
  if (!mon) {

    /* No monitoring. itask can be either CV_ONE_STEP or CV_NORMAL */

    for (itout=0; itout<Ntout; itout++) {

      if (!asa) cv_status = CVode(cvode_mem, tout[itout], y, &tret, itask);
      else      cv_status = CVodeF(cvode_mem, tout[itout], y, &tret, itask, &Nc);
      if (cv_status < 0) goto error_return;

      tdata[itout] = tret;

      GetData(y, &ydata[itout*N], N);

      if (quadr) {
        status = CVodeGetQuad(cvode_mem, &tret, yQ);
        if (status != CV_SUCCESS) goto error_return;
        GetData(yQ, &yQdata[itout*Nq], Nq);
      }

      if (fsa) {
        status = CVodeGetSens(cvode_mem, &tret, yS);
        if (status != CV_SUCCESS) goto error_return;
        for (is=0; is<Ns; is++)
          GetData(yS[is], &ySdata[itout*Ns*N+is*N], N);
      }

    }

  } else {

    /* Monitoring. itask = CV_NORMAL */

    for (itout=0; itout<Ntout; itout++) {
      
      /* In ONE_STEP mode, CVODE reads tout only at the first step.
       * We must therefore check here whether we need to take additional steps,
       * or simply return interpolated solution at tout. */

      status = CVodeGetNumSteps(cvode_mem, &nst);
      if (status != CV_SUCCESS) goto error_return;
      status = CVodeGetCurrentTime(cvode_mem, &tret);
      if (status != CV_SUCCESS) goto error_return;

      if ( (nst>0) && ((tret - tout[itout])*h >= 0.0) ) {

        /* No need to take an additional step */
        cv_status = CV_SUCCESS;

      } else {

        /* Take additional steps */
        while(1) {

          if (!asa) cv_status = CVode(cvode_mem, tout[itout], y, &tret, CV_ONE_STEP);
          else      cv_status = CVodeF(cvode_mem, tout[itout], y, &tret, CV_ONE_STEP, &Nc);
          if (cv_status < 0) goto error_return;

          /* Call the monitoring function */
          if (quadr) {
            status = CVodeGetQuad(cvode_mem, &tret, yQ);
            if (status != CV_SUCCESS) goto error_return;
          }
          if (fsa) {
            status = CVodeGetSens(cvode_mem, &tret, yS);
            if (status != CV_SUCCESS) goto error_return;
          }
          mxW_CVodeMonitor(1, tret, y, yQ, yS, fwdPb);

          /* If a root was found or tstop was reached, break out of while loop */
          if (cv_status == CV_TSTOP_RETURN || cv_status == CV_ROOT_RETURN) break;

          /* If current tout was reached break out of while loop */
          if ( (tret - tout[itout])*h >= 0.0 )  break;

        }

      }
      
      /* On a tstop or root return, return solution at tret.
       * Otherwise (cv_status=CV_SUCCESS), return solution at tout[itout]. */

      if (cv_status == CV_TSTOP_RETURN || cv_status == CV_ROOT_RETURN) {

        if (quadr) {
          status = CVodeGetQuad(cvode_mem, &tret, yQ);
          if (status != CV_SUCCESS) goto error_return;
        }

        if (fsa) {
          status = CVodeGetSens(cvode_mem, &tret, yS);
          if (status != CV_SUCCESS) goto error_return;
        }

      } else {

        tret = tout[itout];
        
        status = CVodeGetDky(cvode_mem, tret, 0, y);
        if (status != CV_SUCCESS) goto error_return;

        if (quadr) {
          status = CVodeGetQuadDky(cvode_mem, tret, 0, yQ);
          if (status != CV_SUCCESS) goto error_return;
        }

        if (fsa) {
          status = CVodeGetSensDky(cvode_mem, tret, 0, yS);
          if (status != CV_SUCCESS) goto error_return;
        }

      }

      tdata[itout] = tret;

      GetData(y, &ydata[itout*N], N);

      if (quadr)  GetData(yQ, &yQdata[itout*Nq], Nq);

      if (fsa)
        for (is=0; is<Ns; is++)
          GetData(yS[is], &ySdata[itout*Ns*N+is*N], N);

    }

  }

  /* CVODE return flag (only non-negative values make it here) */

  plhs[0] = mxCreateDoubleScalar((double)cv_status);
  return(0);

  /* Error return */

 error_return:
  status = -1;
  plhs[0] = mxCreateDoubleScalar((double)status);
  plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(0,0,mxREAL);
  if (quadr) {
    plhs[3] = mxCreateDoubleMatrix(0,0,mxREAL);
  }
  if (fsa) {
    s_idx = quadr ? 4 : 3;
    plhs[s_idx] = mxCreateDoubleMatrix(0,0,mxREAL);
  }
  return(-1);

}


static int CVM_SolveB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData bckPb;

  int buflen;
  char *bufval;

  int nlhs_bad;

  int itaskB, NtoutB;
  double *toutB;

  double tret, h;

  booleantype any_quadrB, any_monB;

  int status, cv_status;

  /*
   * -------------------------------------------------------
   * Check whether quadratures and/or monitoring are enabled
   * for at least one backward problem
   * -------------------------------------------------------
   */

  any_quadrB = SUNFALSE;
  any_monB   = SUNFALSE;
  bckPb = cvmData->bckPb;
  while(bckPb != NULL) {
    if (quadrB) any_quadrB = SUNTRUE;
    if (monB)   any_monB   = SUNTRUE;
    bckPb = bckPb->next;
  }
  
  /*
   * ----------------------------------------------------------------
   * Verify if number of output arguments agrees with current options
   * ----------------------------------------------------------------
   */

  nlhs_bad = 0;

  if (nlhs < 3) nlhs_bad = -1;
  if (nlhs > 4) nlhs_bad = 1;
  if ( (nlhs == 3) && any_quadrB ) nlhs_bad = -1;

  if (nlhs_bad < 0) {
    cvmErrHandler(-999, "CVODES", "CVodeB",
                  "Too few output arguments.", NULL);
    goto error_return;
  }

  if (nlhs_bad > 0) {
    cvmErrHandler(-999, "CVODES", "CVodeB",
                  "Too many output arguments.", NULL);
    goto error_return;
  }

  /*
   * ----------------------------------------------------------------
   * Extract input arguments
   * ----------------------------------------------------------------
   */

  /* Extract tout */

  NtoutB = mxGetM(prhs[0]) * mxGetN(prhs[0]);
  toutB = mxGetPr(prhs[0]);

  /* Check if first tout value is in the right direction */

  status = CVodeGetLastStep(cvode_mem, &h);
  if (status != CV_SUCCESS) goto error_return;
  status = CVodeGetCurrentTime(cvode_mem, &tret);
  if (status != CV_SUCCESS) goto error_return;
 
  /* The stepsize of the forward problem is used to indicate the integration direction */
  if ( (tret - toutB[0])*h < 0.0 ) {
    cvmErrHandler(-999, "CVODES", "CVodeB",
                  "tout value in wrong direction.", NULL);
    goto error_return;
  }

  /* Extract itaskB */

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Normal")) {
    itaskB = CV_NORMAL;
  } else if(!strcmp(bufval,"OneStep")) {
    itaskB = CV_ONE_STEP;
  } else {
    cvmErrHandler(-999, "CVODES", "CVodeB",
                  "Illegal value for itask.", NULL); 
    goto error_return;
  }

  /* If itask == CV_ONE_STEP, then
   * - we do not allow multiple output times
   * - we disable monitoring 
   */

  if ( itaskB == CV_ONE_STEP ) {
    
    if (NtoutB > 1) {
      cvmErrHandler(-999, "CVODES", "CVodeB",
                    "More than one tout value prohibited in ONE_STEP mode.", NULL); 
      goto error_return;
    }

    if (any_monB) {
      cvmErrHandler(+999, "CVODES", "CVodeB",
                    "Monitoring disabled in itask=CV_ONE_STEP", NULL);
      bckPb = cvmData->bckPb;
      while(bckPb != NULL) {
        monB = SUNFALSE;
        bckPb = bckPb->next;
      }
      any_monB = SUNFALSE;
    }

  }

  /* Call the appropriate function to do all the work.
   * Note: if we made it here, we rely on the functions cvmSolveB_one and cvmSolveB_more
   *       to set the output arrays in plhs appropriately. */

  if (NbckPb == 1) cv_status = cvmSolveB_one(plhs, NtoutB, toutB, itaskB);
  else             cv_status = cvmSolveB_more(plhs, NtoutB, toutB, itaskB, any_quadrB, any_monB);

  if (cv_status < 0) return(-1);
  else               return(0);

  /* Error return */

 error_return:
  status = -1;
  plhs[0] = mxCreateDoubleScalar((double)status);
  plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(0,0,mxREAL);
  if (quadrB) {
    plhs[3] = mxCreateDoubleMatrix(0,0,mxREAL);
  }
  return(-1);

}



static int cvmSolveB_one(mxArray *plhs[], int NtoutB, double *toutB, int itaskB)
{
  cvmPbData bckPb;
  
  void *cvode_memB;

  double tretB, hB;
  double *tdata, *ydata, *yQdata;
  int itout;
  sunindextype nstB;

  int status, cv_status;

  bckPb = cvmData->bckPb;

  cvode_memB = CVodeGetAdjCVodeBmem(cvode_mem, indexB);

  /* Check if tout values are legal */

  CVodeGetCurrentTime(cvode_memB, &tretB);
  CVodeGetNumSteps(cvode_memB, &nstB);

  /* hB is used throughout this function as integration direction only */
  if (nstB == 0) {
    hB = toutB[0] - tretB;
  } else {
    CVodeGetLastStep(cvode_memB, &hB);
    if ( (toutB[0] - tretB + hB)*hB < 0.0 ) {
      cvmErrHandler(-999, "CVODES", "CVodeB",
                    "Illegal value of tout.", NULL);
      goto error_return;
    }
  }

  for (itout=1; itout<NtoutB; itout++) {
    if ( (toutB[itout] - toutB[itout-1])*hB < 0.0 ) {
      cvmErrHandler(-999, "CVODES", "CVodeB",
                    "tout values are not monotonic.", NULL);
      goto error_return;
    }
  }

  /*
   * ----------------------------------------------------------------
   * Prepare the output arrays
   * ----------------------------------------------------------------
   */

  /* Return time(s) */

  plhs[1] = mxCreateDoubleMatrix(1,NtoutB,mxREAL);
  tdata = mxGetPr(plhs[1]);
  
  /* Solution vector(s) */
  
  plhs[2] = mxCreateDoubleMatrix(NB,NtoutB,mxREAL);
  ydata = mxGetPr(plhs[2]);
  
  /* Quadrature vector(s) */
  
  if (quadrB) {
    plhs[3] = mxCreateDoubleMatrix(NqB,NtoutB,mxREAL);
    yQdata = mxGetPr(plhs[3]);
  }
  
  /*
   * ----------------------------------------------------------------
   * Call the CVodeB main solver function
   * ----------------------------------------------------------------
   */

  if (!monB) {

    /* No monitoring. itaskB can be either CV_ONE_STEP or CV_NORMAL */

    for (itout=0; itout<NtoutB; itout++) {
      
      cv_status = CVodeB(cvode_mem, toutB[itout], itaskB);
      if (cv_status < 0) goto error_return;

      status = CVodeGetB(cvode_mem, indexB, &tretB, yB);
      if (status != CV_SUCCESS) goto error_return;

      tdata[itout] = tretB;
      
      GetData(yB, &ydata[itout*NB], NB);
      
      if (quadrB) {
        status = CVodeGetQuadB(cvode_mem, indexB, &tretB, yQB);
        if (status != CV_SUCCESS) goto error_return;
        GetData(yQB, &yQdata[itout*NqB], NqB);
      }
      
    }
    
    
  } else {

    /* Monitoring. itask = CV_NORMAL */

    for (itout=0; itout<NtoutB; itout++) {

      /* In ONE_STEP mode, CVODE reads tout only at the first step.
       * We must therefore check here whether we need to take additional steps,
       * or simply return interpolated solution at tout. */

      status = CVodeGetNumSteps(cvode_memB, &nstB);
      if (status != CV_SUCCESS) goto error_return;

      status = CVodeGetCurrentTime(cvode_memB, &tretB);
      if (status != CV_SUCCESS) goto error_return;

      if ( (nstB>0) && ((tretB - toutB[itout])*hB >= 0.0) ) {

        /* No need to take an additional step */
        cv_status = CV_SUCCESS;

      } else {

        /* Take additional steps */
        while(1) {
        
          cv_status = CVodeB(cvode_mem, toutB[itout], CV_ONE_STEP);
          if (cv_status < 0) goto error_return;          

          /* Call the monitoring function */          

          status = CVodeGetB(cvode_mem, indexB, &tretB, yB);
          if (status != CV_SUCCESS) goto error_return;

          if (quadrB) {
            CVodeGetQuadB(cvode_mem, indexB, &tretB, yQB);
            if (status != CV_SUCCESS) goto error_return;
          }

          mxW_CVodeMonitorB(1, indexB, tretB, yB, yQB, bckPb);
          
          /* If current tout was reached break out of while loop */
          if ( (tretB - toutB[itout])*hB >= 0.0 ) break;
          
        }

      }

      tretB = toutB[itout];

      tdata[itout] = tretB;

      status = CVodeGetDky(cvode_memB, tretB, 0, yB);
      if (status != CV_SUCCESS) goto error_return;

      GetData(yB, &ydata[itout*NB], NB);
      
      if (quadrB) {
        status = CVodeGetQuadDky(cvode_memB, tretB, 0, yQB);
        if (status != CV_SUCCESS) goto error_return;
        GetData(yQB, &yQdata[itout*NqB], NqB);
      }

    }

  }

  /* CVodeB return flag (only non-negative values make it here) */

  plhs[0] = mxCreateDoubleScalar((double)cv_status);
  return(0);

 error_return:
  status = -1;
  plhs[0] = mxCreateDoubleScalar((double)status);
  plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(0,0,mxREAL);
  if (quadrB) {
    plhs[3] = mxCreateDoubleMatrix(0,0,mxREAL);
  }
  return(-1);

}


static int cvmSolveB_more(mxArray *plhs[], int NtoutB, double *toutB, int itaskB,
                          booleantype any_quadrB, booleantype any_monB)
{
  cvmPbData bckPb;
  mxArray *cell;
  void *cvode_memB;

  double tretB, h, hB;
  double **tdata, **ydata, **yQdata;
  int itout;
  sunindextype nstB;

  char err_msg[80];

  int status, cv_status;

  /* Check if tout values are legal */

  CVodeGetLastStep(cvode_mem, &h);

  bckPb = cvmData->bckPb;
  while (bckPb != NULL) {

    cvode_memB = CVodeGetAdjCVodeBmem(cvode_mem, indexB);

    status = CVodeGetCurrentTime(cvode_memB, &tretB);
    if (status != CV_SUCCESS) goto error_return;

    status = CVodeGetNumSteps(cvode_memB, &nstB);
    if (status != CV_SUCCESS) goto error_return;

    if (nstB > 0) {
      status = CVodeGetLastStep(cvode_memB, &hB);
      if (status != CV_SUCCESS) goto error_return;
      if ( (toutB[0] - tretB + hB)*hB < 0.0 ) {
        sprintf(err_msg, "CVodeB:: illegal value of tout (pb. %d).",indexB+1);
        cvmErrHandler(-999, "CVODES", "CVodeB", err_msg, NULL);
        goto error_return;
      }
    }
    
    for (itout=1; itout<NtoutB; itout++) {
      if ( (toutB[itout] - toutB[itout-1]) * h > 0.0 ) {
        sprintf(err_msg, "CVodeB:: tout values are not monotonic (pb. %d).", indexB+1);
        cvmErrHandler(-999, "CVODES", "CVodeB", err_msg, NULL);
        goto error_return;
      }
    }

    bckPb = bckPb->next;
    
  }

  /*
   * ----------------------------------------------------------------
   * Prepare the output arrays
   * ----------------------------------------------------------------
   */

  plhs[1] = mxCreateCellMatrix(NbckPb,1);
  tdata = (double **) mxMalloc(NbckPb*sizeof(double *));

  plhs[2] = mxCreateCellMatrix(NbckPb,1);
  ydata = (double **) mxMalloc(NbckPb*sizeof(double *));

  if (any_quadrB) {
    plhs[3] = mxCreateCellMatrix(NbckPb,1);
    yQdata = (double **) mxMalloc(NbckPb*sizeof(double *));
  }


  bckPb = cvmData->bckPb;
  while (bckPb != NULL) {

    /* Return time(s) */
    cell = mxCreateDoubleMatrix(1,NtoutB,mxREAL);
    tdata[indexB] = mxGetPr(cell);
    mxSetCell(plhs[1], indexB, cell);

    /* Solution vector(s) */
    cell = mxCreateDoubleMatrix(NB,NtoutB,mxREAL);
    ydata[indexB] = mxGetPr(cell);
    mxSetCell(plhs[2], indexB, cell);

    /* Quadrature vector(s) */
    if (any_quadrB) {
      if (quadrB) {
        cell = mxCreateDoubleMatrix(NqB,NtoutB,mxREAL);
        yQdata[indexB] = mxGetPr(cell);
      } else {
        cell = mxCreateDoubleMatrix(0,0,mxREAL);
      }
      mxSetCell(plhs[3], indexB, cell);
    }

    bckPb = bckPb->next;

  }

  
  /*
   * ----------------------------------------------------------------
   * Call the CVodeB main solver function
   * ----------------------------------------------------------------
   */

  if (!any_monB) {

    /* No monitoring. itaskB can be either CV_ONE_STEP or CV_NORMAL */

    for (itout=0; itout<NtoutB; itout++) {

      cv_status = CVodeB(cvode_mem, toutB[itout], itaskB);
      if (cv_status < 0) goto error_return;

      bckPb = cvmData->bckPb;
      while (bckPb != NULL) {

        status = CVodeGetB(cvode_mem, indexB, &tretB, yB);
        if (status != CV_SUCCESS) goto error_return;

        tdata[indexB][itout] = tretB;

        GetData(yB, &ydata[indexB][itout*NB], NB);
      
        if (quadrB) {
          status = CVodeGetQuadB(cvode_mem, indexB, &tretB, yQB);
          if (status != CV_SUCCESS) goto error_return;
          GetData(yQB, &yQdata[indexB][itout*NqB], NqB);
        }

        bckPb = bckPb->next;

      }

    }


  } else {

    /* Monitoring for at least one backward problem. itask = CV_NORMAL */

    cvmErrHandler(-999, "CVODES", "CVodeB",
                  "Monitoring currently prohibited with more than one backward problem defined .", NULL); 
    goto error_return;

  }

  /* CVODE return flag (only non-negative values make it here) */

  plhs[0] = mxCreateDoubleScalar((double)cv_status);
  return(0);

 error_return:
  status = -1;
  plhs[0] = mxCreateDoubleScalar((double)status);
  plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(0,0,mxREAL);
  if (quadrB) {
    plhs[3] = mxCreateDoubleMatrix(0,0,mxREAL);
  }
  return(-1);

}


static int CVM_Stats(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData fwdPb;

  const char *fnames_intgr[]={
    "nst",
    "nfe",
    "nge",
    "nsetups",
    "netf",
    "nni",
    "ncfn",
    "qlast",
    "qcur",
    "h0used",
    "hlast",
    "hcur",
    "tcur",
    "RootInfo",
    "QuadInfo",
    "LSInfo",
    "FSAInfo"
  };
  const char *fnames_root[]={
    "nge",
    "roots"
  };
  const char *fnames_dense[]={
    "name",
    "njeD",
    "nfeD"
  };
  const char *fnames_diag[]={
    "name",
    "nfeDI"
  };
  const char *fnames_band[]={
    "name",
    "njeB",
    "nfeB"
  };
  const char *fnames_spils[]={
    "name",
    "nli",
    "npe",
    "nps",
    "ncfl",
    "njeSG",
    "nfeSG"
  };
  const char *fnames_quad[]={
    "nfQe",
    "netfQ"
  };
  const char *fnames_sens[]={
    "nfSe",
    "nfeS",
    "nsetupsS",
    "netfS",
    "nniS",
    "ncfnS",
  };

  sunindextype nst, nfe, nsetups, nni, ncfn, netf, nge;
  int qlast, qcur;
  double h0used, hlast, hcur, tcur;
  int *rootsfound;

  sunindextype njeD, nfeD;
  sunindextype nfeDI;
  sunindextype njeB, nfeB;
  sunindextype nli, npe, nps, ncfl, njeSG, nfeSG;

  sunindextype nfQe, netfQ;

  sunindextype nfSe, nfeS, netfS, nsetupsS;
  sunindextype nniS, ncfnS;

  int i, status;
  mxArray *mxS_root, *mxS_ls, *mxS_quad, *mxS_fsa;
  mxArray *mxS_rootsfound;
  double *tmp;
  int nfields;
  

  if (cvmData == NULL) return(0);

  fwdPb = cvmData->fwdPb;

  status = CVodeGetIntegratorStats(cvode_mem, &nst, &nfe, &nsetups, 
                                   &netf, &qlast, &qcur, &h0used, &hlast, &hcur, &tcur);
  if (status != CV_SUCCESS) goto error_return;

  status = CVodeGetNonlinSolvStats(cvode_mem, &nni, &ncfn);
  if (status != CV_SUCCESS) goto error_return;

  nfields = sizeof(fnames_intgr)/sizeof(*fnames_intgr);
  plhs[0] = mxCreateStructMatrix(1, 1, nfields, fnames_intgr);
 
  mxSetField(plhs[0], 0, "nst",     mxCreateDoubleScalar((double)nst));
  mxSetField(plhs[0], 0, "nfe",     mxCreateDoubleScalar((double)nfe));
  mxSetField(plhs[0], 0, "nsetups", mxCreateDoubleScalar((double)nsetups));
  mxSetField(plhs[0], 0, "netf",    mxCreateDoubleScalar((double)netf));
  mxSetField(plhs[0], 0, "nni",     mxCreateDoubleScalar((double)nni));
  mxSetField(plhs[0], 0, "ncfn",    mxCreateDoubleScalar((double)ncfn));
  mxSetField(plhs[0], 0, "qlast",   mxCreateDoubleScalar((double)qlast));
  mxSetField(plhs[0], 0, "qcur",    mxCreateDoubleScalar((double)qcur));
  mxSetField(plhs[0], 0, "h0used",  mxCreateDoubleScalar(h0used));
  mxSetField(plhs[0], 0, "hlast",   mxCreateDoubleScalar(hlast));
  mxSetField(plhs[0], 0, "hcur",    mxCreateDoubleScalar(hcur));
  mxSetField(plhs[0], 0, "tcur",    mxCreateDoubleScalar(tcur));


  /* Root Finding Statistics */

  if (Ng > 0) {

    status = CVodeGetNumGEvals(cvode_mem, &nge);
    if (status != CV_SUCCESS) goto error_return;

    nfields = sizeof(fnames_root)/sizeof(*fnames_root);
    mxS_root = mxCreateStructMatrix(1, 1, nfields, fnames_root);

    mxSetField(mxS_root, 0, "nge", mxCreateDoubleScalar((double)nge));

    rootsfound = (int *) malloc(Ng*sizeof(int));
    status = CVodeGetRootInfo(cvode_mem, rootsfound);
    if (status != CV_SUCCESS) goto error_return;
    mxS_rootsfound = mxCreateDoubleMatrix(Ng,1,mxREAL);
    tmp = mxGetPr(mxS_rootsfound);
    for (i=0;i<Ng;i++)
      tmp[i] = (double)rootsfound[i];
    mxSetField(mxS_root, 0, "roots", mxS_rootsfound);
    free(rootsfound);

  } else {

    mxS_root = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "RootInfo", mxS_root);

  /* Quadrature Statistics */

  if (quadr) {

    status = CVodeGetQuadStats(cvode_mem, &nfQe, &netfQ);
    if (status != CV_SUCCESS) goto error_return;

    nfields = sizeof(fnames_quad)/sizeof(*fnames_quad);
    mxS_quad = mxCreateStructMatrix(1, 1, nfields, fnames_quad);

    mxSetField(mxS_quad, 0, "nfQe",  mxCreateDoubleScalar((double)nfQe));
    mxSetField(mxS_quad, 0, "netfQ", mxCreateDoubleScalar((double)netfQ));

  } else {

    mxS_quad = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "QuadInfo", mxS_quad);

  /* Linear Solver Statistics */

  switch(ls){

  case LS_NONE:

    mxS_ls = mxCreateDoubleMatrix(0,0,mxREAL);
    break;

  case LS_DENSE:
    
    status = CVDlsGetNumJacEvals(cvode_mem, &njeD);
    if (status != CV_SUCCESS) goto error_return;
    status = CVDlsGetNumRhsEvals(cvode_mem, &nfeD);
    if (status != CV_SUCCESS) goto error_return;

    nfields = sizeof(fnames_dense)/sizeof(*fnames_dense);
    mxS_ls = mxCreateStructMatrix(1, 1, nfields, fnames_dense);
    
    mxSetField(mxS_ls, 0, "name", mxCreateString("Dense"));
    mxSetField(mxS_ls, 0, "njeD", mxCreateDoubleScalar((double)njeD));
    mxSetField(mxS_ls, 0, "nfeD", mxCreateDoubleScalar((double)nfeD));
    
    break;

  case LS_DIAG:

    status = CVDiagGetNumRhsEvals(cvode_mem, &nfeDI);
    if (status != CV_SUCCESS) goto error_return;

    nfields = sizeof(fnames_diag)/sizeof(*fnames_diag);
    mxS_ls = mxCreateStructMatrix(1, 1, nfields, fnames_diag);
 
    mxSetField(mxS_ls, 0, "name", mxCreateString("Diag"));
    mxSetField(mxS_ls, 0, "nfeDI", mxCreateDoubleScalar((double)nfeDI));
      
    break;

  case LS_BAND:
      
    status = CVDlsGetNumJacEvals(cvode_mem, &njeB);
    if (status != CV_SUCCESS) goto error_return;
    status = CVDlsGetNumRhsEvals(cvode_mem, &nfeB);
    if (status != CV_SUCCESS) goto error_return;
  
    nfields = sizeof(fnames_band)/sizeof(*fnames_band);
    mxS_ls = mxCreateStructMatrix(1, 1, nfields, fnames_band);
 
    mxSetField(mxS_ls, 0, "name", mxCreateString("Band"));
    mxSetField(mxS_ls, 0, "njeB", mxCreateDoubleScalar((double)njeB));
    mxSetField(mxS_ls, 0, "nfeB", mxCreateDoubleScalar((double)nfeB));
      
    break;

  case LS_SPGMR:
  case LS_SPBCG:
  case LS_SPTFQMR:

    status = CVSpilsGetNumLinIters(cvode_mem, &nli);
    if (status != CV_SUCCESS) goto error_return;
    status = CVSpilsGetNumPrecEvals(cvode_mem, &npe);
    if (status != CV_SUCCESS) goto error_return;
    status = CVSpilsGetNumPrecSolves(cvode_mem, &nps);
    if (status != CV_SUCCESS) goto error_return;
    status = CVSpilsGetNumConvFails(cvode_mem, &ncfl);
    if (status != CV_SUCCESS) goto error_return;
    status = CVSpilsGetNumJtimesEvals(cvode_mem, &njeSG);
    if (status != CV_SUCCESS) goto error_return;
    status = CVSpilsGetNumRhsEvals(cvode_mem, &nfeSG);
    if (status != CV_SUCCESS) goto error_return;
    
    nfields = sizeof(fnames_spils)/sizeof(*fnames_spils);
    mxS_ls = mxCreateStructMatrix(1, 1, nfields, fnames_spils);
    
    if (ls == LS_SPGMR)
      mxSetField(mxS_ls, 0, "name",  mxCreateString("GMRES"));
    else if (ls == LS_SPBCG)
      mxSetField(mxS_ls, 0, "name",  mxCreateString("BiCGStab"));
    else
      mxSetField(mxS_ls, 0, "name",  mxCreateString("TFQMR"));

    mxSetField(mxS_ls, 0, "nli",   mxCreateDoubleScalar((double)nli));
    mxSetField(mxS_ls, 0, "npe",   mxCreateDoubleScalar((double)npe));
    mxSetField(mxS_ls, 0, "nps",   mxCreateDoubleScalar((double)nps));
    mxSetField(mxS_ls, 0, "ncfl",  mxCreateDoubleScalar((double)ncfl));
    mxSetField(mxS_ls, 0, "njeSG", mxCreateDoubleScalar((double)njeSG));
    mxSetField(mxS_ls, 0, "nfeSG", mxCreateDoubleScalar((double)nfeSG));
    
    break;
    
  }

  mxSetField(plhs[0], 0, "LSInfo", mxS_ls);

  /* Forward Sensitivity Statistics */

  if (fsa) {

    status = CVodeGetSensStats(cvode_mem, &nfSe, &nfeS, &netfS, &nsetupsS); 
    if (status != CV_SUCCESS) goto error_return;

    status = CVodeGetSensNonlinSolvStats(cvode_mem, &nniS, &ncfnS);
    if (status != CV_SUCCESS) goto error_return;

    nfields = sizeof(fnames_sens)/sizeof(*fnames_sens);
    mxS_fsa = mxCreateStructMatrix(1, 1, nfields, fnames_sens);
    
    mxSetField(mxS_fsa, 0, "nfSe",     mxCreateDoubleScalar((double)nfSe));
    mxSetField(mxS_fsa, 0, "nfeS",     mxCreateDoubleScalar((double)nfeS));
    mxSetField(mxS_fsa, 0, "nsetupsS", mxCreateDoubleScalar((double)nsetupsS));
    mxSetField(mxS_fsa, 0, "netfS",    mxCreateDoubleScalar((double)netfS));
    mxSetField(mxS_fsa, 0, "nniS",     mxCreateDoubleScalar((double)nniS));
    mxSetField(mxS_fsa, 0, "ncfnS",    mxCreateDoubleScalar((double)ncfnS));
    
  } else {

    mxS_fsa = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "FSAInfo", mxS_fsa);

  /* Successful return */

  status = 0;
  plhs[1] = mxCreateDoubleScalar((double)status);
  return(0);

  /* Error return */

 error_return:
  status = -1;
  plhs[1] = mxCreateDoubleScalar((double)status);
  return(-1);

}

static int CVM_StatsB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData bckPb;
  int idxB;

  const char *fnames_intgr[]={
    "nst",
    "nfe",
    "nge",
    "nsetups",
    "netf",
    "nni",
    "ncfn",
    "qlast",
    "qcur",
    "h0used",
    "hlast",
    "hcur",
    "tcur",
    "QuadInfo",
    "LSInfo",
  };
  const char *fnames_dense[]={
    "name",
    "njeD",
    "nfeD"
  };
  const char *fnames_diag[]={
    "name",
    "nfeDI"
  };
  const char *fnames_band[]={
    "name",
    "njeB",
    "nfeB"
  };
  const char *fnames_spils[]={
    "name",
    "nli",
    "npe",
    "nps",
    "ncfl",
    "njeSG",
    "nfeSG"
  };
  const char *fnames_quad[]={
    "nfQe",
    "netfQ"
  };

  void *cvode_memB;

  sunindextype nst, nfe, nsetups, nni, ncfn, netf;
  int qlast, qcur;
  double h0used, hlast, hcur, tcur;

  sunindextype njeD, nfeD;
  sunindextype nfeDI;
  sunindextype njeB, nfeB;
  sunindextype nli, npe, nps, ncfl, njeSG, nfeSG;

  sunindextype nfQe, netfQ;

  mxArray *mxS_ls, *mxS_quad;
  int nfields;

  booleantype found_bck;

  int status;

  /* Extract index of current backward problem */
    
  idxB = (int)mxGetScalar(prhs[0]);

  /* Find current backward problem */

  found_bck = SUNFALSE;
  bckPb = cvmData->bckPb;
  while (bckPb != NULL) {
    if (indexB == idxB) {
      found_bck = SUNTRUE;
      break;
    }
    bckPb = bckPb->next;
  }

  if (!found_bck) cvmErrHandler(-999, "CVODES", "CVodeGetStatsB",
                                "idxB has an illegal value.", NULL);

  cvode_memB = CVodeGetAdjCVodeBmem(cvode_mem, indexB);

  status = CVodeGetIntegratorStats(cvode_memB, &nst, &nfe, &nsetups, 
                                 &netf, &qlast, &qcur, &h0used, &hlast, &hcur, &tcur);
  if (status != CV_SUCCESS) goto error_return;

  status = CVodeGetNonlinSolvStats(cvode_memB, &nni, &ncfn);
  if (status != CV_SUCCESS) goto error_return;

  nfields = sizeof(fnames_intgr)/sizeof(*fnames_intgr);
  plhs[0] = mxCreateStructMatrix(1, 1, nfields, fnames_intgr);

  mxSetField(plhs[0], 0, "nst",     mxCreateDoubleScalar((double)nst));
  mxSetField(plhs[0], 0, "nfe",     mxCreateDoubleScalar((double)nfe));
  mxSetField(plhs[0], 0, "nsetups", mxCreateDoubleScalar((double)nsetups));
  mxSetField(plhs[0], 0, "netf",    mxCreateDoubleScalar((double)netf));
  mxSetField(plhs[0], 0, "nni",     mxCreateDoubleScalar((double)nni));
  mxSetField(plhs[0], 0, "ncfn",    mxCreateDoubleScalar((double)ncfn));
  mxSetField(plhs[0], 0, "qlast",   mxCreateDoubleScalar((double)qlast));
  mxSetField(plhs[0], 0, "qcur",    mxCreateDoubleScalar((double)qcur));
  mxSetField(plhs[0], 0, "h0used",  mxCreateDoubleScalar(h0used));
  mxSetField(plhs[0], 0, "hlast",   mxCreateDoubleScalar(hlast));
  mxSetField(plhs[0], 0, "hcur",    mxCreateDoubleScalar(hcur));
  mxSetField(plhs[0], 0, "tcur",    mxCreateDoubleScalar(tcur));


  /* Quadrature Statistics */

  if (quadrB) {

    status = CVodeGetQuadStats(cvode_memB, &nfQe, &netfQ);
    if (status != CV_SUCCESS) goto error_return;

    nfields = sizeof(fnames_quad)/sizeof(*fnames_quad);
    mxS_quad = mxCreateStructMatrix(1, 1, nfields, fnames_quad);

    mxSetField(mxS_quad, 0, "nfQe",  mxCreateDoubleScalar((double)nfQe));
    mxSetField(mxS_quad, 0, "netfQ", mxCreateDoubleScalar((double)netfQ));

  } else {

    mxS_quad = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "QuadInfo", mxS_quad);

  /* Linear Solver Statistics */

  switch(lsB){

  case LS_NONE:

    mxS_ls = mxCreateDoubleMatrix(0,0,mxREAL);
    break;

  case LS_DENSE:
    
    status = CVDlsGetNumJacEvals(cvode_memB, &njeD);
    if (status != CV_SUCCESS) goto error_return;
    status = CVDlsGetNumRhsEvals(cvode_memB, &nfeD);
    if (status != CV_SUCCESS) goto error_return;
    
    nfields = sizeof(fnames_dense)/sizeof(*fnames_dense);
    mxS_ls = mxCreateStructMatrix(1, 1, nfields, fnames_dense);

    mxSetField(mxS_ls, 0, "name", mxCreateString("Dense"));
    mxSetField(mxS_ls, 0, "njeD", mxCreateDoubleScalar((double)njeD));
    mxSetField(mxS_ls, 0, "nfeD", mxCreateDoubleScalar((double)nfeD));
    
    break;

  case LS_DIAG:

    status = CVDiagGetNumRhsEvals(cvode_memB, &nfeDI);
    if (status != CV_SUCCESS) goto error_return;
      
    nfields = sizeof(fnames_diag)/sizeof(*fnames_diag);
    mxS_ls = mxCreateStructMatrix(1, 1, nfields, fnames_diag);
 
    mxSetField(mxS_ls, 0, "name", mxCreateString("Diag"));
    mxSetField(mxS_ls, 0, "nfeDI", mxCreateDoubleScalar((double)nfeDI));
      
    break;

  case LS_BAND:
      
    status = CVDlsGetNumJacEvals(cvode_memB, &njeB);
    if (status != CV_SUCCESS) goto error_return;
    status = CVDlsGetNumRhsEvals(cvode_memB, &nfeB);
    if (status != CV_SUCCESS) goto error_return;
      
    nfields = sizeof(fnames_band)/sizeof(*fnames_band);
    mxS_ls = mxCreateStructMatrix(1, 1, nfields, fnames_band);
 
    mxSetField(mxS_ls, 0, "name", mxCreateString("Band"));
    mxSetField(mxS_ls, 0, "njeB", mxCreateDoubleScalar((double)njeB));
    mxSetField(mxS_ls, 0, "nfeB", mxCreateDoubleScalar((double)nfeB));
      
    break;

  case LS_SPGMR:
  case LS_SPBCG:
  case LS_SPTFQMR:

    status = CVSpilsGetNumLinIters(cvode_memB, &nli);
    if (status != CV_SUCCESS) goto error_return;
    status = CVSpilsGetNumPrecEvals(cvode_memB, &npe);
    if (status != CV_SUCCESS) goto error_return;
    status = CVSpilsGetNumPrecSolves(cvode_memB, &nps);
    if (status != CV_SUCCESS) goto error_return;
    status = CVSpilsGetNumConvFails(cvode_memB, &ncfl);
    if (status != CV_SUCCESS) goto error_return;
    status = CVSpilsGetNumJtimesEvals(cvode_memB, &njeSG);
    if (status != CV_SUCCESS) goto error_return;
    status = CVSpilsGetNumRhsEvals(cvode_memB, &nfeSG);
    if (status != CV_SUCCESS) goto error_return;
    
    nfields = sizeof(fnames_spils)/sizeof(*fnames_spils);
    mxS_ls = mxCreateStructMatrix(1, 1, nfields, fnames_spils);
 
    if (lsB == LS_SPGMR)
      mxSetField(mxS_ls, 0, "name",  mxCreateString("GMRES"));
    else if (lsB == LS_SPBCG)
      mxSetField(mxS_ls, 0, "name",  mxCreateString("BiCGStab"));
    else
      mxSetField(mxS_ls, 0, "name",  mxCreateString("TFQMR"));

    mxSetField(mxS_ls, 0, "nli",   mxCreateDoubleScalar((double)nli));
    mxSetField(mxS_ls, 0, "npe",   mxCreateDoubleScalar((double)npe));
    mxSetField(mxS_ls, 0, "nps",   mxCreateDoubleScalar((double)nps));
    mxSetField(mxS_ls, 0, "ncfl",  mxCreateDoubleScalar((double)ncfl));
    mxSetField(mxS_ls, 0, "njeSG", mxCreateDoubleScalar((double)njeSG));
    mxSetField(mxS_ls, 0, "nfeSG", mxCreateDoubleScalar((double)nfeSG));
    
    break;
  }

  mxSetField(plhs[0], 0, "LSInfo", mxS_ls);

  /* Successful return */

  status = 0;
  plhs[1] = mxCreateDoubleScalar((double)status);
  return(0);

  /* Error return */

 error_return:
  status = -1;
  plhs[1] = mxCreateDoubleScalar((double)status);
  return(-1);

}


static int CVM_Set(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData fwdPb;

  const mxArray *options;
  mxArray *opt;

  double tstop;

  int status;

  fwdPb = cvmData->fwdPb;
  options = prhs[0];

  /* Return now if options was empty */

  if (mxIsEmpty(options)) return(0);

  /* User data */

  opt = mxGetField(options,0,"UserData");
  if ( !mxIsEmpty(opt) ) {
    mxDestroyArray(mtlb_data);
    mtlb_data = mxDuplicateArray(opt);
  }

  /* Stopping time */

  opt = mxGetField(options,0,"StopTime");
  if ( !mxIsEmpty(opt) ) {
    tstop = (double)mxGetScalar(opt);
    status = CVodeSetStopTime(cvode_mem, tstop);
    if (status != CV_SUCCESS) goto error_return;
  }

  /* Successful return */

  status = 0;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(0);

  /* Error return */

 error_return:
  status = -1;
  plhs[0] = mxCreateDoubleScalar((double)status);
  return(-1);

}

static int CVM_SetB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}

static int CVM_Get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData fwdPb;

  double t;
  N_Vector ewt;
  double *this, *next;
  int key, k, i, nfields;

  CVadjCheckPointRec *ckpnt;
  const char *fnames_ckpnt[]={
    "t0",
    "t1",
    "nstep",
    "order",
    "step"
  };

  int status;


  fwdPb = cvmData->fwdPb;

  key = (int) (*mxGetPr(prhs[0]));

  switch (key) {
  case 1:    /* DerivSolution */

    t = *mxGetPr(prhs[1]);
    k = (int) (*mxGetPr(prhs[2]));

    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    status = CVodeGetDky(cvode_mem, t, k, y);
    if (status != CV_SUCCESS) goto error_return;
    GetData(y, mxGetPr(plhs[0]), N);

    break;

  case 2:    /* ErrorWeights */

    ewt = N_VClone(y);

    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    status = CVodeGetErrWeights(cvode_mem, ewt);
    if (status != CV_SUCCESS) goto error_return;
    GetData(ewt, mxGetPr(plhs[0]), N);

    N_VDestroy(ewt);

    break;

  case 3:    /* not used */

    break;

  case 4:    /* CheckPointsInfo */

    ckpnt = (CVadjCheckPointRec *) malloc ( (Nc+1)*sizeof(CVadjCheckPointRec));
    status = CVodeGetAdjCheckPointsInfo(cvode_mem, ckpnt);
    if (status != CV_SUCCESS) {
      free(ckpnt);
      goto error_return;
    }
    nfields = sizeof(fnames_ckpnt)/sizeof(*fnames_ckpnt);
    plhs[0] = mxCreateStructMatrix(Nc+1, 1, nfields, fnames_ckpnt);
    for (i=0; i<=Nc; i++) {
      this = (double *)(ckpnt[Nc-i].my_addr);
      next = (double *)(ckpnt[Nc-i].next_addr);
      mxSetField(plhs[0], i, "t0",    mxCreateDoubleScalar((double)(ckpnt[Nc-i].t0)));
      mxSetField(plhs[0], i, "t1",    mxCreateDoubleScalar((double)(ckpnt[Nc-i].t1)));
      mxSetField(plhs[0], i, "nstep", mxCreateDoubleScalar((double)(ckpnt[Nc-i].nstep)));
      mxSetField(plhs[0], i, "order", mxCreateDoubleScalar((double)(ckpnt[Nc-i].order)));
      mxSetField(plhs[0], i, "step",  mxCreateDoubleScalar((double)(ckpnt[Nc-i].step)));
    }
    free(ckpnt);
    break;

  }

  /* Successful return */

  status = 0;
  plhs[1] = mxCreateDoubleScalar((double)status);
  return(0);

  /* Error return */

 error_return:
  status = -1;
  plhs[1] = mxCreateDoubleScalar((double)status);
  return(-1);

}

static int CVM_Free(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData fwdPb, bckPb;

  if (cvmData == NULL) return(0);

  fwdPb = cvmData->fwdPb;
  if (mon) mxW_CVodeMonitor(2, 0.0, NULL, NULL, NULL, fwdPb);

  bckPb = cvmData->bckPb;
  while (bckPb != NULL) {
    if (monB) mxW_CVodeMonitorB(2, indexB, 0.0, NULL, NULL, bckPb);   
    bckPb = bckPb->next;
  }

  CVodeFree(&cvode_mem);

  return(0);
}


