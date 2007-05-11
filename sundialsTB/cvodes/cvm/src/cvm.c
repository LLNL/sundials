/*
 * -----------------------------------------------------------------
 * $Revision: 1.17 $
 * $Date: 2007-05-11 21:42:53 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials-x.y.z/src/cvodes/LICENSE.
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

static void cvmSolveB_one(mxArray *plhs[], int NtoutB, double *toutB, int itaskB);
static void cvmSolveB_more(mxArray *plhs[], int NtoutB, double *toutB, int itaskB,
                           booleantype any_quadrB, booleantype any_monB);


static int CVM_Stats(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_StatsB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

static int CVM_Get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_Set(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

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
    33 - set one optional input at a time (NYI)

    40 - finalize
  */

  mode = (int)mxGetScalar(prhs[0]);

  mexUnlock();


  if ( (mode != 1) && (cvmData == NULL) ) {
    cvmErrHandler(-999, "CVODES", "-",
                  "Illegal attempt to call before CVodeInit.", NULL);
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

    /* Memory deallocation function */

  case 40:
    CVM_Free(nlhs, plhs, nrhs-1, &prhs[1]);
    cvmFinalCVODESdata();
    return;

  }

  /* Unless this was the CvodeFree call,
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
  mxArray *empty;

  empty = mxCreateDoubleMatrix(0,0,mxREAL);

  /* Allocate space for global CVODES data structure */
  cvmData = (cvmInterfaceData) mxMalloc(sizeof(struct cvmInterfaceData_));

  /* Initialize global CVODES data */

  cvmData->cvode_mem = NULL;

  cvmData->fwdPb     = NULL;
  cvmData->bckPb     = NULL;

  cvmData->NbckPb    = 0;

  cvmData->mx_data   = mxDuplicateArray(empty);

  cvmData->Nd        = 0;
  cvmData->Nc        = 0;
  cvmData->asa       = FALSE;

  mxDestroyArray(empty);

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

  pb->Quadr = FALSE;
  pb->Fsa   = FALSE;
  pb->Mon   = FALSE;

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

  pb->index = 0;
  pb->next  = NULL;

  mxDestroyArray(empty);
}


static void cvmPersistCVODESdata()
{
  cvmPbData tmpPb;

  /* Make global memory persistent */

  mexMakeArrayPersistent(cvmData->mx_data);

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

  CVodeFree(&(cvmData->cvode_mem));

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

  mxDestroyArray(cvmData->mx_data);

  mxFree(cvmData);
  cvmData = NULL;

  return;
}


static void cvmFinalPbData(cvmPbData pb)
{
  if (pb->Y != NULL) N_VDestroy(pb->Y);
  if (pb->YQ != NULL) N_VDestroy(pb->YQ);
  if (pb->YS != NULL) N_VDestroyVectorArray(pb->YS, pb->ns);

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
                   char *msg, void *f_data)
{
  char err_msg[256];

  if (error_code > 0) {
    sprintf(err_msg,"Warning in ==> %s\n%s",function,msg);
    mexWarnMsgTxt(err_msg);    
  } else if (error_code < 0) {
    cvmFinalCVODESdata();
    mexUnlock();
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

#define mx_data     (cvmData->mx_data)


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

#define mx_RHSfct   (fwdPb->RHSfct)
#define mx_QUADfct  (fwdPb->QUADfct)
#define mx_JACfct   (fwdPb->JACfct)
#define mx_PSETfct  (fwdPb->PSETfct)
#define mx_PSOLfct  (fwdPb->PSOLfct)
#define mx_GLOCfct  (fwdPb->GLOCfct)
#define mx_GCOMfct  (fwdPb->GCOMfct)
#define mx_Gfct     (fwdPb->Gfct)
#define mx_SRHSfct  (fwdPb->SRHSfct)

#define mx_MONfct   (fwdPb->MONfct)
#define mx_MONdata  (fwdPb->MONdata)



#define indexB      (bckPb->index)

#define quadrB      (bckPb->Quadr)
#define monB        (bckPb->Mon)

#define yB          (bckPb->Y) 
#define yQB         (bckPb->YQ) 
#define NB          (bckPb->n) 
#define NqB         (bckPb->nq) 
#define lsB         (bckPb->LS) 
#define pmB         (bckPb->PM) 

#define mx_RHSfctB  (bckPb->RHSfct)
#define mx_QUADfctB (bckPb->QUADfct)
#define mx_JACfctB  (bckPb->JACfct)
#define mx_PSETfctB (bckPb->PSETfct)
#define mx_PSOLfctB (bckPb->PSOLfct)
#define mx_GLOCfctB (bckPb->GLOCfct)
#define mx_GCOMfctB (bckPb->GCOMfct)

#define mx_MONfctB  (bckPb->MONfct)
#define mx_MONdataB (bckPb->MONdata)



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

  long int mxsteps;

  int itol;
  realtype reltol, Sabstol, *Vabstol;
  N_Vector NV_abstol;

  double *yQ0;

  double hin, hmax, hmin;

  double tstop;

  booleantype sld;

  booleantype rhs_s; /* ignored */

  int mupper, mlower;
  int ptype, gstype, maxl;
  int mudq, mldq;
  double dqrely;

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

    mxDestroyArray(mx_RHSfct);
    mx_RHSfct = mxDuplicateArray(prhs[0]);

    /* Extract initial time */

    t0 = (double)mxGetScalar(prhs[1]);

    /* Extract initial conditions */

    y0 = mxGetPr(prhs[2]);
    N = mxGetM(prhs[2]);

    /* Create the solution N_Vector */

    y = NewVector(N);

    /* Load initial conditions */

    PutData(y, y0, N);

    /* Extract options structure */
    
    options = prhs[3];

    /* Extract user-data -- optional argument */

    mxDestroyArray(mx_data);
    mx_data = mxDuplicateArray(prhs[4]);

    break;

  case 1:    /* SOLVER RE-INITIALIZATION */

    fwdPb = cvmData->fwdPb;

    /* If monitoring was enabled, finalize it now. */

    if (mon) mtlb_CVodeMonitor(2, 0.0, NULL, NULL, NULL, cvmData);

    /* Extract initial time */

    t0 = (double)mxGetScalar(prhs[0]);

    /* Extract initial conditions */

    y0 = mxGetPr(prhs[1]);
    N = mxGetM(prhs[1]);

    /* Load initial conditions */

    PutData(y, y0, N);

    /* Extract options structure */
    
    options = prhs[2];

    break;

  }

  /* Process the options structure */

  get_IntgrOptions(options, fwdPb, TRUE,
                   &lmm, &iter, &maxord, &sld, &mxsteps,
                   &itol, &reltol, &Sabstol, &Vabstol,
                   &hin, &hmax, &hmin, &tstop, &rhs_s);

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
    /* Attach the global CVODES data as 'user-data' */
    CVodeSetUserData(cvode_mem, cvmData);
    /* Attach error handler function */
    CVodeSetErrHandlerFn(cvode_mem, cvmErrHandler);
    /* Call CVodeInit */
    CVodeInit(cvode_mem, mtlb_CVodeRhs, t0, y);
    /* Redirect output */
    CVodeSetErrFile(cvode_mem, stdout);

    break;

  case 1:

    /* Reinitialize solver */
    CVodeReInit(cvode_mem, t0, y);

    break;

  }

  /*
   * ----------------------------------------
   * Set tolerances
   * ----------------------------------------
   */

  switch (itol) {
    case CV_SS:
      CVodeSStolerances(cvode_mem, reltol, Sabstol);
      break;
    case CV_SV:
      NV_abstol = N_VClone(y);
      PutData(NV_abstol, Vabstol, N);
      CVodeSVtolerances(cvode_mem, reltol, NV_abstol);
      N_VDestroy(NV_abstol);
      break;
    }
  
  /*
   * --------------------------------
   * Set various optional inputs
   * --------------------------------
   */

  /* set maxorder (default is consistent with LMM) */
  CVodeSetMaxOrd(cvode_mem, maxord);

  /* set initial step size (the default value of 0.0 is ignored by CVODES) */
  CVodeSetInitStep(cvode_mem, hin);

  /* set max step (default is infinity) */
  CVodeSetMaxStep(cvode_mem, hmax);

  /* set min step (default is 0) */
  CVodeSetMinStep(cvode_mem, hmin);
 
  /* set number of max steps */
  CVodeSetMaxNumSteps(cvode_mem, mxsteps);

  /* set tstop? */
  if (tstopSet) {
    CVodeSetStopTime(cvode_mem, tstop);
  }

  /* set stability limit detection (default is FALSE) */
  CVodeSetStabLimDet(cvode_mem, sld);
 
  /* Rootfinding? */
  if ( !mxIsEmpty(mx_Gfct) && (Ng > 0) ) {
    CVodeRootInit(cvode_mem, Ng, mtlb_CVodeGfct);
    rootSet = TRUE;
  } else {
    rootSet = FALSE;
  }

  /*
   * ----------------------------------------
   * Need a linear solver?
   * ----------------------------------------
   */

  if (iter == CV_NEWTON) {

    get_LinSolvOptions(options, fwdPb, TRUE,
                       &mupper, &mlower,
                       &mudq, &mldq, &dqrely,
                       &ptype, &gstype, &maxl);

    switch (ls) {

    case LS_DENSE:

      CVDense(cvode_mem, N);
      if (!mxIsEmpty(mx_JACfct)) CVDlsSetDenseJacFn(cvode_mem, mtlb_CVodeDenseJac);

      break;

    case LS_DIAG:

      CVDiag(cvode_mem);

      break;

    case LS_BAND:

      CVBand(cvode_mem, N, mupper, mlower);
      if (!mxIsEmpty(mx_JACfct)) CVDlsSetBandJacFn(cvode_mem, mtlb_CVodeBandJac);

      break;

    case LS_SPGMR:

      CVSpgmr(cvode_mem, ptype, maxl);
      CVSpilsSetGSType(cvode_mem, gstype);

      break;

    case LS_SPBCG:

      CVSpbcg(cvode_mem, ptype, maxl);

      break;

    case LS_SPTFQMR:

      CVSptfqmr(cvode_mem, ptype, maxl);

      break;

    }

    /* Jacobian * vector and preconditioner for SPILS linear solvers */

    if ( (ls==LS_SPGMR) || (ls==LS_SPBCG) || (ls==LS_SPTFQMR) ) {

      if (!mxIsEmpty(mx_JACfct)) CVSpilsSetJacTimesVecFn(cvode_mem, mtlb_CVodeSpilsJac);

      switch (pm) {

      case PM_NONE:

        if (!mxIsEmpty(mx_PSOLfct)) {
          if (!mxIsEmpty(mx_PSETfct)) CVSpilsSetPreconditioner(cvode_mem, mtlb_CVodeSpilsPset, mtlb_CVodeSpilsPsol);
          else                        CVSpilsSetPreconditioner(cvode_mem, NULL, mtlb_CVodeSpilsPsol);
        }

        break;

      case PM_BANDPRE:

        CVBandPrecInit(cvode_mem, N, mupper, mlower);

        break;

      case PM_BBDPRE:

        if (!mxIsEmpty(mx_GCOMfct)) CVBBDPrecInit(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_CVodeBBDgloc, mtlb_CVodeBBDgcom);
        else                        CVBBDPrecInit(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_CVodeBBDgloc, NULL);

        break;
      }

    }

  } else {

    ls = LS_NONE;

  }

  /* Do we monitor? */
  
  if (mon) mtlb_CVodeMonitor(0, t0, NULL, NULL, NULL, cvmData);

  return(0);

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

  fwdPb = cvmData->fwdPb;

  /* 
   * ------------------------------------
   * Process inputs based on action
   * ------------------------------------
   */

  switch (action) {

  case 0:     /* QUADRATURE INITIALIZATION */

    /* Extract user-provided quadrature RHS function */

    mxDestroyArray(mx_QUADfct);
    mx_QUADfct = mxDuplicateArray(prhs[0]);
  
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
    Nq = mxGetM(prhs[0]);

    /* Load quadrature initial conditions */
    
    PutData(yQ, yQ0, Nq);

    /* Extract quadrature options structure */

    options = prhs[1];

    break;

  }

  /* Process the options structure */

  get_QuadOptions(options, fwdPb, TRUE,
                  Nq, &rhs_s,
                  &errconQ, 
                  &itolQ, &reltolQ, &SabstolQ, &VabstolQ);

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
    CVodeQuadInit(cvode_mem, mtlb_CVodeQUADfct, yQ);
    break;
  case 1:
    CVodeQuadReInit(cvode_mem, yQ);
    break;
  }

  /*
   * ----------------------------------------
   * Set tolerances for quadrature variables
   * ----------------------------------------
   */

  CVodeSetQuadErrCon(cvode_mem, errconQ);

  if (errconQ) {
    
    switch (itolQ) {
    case CV_SS:
      CVodeQuadSStolerances(cvode_mem, reltolQ, SabstolQ);
      break;
    case CV_SV:
      NV_abstolQ = N_VClone(yQ);
      PutData(NV_abstolQ, VabstolQ, Nq);
      CVodeQuadSVtolerances(cvode_mem, reltolQ, NV_abstolQ);
      N_VDestroy(NV_abstolQ);
      break;
    }
    
  }

  /* Quadratures will be integrated */

  quadr = TRUE;

  return(0);

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

  int is;

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

    mxDestroyArray(mx_SRHSfct);

    if ( mxIsEmpty(prhs[1]) ) {
      fS_DQ = TRUE;
    } else {
      mx_SRHSfct = mxDuplicateArray(prhs[1]);
      fS_DQ = FALSE;
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

    /* Load sensitivity initial conditions */

    for (is=0;is<Ns;is++)
      PutData(yS[is], &yS0[is*N], N);

    /* Extract qFSA options structure */

    options = prhs[1];

    break;

  }

  /* Process the options structure */

  get_FSAOptions(options, fwdPb, 
                 &ism,
                 &pfield_name, &plist, &pbar,
                 &dqtype, &rho,
                 &errconS, &itolS, &reltolS, &SabstolS, &VabstolS);

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

    if (fS_DQ) {

      if (pfield_name == NULL) cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                                             "CVode FSA Initialization:: pfield required but was not provided.", NULL);

      pfield = mxGetField(mx_data,0,pfield_name);

      if (pfield == NULL) cvmErrHandler(-999, "CVODES", "CVodeSensInit/CVodeSensReInit",
                                        "CVode FSA Initialization:: illegal pfield input.", NULL);

      p = mxGetPr(pfield);

    }

    CVodeSensInit(cvode_mem, Ns, ism, mtlb_CVodeSensRhs, yS);

    break;

  case 1:

    CVodeSensReInit(cvode_mem, ism, yS);

    break;

  }

  /*
   * ----------------------------------------
   * Set tolerances for sensitivity variables
   * ----------------------------------------
   */

  switch (itolS) {
  case CV_SS:
    CVodeSensSStolerances(cvode_mem, reltolS, SabstolS);
    break;
  case CV_SV:
    NV_abstolS = N_VCloneVectorArray(Ns, y);
    for (is=0;is<Ns;is++)
      PutData(NV_abstolS[is], &VabstolS[is*N], N);
    CVodeSensSVtolerances(cvode_mem, reltolS, NV_abstolS);
    N_VDestroyVectorArray(NV_abstolS, Ns);
    break;
  case CV_EE:
    CVodeSensEEtolerances(cvode_mem);
    break;
  }

  /*
   * --------------------------------
   * Set various optional inputs
   * --------------------------------
   */
  
  CVodeSetSensParams(cvode_mem, p, pbar, plist);

  CVodeSetSensDQMethod(cvode_mem, dqtype, rho);
  
  CVodeSetSensErrCon(cvode_mem, errconS);

  fsa = TRUE;

  return(0);
}

/*
 * CVM_SensToggleOff
 *
 * deactivates FSA
 */

static int CVM_SensToggleOff(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData fwdPb;

  fwdPb = cvmData->fwdPb;

  CVodeSensToggleOff(cvode_mem);
  
  fsa = FALSE;

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
    if(status != 0) cvmErrHandler(-999, "CVODES", "CVodeAdjInit", 
                                  "Could not parse InterpType.", NULL);

    if(!strcmp(bufval,"Hermite"))         interp = CV_HERMITE;
    else if(!strcmp(bufval,"Polynomial")) interp = CV_POLYNOMIAL;
    else cvmErrHandler(-999, "CVODES", "CvoeAdjInit",
                       "Interp. type has an illegal value.", NULL);

    CVodeAdjInit(cvode_mem, Nd, interp);

    break;

  case 1:

    CVodeAdjReInit(cvode_mem);

    break;

  }

  asa = TRUE;
  
  return(0);
}


/* CVM_InitializationB
 *
 * action = 0   -> CVodeCreateB + CVodeInitB
 * prhs contains:
 *   fctB
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
  long int mxstepsB;

  int itolB;
  realtype reltolB, SabstolB, *VabstolB;
  N_Vector NV_abstolB;

  double hinB, hmaxB, hminB;

  double tstopB;            /* ignored */
  booleantype sldB;         /* ignored */

  booleantype rhs_s;

  int mupperB, mlowerB;
  int ptypeB, gstypeB, maxlB;
  int mudqB, mldqB;
  double dqrelyB;


  /* 
   * -----------------------------
   * Finalize Forward monitoring
   * -----------------------------
   */

  if (cvmData->fwdPb->Mon) {
    mtlb_CVodeMonitor(2, 0.0, NULL, NULL, NULL, cvmData);
    cvmData->fwdPb->Mon = FALSE;
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

    mxDestroyArray(mx_RHSfctB);
    mx_RHSfctB = mxDuplicateArray(prhs[0]);

    /* Extract final time */

    tB0 = (double)mxGetScalar(prhs[1]);

    /* Extract final conditions */

    yB0 = mxGetPr(prhs[2]);
    NB = mxGetM(prhs[2]);

    /* Create the solution N_Vector */

    yB = NewVector(NB);

    /* Load final conditions */

    PutData(yB, yB0, NB);

    /* Extract options structure */
    
    options = prhs[3];

    break;

  case 1:     /* BACKWARD SOLVER RE-INITIALIZATION */

    bckPb = cvmData->bckPb;

    /* If backward monitoring was enabled, finalize it now. */

    if (monB) mtlb_CVodeMonitor(2, 0.0, NULL, NULL, NULL, cvmData);

    /* Extract index of current backward problem */
    
    idxB = (int)mxGetScalar(prhs[0]);

    /* Extract final time */

    tB0 = (double)mxGetScalar(prhs[1]);

    /* Extract final conditions */

    yB0 = mxGetPr(prhs[2]);
    NB = mxGetM(prhs[2]);

    /* Load final conditions */

    PutData(yB, yB0, NB);

    /* Extract options structure */
    
    options = prhs[3];

    break;

  }

  /* Process the options structure */

  get_IntgrOptions(options, bckPb, FALSE,
                   &lmmB, &iterB, &maxordB, &sldB, &mxstepsB,
                   &itolB, &reltolB, &SabstolB, &VabstolB,
                   &hinB, &hmaxB, &hminB, &tstopB, &rhs_s);

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

    CVodeCreateB(cvode_mem, lmmB, iterB, &idxB);
    CVodeSetUserDataB(cvode_mem, idxB, cvmData);
    CVodeInitB(cvode_mem, idxB, mtlb_CVodeRhsB, tB0, yB);

    /* Return idxB */

    plhs[0] = mxCreateScalarDouble((double)idxB);

    indexB = idxB;

    NbckPb++;

    break;

  case 1:

    CVodeReInitB(cvode_mem, idxB, tB0, yB);

    break;

  }

  /*
   * ----------------------------------------
   * Set tolerances
   * ----------------------------------------
   */

  switch (itolB) {
  case CV_SS:
    CVodeSStolerancesB(cvode_mem, idxB, reltolB, SabstolB);
    break;
  case CV_SV:
    NV_abstolB = N_VClone(yB);
    PutData(NV_abstolB, VabstolB, NB);
    CVodeSVtolerancesB(cvode_mem, idxB, reltolB, NV_abstolB);
    N_VDestroy(NV_abstolB);
    break;
  }

  /*
   * --------------------------------
   * Set various optional inputs
   * --------------------------------
   */

  /* set maxorder (default is consistent with LMM) */
  CVodeSetMaxOrdB(cvode_mem, idxB, maxordB);

  /* set initial step size (the default value of 0.0 is ignored by CVODES) */
  CVodeSetInitStepB(cvode_mem, idxB, hinB);

  /* set max step (default is infinity) */
  CVodeSetMaxStepB(cvode_mem, idxB, hmaxB);

  /* set min step (default is 0) */
  CVodeSetMinStepB(cvode_mem, idxB, hminB);
 
  /* set number of max steps */
  CVodeSetMaxNumStepsB(cvode_mem, idxB, mxstepsB);

  /*
   * ----------------------------------------
   * Need a linear solver?
   * ----------------------------------------
   */

  if (iterB == CV_NEWTON) {

    get_LinSolvOptions(prhs[3], bckPb, FALSE,
                       &mupperB, &mlowerB,
                       &mudqB, &mldqB, &dqrelyB,
                       &ptypeB, &gstypeB, &maxlB);

    switch(lsB) {

    case LS_DENSE:

      CVDenseB(cvode_mem, idxB, NB);
      if (!mxIsEmpty(mx_JACfctB)) CVDlsSetDenseJacFnB(cvode_mem, idxB, mtlb_CVodeDenseJacB);

      break;

    case LS_DIAG:

      CVDiagB(cvode_mem, idxB);

      break;

    case LS_BAND:

      CVBandB(cvode_mem, idxB, NB, mupperB, mlowerB);
      if (!mxIsEmpty(mx_JACfctB)) CVDlsSetBandJacFnB(cvode_mem, idxB, mtlb_CVodeBandJacB);

      break;

    case LS_SPGMR:

      CVSpgmrB(cvode_mem, idxB, ptypeB, maxlB);
      CVSpilsSetGSTypeB(cvode_mem, idxB, gstypeB);

      break;

    case LS_SPBCG:

      CVSpbcgB(cvode_mem, idxB, ptypeB, maxlB);

      break;

    case LS_SPTFQMR:

      CVSptfqmrB(cvode_mem, idxB, ptypeB, maxlB);

      break;

    }

    /* Jacobian * vector and preconditioner for SPILS linear solvers */

    if ( (lsB==LS_SPGMR) || (lsB==LS_SPBCG) || (lsB==LS_SPTFQMR) ) {

      if (!mxIsEmpty(mx_JACfctB)) CVSpilsSetJacTimesVecFnB(cvode_mem, idxB, mtlb_CVodeSpilsJacB);

      switch (pmB) {

      case PM_NONE:

        if (!mxIsEmpty(mx_PSOLfctB)) {
          if (!mxIsEmpty(mx_PSETfctB)) CVSpilsSetPreconditionerB(cvode_mem, idxB, mtlb_CVodeSpilsPsetB, mtlb_CVodeSpilsPsolB);
          else                         CVSpilsSetPreconditionerB(cvode_mem, idxB, NULL, mtlb_CVodeSpilsPsolB);
        }

        break;

      case PM_BANDPRE:

        CVBandPrecInitB(cvode_mem, idxB, NB, mupperB, mlowerB);

        break;

      case PM_BBDPRE:

        if (!mxIsEmpty(mx_GCOMfctB)) CVBBDPrecInitB(cvode_mem, idxB, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_CVodeBBDglocB, mtlb_CVodeBBDgcomB);
        else                         CVBBDPrecInitB(cvode_mem, idxB, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_CVodeBBDglocB, NULL);

        break;

      }

    }

  } else {

    lsB = LS_NONE;

  }

  /* Do we monitor? */

  if (monB) mtlb_CVodeMonitorB(0, idxB, tB0, NULL, NULL, cvmData);

  return(0);
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


  bckPb = cvmData->bckPb;

  /* 
   * ------------------------------------
   * Process inputs based on action
   * ------------------------------------
   */

  switch (action) {

  case 0:     /* BACKWARD QUADRATURE INITIALIZATION */

    /* Extract index of current backward problem */
    
    idxB = (int)mxGetScalar(prhs[0]);

    /* Extract user-provided quadrature RHS function */

    mxDestroyArray(mx_QUADfctB);
    mx_QUADfctB = mxDuplicateArray(prhs[1]);

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

    /* Extract index of current backward problem */
    
    idxB = (int)mxGetScalar(prhs[0]);

    /* Extract quadrature final conditions */

    yQB0 = mxGetPr(prhs[1]);
    NqB = mxGetM(prhs[1]);

    /* Load quadrature final conditions */

    PutData(yQB, yQB0, NqB);

    /* Extract quadrature options structure */

    options = prhs[2];

    break;

  }

  /* Process the options structure */

  get_QuadOptions(options, bckPb, FALSE,
                  NqB, &rhs_s,
                  &errconQB, 
                  &itolQB, &reltolQB, &SabstolQB, &VabstolQB);

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
    CVodeQuadInitB(cvode_mem, idxB, mtlb_CVodeQUADfctB, yQB);
    break;
  case 1:
    CVodeQuadReInitB(cvode_mem, idxB, yQB);
    break;
  }

  /*
   * ----------------------------------------
   * Set tolerances for quadrature variables
   * ----------------------------------------
   */
  
  CVodeSetQuadErrConB(cvode_mem, idxB, errconQB);

  if (errconQB) {

    switch (itolQB) {
    case CV_SS:
      CVodeQuadSStolerancesB(cvode_mem, idxB, reltolQB, SabstolQB);
      break;
    case CV_SV:
      NV_abstolQB = N_VClone(yQB);
      PutData(NV_abstolQB, VabstolQB, NqB);
      CVodeQuadSVtolerancesB(cvode_mem, idxB, reltolQB, NV_abstolQB);
      N_VDestroy(NV_abstolQB);
      break;
    }
    
  }

  quadrB = TRUE;

  return(0);
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

  int itask, is, Ntout, itout, status, s_idx;
  double *tout, tret, h;
  double *tdata, *ydata, *yQdata, *ySdata;
  long int nst;

  fwdPb = cvmData->fwdPb;

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

  if (nlhs_bad < 0) cvmErrHandler(-999, "CVODES", "CVode",
                                  "Too few output arguments.", NULL);
  if (nlhs_bad > 0) cvmErrHandler(-999, "CVODES", "CVode",
                                  "Too many output arguments.", NULL);

  /*
   * ----------------------------------------------------------------
   * Extract input arguments
   * ----------------------------------------------------------------
   */

  /* Extract tout */

  Ntout = mxGetM(prhs[0]) * mxGetN(prhs[0]);
  tout = mxGetPr(prhs[0]);

  /* Check if tout values are legal */

  CVodeGetCurrentTime(cvode_mem, &tret);
  CVodeGetNumSteps(cvode_mem, &nst);

  if (nst == 0) {
    h = tout[0] - tret;
  } else {
    CVodeGetLastStep(cvode_mem, &h);
    if ( (tout[0] - tret + h)*h < 0.0 ) cvmErrHandler(-999, "CVODES", "CVode",
                                                      "Illegal value of tout.", NULL);
  }

  for (itout=1; itout<Ntout; itout++) 
    if ( (tout[itout] - tout[itout-1])*h < 0.0 ) cvmErrHandler(-999, "CVODES", "CVode",
                                                               "tout values are not monotonic.", NULL);

  /* If rootfinding or tstop are enabled, we do not allow multiple output times */

  if (rootSet && (Ntout>1)) cvmErrHandler(-999, "CVODES", "CVode",
                                          "More than one tout value prohibited with rootfinding enabled.", NULL);

  if (tstopSet && (Ntout>1)) cvmErrHandler(-999, "CVODES", "CVode",
                                           "More than one tout value prohibited with tstop enabled.", NULL);

  /* Extract itask */

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Normal")) itask = CV_NORMAL;
  else if(!strcmp(bufval,"OneStep")) itask = CV_ONE_STEP;
  else cvmErrHandler(-999, "CVODES", "CVode",
                     "Illegal value for itask.", NULL); 

  /* If itask==CV_ONE_STEP, we do not allow multiple output times and we do not monitor */

  if (itask==CV_ONE_STEP) {

    if (Ntout>1) cvmErrHandler(-999, "CVODES", "CVode",
                               "More than one tout value prohibited in ONE_STEP mode.", NULL); 

    if (mon) {
      cvmErrHandler(+999, "CVODES", "CVode",
                    "Monitoring disabled in itask=CV_ONE_STEP", NULL);
      mon = FALSE;
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
    s_idx = quadr ? 4 : 3;
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

      if (!asa) status = CVode(cvode_mem, tout[itout], y, &tret, itask);
      else      status = CVodeF(cvode_mem, tout[itout], y, &tret, itask, &Nc);

      tdata[itout] = tret;

      GetData(y, &ydata[itout*N], N);

      if (quadr) {
        CVodeGetQuad(cvode_mem, &tret, yQ);
        GetData(yQ, &yQdata[itout*Nq], Nq);
      }

      if (fsa) {
        CVodeGetSens(cvode_mem, &tret, yS);
        for (is=0; is<Ns; is++)
          GetData(yS[is], &ySdata[itout*Ns*N+is*N], N);
      }

    }

  } else {

    /* Monitoring. itask = CV_NORMAL */

    for (itout=0; itout<Ntout; itout++) {
      
      while(1) {
      
        if (!asa) status = CVode(cvode_mem, tout[itout], y, &tret, CV_ONE_STEP);
        else      status = CVodeF(cvode_mem, tout[itout], y, &tret, CV_ONE_STEP, &Nc);

        /* Call the monitoring function */

        if (quadr)
          CVodeGetQuad(cvode_mem, &tret, yQ);

        if (fsa)
          CVodeGetSens(cvode_mem, &tret, yS);

        mtlb_CVodeMonitor(1, tret, y, yQ, yS, cvmData);

        /* If a root was found or tstop was reached, break out of while loop */

        if (status == CV_TSTOP_RETURN || status == CV_ROOT_RETURN) break;

        /* If current tout was reached break out of while loop */

        CVodeGetCurrentStep(cvode_mem, &h);
        if ( (tret - tout[itout])*h >= 0.0 ) break;
          
      };
      
      /* On a tstop or root return, return solution at tret.
       * Otherwise (status=CV_SUCCESS), return solution at tout[itout]. */

      if (status == CV_TSTOP_RETURN || status == CV_ROOT_RETURN) {

        if (quadr)
          CVodeGetQuad(cvode_mem, &tret, yQ);

        if (fsa)
          CVodeGetSens(cvode_mem, &tret, yS);

      } else {

        tret = tout[itout];
        
        CVodeGetDky(cvode_mem, tret, 0, y);

        if (quadr)
          CVodeGetQuadDky(cvode_mem, tret, 0, yQ);

        if (fsa)
          CVodeGetSensDky(cvode_mem, tret, 0, yS);

      }

      tdata[itout] = tret;

      GetData(y, &ydata[itout*N], N);

      if (quadr)
        GetData(yQ, &yQdata[itout*Nq], Nq);

      if (fsa)
        for (is=0; is<Ns; is++)
          GetData(yS[is], &ySdata[itout*Ns*N+is*N], N);

    }

  }

  /* CVODE return flag */

  plhs[0] = mxCreateScalarDouble((double)status);

  return(0);
}


static int CVM_SolveB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData bckPb;

  int buflen;
  char *bufval;

  int nlhs_bad;

  int itaskB, NtoutB, status;
  double *toutB;

  double tret, h;

  booleantype any_quadrB, any_monB;

  /*
   * -------------------------------------------------------
   * Check whether quadratures and/or monitoring are enabled
   * for at least one backward problem
   * -------------------------------------------------------
   */

  any_quadrB = FALSE;
  any_monB   = FALSE;
  bckPb = cvmData->bckPb;
  while(bckPb != NULL) {
    if (quadrB) any_quadrB = TRUE;
    if (monB)   any_monB   = TRUE;
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

  if (nlhs_bad < 0) cvmErrHandler(-999, "CVODES", "CVodeB",
                                  "Too few output arguments.", NULL);
  if (nlhs_bad > 0) cvmErrHandler(-999, "CVODES", "CVodeB",
                                  "Too many output arguments.", NULL);

  /*
   * ----------------------------------------------------------------
   * Extract input arguments
   * ----------------------------------------------------------------
   */

  /* Extract tout */

  NtoutB = mxGetM(prhs[0]) * mxGetN(prhs[0]);
  toutB = mxGetPr(prhs[0]);

  /* Check if first tout value is in the right direction */

  CVodeGetLastStep(cvode_mem, &h);
  CVodeGetCurrentTime(cvode_mem, &tret);

  if ( (tret - toutB[0])*h < 0.0 ) cvmErrHandler(-999, "CVODES", "CVodeB",
                                                 "tout value in wrong direction.", NULL);
  
  /* Extract itaskB */

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Normal"))       itaskB = CV_NORMAL;
  else if(!strcmp(bufval,"OneStep")) itaskB = CV_ONE_STEP;
  else cvmErrHandler(-999, "CVODES", "CVodeB",
                     "Illegal value for itask.", NULL); 

  /* If itask == CV_ONE_STEP, then
   * - we do not allow multiple output times
   * - we disable monitoring 
   */

  if ( itaskB == CV_ONE_STEP ) {
    
    if (NtoutB > 1) cvmErrHandler(-999, "CVODES", "CVodeB",
                                  "More than one tout value prohibited in ONE_STEP mode.", NULL); 

    if (any_monB) {
      cvmErrHandler(+999, "CVODES", "CVodeB",
                    "Monitoring disabled in itask=CV_ONE_STEP", NULL);
      bckPb = cvmData->bckPb;
      while(bckPb != NULL) {
        monB = FALSE;
        bckPb = bckPb->next;
      }
      any_monB = FALSE;
    }

  }

  if (NbckPb == 1) cvmSolveB_one(plhs, NtoutB, toutB, itaskB);
  else             cvmSolveB_more(plhs, NtoutB, toutB, itaskB, any_quadrB, any_monB);

  return(0);
}



static void cvmSolveB_one(mxArray *plhs[], int NtoutB, double *toutB, int itaskB)
{
  cvmPbData bckPb;
  
  void *cvode_memB;

  double tretB, hB;
  double *tdata, *ydata, *yQdata;
  int itout, status;
  long int nstB;

  bckPb = cvmData->bckPb;

  cvode_memB = CVodeGetAdjCVodeBmem(cvode_mem, indexB);

  /* Check if tout values are legal */

  CVodeGetCurrentTime(cvode_memB, &tretB);
  CVodeGetNumSteps(cvode_memB, &nstB);

  if (nstB == 0) {
    hB = toutB[0] - tretB;
  } else {
    CVodeGetLastStep(cvode_memB, &hB);
    if ( (toutB[0] - tretB + hB)*hB < 0.0 ) cvmErrHandler(-999, "CVODES", "CVodeB",
                                                          "Illegal value of tout.", NULL);
  }

  for (itout=1; itout<NtoutB; itout++) 
    if ( (toutB[itout] - toutB[itout-1])*hB < 0.0 ) cvmErrHandler(-999, "CVODES", "CVodeB",
                                                                  "tout values are not monotonic.", NULL);

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
      
      status = CVodeB(cvode_mem, toutB[itout], itaskB);
      
      CVodeGetB(cvode_mem, indexB, &tretB, yB);
      
      tdata[itout] = tretB;
      
      GetData(yB, &ydata[itout*NB], NB);
      
      if (quadrB) {
        CVodeGetQuadB(cvode_mem, indexB, &tretB, yQB);
        GetData(yQB, &yQdata[itout*NqB], NqB);
      }
      
    }
    
    
  } else {

    /* Monitoring. itask = CV_NORMAL */

    for (itout=0; itout<NtoutB; itout++) {
      
      while(1) {
        
        status = CVodeB(cvode_mem, toutB[itout], CV_ONE_STEP);
        
        /* Call the monitoring function */

        CVodeGetB(cvode_mem, indexB, &tretB, yB);
        
        if (quadrB)
          CVodeGetQuadB(cvode_mem, indexB, &tretB, yQB);
        
        mtlb_CVodeMonitorB(1, indexB, tretB, yB, yQB, cvmData);
        
        /* If current tout was reached break out of while loop */

        CVodeGetCurrentStep(cvode_memB, &hB);
        if ( (tretB - toutB[itout])*hB >= 0.0 ) break;
          
      };


      tretB = toutB[itout];

      tdata[itout] = tretB;

      CVodeGetDky(cvode_memB, tretB, 0, yB);

      GetData(yB, &ydata[itout*NB], NB);
      
      if (quadrB) {
        CVodeGetQuadDky(cvode_memB, tretB, 0, yQB);
        GetData(yQB, &yQdata[itout*NqB], NqB);
      }

    }

  }

  /* CVODE return flag */

  plhs[0] = mxCreateScalarDouble((double)status);
    
  return;
}


static void cvmSolveB_more(mxArray *plhs[], int NtoutB, double *toutB, int itaskB,
                           booleantype any_quadrB, booleantype any_monB)
{
  cvmPbData bckPb;
  mxArray *cell;
  void *cvode_memB;

  double tretB, h, hB;
  double **tdata, **ydata, **yQdata;
  int itout, status;
  long int nstB;

  char err_msg[80];

  /* Check if tout values are legal */

  CVodeGetLastStep(cvode_mem, &h);

  bckPb = cvmData->bckPb;
  while (bckPb != NULL) {

    cvode_memB = CVodeGetAdjCVodeBmem(cvode_mem, indexB);

    CVodeGetCurrentTime(cvode_memB, &tretB);
    CVodeGetNumSteps(cvode_memB, &nstB);

    if (nstB > 0) {
      CVodeGetLastStep(cvode_memB, &hB);
      if ( (toutB[0] - tretB + hB)*hB < 0.0 ) {
        sprintf(err_msg, "CVodeB:: illegal value of tout (pb. %d).",indexB+1);
        cvmErrHandler(-999, "CVODES", "CVodeB", err_msg, NULL);
      }
    }
    
    for (itout=1; itout<NtoutB; itout++) 
      if ( (toutB[itout] - toutB[itout-1]) * h > 0.0 ) {
        sprintf(err_msg, "CVodeB:: tout values are not monotonic (pb. %d).", indexB+1);
        cvmErrHandler(-999, "CVODES", "CVodeB", err_msg, NULL);
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

      status = CVodeB(cvode_mem, toutB[itout], itaskB);
      
      bckPb = cvmData->bckPb;
      while (bckPb != NULL) {

        CVodeGetB(cvode_mem, indexB, &tretB, yB);

        tdata[indexB][itout] = tretB;

        GetData(yB, &ydata[indexB][itout*NB], NB);
      
        if (quadrB) {
          CVodeGetQuadB(cvode_mem, indexB, &tretB, yQB);
          GetData(yQB, &yQdata[indexB][itout*NqB], NqB);
        }

        bckPb = bckPb->next;

      }

    }


  } else {

    /* Monitoring for at least one backward problem. itask = CV_NORMAL */


  }

  /* CVODE return flag */

  plhs[0] = mxCreateScalarDouble((double)status);


  return;
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

  long int nst, nfe, nsetups, nni, ncfn, netf, nge;
  int qlast, qcur;
  double h0used, hlast, hcur, tcur;
  int *rootsfound;

  long int njeD, nfeD;
  long int nfeDI;
  long int njeB, nfeB;
  long int nli, npe, nps, ncfl, njeSG, nfeSG;

  long int nfQe, netfQ;

  long int nfSe, nfeS, netfS, nsetupsS;
  long int nniS, ncfnS;

  int i, flag;
  mxArray *mx_root, *mx_ls, *mx_quad, *mx_fsa;
  mxArray *mx_rootsfound;
  double *tmp;
  int nfields;
  

  if (cvmData == NULL) return(0);

  fwdPb = cvmData->fwdPb;

  flag = CVodeGetIntegratorStats(cvode_mem, &nst, &nfe, &nsetups, 
                                 &netf, &qlast, &qcur, &h0used, &hlast, &hcur, &tcur);

  flag = CVodeGetNonlinSolvStats(cvode_mem, &nni, &ncfn);

  nfields = sizeof(fnames_intgr)/sizeof(*fnames_intgr);
  plhs[0] = mxCreateStructMatrix(1, 1, nfields, fnames_intgr);
 
  mxSetField(plhs[0], 0, "nst",     mxCreateScalarDouble((double)nst));
  mxSetField(plhs[0], 0, "nfe",     mxCreateScalarDouble((double)nfe));
  mxSetField(plhs[0], 0, "nsetups", mxCreateScalarDouble((double)nsetups));
  mxSetField(plhs[0], 0, "netf",    mxCreateScalarDouble((double)netf));
  mxSetField(plhs[0], 0, "nni",     mxCreateScalarDouble((double)nni));
  mxSetField(plhs[0], 0, "ncfn",    mxCreateScalarDouble((double)ncfn));
  mxSetField(plhs[0], 0, "qlast",   mxCreateScalarDouble((double)qlast));
  mxSetField(plhs[0], 0, "qcur",    mxCreateScalarDouble((double)qcur));
  mxSetField(plhs[0], 0, "h0used",  mxCreateScalarDouble(h0used));
  mxSetField(plhs[0], 0, "hlast",   mxCreateScalarDouble(hlast));
  mxSetField(plhs[0], 0, "hcur",    mxCreateScalarDouble(hcur));
  mxSetField(plhs[0], 0, "tcur",    mxCreateScalarDouble(tcur));


  /* Root Finding Statistics */

  if (Ng > 0) {

    flag = CVodeGetNumGEvals(cvode_mem, &nge);

    nfields = sizeof(fnames_root)/sizeof(*fnames_root);
    mx_root = mxCreateStructMatrix(1, 1, nfields, fnames_root);

    mxSetField(mx_root, 0, "nge", mxCreateScalarDouble((double)nge));

    rootsfound = (int *) malloc(nge*sizeof(int));
    flag = CVodeGetRootInfo(cvode_mem, rootsfound);
    mx_rootsfound = mxCreateDoubleMatrix(Ng,1,mxREAL);
    tmp = mxGetPr(mx_rootsfound);
    for (i=0;i<Ng;i++)
      tmp[i] = (double)rootsfound[i];
    mxSetField(mx_root, 0, "roots", mx_rootsfound);

  } else {

    mx_root = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "RootInfo", mx_root);

  /* Quadrature Statistics */

  if (quadr) {

    flag = CVodeGetQuadStats(cvode_mem, &nfQe, &netfQ);

    nfields = sizeof(fnames_quad)/sizeof(*fnames_quad);
    mx_quad = mxCreateStructMatrix(1, 1, nfields, fnames_quad);

    mxSetField(mx_quad, 0, "nfQe",  mxCreateScalarDouble((double)nfQe));
    mxSetField(mx_quad, 0, "netfQ", mxCreateScalarDouble((double)netfQ));

  } else {

    mx_quad = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "QuadInfo", mx_quad);

  /* Linear Solver Statistics */

  switch(ls){

  case LS_NONE:

    mx_ls = mxCreateDoubleMatrix(0,0,mxREAL);
    break;

  case LS_DENSE:
    
    flag = CVDlsGetNumJacEvals(cvode_mem, &njeD);
    flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeD);
    
    nfields = sizeof(fnames_dense)/sizeof(*fnames_dense);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_dense);
    
    mxSetField(mx_ls, 0, "name", mxCreateString("Dense"));
    mxSetField(mx_ls, 0, "njeD", mxCreateScalarDouble((double)njeD));
    mxSetField(mx_ls, 0, "nfeD", mxCreateScalarDouble((double)nfeD));
    
    break;

  case LS_DIAG:

    flag = CVDiagGetNumRhsEvals(cvode_mem, &nfeDI);
      
    nfields = sizeof(fnames_diag)/sizeof(*fnames_diag);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_diag);
 
    mxSetField(mx_ls, 0, "name", mxCreateString("Diag"));
    mxSetField(mx_ls, 0, "nfeDI", mxCreateScalarDouble((double)nfeDI));
      
    break;

  case LS_BAND:
      
    flag = CVDlsGetNumJacEvals(cvode_mem, &njeB);
    flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeB);
      
    nfields = sizeof(fnames_band)/sizeof(*fnames_band);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_band);
 
    mxSetField(mx_ls, 0, "name", mxCreateString("Band"));
    mxSetField(mx_ls, 0, "njeB", mxCreateScalarDouble((double)njeB));
    mxSetField(mx_ls, 0, "nfeB", mxCreateScalarDouble((double)nfeB));
      
    break;

  case LS_SPGMR:
  case LS_SPBCG:
  case LS_SPTFQMR:

    flag = CVSpilsGetNumLinIters(cvode_mem, &nli);
    flag = CVSpilsGetNumPrecEvals(cvode_mem, &npe);
    flag = CVSpilsGetNumPrecSolves(cvode_mem, &nps);
    flag = CVSpilsGetNumConvFails(cvode_mem, &ncfl);
    flag = CVSpilsGetNumJtimesEvals(cvode_mem, &njeSG);
    flag = CVSpilsGetNumRhsEvals(cvode_mem, &nfeSG);
    
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
    mxSetField(mx_ls, 0, "njeSG", mxCreateScalarDouble((double)njeSG));
    mxSetField(mx_ls, 0, "nfeSG", mxCreateScalarDouble((double)nfeSG));
    
    break;
    
  }

  mxSetField(plhs[0], 0, "LSInfo", mx_ls);

  /* Forward Sensitivity Statistics */

  if (fsa) {

    flag = CVodeGetSensStats(cvode_mem, &nfSe, &nfeS, &netfS, &nsetupsS); 

    flag = CVodeGetSensNonlinSolvStats(cvode_mem, &nniS, &ncfnS);

    nfields = sizeof(fnames_sens)/sizeof(*fnames_sens);
    mx_fsa = mxCreateStructMatrix(1, 1, nfields, fnames_sens);
    
    mxSetField(mx_fsa, 0, "nfSe",     mxCreateScalarDouble((double)nfSe));
    mxSetField(mx_fsa, 0, "nfeS",     mxCreateScalarDouble((double)nfeS));
    mxSetField(mx_fsa, 0, "nsetupsS", mxCreateScalarDouble((double)nsetupsS));
    mxSetField(mx_fsa, 0, "netfS",    mxCreateScalarDouble((double)netfS));
    mxSetField(mx_fsa, 0, "nniS",     mxCreateScalarDouble((double)nniS));
    mxSetField(mx_fsa, 0, "ncfnS",    mxCreateScalarDouble((double)ncfnS));
    
  } else {

    mx_fsa = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "FSAInfo", mx_fsa);

  return(0);
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

  long int nst, nfe, nsetups, nni, ncfn, netf;
  int qlast, qcur;
  double h0used, hlast, hcur, tcur;

  long int njeD, nfeD;
  long int nfeDI;
  long int njeB, nfeB;
  long int nli, npe, nps, ncfl, njeSG, nfeSG;

  long int nfQe, netfQ;

  int flag;
  mxArray *mx_ls, *mx_quad;
  int nfields;


  /* Extract index of current backward problem */
    
  idxB = (int)mxGetScalar(prhs[0]);

  /* Find current backward problem */

  bckPb = cvmData->bckPb;
  while (bckPb != NULL) {
    if (indexB == idxB) break;
    bckPb = bckPb->next;
  }

  cvode_memB = CVodeGetAdjCVodeBmem(cvode_mem, indexB);

  flag = CVodeGetIntegratorStats(cvode_memB, &nst, &nfe, &nsetups, 
                                 &netf, &qlast, &qcur, &h0used, &hlast, &hcur, &tcur);

  flag = CVodeGetNonlinSolvStats(cvode_memB, &nni, &ncfn);

  nfields = sizeof(fnames_intgr)/sizeof(*fnames_intgr);
  plhs[0] = mxCreateStructMatrix(1, 1, nfields, fnames_intgr);

  mxSetField(plhs[0], 0, "nst",     mxCreateScalarDouble((double)nst));
  mxSetField(plhs[0], 0, "nfe",     mxCreateScalarDouble((double)nfe));
  mxSetField(plhs[0], 0, "nsetups", mxCreateScalarDouble((double)nsetups));
  mxSetField(plhs[0], 0, "netf",    mxCreateScalarDouble((double)netf));
  mxSetField(plhs[0], 0, "nni",     mxCreateScalarDouble((double)nni));
  mxSetField(plhs[0], 0, "ncfn",    mxCreateScalarDouble((double)ncfn));
  mxSetField(plhs[0], 0, "qlast",   mxCreateScalarDouble((double)qlast));
  mxSetField(plhs[0], 0, "qcur",    mxCreateScalarDouble((double)qcur));
  mxSetField(plhs[0], 0, "h0used",  mxCreateScalarDouble(h0used));
  mxSetField(plhs[0], 0, "hlast",   mxCreateScalarDouble(hlast));
  mxSetField(plhs[0], 0, "hcur",    mxCreateScalarDouble(hcur));
  mxSetField(plhs[0], 0, "tcur",    mxCreateScalarDouble(tcur));


  /* Quadrature Statistics */

  if (quadrB) {

    flag = CVodeGetQuadStats(cvode_memB, &nfQe, &netfQ);

    nfields = sizeof(fnames_quad)/sizeof(*fnames_quad);
    mx_quad = mxCreateStructMatrix(1, 1, nfields, fnames_quad);

    mxSetField(mx_quad, 0, "nfQe",  mxCreateScalarDouble((double)nfQe));
    mxSetField(mx_quad, 0, "netfQ", mxCreateScalarDouble((double)netfQ));

  } else {

    mx_quad = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "QuadInfo", mx_quad);

  /* Linear Solver Statistics */

  switch(lsB){

  case LS_NONE:

    mx_ls = mxCreateDoubleMatrix(0,0,mxREAL);
    break;

  case LS_DENSE:
    
    flag = CVDlsGetNumJacEvals(cvode_memB, &njeD);
    flag = CVDlsGetNumRhsEvals(cvode_memB, &nfeD);
    
    nfields = sizeof(fnames_dense)/sizeof(*fnames_dense);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_dense);

    mxSetField(mx_ls, 0, "name", mxCreateString("Dense"));
    mxSetField(mx_ls, 0, "njeD", mxCreateScalarDouble((double)njeD));
    mxSetField(mx_ls, 0, "nfeD", mxCreateScalarDouble((double)nfeD));
    
    break;

  case LS_DIAG:

    flag = CVDiagGetNumRhsEvals(cvode_memB, &nfeDI);
      
    nfields = sizeof(fnames_diag)/sizeof(*fnames_diag);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_diag);
 
    mxSetField(mx_ls, 0, "name", mxCreateString("Diag"));
    mxSetField(mx_ls, 0, "nfeDI", mxCreateScalarDouble((double)nfeDI));
      
    break;

  case LS_BAND:
      
    flag = CVDlsGetNumJacEvals(cvode_memB, &njeB);
    flag = CVDlsGetNumRhsEvals(cvode_memB, &nfeB);
      
    nfields = sizeof(fnames_band)/sizeof(*fnames_band);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_band);
 
    mxSetField(mx_ls, 0, "name", mxCreateString("Band"));
    mxSetField(mx_ls, 0, "njeB", mxCreateScalarDouble((double)njeB));
    mxSetField(mx_ls, 0, "nfeB", mxCreateScalarDouble((double)nfeB));
      
    break;

  case LS_SPGMR:
  case LS_SPBCG:
  case LS_SPTFQMR:

    flag = CVSpilsGetNumLinIters(cvode_memB, &nli);
    flag = CVSpilsGetNumPrecEvals(cvode_memB, &npe);
    flag = CVSpilsGetNumPrecSolves(cvode_memB, &nps);
    flag = CVSpilsGetNumConvFails(cvode_memB, &ncfl);
    flag = CVSpilsGetNumJtimesEvals(cvode_memB, &njeSG);
    flag = CVSpilsGetNumRhsEvals(cvode_memB, &nfeSG);
    
    nfields = sizeof(fnames_spils)/sizeof(*fnames_spils);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_spils);
 
    if (lsB == LS_SPGMR)
      mxSetField(mx_ls, 0, "name",  mxCreateString("GMRES"));
    else if (lsB == LS_SPBCG)
      mxSetField(mx_ls, 0, "name",  mxCreateString("BiCGStab"));
    else
      mxSetField(mx_ls, 0, "name",  mxCreateString("TFQMR"));

    mxSetField(mx_ls, 0, "nli",   mxCreateScalarDouble((double)nli));
    mxSetField(mx_ls, 0, "npe",   mxCreateScalarDouble((double)npe));
    mxSetField(mx_ls, 0, "nps",   mxCreateScalarDouble((double)nps));
    mxSetField(mx_ls, 0, "ncfl",  mxCreateScalarDouble((double)ncfl));
    mxSetField(mx_ls, 0, "njeSG", mxCreateScalarDouble((double)njeSG));
    mxSetField(mx_ls, 0, "nfeSG", mxCreateScalarDouble((double)nfeSG));
    
    break;
  }

  mxSetField(plhs[0], 0, "LSInfo", mx_ls);

  return(0);
}


static int CVM_Set(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}

static int CVM_Get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData fwdPb;

  double t;
  N_Vector yd, ewt;
  double *this, *next;
  int key, k, which, i, nfields;
  mxArray *mx_y, *mx_yd;

  CVadjCheckPointRec *ckpnt;
  const char *fnames_ckpnt[]={
    "t0",
    "t1",
    "nstep",
    "order",
    "step"
  };


  fwdPb = cvmData->fwdPb;

  key = (int) (*mxGetPr(prhs[0]));

  switch (key) {
  case 1:    /* DerivSolution */

    t = *mxGetPr(prhs[1]);
    k = (int) (*mxGetPr(prhs[2]));

    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    CVodeGetDky(cvode_mem, t, k, y);
    GetData(y, mxGetPr(plhs[0]), N);

    break;

  case 2:    /* ErrorWeights */

    ewt = N_VClone(y);

    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    CVodeGetErrWeights(cvode_mem, ewt);
    GetData(ewt, mxGetPr(plhs[0]), N);

    N_VDestroy(ewt);

    break;

  case 3:    /* not used */

    break;

  case 4:    /* CheckPointsInfo */

    ckpnt = (CVadjCheckPointRec *) malloc ( (Nc+1)*sizeof(CVadjCheckPointRec));
    CVodeGetAdjCheckPointsInfo(cvode_mem, ckpnt);
    nfields = sizeof(fnames_ckpnt)/sizeof(*fnames_ckpnt);
    plhs[0] = mxCreateStructMatrix(Nc+1, 1, nfields, fnames_ckpnt);
    for (i=0; i<=Nc; i++) {
      this = (double *)(ckpnt[Nc-i].my_addr);
      next = (double *)(ckpnt[Nc-i].next_addr);
      mxSetField(plhs[0], i, "t0",    mxCreateScalarDouble((double)(ckpnt[Nc-i].t0)));
      mxSetField(plhs[0], i, "t1",    mxCreateScalarDouble((double)(ckpnt[Nc-i].t1)));
      mxSetField(plhs[0], i, "nstep", mxCreateScalarDouble((double)(ckpnt[Nc-i].nstep)));
      mxSetField(plhs[0], i, "order", mxCreateScalarDouble((double)(ckpnt[Nc-i].order)));
      mxSetField(plhs[0], i, "step",  mxCreateScalarDouble((double)(ckpnt[Nc-i].step)));
    }

    break;

  case 5:    /* CurrentCheckPoint */

    CVodeGetAdjCurrentCheckPoint(cvode_mem, (void **)(&this));

    plhs[0] = mxCreateScalarDouble(*this);

    break;

    break;

  }

  return(0);
}

static int CVM_Free(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cvmPbData fwdPb, bckPb;

  if (cvmData == NULL) return(0);

  fwdPb = cvmData->fwdPb;
  if (mon) mtlb_CVodeMonitor(2, 0.0, NULL, NULL, NULL, cvmData);

  bckPb = cvmData->bckPb;
  while (bckPb != NULL) {
    if (monB) mtlb_CVodeMonitorB(2, indexB, 0.0, NULL, NULL, cvmData);   
    bckPb = bckPb->next;
  }

  return(0);
}


