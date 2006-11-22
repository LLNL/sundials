/*
 * -----------------------------------------------------------------
 * $Revision: 1.10 $
 * $Date: 2006-11-22 00:12:51 $
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
 * Definitions for global variables (declared in cvm.h)
 * ---------------------------------------------------------------------------------
 */

booleantype cvm_quad;      /* Forward quadratures? */
booleantype cvm_quadB;     /* Backward quadratures? */
booleantype cvm_asa;       /* Adjoint sensitivity? */
booleantype cvm_fsa;       /* Forward sensitivity? */
booleantype cvm_mon;       /* Forward monitoring? */ 
booleantype cvm_monB;      /* Backward monitoring? */ 

cvm_CVODESdata cvm_Cdata = NULL;  /* CVODES data */
cvm_MATLABdata cvm_Mdata = NULL;  /* MATLAB data */

/*
 * ---------------------------------------------------------------------------------
 * Static function prototypes
 * ---------------------------------------------------------------------------------
 */

static void CVM_init();
static void CVM_makePersistent();
static void CVM_final();

static int CVM_Initialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_SensInitialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_SensToggleOff(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_AdjInitialization(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_InitializationB(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_Solve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_SolveB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
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
     2 - initialize forward sensitivity calculations
     3 - initialize adjoint sensitivity calculations
     4 - initialize backward solver

    11 - reinitialize CVODES solver
    12 - reinitialize forward sensitivity calculations
    13 - toggle FSA off
    14 - reinitialize backward solver

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

  switch(mode) {
  case 1:
    if (cvm_Cdata != NULL) {
      /* a previous pb ws initialized, we must clear  memory */
      CVM_Free(nlhs, plhs, nrhs-1, &prhs[1]);
      CVM_final();
    }
    CVM_init();
    CVM_Initialization(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 2:
    CVM_SensInitialization(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 3:
    CVM_AdjInitialization(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 4:
    CVM_InitializationB(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 11:
    CVM_Initialization(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 12:
    CVM_SensInitialization(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 13:
    CVM_SensToggleOff(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 14:
    CVM_InitializationB(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 20:
    CVM_Solve(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 21:
    CVM_SolveB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
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
  case 40:
    CVM_Free(nlhs, plhs, nrhs-1, &prhs[1]);
    CVM_final();
    return;
  }

  /* do not call CVM_makePersistent after free */
  if (mode != 40) CVM_makePersistent();
  mexLock();

  return;
}

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */

#define cvode_mem   (cvm_Cdata->cvode_mem)
#define bp_data     (cvm_Cdata->bp_data) 
#define y           (cvm_Cdata->y) 
#define yQ          (cvm_Cdata->yQ) 
#define yS          (cvm_Cdata->yS) 
#define N           (cvm_Cdata->N) 
#define Nq          (cvm_Cdata->Nq) 
#define Ng          (cvm_Cdata->Ng) 
#define Ns          (cvm_Cdata->Ns) 
#define Nd          (cvm_Cdata->Nd) 
#define Nc          (cvm_Cdata->Nc) 
#define ls          (cvm_Cdata->ls) 
#define pm          (cvm_Cdata->pm) 
#define ism         (cvm_Cdata->ism) 
#define cvadj_mem   (cvm_Cdata->cvadj_mem) 
#define interp      (cvm_Cdata->interp) 
#define yB          (cvm_Cdata->yB) 
#define yQB         (cvm_Cdata->yQB) 
#define NB          (cvm_Cdata->NB) 
#define NqB         (cvm_Cdata->NqB) 
#define lsB         (cvm_Cdata->lsB) 
#define pmB         (cvm_Cdata->pmB) 
#define errmsg      (cvm_Cdata->errmsg)

#define mx_data     (cvm_Mdata->mx_data)

#define mx_RHSfct   (cvm_Mdata->mx_RHSfct)
#define mx_QUADfct  (cvm_Mdata->mx_QUADfct)
#define mx_JACfct   (cvm_Mdata->mx_JACfct)
#define mx_PSETfct  (cvm_Mdata->mx_PSETfct)
#define mx_PSOLfct  (cvm_Mdata->mx_PSOLfct)
#define mx_GLOCfct  (cvm_Mdata->mx_GLOCfct)
#define mx_GCOMfct  (cvm_Mdata->mx_GCOMfct)
#define mx_Gfct     (cvm_Mdata->mx_Gfct)
#define mx_SRHSfct  (cvm_Mdata->mx_SRHSfct)

#define mx_RHSfctB  (cvm_Mdata->mx_RHSfctB)
#define mx_QUADfctB (cvm_Mdata->mx_QUADfctB)
#define mx_JACfctB  (cvm_Mdata->mx_JACfctB)
#define mx_PSETfctB (cvm_Mdata->mx_PSETfctB)
#define mx_PSOLfctB (cvm_Mdata->mx_PSOLfctB)
#define mx_GLOCfctB (cvm_Mdata->mx_GLOCfctB)
#define mx_GCOMfctB (cvm_Mdata->mx_GCOMfctB)

#define mx_MONfct   (cvm_Mdata->mx_MONfct)
#define mx_MONdata  (cvm_Mdata->mx_MONdata)

#define mx_MONfctB  (cvm_Mdata->mx_MONfctB)
#define mx_MONdataB (cvm_Mdata->mx_MONdataB)

/*
 * ---------------------------------------------------------------------------------
 * Exported procedures
 * ---------------------------------------------------------------------------------
 */

/* CVM_Initialization
 *
 * action = 0   -> CVodeCreate + CVodeMalloc
 * action = 1   -> CVodeReInit
 *
 * prhs contains:
 *   fct
 *   t0
 *   y0
 *   options
 *   data
 *
 * plhs contains:
 *   status
 *
 */

static int CVM_Initialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double t0, *y0;
  int status;

  int lmm, iter, maxord;

  long int mxsteps;

  int itol;
  realtype reltol, Sabstol, *Vabstol;
  N_Vector NV_abstol;
  void *abstol;

  double *yQ0;

  booleantype errconQ;
  int itolQ;
  realtype reltolQ, SabstolQ, *VabstolQ;
  N_Vector NV_abstolQ;
  void *abstolQ;

  double hin, hmax, hmin;

  double tstop;
  booleantype tstopSet;

  booleantype sld;

  int mupper, mlower;
  int ptype, gstype, maxl;
  int mudq, mldq;
  double dqrely;

  /* ------------------------------------
   * Initialize appropriate vector module
   * ------------------------------------
   */

  if (action == 0) InitVectors();

  /* 
   * -----------------------------
   * Extract stuff from arguments:
   * - RHS function
   * - initial time
   * - initial conditions
   * - integration options
   * - user data
   * -----------------------------
   */

  /* Matlab user-provided function */

  mxDestroyArray(mx_RHSfct);
  mx_RHSfct = mxDuplicateArray(prhs[0]);
  
  /* Initial time */

  t0 = (double)mxGetScalar(prhs[1]);

  /* Initial conditions */

  y0 = mxGetPr(prhs[2]);
  N = mxGetM(prhs[2]);

  /* Integrator Options -- optional argument */

  status = get_IntgrOptions(prhs[3], TRUE,
                            &lmm, &iter, &maxord, &sld, &mxsteps,
                            &itol, &reltol, &Sabstol, &Vabstol,
                            &hin, &hmax, &hmin, &tstop, &tstopSet);

  /* User data -- optional argument */

  mxDestroyArray(mx_data);
  mx_data = mxDuplicateArray(prhs[4]);

  /* 
   * ---------------------------------------
   * Set initial conditions and tolerances
   * ---------------------------------------
   */
  
  if (action == 0)  y = NewVector(N);

  PutData(y, y0, N);

  switch (itol) {
  case CV_SS:
    abstol = (void *) &Sabstol;
    break;
  case CV_SV:
    NV_abstol = N_VClone(y);
    PutData(NV_abstol, Vabstol, N);
    abstol = (void *) NV_abstol;
    break;
  }

  /* 
   * ----------------------------------------
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
    /* attach error handler function */
    status = CVodeSetErrHandlerFn(cvode_mem, mtlb_CVodeErrHandler, NULL);
    /* Call CVodeMalloc */
    status = CVodeMalloc(cvode_mem, mtlb_CVodeRhs, t0, y, itol, reltol, abstol);
    /* Redirect output */
    status = CVodeSetErrFile(cvode_mem, stdout);

    break;

  case 1:

    /* Reinitialize solver */
    status = CVodeReInit(cvode_mem, mtlb_CVodeRhs, t0, y, itol, reltol, abstol);

    break;

  }

  /* free NV_abstol if allocated */
  if (itol == CV_SV)  N_VDestroy(NV_abstol);

  /*
   * --------------------------------
   * Set various optional inputs
   * --------------------------------
   */

  /* set maxorder (default is consistent with LMM) */
  status = CVodeSetMaxOrd(cvode_mem, maxord);

  /* set initial step size (the default value of 0.0 is ignored by CVODES) */
  status = CVodeSetInitStep(cvode_mem, hin);

  /* set max step (default is infinity) */
  status = CVodeSetMaxStep(cvode_mem, hmax);

  /* set min step (default is 0) */
  status = CVodeSetMinStep(cvode_mem, hmin);
 
  /* set number of max steps */
  status = CVodeSetMaxNumSteps(cvode_mem, mxsteps);

  /* set tstop? */
  if (tstopSet)
    status = CVodeSetStopTime(cvode_mem, tstop);
  
  /* set stability limit detection (default is FALSE) */
  status = CVodeSetStabLimDet(cvode_mem, sld);
 
  /* Rootfinding? */
  if ( !mxIsEmpty(mx_Gfct) && (Ng > 0) ) {
    status = CVodeRootInit(cvode_mem, Ng, mtlb_CVodeGfct, NULL);
  }

  /* Quadratures? */
  if ( cvm_quad ) {

    status = get_QuadOptions(prhs[3], TRUE,
                               &yQ0, &errconQ, 
                               &itolQ, &reltolQ, &SabstolQ, &VabstolQ);

    if(status) mexErrMsgTxt("CVode Initialization:: illegal quadrature input.");      

    if (action == 0) yQ = NewVector(Nq);

    PutData(yQ, yQ0, Nq);

    switch (action) {
    case 0:
      status = CVodeQuadMalloc(cvode_mem, mtlb_CVodeQUADfct, yQ);
      break;
    case 1:
      status = CVodeQuadReInit(cvode_mem, mtlb_CVodeQUADfct, yQ);
      break;
    }

    if (errconQ) {
    
      switch (itolQ) {
      case CV_SS:
        abstolQ = (void *) &SabstolQ;
        break;
      case CV_SV:
        NV_abstolQ = N_VClone(yQ);
        PutData(NV_abstolQ, VabstolQ, Nq);
        abstolQ = (void *) NV_abstolQ;
        break;
      }
      
      status = CVodeSetQuadErrCon(cvode_mem, errconQ, itolQ, reltolQ, abstolQ);
      
      if (itolQ == CV_SV) N_VDestroy(NV_abstolQ);

    }

  }

  /* Need a linear solver? */
  if (iter == CV_NEWTON) {

    status = get_LinSolvOptions(prhs[3], TRUE,
                                &mupper, &mlower,
                                &mudq, &mldq, &dqrely,
                                &ptype, &gstype, &maxl);

    switch (ls) {
    case LS_DENSE:
      status = CVDense(cvode_mem, N);
      if (!mxIsEmpty(mx_JACfct))
        status = CVDlsSetJacFn(cvode_mem, mtlb_CVodeDenseJac, NULL);
      break;
    case LS_DIAG:
      status = CVDiag(cvode_mem);
      break;
    case LS_BAND:
      status = CVBand(cvode_mem, N, mupper, mlower);
      if (!mxIsEmpty(mx_JACfct))
        status = CVDlsSetJacFn(cvode_mem, mtlb_CVodeBandJac, NULL);
      break;
    case LS_SPGMR:
      switch (pm) {
      case PM_NONE:
        status = CVSpgmr(cvode_mem, ptype, maxl);
        if (!mxIsEmpty(mx_PSOLfct)) {
          if (!mxIsEmpty(mx_PSETfct))
            status = CVSpilsSetPreconditioner(cvode_mem, mtlb_CVodeSpilsPset, mtlb_CVodeSpilsPsol, NULL);
          else
            status = CVSpilsSetPreconditioner(cvode_mem, NULL, mtlb_CVodeSpilsPsol, NULL);
        }
        break;
      case PM_BANDPRE:
        bp_data = CVBandPrecAlloc(cvode_mem, N, mupper, mlower);
        status = CVBPSpgmr(cvode_mem, ptype, maxl, bp_data);
        break;
      case PM_BBDPRE:
        if (!mxIsEmpty(mx_GCOMfct))
          bp_data = CVBBDPrecAlloc(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_CVodeBBDgloc, mtlb_CVodeBBDgcom);
        else
          bp_data = CVBBDPrecAlloc(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_CVodeBBDgloc, NULL);
        status = CVBBDSpgmr(cvode_mem, ptype, maxl, bp_data);
        break;
      }
      status = CVSpilsSetGSType(cvode_mem, gstype);
      if (!mxIsEmpty(mx_JACfct))
        status = CVSpilsSetJacTimesVecFn(cvode_mem, mtlb_CVodeSpilsJac, NULL);
      break;
    case LS_SPBCG:
      switch (pm) {
      case PM_NONE:
        status = CVSpbcg(cvode_mem, ptype, maxl);
        if (!mxIsEmpty(mx_PSOLfct)) {
          if (!mxIsEmpty(mx_PSETfct))
            status = CVSpilsSetPreconditioner(cvode_mem, mtlb_CVodeSpilsPset, mtlb_CVodeSpilsPsol, NULL);
          else
            status = CVSpilsSetPreconditioner(cvode_mem, NULL, mtlb_CVodeSpilsPsol, NULL);
        }
        break;
      case PM_BANDPRE:
        bp_data = CVBandPrecAlloc(cvode_mem, N, mupper, mlower);
        status = CVBPSpbcg(cvode_mem, ptype, maxl, bp_data);
        break;
      case PM_BBDPRE:
        if (!mxIsEmpty(mx_GCOMfct))
          bp_data = CVBBDPrecAlloc(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_CVodeBBDgloc, mtlb_CVodeBBDgcom);
        else
          bp_data = CVBBDPrecAlloc(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_CVodeBBDgloc, NULL);
        CVBBDSpbcg(cvode_mem, ptype, maxl, bp_data);
        break;
      }
      if (!mxIsEmpty(mx_JACfct))
        status = CVSpilsSetJacTimesVecFn(cvode_mem, mtlb_CVodeSpilsJac, NULL);
      break;
    case LS_SPTFQMR:
      switch (pm) {
      case PM_NONE:
        status = CVSptfqmr(cvode_mem, ptype, maxl);
        if (!mxIsEmpty(mx_PSOLfct)) {
          if (!mxIsEmpty(mx_PSETfct))
            status = CVSpilsSetPreconditioner(cvode_mem, mtlb_CVodeSpilsPset, mtlb_CVodeSpilsPsol, NULL);
          else
            status = CVSpilsSetPreconditioner(cvode_mem, NULL, mtlb_CVodeSpilsPsol, NULL);
        }
        break;
      case PM_BANDPRE:
        bp_data = CVBandPrecAlloc(cvode_mem, N, mupper, mlower);
        status = CVBPSptfqmr(cvode_mem, ptype, maxl, bp_data);
        break;
      case PM_BBDPRE:
        if (!mxIsEmpty(mx_GCOMfct))
          bp_data = CVBBDPrecAlloc(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_CVodeBBDgloc, mtlb_CVodeBBDgcom);
        else
          bp_data = CVBBDPrecAlloc(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_CVodeBBDgloc, NULL);
        CVBBDSptfqmr(cvode_mem, ptype, maxl, bp_data);
        break;
      }
      if (!mxIsEmpty(mx_JACfct))
        status = CVSpilsSetJacTimesVecFn(cvode_mem, mtlb_CVodeSpilsJac, NULL);
      break;
    }

  } else {

    ls = LS_NONE;

  }

  /* Do we monitor? */
  
  if (cvm_mon) {
    mtlb_CVodeMonitor(0, t0, NULL, NULL, NULL);
  }

  return(0);

}

/* CVM_SensInitialization
 *
 * action = 0 -> CVodeSensMalloc
 * action = 1 -> CVodeSensReInit
 *
 * prhs contains:
 *   Ns (should be 0, if action=1)
 *   sensi_meth
 *   yS0
 *   options
 *
 * plhs contains:
 *   status
 *
 */

static int CVM_SensInitialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *yS0;
  int buflen, status;
  char *bufval;

  mxArray *pfield;
  char *pfield_name;

  booleantype userSRHS, errconS;
  int itolS;
  realtype reltolS;
  realtype *SabstolS, *VabstolS;
  N_Vector *NV_abstolS;
  void *abstolS;

  int *plist, dqtype;
  double *p, *pbar, rho;

  int is;

  p = NULL;
  plist = NULL;
  pbar = NULL;

  /* Number of sensitivities */

  if (action==0) Ns = (int)mxGetScalar(prhs[0]);

  /* Sensitivity method */

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Simultaneous")) ism = CV_SIMULTANEOUS;
  if(!strcmp(bufval,"Staggered"))    ism = CV_STAGGERED;

  /* Sensitivity initial conditions */

  yS0 = mxGetPr(prhs[2]);

  /* Extract Options */

  status = get_FSAOptions(prhs[3], 
                          &pfield_name, &plist, &pbar,
                          &userSRHS, &dqtype, &rho,
                          &errconS, &itolS, &reltolS, &SabstolS, &VabstolS);

  if(status) mexErrMsgTxt("CVode FSA Initialization:: illegal forward sensitivity input.");   

  /* Prepare arguments for CVODES functions */

  if (mxIsEmpty(mx_SRHSfct)) {
    if (pfield_name == NULL)
      mexErrMsgTxt("CVode FSA Initialization:: pfield required but was not provided.");
    pfield = mxGetField(mx_data,0,pfield_name);
    if (pfield == NULL)
      mexErrMsgTxt("CVode FSA Initialization:: illegal pfield input.");
    p = mxGetPr(pfield);
    mxFree(pfield_name);
  }

  if (action == 0) yS = N_VCloneVectorArray(Ns, y);

  for (is=0;is<Ns;is++) {
    PutData(yS[is], &yS0[is*N], N);
  }

  switch (action) {
  case 0:
    status = CVodeSensMalloc(cvode_mem, Ns, ism, yS);
    break;
  case 1:
    status = CVodeSensReInit(cvode_mem, ism, yS);
    break;
  }

  switch (itolS) {
  case CV_SS:
    abstolS = (void *) SabstolS;
    break;
  case CV_SV:
    NV_abstolS = N_VCloneVectorArray(Ns, y);
    for (is=0;is<Ns;is++)
      PutData(NV_abstolS[is], &VabstolS[is*N], N);
    abstolS = (void *) NV_abstolS;
    break;
  case CV_EE:
    abstolS = NULL;
    break;
  }
  
  status = CVodeSetSensTolerances(cvode_mem, itolS, reltolS, abstolS);
  
  if (itolS == CV_SS)      free(SabstolS);
  else if (itolS == CV_SV) N_VDestroyVectorArray(NV_abstolS, Ns);
  
  status = CVodeSetSensParams(cvode_mem, p, pbar, plist);
  
  if (plist != NULL) free(plist);
  if (pbar != NULL)  free(pbar);
  
  status = CVodeSetSensDQMethod(cvode_mem, dqtype, rho);
  
  status = CVodeSetSensErrCon(cvode_mem, errconS);

  if (userSRHS) {
    status = CVodeSetSensRhsFn(cvode_mem, mtlb_CVodeSensRhs, NULL);
  }

  cvm_fsa = TRUE;

  return(0);
}

/*
 * CVM_SenstoggleOff
 *
 * deactivates FSA
 */

static int CVM_SensToggleOff(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int status;

  status = CVodeSensToggleOff(cvode_mem);
  
  cvm_fsa = FALSE;

  return(0);
}



/* CVM_AdjInitialization
 *
 * prhs contains:
 *
 * plhs contains:
 *   status
 */

static int CVM_AdjInitialization(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int buflen, status;
  char *bufval;

  /* Number of steps */

  Nd = (int)mxGetScalar(prhs[0]);

  /* Interpolation method */

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Hermite"))    interp = CV_HERMITE;
  if(!strcmp(bufval,"Polynomial")) interp = CV_POLYNOMIAL;

  cvadj_mem = CVadjMalloc(cvode_mem, Nd, interp);

  cvm_asa = TRUE;
  
  return(0);
}

/* CVM_InitializationB
 *
 * action = 0   -> CVodeCreateB + CVodeMallocB
 * action = 1   -> CVodeReInitB
 *
 * prhs contains:
 *   fctB
 *   tF
 *   yB0
 *   options
 *   data
 *
 * plhs contains:
 *   status
 *
 */

static int CVM_InitializationB(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double tB0, *yB0;
  int status;

  int lmmB, iterB, maxordB;

  long int mxstepsB;

  int itolB;
  realtype reltolB, SabstolB, *VabstolB;
  N_Vector NV_abstolB;
  void *abstolB;

  double *yQB0;

  booleantype errconQB;
  int itolQB;
  realtype reltolQB, SabstolQB, *VabstolQB;
  N_Vector NV_abstolQB;
  void *abstolQB;

  double hinB, hmaxB, hminB;

  double tstopB;            /* ignored */
  booleantype tstopSetB;    /* ignored */

  booleantype sldB;         /* ignored */

  int mupperB, mlowerB;
  int ptypeB, gstypeB, maxlB;
  int mudqB, mldqB;
  double dqrelyB;

  /* 
   * -----------------------------
   * Finalize Forward monitoring
   * -----------------------------
   */


  if (cvm_mon) {
    mtlb_CVodeMonitor(2, 0.0, NULL, NULL, NULL);
    cvm_mon = FALSE;
  }

  /* 
   * -----------------------------
   * Extract stuff from arguments:
   * - Backward RHS function
   * - final time
   * - final conditions
   * - Backward integration options
   * -----------------------------
   */

  /* Matlab user-provided function */

  mxDestroyArray(mx_RHSfctB);
  mx_RHSfctB = mxDuplicateArray(prhs[0]);
  
  /* Final time */

  tB0 = (double)mxGetScalar(prhs[1]);

  /* Final conditions */

  yB0 = mxGetPr(prhs[2]);
  NB = mxGetM(prhs[2]);

  /* Integrator Options -- optional argument */

  status = get_IntgrOptions(prhs[3], FALSE,
                            &lmmB, &iterB, &maxordB, &sldB, &mxstepsB,
                            &itolB, &reltolB, &SabstolB, &VabstolB,
                            &hinB, &hmaxB, &hminB, &tstopB, &tstopSetB);

  /* 
   * ---------------------------------------
   * Set final conditions and tolerances
   * ---------------------------------------
   */

  if (action == 0) yB = NewVector(NB);

  PutData(yB, yB0, NB);

  switch (itolB) {
  case CV_SS:
    abstolB = (void *) &SabstolB;
    break;
  case CV_SV:
    NV_abstolB = N_VClone(yB);
    PutData(NV_abstolB, VabstolB, NB);
    abstolB = (void *) VabstolB;
    break;
  }

  /* 
   * -----------------------------------------------------------
   * If action = 0
   *    Create CVODES object for backward phase
   *    Allocate memory
   *    Redirect output
   * If action = 1
   *    Reinitialize solver
   * -----------------------------------------------------------
   */

  switch (action) {

  case 0:

    /* Create CVODES object */
    status = CVodeCreateB(cvadj_mem, lmmB, iterB);
    /* Call CVodeMallocB */
    status = CVodeMallocB(cvadj_mem, mtlb_CVodeRhsB, tB0, yB, itolB, reltolB, abstolB);
    /* Redirect output */
    status = CVodeSetErrFileB(cvadj_mem, stdout);

    break;

  case 1:

    /* Reinitialize solver */
    status = CVodeReInitB(cvadj_mem, mtlb_CVodeRhsB, tB0, yB, itolB, reltolB, abstolB);

    break;

  }

  if (itolB == CV_SV) {
    N_VDestroy(NV_abstolB);
  }

  /* set maxorder (default is consistent with LMM) */
  status = CVodeSetMaxOrdB(cvadj_mem, maxordB);

  /* set initial step size (the default value of 0.0 is ignored by CVODES) */
  status = CVodeSetInitStepB(cvadj_mem, hinB);

  /* set max step (default is infinity) */
  status = CVodeSetMaxStepB(cvadj_mem, hmaxB);

  /* set min step (default is 0) */
  status = CVodeSetMinStepB(cvadj_mem, hminB);
 
  /* set number of max steps */
  status = CVodeSetMaxNumStepsB(cvadj_mem, mxstepsB);

  /* Quadratures? */
  if ( cvm_quadB ) {
    
    status = get_QuadOptions(prhs[3], FALSE,
                             &yQB0, &errconQB, 
                             &itolQB, &reltolQB, &SabstolQB, &VabstolQB);

    if(status)
      mexErrMsgTxt("CVode Backward Initialization:: illegal quadrature input.");      

    if (action == 0) yQB = NewVector(NqB);

    PutData(yQB, yQB0, NqB);

    switch (action) {
    case 0:
      status = CVodeQuadMallocB(cvadj_mem, mtlb_CVodeQUADfctB, yQB);
      break;
    case 1:
      status = CVodeQuadReInitB(cvadj_mem, mtlb_CVodeQUADfctB, yQB);
      break;
    }

    switch (itolQB) {
    case CV_SS:
      abstolQB = (void *) &SabstolQB;
      break;
    case CV_SV:
      NV_abstolQB = N_VClone(yQB);
      PutData(NV_abstolQB, VabstolQB, NqB);
      abstolQB = (void *) NV_abstolQB;
      break;
    }

    status = CVodeSetQuadErrConB(cvadj_mem, errconQB, itolQB, reltolQB, abstolQB);

    if (itolQB == CV_SV) {
      N_VDestroy(NV_abstolQB);
    }

  }

  /* Need a linear solver? */

  if (iterB == CV_NEWTON) {

    status = get_LinSolvOptions(prhs[3], FALSE,
                                &mupperB, &mlowerB,
                                &mudqB, &mldqB, &dqrelyB,
                                &ptypeB, &gstypeB, &maxlB);

    switch(lsB) {
    case LS_DENSE:
      status = CVDenseB(cvadj_mem, NB);
      if (!mxIsEmpty(mx_JACfctB))
        status = CVDlsSetJacFnB(cvadj_mem, mtlb_CVodeDenseJacB, NULL);
      break;
    case LS_DIAG:
      status = CVDiagB(cvadj_mem);
      break;
    case LS_BAND:
      status = CVBandB(cvadj_mem, NB, mupperB, mlowerB);
      if (!mxIsEmpty(mx_JACfctB))
        status = CVDlsSetJacFnB(cvadj_mem, mtlb_CVodeBandJacB, NULL);
      break;
    case LS_SPGMR:
      switch (pmB) {
      case PM_NONE:
        status = CVSpgmrB(cvadj_mem, ptypeB, maxlB);
        if (!mxIsEmpty(mx_PSOLfctB)) {
          if (!mxIsEmpty(mx_PSETfctB))
            status = CVSpilsSetPreconditionerB(cvadj_mem, mtlb_CVodeSpilsPsetB, mtlb_CVodeSpilsPsolB, NULL);
          else
            status = CVSpilsSetPreconditionerB(cvadj_mem, NULL, mtlb_CVodeSpilsPsolB, NULL);
        }
        break;
      case PM_BANDPRE:
        status = CVBandPrecAllocB(cvadj_mem, NB, mupperB, mlowerB);
        status = CVBPSpgmrB(cvadj_mem, ptypeB, maxlB);
        break;
      case PM_BBDPRE:
        if (!mxIsEmpty(mx_GCOMfctB)) {
          status = CVBBDPrecAllocB(cvadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_CVodeBBDglocB, mtlb_CVodeBBDgcomB);
        } else {
          status = CVBBDPrecAllocB(cvadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_CVodeBBDglocB, NULL);
        }
        CVBBDSpgmrB(cvadj_mem, ptypeB, maxlB);
        break;
      }
      status = CVSpilsSetGSTypeB(cvadj_mem, gstypeB);
      if (!mxIsEmpty(mx_JACfctB))
        status = CVSpilsSetJacTimesVecFnB(cvadj_mem, mtlb_CVodeSpilsJacB, NULL);
      break;
    case LS_SPBCG:
      switch (pmB) {
      case PM_NONE:
        status = CVSpbcgB(cvadj_mem, ptypeB, maxlB);
        if (!mxIsEmpty(mx_PSOLfctB)) {
          if (!mxIsEmpty(mx_PSETfctB))
            status = CVSpilsSetPreconditionerB(cvadj_mem, mtlb_CVodeSpilsPsetB, mtlb_CVodeSpilsPsolB, NULL);
          else
            status = CVSpilsSetPreconditionerB(cvadj_mem, NULL, mtlb_CVodeSpilsPsolB, NULL);
        }
        break;
      case PM_BANDPRE:
        status = CVBandPrecAllocB(cvadj_mem, NB, mupperB, mlowerB);
        status = CVBPSpbcgB(cvadj_mem, ptypeB, maxlB);
        break;
      case PM_BBDPRE:
        if (!mxIsEmpty(mx_GCOMfctB)) {
          status = CVBBDPrecAllocB(cvadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_CVodeBBDglocB, mtlb_CVodeBBDgcomB);
        } else {
          status = CVBBDPrecAllocB(cvadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_CVodeBBDglocB, NULL);
        }
        CVBBDSpbcgB(cvadj_mem, ptypeB, maxlB);
        break;
      }
      if (!mxIsEmpty(mx_JACfctB))
        status = CVSpilsSetJacTimesVecFnB(cvadj_mem, mtlb_CVodeSpilsJacB, NULL);
      break;
    case LS_SPTFQMR:
      switch (pmB) {
      case PM_NONE:
        status = CVSptfqmrB(cvadj_mem, ptypeB, maxlB);
        if (!mxIsEmpty(mx_PSOLfctB)) {
          if (!mxIsEmpty(mx_PSETfctB))
            status = CVSpilsSetPreconditionerB(cvadj_mem, mtlb_CVodeSpilsPsetB, mtlb_CVodeSpilsPsolB, NULL);
          else
            status = CVSpilsSetPreconditionerB(cvadj_mem, NULL, mtlb_CVodeSpilsPsolB, NULL);
        }
        break;
      case PM_BANDPRE:
        status = CVBandPrecAllocB(cvadj_mem, NB, mupperB, mlowerB);
        status = CVBPSptfqmrB(cvadj_mem, ptypeB, maxlB);
        break;
      case PM_BBDPRE:
        if (!mxIsEmpty(mx_GCOMfctB)) {
          status = CVBBDPrecAllocB(cvadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_CVodeBBDglocB, mtlb_CVodeBBDgcomB);
        } else {
          status = CVBBDPrecAllocB(cvadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_CVodeBBDglocB, NULL);
        }
        CVBBDSptfqmrB(cvadj_mem, ptypeB, maxlB);
        break;
      }
      if (!mxIsEmpty(mx_JACfctB))
        status = CVSpilsSetJacTimesVecFnB(cvadj_mem, mtlb_CVodeSpilsJacB, NULL);
      break;
    }

  } else {

    lsB = LS_NONE;

  }

  /* Do we monitor? */
  
  if (cvm_monB) {
    mtlb_CVodeMonitorB(0, tB0, NULL, NULL);
  }

  return(0);
}


/*
 * CVM_Solve  - Main solution function
 *
 */

static int CVM_Solve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double tout, tret, *ydata;
  int buflen, status, itask, is;
  char *bufval;

  int itask1;
  booleantype iret;
  double h;

  /* Exract tout */

  tout = (double)mxGetScalar(prhs[0]);

  /* Extract itask */

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Normal")) itask = CV_NORMAL;
  if(!strcmp(bufval,"OneStep")) itask = CV_ONE_STEP;
  if(!strcmp(bufval,"NormalTstop")) itask = CV_NORMAL_TSTOP;
  if(!strcmp(bufval,"OneStepTstop")) itask = CV_ONE_STEP_TSTOP;

  /* Call CVode */

  if (!cvm_mon) {

    if (!cvm_asa)
      status = CVode(cvode_mem, tout, y, &tret, itask);
    else
      status = CVodeF(cvadj_mem, tout, y, &tret, itask, &Nc);

  } else {

    if      (itask == CV_NORMAL)         {iret = FALSE; itask1 = CV_ONE_STEP;}
    else if (itask == CV_ONE_STEP)       {iret = TRUE;  itask1 = CV_ONE_STEP;}
    else if (itask == CV_NORMAL_TSTOP)   {iret = FALSE; itask1 = CV_ONE_STEP_TSTOP;}
    else if (itask == CV_ONE_STEP_TSTOP) {iret = TRUE;  itask1 = CV_ONE_STEP_TSTOP;}
    
    while(1) {
      
      if (!cvm_asa)
        status = CVode(cvode_mem, tout, y, &tret, itask1);
      else
        status = CVodeF(cvadj_mem, tout, y, &tret, itask1, &Nc);

      /* break on CVode error */
      if (status < 0) break;   
      
      /* In NORMAL_MODE test if tout was reached */
      if (!iret) {
        CVodeGetCurrentStep(cvode_mem, &h);
        if ( (tret - tout)*h >= 0.0 ) {
          tret = tout;
          CVodeGetDky(cvode_mem, tout, 0, y);
          iret = TRUE;
        }
      }

      /* If root or tstop return, we'll need to break */
      if (status != CV_SUCCESS) iret = TRUE; 

      if (cvm_quad)
        status = CVodeGetQuad(cvode_mem, tret, yQ);

      if (cvm_fsa)
        status = CVodeGetSens(cvode_mem, tret, yS);

      /* Call the monitoring function */
      mtlb_CVodeMonitor(1, tret, y, yQ, yS);

      /* break if we need to */
      if(iret)  break;
      
    };

  }

  /* CVODE return flag */
  plhs[0] = mxCreateScalarDouble((double)status);

  /* Return time */
  plhs[1] = mxCreateScalarDouble(tret);

  /* Solution vector */
  plhs[2] = mxCreateDoubleMatrix(N,1,mxREAL);
  GetData(y, mxGetPr(plhs[2]), N);


  if (nlhs == 3)
    return(0);

  if (nlhs == 4) {

    if (cvm_quad) {
      plhs[3] = mxCreateDoubleMatrix(Nq,1,mxREAL);
      status = CVodeGetQuad(cvode_mem, tret, yQ);
      GetData(yQ, mxGetPr(plhs[3]), Nq);
      return(0);
    } else if (cvm_fsa) {
      plhs[3] = mxCreateDoubleMatrix(N,Ns,mxREAL);
      ydata = mxGetPr(plhs[3]);
      status = CVodeGetSens(cvode_mem, tret, yS);
      for (is=0; is<Ns; is++)
        GetData(yS[is], &ydata[is*N], N);
      return(0);
    } else {
      mexErrMsgTxt("CVode:: too many output arguments (4).");
      return(-1);
    }

  }

  if ( (!cvm_quad) || (!cvm_fsa) ) {
    mexErrMsgTxt("CVode:: too many output arguments (5).");
    return(-1);
  }

  /* Quadratures */
  plhs[3] = mxCreateDoubleMatrix(Nq,1,mxREAL);
  status = CVodeGetQuad(cvode_mem, tret, yQ);
  GetData(yQ, mxGetPr(plhs[3]), Nq);

  /* Sensitivities */
  plhs[4] = mxCreateDoubleMatrix(N,Ns,mxREAL);
  ydata = mxGetPr(plhs[4]);
  status = CVodeGetSens(cvode_mem, tret, yS);
  for (is=0; is<Ns; is++)
    GetData(yS[is], &ydata[is*N], N);

  return(0);
}


static int CVM_SolveB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double toutB, tretB;
  int itaskB;
  int buflen, status;
  char *bufval;

  void *cvode_memB;
  booleantype iretB;
  double hB;

  /* Exract toutB */

  toutB = (double)mxGetScalar(prhs[0]);

  /* Extract itaskB */

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Normal"))       itaskB = CV_NORMAL;
  else if(!strcmp(bufval,"OneStep")) itaskB = CV_ONE_STEP;
  else status = -1;

  /* Call CVodeB */

  if (!cvm_monB) {

    status = CVodeB(cvadj_mem, toutB, yB, &tretB, itaskB);

  } else {
    
    if      (itaskB == CV_NORMAL)   iretB = FALSE;
    else if (itaskB == CV_ONE_STEP) iretB = TRUE;

    while(1) {

      status = CVodeB(cvadj_mem, toutB, yB, &tretB, CV_ONE_STEP);

      /* break on CVode error */
      if (status < 0) break;   
      
      /* In NORMAL_MODE test if tout was reached */
      if (!iretB) {
        cvode_memB = CVadjGetCVodeBmem(cvadj_mem);
        CVodeGetCurrentStep(cvode_memB, &hB);
        if ( (tretB - toutB)*hB >= 0.0 ) {
          tretB = toutB;
          CVodeGetDky(cvode_memB, toutB, 0, yB);
          iretB = TRUE;
        }
      }

      if (cvm_quadB)
        status = CVodeGetQuadB(cvadj_mem, yQB);

      /* Call the monitoring function */
      mtlb_CVodeMonitorB(1, tretB, yB, yQB);

      /* break if we need to */
      if(iretB)  break;

    };

  }

  /* CVodeB return flag */
  plhs[0] = mxCreateScalarDouble((double)status);

  /* Return time */
  plhs[1] = mxCreateScalarDouble(tretB);

  /* Solution vector */
  plhs[2] = mxCreateDoubleMatrix(NB,1,mxREAL);
  GetData(yB, mxGetPr(plhs[2]), NB);

  if (nlhs == 3)
    return(0);

  if(!cvm_quadB) {
    mexErrMsgTxt("CVodeB:: too many output arguments.");
    return(-1);
  }

  plhs[3] = mxCreateDoubleMatrix(NqB,1,mxREAL);
  status = CVodeGetQuadB(cvadj_mem, yQB);
  GetData(yQB, mxGetPr(plhs[3]), NqB);

  return(0);
}


static int CVM_Stats(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
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

  if (cvm_Cdata == NULL) return(0);

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

  if (cvm_quad) {

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

  if (cvm_fsa) {

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

  cvode_memB = CVadjGetCVodeBmem(cvadj_mem);

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

  if (cvm_quadB) {

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
  double t;
  N_Vector yd, ewt;
  double *this, *next;
  int key, k, which, status, i, nfields;
  mxArray *mx_y, *mx_yd;

  CVadjCheckPointRec *ckpnt;
  const char *fnames_ckpnt[]={
    "t0",
    "t1",
    "nstep",
    "order",
    "step"
  };

  key = (int) (*mxGetPr(prhs[0]));

  switch (key) {
  case 1:    /* DerivSolution */

    t = *mxGetPr(prhs[1]);
    k = (int) (*mxGetPr(prhs[2]));

    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    status = CVodeGetDky(cvode_mem, t, k, y);
    GetData(y, mxGetPr(plhs[0]), N);

    break;

  case 2:    /* ErrorWeights */

    ewt = N_VClone(y);

    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    status = CVodeGetErrWeights(cvode_mem, ewt);
    GetData(ewt, mxGetPr(plhs[0]), N);

    N_VDestroy(ewt);

    break;

  case 3:    /* not used */

    break;

  case 4:    /* CheckPointsInfo */

    ckpnt = (CVadjCheckPointRec *) malloc ( (Nc+1)*sizeof(CVadjCheckPointRec));
    CVadjGetCheckPointsInfo(cvadj_mem, ckpnt);
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

    status = CVadjGetCurrentCheckPoint(cvadj_mem, (void **)(&this));

    plhs[0] = mxCreateScalarDouble(*this);

    break;

  case 6:    /* DataPointInfo */
    
    which = (int) (*mxGetPr(prhs[1]));

    if (interp == CV_HERMITE) {
    
      yd = N_VClone(y);

      status = CVadjGetDataPointHermite(cvadj_mem, which, &t, y, yd);

      plhs[0] = mxCreateCellMatrix(1, 3);

      mxSetCell(plhs[0],1,mxCreateScalarDouble(t));

      mx_y  = mxCreateDoubleMatrix(N,1,mxREAL);
      GetData(y, mxGetPr(mx_y), N);
      mxSetCell(plhs[0],2,mx_y);

      mx_yd = mxCreateDoubleMatrix(N,1,mxREAL);
      GetData(yd, mxGetPr(mx_yd), N);
      mxSetCell(plhs[0],3,mx_y);


      N_VDestroy(yd);
    }

    break;

  }

  return(0);
}

static int CVM_Free(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  if (cvm_Cdata == NULL) return(0);

  if (cvm_mon) {
    mtlb_CVodeMonitor(2, 0.0, NULL, NULL, NULL);
  }

  N_VDestroy(y);
  if (cvm_quad) N_VDestroy(yQ);
  if (pm == PM_BANDPRE) CVBandPrecFree(&bp_data);
  if (pm == PM_BBDPRE)  CVBBDPrecFree(&bp_data);

  if (cvm_fsa) N_VDestroyVectorArray(yS, Ns);

  CVodeFree(&cvode_mem);

  if (cvm_asa) {

    if (cvm_monB) {
      mtlb_CVodeMonitorB(2, 0.0, NULL, NULL);
    }

    N_VDestroy(yB);

    if (cvm_quadB) N_VDestroy(yQB);
   
    CVadjFree(&cvadj_mem);

  }

  return(0);
}


static void CVM_init()
{
  mxArray *empty;

  /* Allocate space for global CVODES and MATLAB data structures */
  
  cvm_Cdata = (cvm_CVODESdata) mxMalloc(sizeof(struct cvm_CVODESdataStruct));
  cvm_Mdata = (cvm_MATLABdata) mxMalloc(sizeof(struct cvm_MATLABdataStruct));

  /* Initialize global CVODES data */

  cvode_mem = NULL;
  cvadj_mem = NULL;
  bp_data   = NULL;

  y         = NULL;
  yQ        = NULL;
  yS        = NULL;
  yB        = NULL;
  yQB       = NULL;

  N   = 0;
  Nq  = 0;
  Ng  = 0;
  Ns  = 0;
  Nd  = 0;
  Nc  = 0;
  NB  = 0;
  NqB = 0;

  ism = CV_STAGGERED;

  interp = CV_POLYNOMIAL;

  ls  = LS_DENSE;
  lsB = LS_DENSE;
  pm  = PM_NONE;
  pmB = PM_NONE;

  errmsg = TRUE;

  /* Initialize global control variables */

  cvm_quad  = FALSE;
  cvm_quadB = FALSE;
  cvm_fsa   = FALSE;
  cvm_asa   = FALSE;
  cvm_mon   = FALSE;
  cvm_monB  = FALSE;

  /* Initialize global MATLAB data */

  empty = mxCreateDoubleMatrix(0,0,mxREAL);

  mx_data     = mxDuplicateArray(empty);

  mx_RHSfct   = mxDuplicateArray(empty);
  mx_Gfct     = mxDuplicateArray(empty);
  mx_QUADfct  = mxDuplicateArray(empty);
  mx_SRHSfct  = mxDuplicateArray(empty);
  mx_JACfct   = mxDuplicateArray(empty);
  mx_PSETfct  = mxDuplicateArray(empty);
  mx_PSOLfct  = mxDuplicateArray(empty);
  mx_GLOCfct  = mxDuplicateArray(empty);
  mx_GCOMfct  = mxDuplicateArray(empty);

  mx_RHSfctB  = mxDuplicateArray(empty);
  mx_QUADfctB = mxDuplicateArray(empty);
  mx_JACfctB  = mxDuplicateArray(empty);
  mx_PSETfctB = mxDuplicateArray(empty);
  mx_PSOLfctB = mxDuplicateArray(empty);
  mx_GLOCfctB = mxDuplicateArray(empty);
  mx_GCOMfctB = mxDuplicateArray(empty);

  mx_MONfct   = mxDuplicateArray(empty);
  mx_MONdata  = mxDuplicateArray(empty);

  mx_MONfctB  = mxDuplicateArray(empty);
  mx_MONdataB = mxDuplicateArray(empty);  

  mxDestroyArray(empty);

  return;
}

static void CVM_makePersistent()
{
  /* Make global memory persistent */

  mexMakeArrayPersistent(mx_data);

  mexMakeArrayPersistent(mx_RHSfct);
  mexMakeArrayPersistent(mx_Gfct);
  mexMakeArrayPersistent(mx_QUADfct);
  mexMakeArrayPersistent(mx_SRHSfct);
  mexMakeArrayPersistent(mx_JACfct);
  mexMakeArrayPersistent(mx_PSETfct);
  mexMakeArrayPersistent(mx_PSOLfct);
  mexMakeArrayPersistent(mx_GLOCfct);
  mexMakeArrayPersistent(mx_GCOMfct);

  mexMakeArrayPersistent(mx_RHSfctB);
  mexMakeArrayPersistent(mx_QUADfctB);
  mexMakeArrayPersistent(mx_JACfctB);
  mexMakeArrayPersistent(mx_PSETfctB);
  mexMakeArrayPersistent(mx_PSOLfctB);
  mexMakeArrayPersistent(mx_GLOCfctB);
  mexMakeArrayPersistent(mx_GCOMfctB);

  mexMakeArrayPersistent(mx_MONfct);
  mexMakeArrayPersistent(mx_MONdata);

  mexMakeArrayPersistent(mx_MONfctB);
  mexMakeArrayPersistent(mx_MONdataB);

  mexMakeMemoryPersistent(cvm_Cdata);
  mexMakeMemoryPersistent(cvm_Mdata);

  return;
}


static void CVM_final()
{  
  
  if (cvm_Cdata == NULL) return;

  mxDestroyArray(mx_data);
  
  mxDestroyArray(mx_RHSfct);
  mxDestroyArray(mx_Gfct);
  mxDestroyArray(mx_QUADfct);
  mxDestroyArray(mx_SRHSfct);
  mxDestroyArray(mx_JACfct);
  mxDestroyArray(mx_PSETfct);
  mxDestroyArray(mx_PSOLfct);
  mxDestroyArray(mx_GLOCfct);
  mxDestroyArray(mx_GCOMfct);
  
  mxDestroyArray(mx_RHSfctB);
  mxDestroyArray(mx_QUADfctB);
  mxDestroyArray(mx_JACfctB);
  mxDestroyArray(mx_PSETfctB);
  mxDestroyArray(mx_PSOLfctB);
  mxDestroyArray(mx_GLOCfctB);
  mxDestroyArray(mx_GCOMfctB);

  mxDestroyArray(mx_MONfct);
  mxDestroyArray(mx_MONdata);
  
  mxDestroyArray(mx_MONfctB);
  mxDestroyArray(mx_MONdataB);

  mxFree(cvm_Cdata);
  cvm_Cdata = NULL;
  mxFree(cvm_Mdata);
  cvm_Mdata = NULL;

  return;
}
