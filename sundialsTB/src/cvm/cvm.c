/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-02-13 23:01:29 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
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

/* CVODE data */

void *cvode_mem;   /* CVODES solver memory */
void *bp_data;     /* Preconditioner memory (BandPre or BBDPre) */
N_Vector y;        /* solution vector */
N_Vector yQ;       /* quadratures vector */
N_Vector *yS;      /* sensitivity vectors */
int N;             /* problem dimension */
int Nq;            /* number of quadratures */
int Ng;            /* number of root functions */
int Ns;            /* number of sensitivities */
int Nd;            /* number of data points */
int Nc;            /* number of check points */
int ls;            /* linear solver type */
int pm;            /* preconditioner module */
int ism;           /* sensitivity method */

void *cvadj_mem;   /* CVODES adjoint memory */
int interp;
N_Vector yB;
N_Vector yQB;
int NB;
int NqB;
int lsB;
int pmB;

/* Matlab data */

mxArray *mx_mtlb_RHSfct;
mxArray *mx_mtlb_QUADfct;
mxArray *mx_mtlb_JACfct;
mxArray *mx_mtlb_PSETfct;
mxArray *mx_mtlb_PSOLfct;
mxArray *mx_mtlb_GLOCfct;
mxArray *mx_mtlb_GCOMfct;

mxArray *mx_mtlb_Gfct;

mxArray *mx_mtlb_SRHSfct;

mxArray *mx_mtlb_RHSfctB;
mxArray *mx_mtlb_QUADfctB;
mxArray *mx_mtlb_JACfctB;
mxArray *mx_mtlb_PSETfctB;
mxArray *mx_mtlb_PSOLfctB;
mxArray *mx_mtlb_GLOCfctB;
mxArray *mx_mtlb_GCOMfctB;

mxArray *mx_mtlb_data;

/* Monitors */

booleantype monitor;
mxArray *mx_mtlb_MONfct;
mxArray *mx_mtlb_MONdata;

booleantype monitorB;
mxArray *mx_mtlb_MONfctB;
mxArray *mx_mtlb_MONdataB;

/*
 * ---------------------------------------------------------------------------------
 * static function prototypes
 * ---------------------------------------------------------------------------------
 */

static int CVM_init();
static int CVM_initB();
static int CVM_final();

static int CVM_Malloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_MallocB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_ReInit(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int CVM_ReInitB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
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
     2 - initialize backward solver
     3 - reinitialize CVODES solver (NYI)
     4 - reinitialize backward solver (NYI)
     5 - solve problem
     6 - solve backward problem
     7 - get integrator stats
     8 - get backward integrator stats
     9 - extract data from cvode_mem
    10 - set one optional input at a time (NYI)
    11 - finalize
  */

  mode = (int)mxGetScalar(prhs[0]);

  switch(mode) {
  case 1:
    CVM_init();
    CVM_Malloc(nlhs, plhs,nrhs-1,&prhs[1]);
    mexLock();
    mexMakeArrayPersistent(mx_mtlb_data);
    mexMakeArrayPersistent(mx_mtlb_RHSfct);
    mexMakeArrayPersistent(mx_mtlb_QUADfct);
    mexMakeArrayPersistent(mx_mtlb_Gfct);
    mexMakeArrayPersistent(mx_mtlb_JACfct);
    mexMakeArrayPersistent(mx_mtlb_PSETfct);
    mexMakeArrayPersistent(mx_mtlb_PSOLfct);
    mexMakeArrayPersistent(mx_mtlb_GLOCfct);
    mexMakeArrayPersistent(mx_mtlb_GCOMfct);
    mexMakeArrayPersistent(mx_mtlb_SRHSfct);
    mexMakeArrayPersistent(mx_mtlb_MONfct);
    mexMakeArrayPersistent(mx_mtlb_MONdata);
    break;
  case 2:
    CVM_initB();
    CVM_MallocB(nlhs, plhs,nrhs-1,&prhs[1]);
    mexMakeArrayPersistent(mx_mtlb_RHSfctB);
    mexMakeArrayPersistent(mx_mtlb_QUADfctB);
    mexMakeArrayPersistent(mx_mtlb_JACfctB);
    mexMakeArrayPersistent(mx_mtlb_PSETfctB);
    mexMakeArrayPersistent(mx_mtlb_PSOLfctB);
    mexMakeArrayPersistent(mx_mtlb_GLOCfctB);
    mexMakeArrayPersistent(mx_mtlb_GCOMfctB);
    mexMakeArrayPersistent(mx_mtlb_MONfctB);
    mexMakeArrayPersistent(mx_mtlb_MONdataB);
    break;
  case 3:
    CVM_ReInit(nlhs, plhs,nrhs-1,&prhs[1]);
    break;
  case 4:
    CVM_ReInitB(nlhs, plhs,nrhs-1,&prhs[1]);
    break;
  case 5:
    CVM_Solve(nlhs, plhs,nrhs-1,&prhs[1]);
    break;
  case 6:
    CVM_SolveB(nlhs, plhs,nrhs-1,&prhs[1]);
    break;
  case 7:
    CVM_Stats(nlhs, plhs,nrhs-1,&prhs[1]);
    break;
  case 8:
    CVM_StatsB(nlhs, plhs,nrhs-1,&prhs[1]);
    break;
  case 9:
    CVM_Get(nlhs, plhs,nrhs-1,&prhs[1]);
    break;
  case 10:
    CVM_Set(nlhs, plhs,nrhs-1,&prhs[1]);
    break;
  case 11:
    CVM_Free(nlhs, plhs,nrhs-1,&prhs[1]);
    CVM_final();
    break;
  }

}

/*
 * ---------------------------------------------------------------------------------
 * Exported procedures
 * ---------------------------------------------------------------------------------
 */

static int CVM_Malloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double t0, *y0;
  int i, is, status;
  double *tmp;
  mxArray *pfield;
  char *pfield_name;

  int vec_type;
  mxArray *mx_in[1], *mx_out[2];
  mxArray *mx_comm;

  int lmm, iter, maxord;
  long int mxsteps;

  booleantype quad, fsa, asa;

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

  double *yS0;

  booleantype errconS;
  int itolS;
  realtype reltolS;
  realtype *SabstolS, *VabstolS;
  N_Vector *NV_abstolS;
  void *abstolS;

  int *plist;
  double *p;
  double *pbar;

  int Srhs;
  double rho;

  double hin, hmax, hmin;

  double tstop;
  booleantype tstopSet;

  booleantype sld;

  int mupper, mlower;
  int ptype, gstype, maxl;
  int mudq, mldq;
  double dqrely;

  mxArray *options;

  N_Vector yQ_tmp, *yS_tmp;

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
   * - RHS function
   * - initial time
   * - initial conditions
   * - integration options
   * - user data
   * -----------------------------
   */

  /* Matlab user-provided function */

  mx_mtlb_RHSfct = mxDuplicateArray(prhs[0]);
  
  /* initial time */

  t0 = (double)mxGetScalar(prhs[1]);

  /* Initial conditions */

  y0 = mxGetPr(prhs[2]);
  N = mxGetM(prhs[2]);

  /* Integrator Options -- optional argument */

  status = get_IntgrOptions(prhs[3],
                            &lmm, &iter, &maxord, &sld, &mxsteps,
                            &itol, &reltol, &Sabstol, &Vabstol,
                            &hin, &hmax, &hmin, &tstop, &tstopSet,
                            &Ng, &mx_mtlb_Gfct,
                            &quad, &fsa, &asa,
                            &monitor, &mx_mtlb_MONfct, &mx_mtlb_MONdata);

  /* User data -- optional argument */

  mx_mtlb_data = mxDuplicateArray(prhs[4]);

  /* 
   * ---------------------------------------
   * Set initial conditions and tolerances
   * ---------------------------------------
   */

  y = NewVector(N);

  N_VSetArrayPointer(y0, y);

  switch (itol) {
  case CV_SS:
    abstol = (void *) &Sabstol;
    break;
  case CV_SV:
    NV_abstol = N_VCloneEmpty(y);
    N_VSetArrayPointer(Vabstol, NV_abstol);
    abstol = (void *) NV_abstol;
    break;
  }

  /* 
   * ---------------------------------------
   * Create cvode object and allocate memory
   * ---------------------------------------
   */

  cvode_mem = CVodeCreate(lmm, iter);

  /* attach error handler function */
  status = CVodeSetErrHandlerFn(cvode_mem, mtlb_CVodeErrHandler, NULL);

  /* Call CVodeMalloc */

  status = CVodeMalloc(cvode_mem, mtlb_CVodeRhs, t0, y, itol, reltol, abstol);

  if (itol == CV_SV) {
    N_VDestroy(NV_abstol);
  }

  /* Redirect output */
  status = CVodeSetErrFile(cvode_mem, stdout);

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
  if ( !mxIsEmpty(mx_mtlb_Gfct) && (Ng > 0) ) {
    status = CVodeRootInit(cvode_mem, Ng, mtlb_CVodeGfct, NULL);
  }

  /* Quadratures? */
  if ( quad ) {

    status = get_QuadOptions(prhs[3],
                             &Nq, &yQ0, &mx_mtlb_QUADfct,
                             &errconQ, &itolQ, &reltolQ, &SabstolQ, &VabstolQ);

    if(status)
      mexErrMsgTxt("CVodeMalloc:: illegal quadrature input.");      

    yQ = NewVector(Nq);

    N_VSetArrayPointer(yQ0, yQ);

    status = CVodeQuadMalloc(cvode_mem, mtlb_CVodeQUADfct, yQ);

    switch (itolQ) {
    case CV_SS:
      abstolQ = (void *) &SabstolQ;
      break;
    case CV_SV:
      NV_abstolQ = N_VCloneEmpty(yQ);
      N_VSetArrayPointer(VabstolQ, NV_abstolQ);
      abstolQ = (void *) NV_abstolQ;
      break;
    }

    status = CVodeSetQuadErrCon(cvode_mem, errconQ, itolQ, reltolQ, abstolQ);

    if (itolQ == CV_SV) {
      N_VDestroy(NV_abstolQ);
    }

  }

  /* Forward Sensitivities? */
  if (fsa ) {

    plist = NULL;
    pbar = NULL;

    status = get_FSAOptions(prhs[3], 
                            &Ns, &yS0, &ism, 
                            &pfield_name, &plist, &pbar,
                            &Srhs, &mx_mtlb_SRHSfct, &rho,
                            &errconS, &itolS, &reltolS, &SabstolS, &VabstolS);

    if(status)
      mexErrMsgTxt("CVodeMalloc:: illegal forward sensitivity input.");   

    if (mxIsEmpty(mx_mtlb_SRHSfct)) {
      if (pfield_name == NULL)
        mexErrMsgTxt("CVodeMalloc:: pfield required but was not provided.");
      pfield = mxGetField(mx_mtlb_data,0,pfield_name);
      if (pfield == NULL)
        mexErrMsgTxt("CVodeMalloc:: illegal pfield input.");
      p = mxGetPr(pfield);
      mxFree(pfield_name);
    }

    yS = N_VCloneEmptyVectorArray(Ns, y);
    for (is=0;is<Ns;is++)
      N_VSetArrayPointer(&yS0[is*N], yS[is]);
    
    status = CVodeSensMalloc(cvode_mem, Ns, ism, yS);

    switch (itolS) {
    case CV_SS:
      abstolS = (void *) SabstolS;
      break;
    case CV_SV:
      NV_abstolS = N_VCloneEmptyVectorArray(Ns, y);
      for (is=0;is<Ns;is++)
        N_VSetArrayPointer(&VabstolS[is*N], NV_abstolS[is]);
      abstolS = (void *) NV_abstolS;
      break;
    case CV_EE:
      abstolS = NULL;
      break;
    }
    
    status = CVodeSetSensTolerances(cvode_mem, itolS, reltolS, abstolS);

    switch (itolS) {
    case CV_SS:
      free(SabstolS);
      break;
    case CV_SV:
      N_VDestroyVectorArray(NV_abstolS, Ns);
      break;
    }

    status = CVodeSetSensParams(cvode_mem, p, pbar, plist);

    if (plist != NULL) free(plist);
    if (pbar != NULL)  free(pbar);

    status = CVodeSetSensRho(cvode_mem, rho);

    status = CVodeSetSensErrCon(cvode_mem, errconS);
    
    if (Srhs == 1) {
      status = CVodeSetSensRhs1Fn(cvode_mem, mtlb_CVodeSensRhs1);
    } else if (Srhs == 2) {
      status = CVodeSetSensRhsFn(cvode_mem, mtlb_CVodeSensRhs);
    }
    
  }

  /* Need a linear solver? */
  if (iter == CV_NEWTON) {

    status = get_LinSolvOptions(prhs[3], &ls, 
                                &mupper, &mlower,
                                &mudq, &mldq, &dqrely,
                                &ptype, &gstype, &maxl, &pm,
                                &mx_mtlb_JACfct, &mx_mtlb_PSETfct, &mx_mtlb_PSOLfct,
                                &mx_mtlb_GLOCfct, &mx_mtlb_GCOMfct);

    switch (ls) {
    case LS_DENSE:
      status = CVDense(cvode_mem, N);
      if (!mxIsEmpty(mx_mtlb_JACfct))
        status = CVDenseSetJacFn(cvode_mem, mtlb_CVodeDenseJac, NULL);
      break;
    case LS_DIAG:
      status = CVDiag(cvode_mem);
      break;
    case LS_BAND:
      status = CVBand(cvode_mem, N, mupper, mlower);
      if (!mxIsEmpty(mx_mtlb_JACfct))
        status = CVBandSetJacFn(cvode_mem, mtlb_CVodeBandJac, NULL);
      break;
    case LS_SPGMR:
      switch (pm) {
      case PM_NONE:
        status = CVSpgmr(cvode_mem, ptype, maxl);
        if (!mxIsEmpty(mx_mtlb_PSOLfct)) {
          if (!mxIsEmpty(mx_mtlb_PSETfct))
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
        if (!mxIsEmpty(mx_mtlb_GCOMfct)) {
          bp_data = CVBBDPrecAlloc(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_CVodeBBDgloc, mtlb_CVodeBBDgcom);
        } else {
          bp_data = CVBBDPrecAlloc(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_CVodeBBDgloc, NULL);
        }
        CVBBDSpgmr(cvode_mem, ptype, maxl, bp_data);
        break;
      }
      status = CVSpilsSetGSType(cvode_mem, gstype);
      if (!mxIsEmpty(mx_mtlb_JACfct))
        status = CVSpilsSetJacTimesVecFn(cvode_mem, mtlb_CVodeSpilsJac, NULL);
      break;
    case LS_SPBCG:
      switch (pm) {
      case PM_NONE:
        status = CVSpbcg(cvode_mem, ptype, maxl);
        if (!mxIsEmpty(mx_mtlb_PSOLfct)) {
          if (!mxIsEmpty(mx_mtlb_PSETfct))
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
        if (!mxIsEmpty(mx_mtlb_GCOMfct)) {
          bp_data = CVBBDPrecAlloc(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_CVodeBBDgloc, mtlb_CVodeBBDgcom);
        } else {
          bp_data = CVBBDPrecAlloc(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_CVodeBBDgloc, NULL);
        }
        CVBBDSpbcg(cvode_mem, ptype, maxl, bp_data);
        break;
      }
      if (!mxIsEmpty(mx_mtlb_JACfct))
        status = CVSpilsSetJacTimesVecFn(cvode_mem, mtlb_CVodeSpilsJac, NULL);
      break;
    case LS_SPTFQMR:
      switch (pm) {
      case PM_NONE:
        status = CVSptfqmr(cvode_mem, ptype, maxl);
        if (!mxIsEmpty(mx_mtlb_PSOLfct)) {
          if (!mxIsEmpty(mx_mtlb_PSETfct))
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
        if (!mxIsEmpty(mx_mtlb_GCOMfct)) {
          bp_data = CVBBDPrecAlloc(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_CVodeBBDgloc, mtlb_CVodeBBDgcom);
        } else {
          bp_data = CVBBDPrecAlloc(cvode_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_CVodeBBDgloc, NULL);
        }
        CVBBDSptfqmr(cvode_mem, ptype, maxl, bp_data);
        break;
      }
      if (!mxIsEmpty(mx_mtlb_JACfct))
        status = CVSpilsSetJacTimesVecFn(cvode_mem, mtlb_CVodeSpilsJac, NULL);
      break;
    }

  } else {

    ls = LS_NONE;

  }

  /* Need to initialize adjoint memeory? */

  if (asa) {

    status = get_ASAOptions(prhs[3], &Nd, &interp);

    if (status)
      mexErrMsgTxt("CVodeMalloc:: illegal ASA input.");      

    cvadj_mem = CVadjMalloc(cvode_mem, Nd, interp);

  } else {

    cvadj_mem = NULL;

  }

  /* Do we monitor? */
  
  if (monitor) {
    if (Nq > 0)
      yQ_tmp = yQ;
    else
      yQ_tmp = NULL;
    if (Ns > 0) 
      yS_tmp = yS;
    else
      yS_tmp = NULL;
    mtlb_CVodeMonitor(0, t0, y, yQ_tmp, yS_tmp);    
  }

  return(0);

}

static int CVM_MallocB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double tB0, *yB0;
  int i, status;

  int lmmB, iterB, maxordB;

  long int mxstepsB;

  booleantype quadB;
  booleantype fsaB;        /* ignored */
  booleantype asaB;        /* ignored */

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
  int NgB;                  /* ignored */
  mxArray *mx_mtlb_GfctB;   /* ignored */

  int mupperB, mlowerB;
  int ptypeB, gstypeB, maxlB;
  int mudqB, mldqB;
  double dqrelyB;

  mxArray *optionsB;

  N_Vector yQB_tmp;

  /* 
   * -----------------------------
   * Finalize Forward monitoring
   * -----------------------------
   */


  if (monitor) {
    mtlb_CVodeMonitor(2, 0.0, NULL, NULL, NULL);
    monitor = FALSE;
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

  mx_mtlb_RHSfctB = mxDuplicateArray(prhs[0]);
  
  /* Final time */

  tB0 = (double)mxGetScalar(prhs[1]);

  /* Final conditions */

  yB0 = mxGetPr(prhs[2]);
  NB = mxGetM(prhs[2]);

  /* Integrator Options -- optional argument */

  status = get_IntgrOptions(prhs[3],
                            &lmmB, &iterB, &maxordB, &sldB, &mxstepsB,
                            &itolB, &reltolB, &SabstolB, &VabstolB,
                            &hinB, &hmaxB, &hminB, &tstopB, &tstopSetB,
                            &NgB, &mx_mtlb_GfctB, 
                            &quadB, &fsaB, &asaB,
                            &monitorB, &mx_mtlb_MONfctB, &mx_mtlb_MONdataB);
  /* 
   * ---------------------------------------
   * Set final conditions and tolerances
   * ---------------------------------------
   */

  yB = NewVector(NB);

  N_VSetArrayPointer(yB0, yB);

  switch (itolB) {
  case CV_SS:
    abstolB = (void *) &SabstolB;
    break;
  case CV_SV:
    NV_abstolB = N_VCloneEmpty(yB);
    N_VSetArrayPointer(VabstolB, NV_abstolB);
    abstolB = (void *) VabstolB;
    break;
  }

  /* 
   * ----------------------------------------------------------
   * Create cvode object for backward phase and allocate memory
   * ----------------------------------------------------------
   */

  status = CVodeCreateB(cvadj_mem, lmmB, iterB);

  /* Call CVodeMallocB */
  status = CVodeMallocB(cvadj_mem, mtlb_CVodeRhsB, tB0, yB, itolB, reltolB, abstolB);

  if (itolB == CV_SV) {
    N_VDestroy(NV_abstolB);
  }

  /* Redirect output */
  status = CVodeSetErrFileB(cvadj_mem, stdout);

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
  if ( quadB ) {
    
    status = get_QuadOptions(prhs[3],
                             &NqB, &yQB0, &mx_mtlb_QUADfctB,
                             &errconQB, &itolQB, &reltolQB, &SabstolQB, &VabstolQB);

    if(status)
      mexErrMsgTxt("CVodeMallocB:: illegal quadrature input.");      

    yQB = NewVector(NqB);

    N_VSetArrayPointer(yQB0, yQB);

    status = CVodeQuadMallocB(cvadj_mem, mtlb_CVodeQUADfctB, yQB);
    
    switch (itolQB) {
    case CV_SS:
      abstolQB = (void *) &SabstolQB;
      break;
    case CV_SV:
      NV_abstolQB = N_VCloneEmpty(yQB);
      N_VSetArrayPointer(VabstolQB, NV_abstolQB);
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

    status = get_LinSolvOptions(prhs[3], &lsB, 
                                &mupperB, &mlowerB,
                                &mudqB, &mldqB, &dqrelyB,
                                &ptypeB, &gstypeB, &maxlB, &pmB,
                                &mx_mtlb_JACfctB, &mx_mtlb_PSETfctB, &mx_mtlb_PSOLfctB,
                                &mx_mtlb_GLOCfctB, &mx_mtlb_GCOMfctB);

    switch(lsB) {
    case LS_DENSE:
      status = CVDenseB(cvadj_mem, NB);
      if (!mxIsEmpty(mx_mtlb_JACfctB))
        status = CVDenseSetJacFnB(cvadj_mem, mtlb_CVodeDenseJacB, NULL);
      break;
    case LS_DIAG:
      status = CVDiagB(cvadj_mem);
      break;
    case LS_BAND:
      status = CVBandB(cvadj_mem, NB, mupperB, mlowerB);
      if (!mxIsEmpty(mx_mtlb_JACfctB))
        status = CVBandSetJacFnB(cvadj_mem, mtlb_CVodeBandJacB, NULL);
      break;
    case LS_SPGMR:
      switch (pmB) {
      case PM_NONE:
        status = CVSpgmrB(cvadj_mem, ptypeB, maxlB);
        if (!mxIsEmpty(mx_mtlb_PSOLfctB)) {
          if (!mxIsEmpty(mx_mtlb_PSETfctB))
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
        if (!mxIsEmpty(mx_mtlb_GCOMfctB)) {
          status = CVBBDPrecAllocB(cvadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_CVodeBBDglocB, mtlb_CVodeBBDgcomB);
        } else {
          status = CVBBDPrecAllocB(cvadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_CVodeBBDglocB, NULL);
        }
        CVBBDSpgmrB(cvadj_mem, ptypeB, maxlB);
        break;
      }
      status = CVSpilsSetGSTypeB(cvadj_mem, gstypeB);
      if (!mxIsEmpty(mx_mtlb_JACfctB))
        status = CVSpilsSetJacTimesVecFnB(cvadj_mem, mtlb_CVodeSpilsJacB, NULL);
      break;
    case LS_SPBCG:
      switch (pmB) {
      case PM_NONE:
        status = CVSpbcgB(cvadj_mem, ptypeB, maxlB);
        if (!mxIsEmpty(mx_mtlb_PSOLfctB)) {
          if (!mxIsEmpty(mx_mtlb_PSETfctB))
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
        if (!mxIsEmpty(mx_mtlb_GCOMfctB)) {
          status = CVBBDPrecAllocB(cvadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_CVodeBBDglocB, mtlb_CVodeBBDgcomB);
        } else {
          status = CVBBDPrecAllocB(cvadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_CVodeBBDglocB, NULL);
        }
        CVBBDSpbcgB(cvadj_mem, ptypeB, maxlB);
        break;
      }
      if (!mxIsEmpty(mx_mtlb_JACfctB))
        status = CVSpilsSetJacTimesVecFnB(cvadj_mem, mtlb_CVodeSpilsJacB, NULL);
      break;
    case LS_SPTFQMR:
      switch (pmB) {
      case PM_NONE:
        status = CVSptfqmrB(cvadj_mem, ptypeB, maxlB);
        if (!mxIsEmpty(mx_mtlb_PSOLfctB)) {
          if (!mxIsEmpty(mx_mtlb_PSETfctB))
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
        if (!mxIsEmpty(mx_mtlb_GCOMfctB)) {
          status = CVBBDPrecAllocB(cvadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_CVodeBBDglocB, mtlb_CVodeBBDgcomB);
        } else {
          status = CVBBDPrecAllocB(cvadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_CVodeBBDglocB, NULL);
        }
        CVBBDSptfqmrB(cvadj_mem, ptypeB, maxlB);
        break;
      }
      if (!mxIsEmpty(mx_mtlb_JACfctB))
        status = CVSpilsSetJacTimesVecFnB(cvadj_mem, mtlb_CVodeSpilsJacB, NULL);
      break;
    }

  } else {

    lsB = LS_NONE;

  }

  /* Do we monitor? */
  
  if (monitorB) {
    if (NqB > 0)
      yQB_tmp = yQB;
    else
      yQB_tmp = NULL;
    mtlb_CVodeMonitorB(0, tB0, yB, yQB_tmp);
  }

  return(0);
}

static int CVM_ReInit(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}

static int CVM_ReInitB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}

static int CVM_Solve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double tout, tret, *ydata;
  int buflen, status, itask, i, is;
  char *bufval;

  int itask1;
  booleantype iret;
  double h;

  N_Vector yQ_tmp, *yS_tmp;

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

  /* Solution vector */
  plhs[2] = mxCreateDoubleMatrix(N,1,mxREAL);
  N_VSetArrayPointer(mxGetPr(plhs[2]), y);

  if (!monitor) {

    if (Nd==0) {
      status = CVode(cvode_mem, tout, y, &tret, itask);
    } else {
      status = CVodeF(cvadj_mem, tout, y, &tret, itask, &Nc);
    }

  } else {

    if (itask == CV_NORMAL)              {iret = FALSE; itask1 = CV_ONE_STEP;}
    else if (itask == CV_ONE_STEP)       {iret = TRUE;  itask1 = CV_ONE_STEP;}
    else if (itask == CV_NORMAL_TSTOP)   {iret = FALSE; itask1 = CV_ONE_STEP_TSTOP;}
    else if (itask == CV_ONE_STEP_TSTOP) {iret = TRUE;  itask1 = CV_ONE_STEP_TSTOP;}
    
    while(1) {
      
      if (Nd==0) {
        status = CVode(cvode_mem, tout, y, &tret, itask1);
      } else {
        status = CVodeF(cvadj_mem, tout, y, &tret, itask1, &Nc);
      }

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

      if (Nq > 0) { 
        status = CVodeGetQuad(cvode_mem, tret, yQ);
        yQ_tmp = yQ;
      } else {
        yQ_tmp == NULL;
      }

      if (Ns > 0) {
        status = CVodeGetSens(cvode_mem, tret, yS);
        yS_tmp = yS;
      } else {
        yS_tmp = NULL;
      }

      mtlb_CVodeMonitor(1, tret, y, yQ_tmp, yS_tmp);
      
      /* break if we need to */
      if(iret)  break;
      
    };

  }

  /* CVODE return flag */
  plhs[0] = mxCreateScalarDouble((double)status);

  /* Return time */
  plhs[1] = mxCreateScalarDouble(tret);

  if (nlhs == 3)
    return;

  if (nlhs == 4) {
  
    if (Nq > 0) {
      plhs[3] = mxCreateDoubleMatrix(Nq,1,mxREAL);
      N_VSetArrayPointer(mxGetPr(plhs[3]), yQ);
      status = CVodeGetQuad(cvode_mem, tret, yQ);
      return;
    } else if (Ns > 0) {
      plhs[3] = mxCreateDoubleMatrix(N,Ns,mxREAL);
      ydata = mxGetPr(plhs[3]);
      for (is=0; is<Ns; is++)
        N_VSetArrayPointer(&ydata[is*N], yS[is]);
      status = CVodeGetSens(cvode_mem, tret, yS);
      return;
    } else {
      mexErrMsgTxt("CVode:: too many output arguments.");
      return;
    }

  }

  if ( (Nq ==0) | (Ns == 0) ) {
    mexErrMsgTxt("CVode:: too many output arguments.");
    return;
  }

  /* Quadratures */
  plhs[3] = mxCreateDoubleMatrix(Nq,1,mxREAL);
  N_VSetArrayPointer(mxGetPr(plhs[3]), yQ);
  status = CVodeGetQuad(cvode_mem, tret, yQ);

  /* Sensitivities */
  plhs[4] = mxCreateDoubleMatrix(N,Ns,mxREAL);
  ydata = mxGetPr(plhs[4]);
  for (is=0; is<Ns; is++)
    N_VSetArrayPointer(&ydata[is*N], yS[is]);
  status = CVodeGetSens(cvode_mem, tret, yS);

  return(0);
}


static int CVM_SolveB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double toutB, tretB;
  int itaskB;
  double *tmp;
  int buflen, status, i;
  char *bufval;

  void *cvode_memB;
  booleantype iretB;
  double hB;

  N_Vector yQB_tmp;

  /* Exract toutB */

  toutB = (double)mxGetScalar(prhs[0]);

  /* Extract itaskB */

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Normal")) itaskB = CV_NORMAL;
  if(!strcmp(bufval,"OneStep")) itaskB = CV_ONE_STEP;

  /* Call CVodeB */

  /* Solution vector */
  plhs[2] = mxCreateDoubleMatrix(NB,1,mxREAL);
  N_VSetArrayPointer(mxGetPr(plhs[2]), yB);

  if (!monitorB) {

    status = CVodeB(cvadj_mem, toutB, yB, &tretB, itaskB);

  } else {
    
    if (itaskB == CV_NORMAL)        iretB = FALSE;
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

      if (NqB > 0) {
        status = CVodeGetQuadB(cvadj_mem, yQB);
        yQB_tmp = yQB;
      } else {
        yQB_tmp = NULL;
      }

      mtlb_CVodeMonitorB(1, tretB, yB, yQB_tmp);

      /* break if we need to */
      if(iretB)  break;

    };

  }

  /* CVodeB return flag */
  plhs[0] = mxCreateScalarDouble((double)status);

  /* Return time */
  plhs[1] = mxCreateScalarDouble(tretB);

  if (nlhs == 3)
    return;

  if (NqB > 0) {
    plhs[3] = mxCreateDoubleMatrix(NqB,1,mxREAL);
    N_VSetArrayPointer(mxGetPr(plhs[3]), yQB);
    status = CVodeGetQuadB(cvadj_mem, yQB);
    return;
  } else {
    mexErrMsgTxt("CVodeB:: too many output arguments.");
    return;
  }

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
    "nniSTGR1",
    "ncfnSTGR1"
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
  long int *nniSTGR1, *ncfnSTGR1;

  int i, flag;
  mxArray *mx_root, *mx_ls, *mx_quad, *mx_fsa, *mx_asa;
  mxArray *mx_rootsfound;
  mxArray *mx_nniSTGR1, *mx_ncfnSTGR1;
  double *tmp;
  int nfields;

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

  if (Nq >0) {

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
    
    flag = CVDenseGetNumJacEvals(cvode_mem, &njeD);
    flag = CVDenseGetNumRhsEvals(cvode_mem, &nfeD);
    
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
      
    flag = CVBandGetNumJacEvals(cvode_mem, &njeB);
    flag = CVBandGetNumRhsEvals(cvode_mem, &nfeB);
      
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

  /* forward Sensitivity Statistics */

  if (Ns > 0) {

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
    
    if ( ism == CV_STAGGERED1 ) {
      
      nniSTGR1 = (long int *) malloc(Ns*sizeof(long int));
      ncfnSTGR1 = (long int *) malloc(Ns*sizeof(long int));
      
      flag = CVodeGetNumStgrSensNonlinSolvIters(cvode_mem, nniSTGR1);
      flag = CVodeGetNumStgrSensNonlinSolvConvFails(cvode_mem, ncfnSTGR1);
      
      mx_nniSTGR1 = mxCreateDoubleMatrix(1,Ns,mxREAL);
      tmp = mxGetPr(mx_nniSTGR1);
      for (i=0;i<Ns;i++) 
        tmp[i] = nniSTGR1[i];
      mx_ncfnSTGR1 = mxCreateDoubleMatrix(1,Ns,mxREAL);
      tmp = mxGetPr(mx_ncfnSTGR1);
      for (i=0;i<Ns;i++) 
        tmp[i] = ncfnSTGR1[i];
      
    } else {

      mx_nniSTGR1 = mxCreateDoubleMatrix(0,0,mxREAL);
      mx_ncfnSTGR1 = mxCreateDoubleMatrix(0,0,mxREAL);

    }
    
    mxSetField(mx_fsa, 0, "nniSTGR1",  mx_nniSTGR1);
    mxSetField(mx_fsa, 0, "ncfnSTGR1", mx_ncfnSTGR1);
    

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

  int i, flag;
  mxArray *mx_ls, *mx_quad;
  double *tmp;
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

  if (NqB >0) {

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
    
    flag = CVDenseGetNumJacEvals(cvode_memB, &njeD);
    flag = CVDenseGetNumRhsEvals(cvode_memB, &nfeD);
    
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
      
    flag = CVBandGetNumJacEvals(cvode_memB, &njeB);
    flag = CVBandGetNumRhsEvals(cvode_memB, &nfeB);
      
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
  unsigned int addr;
  int key, k, which, status, i, nfields;
  mxArray *mx_y, *mx_yd;

  CheckPointRec *ckpnt;
  const char *fnames_ckpnt[]={
    "addr",
    "next",
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
    N_VSetArrayPointer(mxGetPr(plhs[0]), y);
    status = CVodeGetDky(cvode_mem, t, k, y);

    break;

  case 2:    /* ErrorWeights */

    ewt = N_VClone(y);

    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    N_VSetArrayPointer(mxGetPr(plhs[0]), ewt);
    status = CVodeGetErrWeights(cvode_mem, ewt);

    N_VDestroy(ewt);

    break;

  case 3:    /* not used */

    break;

  case 4:    /* CheckPointsInfo */

    ckpnt = (CheckPointRec *) malloc ( (Nc+1)*sizeof(CheckPointRec));
    CVadjGetCheckPointsInfo(cvadj_mem, ckpnt);
    nfields = sizeof(fnames_ckpnt)/sizeof(*fnames_ckpnt);
    plhs[0] = mxCreateStructMatrix(Nc+1, 1, nfields, fnames_ckpnt);
    for (i=0; i<=Nc; i++) {
      mxSetField(plhs[0], i, "addr",  mxCreateScalarDouble((double)(ckpnt[Nc-i].my_addr)));
      mxSetField(plhs[0], i, "next",  mxCreateScalarDouble((double)(ckpnt[Nc-i].next_addr)));
      mxSetField(plhs[0], i, "t0",    mxCreateScalarDouble((double)(ckpnt[Nc-i].t0)));
      mxSetField(plhs[0], i, "t1",    mxCreateScalarDouble((double)(ckpnt[Nc-i].t1)));
      mxSetField(plhs[0], i, "nstep", mxCreateScalarDouble((double)(ckpnt[Nc-i].nstep)));
      mxSetField(plhs[0], i, "order", mxCreateScalarDouble((double)(ckpnt[Nc-i].order)));
      mxSetField(plhs[0], i, "step",  mxCreateScalarDouble((double)(ckpnt[Nc-i].step)));
    }

    break;

  case 5:    /* CurrentCheckPoint */

    status = CVadjGetCurrentCheckPoint(cvadj_mem, &addr);

    plhs[0] = mxCreateScalarDouble((double)addr);

    break;

  case 6:    /* DataPointInfo */
    
    which = (int) (*mxGetPr(prhs[1]));

    if (interp == CV_HERMITE) {
    
      yd = N_VClone(y);

      plhs[0] = mxCreateCellMatrix(1, 3);
      mx_y  = mxCreateDoubleMatrix(N,1,mxREAL);
      N_VSetArrayPointer(mxGetPr(mx_y), y);
      mxSetCell(plhs[0],2,mx_y);
      mx_yd = mxCreateDoubleMatrix(N,1,mxREAL);
      N_VSetArrayPointer(mxGetPr(mx_yd), yd);
      mxSetCell(plhs[0],3,mx_y);

      status = CVadjGetDataPointHermite(cvadj_mem, which, &t, y, yd);

      mxSetCell(plhs[0],2,mxCreateScalarDouble(t));

      N_VDestroy(yd);
    }

    break;

  }

  return(0);
}

static int CVM_Free(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  if (monitor) {
    mtlb_CVodeMonitor(2, 0.0, NULL, NULL, NULL);
  }

  N_VDestroy(y);
  if (Nq > 0) N_VDestroy(yQ);
  if (pm == PM_BANDPRE) CVBandPrecFree(&bp_data);
  if (pm == PM_BBDPRE)  CVBBDPrecFree(&bp_data);

  if (Ns > 0) N_VDestroyVectorArray(yS, Ns);

  CVodeFree(&cvode_mem);

  if (Nd > 0) {

    if (NB > 0) {
      N_VDestroy(yB);
      if (monitorB) 
        mtlb_CVodeMonitorB(2, 0.0, NULL, NULL);
    }

    if (NqB > 0) N_VDestroy(yQB);
   
    CVadjFree(&cvadj_mem);
  }

  return(0);
}


static int CVM_init()
{
  mxArray *empty;

  /* Initialize global dimensions to zero */
  N   = 0;
  Nq  = 0;
  Ng  = 0;
  Ns  = 0;
  Nd  = 0;
  Nc  = 0;
  NB  = 0;
  NqB = 0;
  monitor = FALSE;

  empty = mxCreateDoubleMatrix(0,0,mxREAL);

  mx_mtlb_data     = mxDuplicateArray(empty);

  mx_mtlb_RHSfct   = mxDuplicateArray(empty);
  mx_mtlb_Gfct     = mxDuplicateArray(empty);
  mx_mtlb_QUADfct  = mxDuplicateArray(empty);
  mx_mtlb_JACfct   = mxDuplicateArray(empty);
  mx_mtlb_PSETfct  = mxDuplicateArray(empty);
  mx_mtlb_PSOLfct  = mxDuplicateArray(empty);
  mx_mtlb_GLOCfct  = mxDuplicateArray(empty);
  mx_mtlb_GCOMfct  = mxDuplicateArray(empty);

  mx_mtlb_MONfct   = mxDuplicateArray(empty);
  mx_mtlb_MONdata  = mxDuplicateArray(empty);
  
  mx_mtlb_SRHSfct  = mxDuplicateArray(empty);

  mxDestroyArray(empty);

}

static int CVM_initB()
{
  mxArray *empty;

  empty = mxCreateDoubleMatrix(0,0,mxREAL);

  mx_mtlb_RHSfctB  = mxDuplicateArray(empty);
  mx_mtlb_QUADfctB = mxDuplicateArray(empty);
  mx_mtlb_JACfctB  = mxDuplicateArray(empty);
  mx_mtlb_PSETfctB = mxDuplicateArray(empty);
  mx_mtlb_PSOLfctB = mxDuplicateArray(empty);
  mx_mtlb_GLOCfctB = mxDuplicateArray(empty);
  mx_mtlb_GCOMfctB = mxDuplicateArray(empty);
  mx_mtlb_MONfctB  = mxDuplicateArray(empty);
  mx_mtlb_MONdataB = mxDuplicateArray(empty);

  mxDestroyArray(empty);

  return(0);
}

static int CVM_final()
{
  mexUnlock();
  
  mxDestroyArray(mx_mtlb_data);
  
  mxDestroyArray(mx_mtlb_RHSfct);
  mxDestroyArray(mx_mtlb_QUADfct);
  mxDestroyArray(mx_mtlb_Gfct);
  mxDestroyArray(mx_mtlb_JACfct);
  mxDestroyArray(mx_mtlb_PSETfct);
  mxDestroyArray(mx_mtlb_PSOLfct);

  mxDestroyArray(mx_mtlb_MONfct);
  mxDestroyArray(mx_mtlb_MONdata);

  mxDestroyArray(mx_mtlb_SRHSfct);
  
  if (NB > 0) {
    mxDestroyArray(mx_mtlb_RHSfctB);
    mxDestroyArray(mx_mtlb_QUADfctB);
    mxDestroyArray(mx_mtlb_JACfctB);
    mxDestroyArray(mx_mtlb_PSETfctB);
    mxDestroyArray(mx_mtlb_PSOLfctB);
    mxDestroyArray(mx_mtlb_MONfctB);
    mxDestroyArray(mx_mtlb_MONdataB);
  }

  return(0);
}
