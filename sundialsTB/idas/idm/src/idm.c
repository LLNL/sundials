/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-07-19 22:10:53 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/idas/LICENSE.
 * -----------------------------------------------------------------
 * MEX implementation for IDAS Matlab interface.
 * -----------------------------------------------------------------
 */

#include <string.h>
#include <stdlib.h>
#include "idm.h"
#include "nvm.h"

/*
 * ---------------------------------------------------------------------------------
 * Definitions for global variables (declared in cvm.h)
 * ---------------------------------------------------------------------------------
 */

idm_IDASdata idm_Cdata;    /* IDA data */
booleantype idm_quad;      /* Forward quadratures? */
booleantype idm_quadB;     /* Backward quadratures? */
booleantype idm_asa;       /* Adjoint sensitivity? */
booleantype idm_fsa;       /* Forward sensitivity? */
booleantype idm_mon;       /* Forward monitoring? */ 
booleantype idm_monB;      /* Backward monitoring? */ 

idm_MATLABdata idm_Mdata;  /* MATLAB data */

/*
 * ---------------------------------------------------------------------------------
 * Static function prototypes
 * ---------------------------------------------------------------------------------
 */

static void IDM_init();
static void IDM_makePersistent();
static void IDM_final();

static int IDM_Malloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_SensMalloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_AdjMalloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_MallocB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_ReInit(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_ReInitB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_CalcIC(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_CalcICB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_Solve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_SolveB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_Stats(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_StatsB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_Get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_Set(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_Free(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/*
 * ---------------------------------------------------------------------------------
 * Main entry point
 * ---------------------------------------------------------------------------------
 */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  int mode;
  /* 
     Modes:
     
     1 - initialize IDAS solver
     2 - initialize forward sensitivity calculations
     3 - initialize adjoint sensitivity calculations
     4 - initialize backward solver
     5 - reinitialize IDAS solver (NYI)
     6 - reinitialize backward solver (NYI)
     7 - calculate consistent IC
     8 - calculate backward consistent IC
    10 - solve problem
    11 - solve backward problem
    20 - get integrator stats
    21 - get backward integrator stats
    22 - extract data from ida_mem
    23 - set one optional input at a time (NYI)
    30 - finalize
  */

  mode = (int)mxGetScalar(prhs[0]);

  mexUnlock();

  switch(mode) {
  case 1:
    IDM_init();
    IDM_Malloc(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 2:
    IDM_SensMalloc(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 3:
    IDM_AdjMalloc(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 4:
    IDM_MallocB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 5:
    IDM_ReInit(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 6:
    IDM_ReInitB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 7:
    IDM_CalcIC(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 8:
    IDM_CalcICB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 10:
    IDM_Solve(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 11:
    IDM_SolveB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 20:
    IDM_Stats(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 21:
    IDM_StatsB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 22:
    IDM_Get(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 23:
    IDM_Set(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
  case 30:
    IDM_Free(nlhs, plhs, nrhs-1, &prhs[1]);
    IDM_final();
    return;
  }

  IDM_makePersistent();
  mexLock();

  return;
}

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */

#define ida_mem     (idm_Cdata->ida_mem)
#define bp_data     (idm_Cdata->bp_data) 
#define yy          (idm_Cdata->yy) 
#define yp          (idm_Cdata->yp) 
#define yQ          (idm_Cdata->yQ) 
#define yyS         (idm_Cdata->yyS) 
#define ypS         (idm_Cdata->ypS) 
#define N           (idm_Cdata->N) 
#define Nq          (idm_Cdata->Nq) 
#define Ng          (idm_Cdata->Ng) 
#define Ns          (idm_Cdata->Ns) 
#define Nd          (idm_Cdata->Nd) 
#define Nc          (idm_Cdata->Nc) 
#define ls          (idm_Cdata->ls) 
#define pm          (idm_Cdata->pm) 
#define ism         (idm_Cdata->ism) 
#define idaadj_mem  (idm_Cdata->idaadj_mem) 
#define interp      (idm_Cdata->interp) 
#define yyB         (idm_Cdata->yyB) 
#define ypB         (idm_Cdata->ypB) 
#define yQB         (idm_Cdata->yQB) 
#define NB          (idm_Cdata->NB) 
#define NqB         (idm_Cdata->NqB) 
#define lsB         (idm_Cdata->lsB) 
#define pmB         (idm_Cdata->pmB) 

#define mx_data     (idm_Mdata->mx_data)

#define mx_RESfct   (idm_Mdata->mx_RESfct)
#define mx_QUADfct  (idm_Mdata->mx_QUADfct)
#define mx_JACfct   (idm_Mdata->mx_JACfct)
#define mx_PSETfct  (idm_Mdata->mx_PSETfct)
#define mx_PSOLfct  (idm_Mdata->mx_PSOLfct)
#define mx_GLOCfct  (idm_Mdata->mx_GLOCfct)
#define mx_GCOMfct  (idm_Mdata->mx_GCOMfct)
#define mx_Gfct     (idm_Mdata->mx_Gfct)
#define mx_SRESfct  (idm_Mdata->mx_SRESfct)

#define mx_RESfctB  (idm_Mdata->mx_RESfctB)
#define mx_QUADfctB (idm_Mdata->mx_QUADfctB)
#define mx_JACfctB  (idm_Mdata->mx_JACfctB)
#define mx_PSETfctB (idm_Mdata->mx_PSETfctB)
#define mx_PSOLfctB (idm_Mdata->mx_PSOLfctB)
#define mx_GLOCfctB (idm_Mdata->mx_GLOCfctB)
#define mx_GCOMfctB (idm_Mdata->mx_GCOMfctB)

#define mx_MONfct   (idm_Mdata->mx_MONfct)
#define mx_MONdata  (idm_Mdata->mx_MONdata)

#define mx_MONfctB  (idm_Mdata->mx_MONfctB)
#define mx_MONdataB (idm_Mdata->mx_MONdataB)

/*
 * ---------------------------------------------------------------------------------
 * Exported procedures
 * ---------------------------------------------------------------------------------
 */

/* IDM_Malloc
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

static int IDM_Malloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double t0, *yy0, *yp0;
  int i, status;

  int maxord;

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

  double hin, hmax;

  double tstop;
  booleantype tstopSet;

  int mupper, mlower;
  int gstype, maxl;
  int mudq, mldq;
  double dqrely;

  double *id, *cnstr;
  N_Vector NV_id, NV_cnstr;

  dqrely = 0.0;

  /* ------------------------------------
   * Initialize appropriate vector module
   * ------------------------------------
   */
  
  InitVectors();

  /* 
   * -----------------------------
   * Extract stuff from arguments:
   * - RES function
   * - initial time
   * - initial conditions
   * - integration options
   * - user data
   * -----------------------------
   */

  /* Matlab user-provided function */

  mx_RESfct = mxDuplicateArray(prhs[0]);
  
  /* Initial time */

  t0 = (double)mxGetScalar(prhs[1]);

  /* Initial conditions */

  yy0 = mxGetPr(prhs[2]);
  yp0 = mxGetPr(prhs[3]);
  N = mxGetM(prhs[2]);

  /* Integrator Options -- optional argument */

  status = get_IntgrOptions(prhs[4], TRUE,
                            &maxord, &mxsteps,
                            &itol, &reltol, &Sabstol, &Vabstol,
                            &hin, &hmax, &tstop, &tstopSet,
                            &id, &cnstr);

  /* User data -- optional argument */

  mx_data = mxDuplicateArray(prhs[5]);

  /* 
   * ---------------------------------------
   * Set initial conditions and tolerances
   * ---------------------------------------
   */

  yy = NewVector(N);
  PutData(yy, yy0, N);

  yp = NewVector(N);
  PutData(yp, yp0, N);

  switch (itol) {
  case IDA_SS:
    abstol = (void *) &Sabstol;
    break;
  case IDA_SV:
    NV_abstol = N_VClone(yy);
    PutData(NV_abstol, Vabstol, N);
    abstol = (void *) NV_abstol;
    break;
  }

  /* 
   * ----------------------------------------
   * Create IDAS object and allocate memory
   * ----------------------------------------
   */

  ida_mem = IDACreate();

  /* Attach error handler function */
  status = IDASetErrHandlerFn(ida_mem, mtlb_IdaErrHandler, NULL);

  /* Call IDAMalloc */
  status = IDAMalloc(ida_mem, mtlb_IdaRes, t0, yy, yp, itol, reltol, abstol);

  if (itol == IDA_SV) {
    N_VDestroy(NV_abstol);
  }

  /* Redirect output */
  status = IDASetErrFile(ida_mem, stdout);

  /* set maxorder (default is 5) */
  status = IDASetMaxOrd(ida_mem, maxord);

  /* set initial step size (the default value of 0.0 is ignored by IDAS) */
  status = IDASetInitStep(ida_mem, hin);

  /* set max step (default is infinity) */
  status = IDASetMaxStep(ida_mem, hmax);

  /* set number of max steps */
  status = IDASetMaxNumSteps(ida_mem, mxsteps);

  /* set tstop? */
  if (tstopSet)
    status = IDASetStopTime(ida_mem, tstop);
  
  /* Rootfinding? */
  if ( !mxIsEmpty(mx_Gfct) && (Ng > 0) ) {
    status = IDARootInit(ida_mem, Ng, mtlb_IdaGfct, NULL);
  }

  /* ID vector specified? */
  if (id != NULL) {
    NV_id = N_VClone(yy);
    PutData(NV_id, id, N);
    status = IDASetId(ida_mem, NV_id);
    N_VDestroy(NV_id);
  }

  /* Constraint vector specified? */
  if (cnstr != NULL) {
    NV_cnstr = N_VClone(yy);
    PutData(NV_cnstr, cnstr, N);
    status = IDASetConstraints(ida_mem, NV_cnstr);
    N_VDestroy(NV_cnstr);
  }

  /* Quadratures? */
  /*
  if ( idm_quad ) {

    status = get_QuadOptions(prhs[4], TRUE,
                             &yQ0, &errconQ, 
                             &itolQ, &reltolQ, &SabstolQ, &VabstolQ);

    if(status) mexErrMsgTxt("IDAMalloc:: illegal quadrature input.");      

    yQ = NewVector(Nq);
    PutData(yQ, yQ0, Nq);

    status = IDAQuadMalloc(ida_mem, mtlb_IdaQuadFct, yQ);

    if (errconQ) {
    
      switch (itolQ) {
      case IDA_SS:
        abstolQ = (void *) &SabstolQ;
        break;
      case IDA_SV:
        NV_abstolQ = N_VClone(yQ);
        PutData(NV_abstolQ, VabstolQ, Nq);
        abstolQ = (void *) NV_abstolQ;
        break;
      }
      
      status = IDASetQuadErrCon(ida_mem, errconQ, itolQ, reltolQ, abstolQ);
      
      if (itolQ == IDA_SV) N_VDestroy(NV_abstolQ);

    }

  }
  */

  /* Linear solver */

  status = get_LinSolvOptions(prhs[4], TRUE,
                              &mupper, &mlower,
                              &mudq, &mldq,
                              &gstype, &maxl);
  
  switch (ls) {
  case LS_DENSE:
    status = IDADense(ida_mem, N);
    if (!mxIsEmpty(mx_JACfct))
      status = IDADenseSetJacFn(ida_mem, mtlb_IdaDenseJac, NULL);
    break;
  case LS_BAND:
    status = IDABand(ida_mem, N, mupper, mlower);
    if (!mxIsEmpty(mx_JACfct))
      status = IDABandSetJacFn(ida_mem, mtlb_IdaBandJac, NULL);
    break;
  case LS_SPGMR:
    switch (pm) {
    case PM_NONE:
      status = IDASpgmr(ida_mem, maxl);
      if (!mxIsEmpty(mx_PSOLfct)) {
        if (!mxIsEmpty(mx_PSETfct))
          status = IDASpilsSetPreconditioner(ida_mem, mtlb_IdaSpilsPset, mtlb_IdaSpilsPsol, NULL);
        else
          status = IDASpilsSetPreconditioner(ida_mem, NULL, mtlb_IdaSpilsPsol, NULL);
      }
      break;
    case PM_BBDPRE:
      if (!mxIsEmpty(mx_GCOMfct))
        bp_data = IDABBDPrecAlloc(ida_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_IdaBBDgloc, mtlb_IdaBBDgcom);
      else
        bp_data = IDABBDPrecAlloc(ida_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_IdaBBDgloc, NULL);
      status = IDABBDSpgmr(ida_mem, maxl, bp_data);
      break;
    }
    status = IDASpilsSetGSType(ida_mem, gstype);
    if (!mxIsEmpty(mx_JACfct))
      status = IDASpilsSetJacTimesVecFn(ida_mem, mtlb_IdaSpilsJac, NULL);
    break;
  case LS_SPBCG:
    switch (pm) {
    case PM_NONE:
      status = IDASpbcg(ida_mem, maxl);
      if (!mxIsEmpty(mx_PSOLfct)) {
        if (!mxIsEmpty(mx_PSETfct))
          status = IDASpilsSetPreconditioner(ida_mem, mtlb_IdaSpilsPset, mtlb_IdaSpilsPsol, NULL);
        else
          status = IDASpilsSetPreconditioner(ida_mem, NULL, mtlb_IdaSpilsPsol, NULL);
      }
      break;
    case PM_BBDPRE:
      if (!mxIsEmpty(mx_GCOMfct))
        bp_data = IDABBDPrecAlloc(ida_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_IdaBBDgloc, mtlb_IdaBBDgcom);
      else
        bp_data = IDABBDPrecAlloc(ida_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_IdaBBDgloc, NULL);
      IDABBDSpbcg(ida_mem, maxl, bp_data);
      break;
    }
    if (!mxIsEmpty(mx_JACfct))
      status = IDASpilsSetJacTimesVecFn(ida_mem, mtlb_IdaSpilsJac, NULL);
    break;
  case LS_SPTFQMR:
    switch (pm) {
    case PM_NONE:
      status = IDASptfqmr(ida_mem, maxl);
      if (!mxIsEmpty(mx_PSOLfct)) {
        if (!mxIsEmpty(mx_PSETfct))
          status = IDASpilsSetPreconditioner(ida_mem, mtlb_IdaSpilsPset, mtlb_IdaSpilsPsol, NULL);
        else
          status = IDASpilsSetPreconditioner(ida_mem, NULL, mtlb_IdaSpilsPsol, NULL);
      }
      break;
    case PM_BBDPRE:
      if (!mxIsEmpty(mx_GCOMfct))
        bp_data = IDABBDPrecAlloc(ida_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_IdaBBDgloc, mtlb_IdaBBDgcom);
      else
        bp_data = IDABBDPrecAlloc(ida_mem, N, mudq, mldq, mupper, mlower, dqrely, mtlb_IdaBBDgloc, NULL);
      IDABBDSptfqmr(ida_mem, maxl, bp_data);
      break;
    }
    if (!mxIsEmpty(mx_JACfct))
      status = IDASpilsSetJacTimesVecFn(ida_mem, mtlb_IdaSpilsJac, NULL);
    break;
  }

  /* Do we monitor? */
  
  if (idm_mon) {
    mtlb_IdaMonitor(0, t0, NULL, NULL, NULL, NULL, NULL);
  }

  return(0);

}

/* IDM_SensMalloc
 *
 * prhs contains:
 *   Ns
 *   sensi_meth
 *   yS0
 *   options
 *
 * plhs contains:
 *   status
 *
 */

static int IDM_SensMalloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}
/*
static int IDM_SensMalloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double t0, *yyS0, *ypS0;
  int buflen, status;
  char *bufval;

  mxArray *pfield;
  char *pfield_name;

  booleantype userSRES, errconS;
  int itolS;
  realtype reltolS;
  realtype *SabstolS, *VabstolS;
  N_Vector *NV_abstolS;
  void *abstolS;

  int *plist;
  double *p, *pbar, rho;

  int is;

  p = NULL;
  plist = NULL;
  pbar = NULL;

  Ns = (int)mxGetScalar(prhs[0]);

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Simultaneous")) ism = IDA_SIMULTANEOUS;
  if(!strcmp(bufval,"Staggered"))    ism = IDA_STAGGERED;

  yyS0 = mxGetPr(prhs[2]);
  ypS0 = mxGetPr(prhs[3]);

  status = get_FSAOptions(prhs[4], 
                          &pfield_name, &plist, &pbar,
                          &userSRES, &rho,
                          &errconS, &itolS, &reltolS, &SabstolS, &VabstolS);

  if(status) mexErrMsgTxt("IDAMalloc:: illegal forward sensitivity input.");   

  if (mxIsEmpty(mx_SRESfct)) {
    if (pfield_name == NULL)
      mexErrMsgTxt("IDAMalloc:: pfield required but was not provided.");
    pfield = mxGetField(mx_data,0,pfield_name);
    if (pfield == NULL)
      mexErrMsgTxt("IDAMalloc:: illegal pfield input.");
    p = mxGetPr(pfield);
    mxFree(pfield_name);
  }

  yyS = N_VCloneVectorArray(Ns, yy);
  ypS = N_VCloneVectorArray(Ns, yy);
  for (is=0;is<Ns;is++) {
    PutData(yyS[is], &yyS0[is*N], N);
    PutData(ypS[is], &ypS0[is*N], N);
  }

  status = IDASensMalloc(ida_mem, Ns, ism, yyS, ypS);

  switch (itolS) {
  case IDA_SS:
    abstolS = (void *) SabstolS;
    break;
  case IDA_SV:
    NV_abstolS = N_VCloneVectorArray(Ns, y);
    for (is=0;is<Ns;is++)
      PutData(NV_abstolS[is], &VabstolS[is*N], N);
    abstolS = (void *) NV_abstolS;
    break;
  case IDA_EE:
    abstolS = NULL;
    break;
  }
  
  status = IDASetSensTolerances(ida_mem, itolS, reltolS, abstolS);
  
  if (itolS == IDA_SS)      free(SabstolS);
  else if (itolS == IDA_SV) N_VDestroyVectorArray(NV_abstolS, Ns);
  
  status = IDASetSensParams(ida_mem, p, pbar, plist);
  
  if (plist != NULL) free(plist);
  if (pbar != NULL)  free(pbar);
  
  status = IDASetSensRho(ida_mem, rho);
  
  status = IDASetSensErrCon(ida_mem, errconS);

  if (userSRES) {
    status = IDASetSensResFn(ida_mem, mtlb_IdaSensRhs, NULL);
  }

  idm_fsa = TRUE;

  return(0);
}
*/

/* IDM_AdjMalloc
 *
 * prhs contains:
 *
 * plhs contains:
 *   status
 */

static int IDM_AdjMalloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}
/*
static int IDM_AdjMalloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int buflen, status;
  char *bufval;

  Nd = (int)mxGetScalar(prhs[0]);

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Hermite"))    interp = IDA_HERMITE;
  if(!strcmp(bufval,"Polynomial")) interp = IDA_POLYNOMIAL;

  idaadj_mem = IDAadjMalloc(ida_mem, Nd, interp);

  idm_asa = TRUE;
  
  return(0);
}
*/
static int IDM_MallocB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}
/*
static int IDM_MallocB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double tB0, *yyB0, *ypB0;
  int i, status;

  int maxordB;

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

  double hinB, hmaxB;

  double tstopB;
  booleantype tstopSetB;

  int mupperB, mlowerB;
  int gstypeB, maxlB;
  int mudqB, mldqB;
  double dqrelyB;

  double *idB, *cnstrB;

  dqrely = 0.0;

  if (idm_mon) {
    mtlb_IdaMonitor(2, 0.0, NULL, NULL, NULL);
    idm_mon = FALSE;
  }

  mx_RESfctB = mxDuplicateArray(prhs[0]);
  
  tB0 = (double)mxGetScalar(prhs[1]);

  yyB0 = mxGetPr(prhs[2]);
  ypB0 = mxGetPr(prhs[3]);
  NB = mxGetM(prhs[2]);

  status = get_IntgrOptions(prhs[4], FALSE,
                            &maxordB, &mxstepsB,
                            &itolB, &reltolB, &SabstolB, &VabstolB,
                            &hinB, &hmaxB, &tstopB, &tstopSetB,
                            &idB, &cnstrB);

  yyB = NewVector(NB);
  PutData(yyB, yyB0, NB);

  ypB = NewVector(NB);
  PutData(ypB, ypB0, NB);

  switch (itolB) {
  case IDA_SS:
    abstolB = (void *) &SabstolB;
    break;
  case IDA_SV:
    NV_abstolB = N_VClone(yB);
    PutData(NV_abstolB, VabstolB, NB);
    abstolB = (void *) VabstolB;
    break;
  }

  status = IDACreateB(idaadj_mem);

  status = IDAMallocB(idaadj_mem, mtlb_IdaResB, tB0, yyB, ypB, itolB, reltolB, abstolB);

  if (itolB == IDA_SV) {
    N_VDestroy(NV_abstolB);
  }

  status = IDASetErrFileB(idaadj_mem, stdout);

  status = IDASetMaxOrdB(idaadj_mem, maxordB);

  status = IDASetInitStepB(idaadj_mem, hinB);

  status = IDASetMaxStepB(idaadj_mem, hmaxB);

  status = IDASetMaxNumStepsB(idaadj_mem, mxstepsB);

  if ( idm_quadB ) {
    
    status = get_QuadOptions(prhs[4], FALSE,
                             &yQB0, &errconQB, 
                             &itolQB, &reltolQB, &SabstolQB, &VabstolQB);

    if(status)
      mexErrMsgTxt("IDAMallocB:: illegal quadrature input.");      

    yQB = NewVector(NqB);
    PutData(yQB, yQB0, NqB);

    status = IDAQuadMallocB(idaadj_mem, mtlb_IdaQuadFctB, yQB);
    
    switch (itolQB) {
    case IDA_SS:
      abstolQB = (void *) &SabstolQB;
      break;
    case IDA_SV:
      NV_abstolQB = N_VClone(yQB);
      PutData(NV_abstolQB, VabstolQB, NqB);
      abstolQB = (void *) NV_abstolQB;
      break;
    }

    status = IDASetQuadErrConB(idaadj_mem, errconQB, itolQB, reltolQB, abstolQB);

    if (itolQB == IDA_SV) {
      N_VDestroy(NV_abstolQB);
    }

  }

  status = get_LinSolvOptions(prhs[4], FALSE,
                              &mupperB, &mlowerB,
                              &mudqB, &mldqB,
                              &gstypeB, &maxlB);
  
  switch(lsB) {
  case LS_DENSE:
    status = IDADenseB(idaadj_mem, NB);
    if (!mxIsEmpty(mx_JACfctB))
      status = IDADenseSetJacFnB(idaadj_mem, mtlb_IdaDenseJacB, NULL);
    break;
  case LS_BAND:
    status = IDABandB(idaadj_mem, NB, mupperB, mlowerB);
    if (!mxIsEmpty(mx_JACfctB))
      status = IDABandSetJacFnB(idaadj_mem, mtlb_IdaBandJacB, NULL);
    break;
  case LS_SPGMR:
    switch (pmB) {
    case PM_NONE:
      status = IDASpgmrB(idaadj_mem, maxlB);
      if (!mxIsEmpty(mx_PSOLfctB)) {
        if (!mxIsEmpty(mx_PSETfctB))
          status = IDASpilsSetPreconditionerB(idaadj_mem, mtlb_IdaSpilsPsetB, mtlb_IdaSpilsPsolB, NULL);
        else
          status = IDASpilsSetPreconditionerB(idaadj_mem, NULL, mtlb_IdaSpilsPsolB, NULL);
      }
      break;
    case PM_BBDPRE:
      if (!mxIsEmpty(mx_GCOMfctB)) {
        status = IDABBDPrecAllocB(idaadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_IdaBBDglocB, mtlb_IdaBBDgcomB);
      } else {
        status = IDABBDPrecAllocB(idaadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_IdaBBDglocB, NULL);
      }
      IDABBDSpgmrB(idaadj_mem, maxlB);
      break;
    }
    status = IDASpilsSetGSTypeB(idaadj_mem, gstypeB);
    if (!mxIsEmpty(mx_JACfctB))
      status = IDASpilsSetJacTimesVecFnB(idaadj_mem, mtlb_IdaSpilsJacB, NULL);
    break;
  case LS_SPBCG:
    switch (pmB) {
    case PM_NONE:
      status = IDASpbcgB(idaadj_mem, maxlB);
      if (!mxIsEmpty(mx_PSOLfctB)) {
        if (!mxIsEmpty(mx_PSETfctB))
          status = IDASpilsSetPreconditionerB(idaadj_mem, mtlb_IdaSpilsPsetB, mtlb_IdaSpilsPsolB, NULL);
        else
          status = IDASpilsSetPreconditionerB(idaadj_mem, NULL, mtlb_IdaSpilsPsolB, NULL);
      }
      break;
    case PM_BBDPRE:
      if (!mxIsEmpty(mx_GCOMfctB)) {
        status = IDABBDPrecAllocB(idaadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_IdaBBDglocB, mtlb_IdaBBDgcomB);
      } else {
        status = IDABBDPrecAllocB(idaadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_IdaBBDglocB, NULL);
      }
      IDABBDSpbcgB(idaadj_mem, maxlB);
      break;
    }
    if (!mxIsEmpty(mx_JACfctB))
      status = IDASpilsSetJacTimesVecFnB(idaadj_mem, mtlb_IdaSpilsJacB, NULL);
    break;
  case LS_SPTFQMR:
    switch (pmB) {
    case PM_NONE:
      status = IDASptfqmrB(idaadj_mem, maxlB);
      if (!mxIsEmpty(mx_PSOLfctB)) {
        if (!mxIsEmpty(mx_PSETfctB))
          status = IDASpilsSetPreconditionerB(idaadj_mem, mtlb_IdaSpilsPsetB, mtlb_IdaSpilsPsolB, NULL);
        else
          status = IDASpilsSetPreconditionerB(idaadj_mem, NULL, mtlb_IdaSpilsPsolB, NULL);
      }
      break;
    case PM_BBDPRE:
      if (!mxIsEmpty(mx_GCOMfctB)) {
        status = IDABBDPrecAllocB(idaadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_IdaBBDglocB, mtlb_IdaBBDgcomB);
      } else {
        status = IDABBDPrecAllocB(idaadj_mem, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mtlb_IdaBBDglocB, NULL);
      }
      IDABBDSptfqmrB(idaadj_mem, maxlB);
      break;
    }
    if (!mxIsEmpty(mx_JACfctB))
      status = IDASpilsSetJacTimesVecFnB(idaadj_mem, mtlb_IdaSpilsJacB, NULL);
    break;
  }

  if (idm_monB) {
    mtlb_IdaMonitorB(0, tB0, NULL, NULL);
  }

  return(0);
}
*/
static int IDM_ReInit(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}

static int IDM_ReInitB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}

static int IDM_CalcIC(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double tout, t0;
  int buflen, status, icopt;
  char *bufval;

  /* Extract tout */
  tout = (double) mxGetScalar(prhs[0]);

  /* Extract icopt */
  icopt = -1;
  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"FindAlgebraic")) icopt = IDA_YA_YDP_INIT;
  else if(!strcmp(bufval,"FindAll"))  icopt = IDA_Y_INIT;
  
  /* Call IDACalcIC */
  status = IDACalcIC(ida_mem, icopt, tout);

  /* IDACalcIC return flag */
  plhs[0] = mxCreateScalarDouble((double)status);

  if (nlhs == 1) return(0);

  /* Extract and return corrected IC */
  IDAGetConsistentIC(ida_mem, yy, yp);
  plhs[1] = mxCreateDoubleMatrix(N,1,mxREAL);
  GetData(yy, mxGetPr(plhs[1]), N);
  plhs[2] = mxCreateDoubleMatrix(N,1,mxREAL);
  GetData(yp, mxGetPr(plhs[2]), N);

  return(0);
}

static int IDM_CalcICB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
}

static int IDM_Solve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double tout, tret, *yySdata, *ypSdata;
  int buflen, status, itask, i, is;
  char *bufval;

  long int nst;
  double t0;

  int itask1;
  booleantype iret;
  double h;

  /* Extract tout */
  tout = (double)mxGetScalar(prhs[0]);

  /* Extract itask */
  itask = -1;
  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Normal")) itask = IDA_NORMAL;
  else if(!strcmp(bufval,"OneStep")) itask = IDA_ONE_STEP;
  else if(!strcmp(bufval,"NormalTstop")) itask = IDA_NORMAL_TSTOP;
  else if(!strcmp(bufval,"OneStepTstop")) itask = IDA_ONE_STEP_TSTOP;

  /* Call IDA */

  if (!idm_mon) {

    status = IDASolve(ida_mem, tout, &tret, yy, yp, itask);
    /*
    if (!idm_asa)
      status = IDASolve(ida_mem, tout, &tret, yy, yp, itask);
    else
      status = IDASolveF(idaadj_mem, tout, &tret, yy, yp, itask, &Nc);
    */

  } else {

    if      (itask == IDA_NORMAL)         {iret = FALSE; itask1 = IDA_ONE_STEP;}
    else if (itask == IDA_ONE_STEP)       {iret = TRUE;  itask1 = IDA_ONE_STEP;}
    else if (itask == IDA_NORMAL_TSTOP)   {iret = FALSE; itask1 = IDA_ONE_STEP_TSTOP;}
    else if (itask == IDA_ONE_STEP_TSTOP) {iret = TRUE;  itask1 = IDA_ONE_STEP_TSTOP;}
    
    while(1) {
      
      status = IDASolve(ida_mem, tout, &tret, yy, yp, itask1);
      /*
      if (!idm_asa)
        status = IDASolve(ida_mem, tout, &tret, yy, yp, itask1);
      else
        status = IDASolveF(idaadj_mem, tout, &tret, yy, yp, itask1, &Nc);
      */

      /* break on IDA error */
      if (status < 0) break;   
      
      /* In NORMAL_MODE test if tout was reached */
      if (!iret) {
        IDAGetCurrentStep(ida_mem, &h);
        if ( (tret - tout)*h >= 0.0 ) {
          tret = tout;
          IDAGetSolution(ida_mem, tout, yy, yp);
          iret = TRUE;
        }
      }

      /* If root or tstop return, we'll need to break */
      if (status != IDA_SUCCESS) iret = TRUE; 

      /*
      if (idm_quad)
        status = IDAGetQuad(ida_mem, tret, yQ);
      */

      /*
      if (idm_fsa)
        status = IDAGetSens(ida_mem, tret, yyS, ypS);
      */

      /* Call the monitoring function */
      mtlb_IdaMonitor(1, tret, yy, yp, yQ, yyS, ypS);

      /* break if we need to */
      if(iret)  break;
      
    };

  }

  /* IDA return flag */
  plhs[0] = mxCreateScalarDouble((double)status);

  /* Return time */
  plhs[1] = mxCreateScalarDouble(tret);

  /* Solution vectors */
  plhs[2] = mxCreateDoubleMatrix(N,1,mxREAL);
  GetData(yy, mxGetPr(plhs[2]), N);
  plhs[3] = mxCreateDoubleMatrix(N,1,mxREAL);
  GetData(yp, mxGetPr(plhs[3]), N);

  if (nlhs == 4) return(0);

  /*
  if (nlhs == 5) {

    if (idm_quad) {
      plhs[4] = mxCreateDoubleMatrix(Nq,1,mxREAL);
      status = IDAGetQuad(ida_mem, tret, yQ);
      GetData(yQ, mxGetPr(plhs[4]), Nq);
      return(0);
    } else {
      mexErrMsgTxt("IDA:: incorrect number of outputs for IDASolve.");
      return(1);
    }

  }

  if (nlhs == 6) {

    if (idm_quad) {
      mexErrMsgTxt("IDA:: incorrect number of outputs for IDASolve.");
      return(1);
    } 
    if (idm_fsa) {
      plhs[4] = mxCreateDoubleMatrix(N,Ns,mxREAL);
      plhs[5] = mxCreateDoubleMatrix(N,Ns,mxREAL);
      yySdata = mxGetPr(plhs[4]);
      ypSdata = mxGetPr(plhs[5]);
      status = IDAGetSens(ida_mem, tret, yyS, ypS);
      for (is=0; is<Ns; is++) {
        GetData(yyS[is], &yySdata[is*N], N);
        GetData(ypS[is], &ypSdata[is*N], N);
      }
      return(0);
    } else {
      mexErrMsgTxt("IDA:: incorrect number of outputs for IDASolve.");
      return(1);
    }

  }

  if (nlhs == 7) {

    if ( (!idm_quad) || (!idm_fsa) ) {
      mexErrMsgTxt("IDA:: incorrect number of outputs for IDASolve.");
      return(1);
    }

    plhs[4] = mxCreateDoubleMatrix(Nq,1,mxREAL);
    status = IDAGetQuad(ida_mem, tret, yQ);
    GetData(yQ, mxGetPr(plhs[4]), Nq);

    plhs[5] = mxCreateDoubleMatrix(N,Ns,mxREAL);
    plhs[6] = mxCreateDoubleMatrix(N,Ns,mxREAL);
    yySdata = mxGetPr(plhs[5]);
    ypSdata = mxGetPr(plhs[6]);
    status = IDAGetSens(ida_mem, tret, yyS, ypS);
    for (is=0; is<Ns; is++) {
      GetData(yyS[is], &yySdata[is*N], N);
      GetData(ypS[is], &ypSdata[is*N], N);
    }
    return(0);
  }
  */

  mexErrMsgTxt("IDA:: incorrect number of outputs for IDASolve.");
  return(1);

}


static int IDM_SolveB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}
/*
static int IDM_SolveB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double toutB, tretB;
  int itaskB;
  double *tmp;
  int buflen, status, i;
  char *bufval;

  void *ida_memB;
  booleantype iretB;
  double hB;


  toutB = (double)mxGetScalar(prhs[0]);

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Normal"))       itaskB = IDA_NORMAL;
  else if(!strcmp(bufval,"OneStep")) itaskB = IDA_ONE_STEP;
  else status = -1;

  if (!idm_monB) {

    status = IDAsolveB(idaadj_mem, toutB, &tretB, yyB, ypB, itaskB);

  } else {
    
    if      (itaskB == IDA_NORMAL)   iretB = FALSE;
    else if (itaskB == IDA_ONE_STEP) iretB = TRUE;

    while(1) {

      status = IDASolveB(idaadj_mem, toutB, &tretB, yyB, ypB, IDA_ONE_STEP);

      if (status < 0) break;   
      
      if (!iretB) {
        ida_memB = IDAadjGetIDABmem(idaadj_mem);
        IDAGetCurrentStep(ida_memB, &hB);
        if ( (tretB - toutB)*hB >= 0.0 ) {
          tretB = toutB;
          IDAGetSolution(ida_memB, toutB, yyB, ypB);
          iretB = TRUE;
        }
      }

      if (idm_quadB)
        status = IDAGetQuadB(idaadj_mem, tretB, yQB);

      mtlb_IdaMonitorB(1, tretB, yyB, ypB, yQB);

      if(iretB)  break;

    };

  }

  plhs[0] = mxCreateScalarDouble((double)status);

  plhs[1] = mxCreateScalarDouble(tretB);

  plhs[2] = mxCreateDoubleMatrix(NB,1,mxREAL);
  GetData(yyB, mxGetPr(plhs[2]), NB);
  plhs[3] = mxCreateDoubleMatrix(NB,1,mxREAL);
  GetData(ypB, mxGetPr(plhs[3]), NB);

  if (nlhs == 4) return(0);

  if (nlhs == 5) {
    if (idm_quadB) {
      plhs[4] = mxCreateDoubleMatrix(NqB,1,mxREAL);
      status = IDAGetQuadB(idaadj_mem, yQB);
      GetData(yQB, mxGetPr(plhs[4]), NqB);
      return;
    } else {
      mexErrMsgTxt("IDA:: incorrect number of outputs for IDASolveB.");
      return;
    }
  }

  mexErrMsgTxt("IDA:: incorrect number of outputs for IDASolveB.");
  return(1);
}

*/
static int IDM_Stats(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const char *fnames_intgr[]={
    "nst",
    "nre",
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
    "nreD"
  };
  const char *fnames_band[]={
    "name",
    "njeB",
    "nreB"
  };
  const char *fnames_spils[]={
    "name",
    "nli",
    "npe",
    "nps",
    "ncfl",
    "njeSG",
    "nreSG"
  };
  const char *fnames_quad[]={
    "nfQe",
    "netfQ"
  };
  const char *fnames_sens[]={
    "nrSe",
    "nreS",
    "nsetupsS",
    "netfS",
    "nniS",
    "ncfnS",
  };

  long int nst, nre, nsetups, nni, ncfn, netf, nge;
  int qlast, qcur;
  double h0used, hlast, hcur, tcur;
  int *rootsfound;

  long int njeD, nreD;
  long int njeB, nreB;
  long int nli, npe, nps, ncfl, njeSG, nreSG;

  long int nfQe, netfQ;

  long int nrSe, nreS, netfS, nsetupsS;
  long int nniS, ncfnS;
  long int *nniSTGR1, *ncfnSTGR1;

  int i, flag;
  mxArray *mx_root, *mx_ls, *mx_quad, *mx_fsa, *mx_asa;
  mxArray *mx_rootsfound;
  mxArray *mx_nniSTGR1, *mx_ncfnSTGR1;
  double *tmp;
  int nfields;

  flag = IDAGetIntegratorStats(ida_mem, &nst, &nre, &nsetups, 
                               &netf, &qlast, &qcur, &h0used, &hlast, &hcur, &tcur);

  flag = IDAGetNonlinSolvStats(ida_mem, &nni, &ncfn);

  nfields = sizeof(fnames_intgr)/sizeof(*fnames_intgr);
  plhs[0] = mxCreateStructMatrix(1, 1, nfields, fnames_intgr);
 
  mxSetField(plhs[0], 0, "nst",     mxCreateScalarDouble((double)nst));
  mxSetField(plhs[0], 0, "nre",     mxCreateScalarDouble((double)nre));
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

    flag = IDAGetNumGEvals(ida_mem, &nge);

    nfields = sizeof(fnames_root)/sizeof(*fnames_root);
    mx_root = mxCreateStructMatrix(1, 1, nfields, fnames_root);

    mxSetField(mx_root, 0, "nge", mxCreateScalarDouble((double)nge));

    rootsfound = (int *) malloc(nge*sizeof(int));
    flag = IDAGetRootInfo(ida_mem, rootsfound);
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

  /*
  if (idm_quad) {

    flag = IDAGetQuadStats(ida_mem, &nfQe, &netfQ);

    nfields = sizeof(fnames_quad)/sizeof(*fnames_quad);
    mx_quad = mxCreateStructMatrix(1, 1, nfields, fnames_quad);

    mxSetField(mx_quad, 0, "nfQe",  mxCreateScalarDouble((double)nfQe));
    mxSetField(mx_quad, 0, "netfQ", mxCreateScalarDouble((double)netfQ));

  } else {

    mx_quad = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "QuadInfo", mx_quad);
  */

  /* Linear Solver Statistics */

  switch(ls){

  case LS_DENSE:
    
    flag = IDADenseGetNumJacEvals(ida_mem, &njeD);
    flag = IDADenseGetNumResEvals(ida_mem, &nreD);
    
    nfields = sizeof(fnames_dense)/sizeof(*fnames_dense);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_dense);
    
    mxSetField(mx_ls, 0, "name", mxCreateString("Dense"));
    mxSetField(mx_ls, 0, "njeD", mxCreateScalarDouble((double)njeD));
    mxSetField(mx_ls, 0, "nreD", mxCreateScalarDouble((double)nreD));
    
    break;

  case LS_BAND:
      
    flag = IDABandGetNumJacEvals(ida_mem, &njeB);
    flag = IDABandGetNumResEvals(ida_mem, &nreB);
      
    nfields = sizeof(fnames_band)/sizeof(*fnames_band);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_band);
 
    mxSetField(mx_ls, 0, "name", mxCreateString("Band"));
    mxSetField(mx_ls, 0, "njeB", mxCreateScalarDouble((double)njeB));
    mxSetField(mx_ls, 0, "nreB", mxCreateScalarDouble((double)nreB));
      
    break;

  case LS_SPGMR:
  case LS_SPBCG:
  case LS_SPTFQMR:

    flag = IDASpilsGetNumLinIters(ida_mem, &nli);
    flag = IDASpilsGetNumPrecEvals(ida_mem, &npe);
    flag = IDASpilsGetNumPrecSolves(ida_mem, &nps);
    flag = IDASpilsGetNumConvFails(ida_mem, &ncfl);
    flag = IDASpilsGetNumJtimesEvals(ida_mem, &njeSG);
    flag = IDASpilsGetNumResEvals(ida_mem, &nreSG);
    
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
    mxSetField(mx_ls, 0, "nreSG", mxCreateScalarDouble((double)nreSG));
    
    break;
    
  }

  mxSetField(plhs[0], 0, "LSInfo", mx_ls);

  /* Forward Sensitivity Statistics */

  /*
  if (idm_fsa) {

    flag = IDAGetSensStats(ida_mem, &nrSe, &nreS, &netfS, &nsetupsS); 

    flag = IDAGetSensNonlinSolvStats(ida_mem, &nniS, &ncfnS);

    nfields = sizeof(fnames_sens)/sizeof(*fnames_sens);
    mx_fsa = mxCreateStructMatrix(1, 1, nfields, fnames_sens);
    
    mxSetField(mx_fsa, 0, "nrSe",     mxCreateScalarDouble((double)nrSe));
    mxSetField(mx_fsa, 0, "nreS",     mxCreateScalarDouble((double)nreS));
    mxSetField(mx_fsa, 0, "nsetupsS", mxCreateScalarDouble((double)nsetupsS));
    mxSetField(mx_fsa, 0, "netfS",    mxCreateScalarDouble((double)netfS));
    mxSetField(mx_fsa, 0, "nniS",     mxCreateScalarDouble((double)nniS));
    mxSetField(mx_fsa, 0, "ncfnS",    mxCreateScalarDouble((double)ncfnS));
    
  } else {

    mx_fsa = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "FSAInfo", mx_fsa);
  */

  return(0);
}

static int IDM_StatsB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}
/*
static int IDM_StatsB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const char *fnames_intgr[]={
    "nst",
    "nre",
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
    "nreD"
  };
  const char *fnames_band[]={
    "name",
    "njeB",
    "nreB"
  };
  const char *fnames_spils[]={
    "name",
    "nli",
    "npe",
    "nps",
    "ncfl",
    "njeSG",
    "nreSG"
  };
  const char *fnames_quad[]={
    "nfQe",
    "netfQ"
  };

  void *ida_memB;

  long int nst, nre, nsetups, nni, ncfn, netf;
  int qlast, qcur;
  double h0used, hlast, hcur, tcur;

  long int njeD, nreD;
  long int njeB, nreB;
  long int nli, npe, nps, ncfl, njeSG, nreSG;

  long int nfQe, netfQ;

  int i, flag;
  mxArray *mx_ls, *mx_quad;
  double *tmp;
  int nfields;

  ida_memB = IDAadjGetIDABmem(idaadj_mem);

  flag = IDAGetIntegratorStats(ida_memB, &nst, &nre, &nsetups, 
                               &netf, &qlast, &qcur, &h0used, &hlast, &hcur, &tcur);

  flag = IDAGetNonlinSolvStats(ida_memB, &nni, &ncfn);

  nfields = sizeof(fnames_intgr)/sizeof(*fnames_intgr);
  plhs[0] = mxCreateStructMatrix(1, 1, nfields, fnames_intgr);
 
  mxSetField(plhs[0], 0, "nst",     mxCreateScalarDouble((double)nst));
  mxSetField(plhs[0], 0, "nre",     mxCreateScalarDouble((double)nre));
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

  if (idm_quadB) {

    flag = IDAGetQuadStats(ida_memB, &nfQe, &netfQ);

    nfields = sizeof(fnames_quad)/sizeof(*fnames_quad);
    mx_quad = mxCreateStructMatrix(1, 1, nfields, fnames_quad);

    mxSetField(mx_quad, 0, "nfQe",  mxCreateScalarDouble((double)nfQe));
    mxSetField(mx_quad, 0, "netfQ", mxCreateScalarDouble((double)netfQ));

  } else {

    mx_quad = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "QuadInfo", mx_quad);

  switch(lsB){

  case LS_DENSE:
    
    flag = IDADenseGetNumJacEvals(ida_memB, &njeD);
    flag = IDADenseGetNumResEvals(ida_memB, &nreD);
    
    nfields = sizeof(fnames_dense)/sizeof(*fnames_dense);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_dense);

    mxSetField(mx_ls, 0, "name", mxCreateString("Dense"));
    mxSetField(mx_ls, 0, "njeD", mxCreateScalarDouble((double)njeD));
    mxSetField(mx_ls, 0, "nreD", mxCreateScalarDouble((double)nreD));
    
    break;

  case LS_BAND:
      
    flag = IDABandGetNumJacEvals(ida_memB, &njeB);
    flag = IDABandGetNumResEvals(ida_memB, &nreB);
      
    nfields = sizeof(fnames_band)/sizeof(*fnames_band);
    mx_ls = mxCreateStructMatrix(1, 1, nfields, fnames_band);
 
    mxSetField(mx_ls, 0, "name", mxCreateString("Band"));
    mxSetField(mx_ls, 0, "njeB", mxCreateScalarDouble((double)njeB));
    mxSetField(mx_ls, 0, "nreB", mxCreateScalarDouble((double)nreB));
      
    break;

  case LS_SPGMR:
  case LS_SPBCG:
  case LS_SPTFQMR:

    flag = IDASpilsGetNumLinIters(ida_memB, &nli);
    flag = IDASpilsGetNumPrecEvals(ida_memB, &npe);
    flag = IDASpilsGetNumPrecSolves(ida_memB, &nps);
    flag = IDASpilsGetNumConvFails(ida_memB, &ncfl);
    flag = IDASpilsGetNumJtimesEvals(ida_memB, &njeSG);
    flag = IDASpilsGetNumResEvals(ida_memB, &nreSG);
    
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
    mxSetField(mx_ls, 0, "nreSG", mxCreateScalarDouble((double)nreSG));
    
    break;
  }

  mxSetField(plhs[0], 0, "LSInfo", mx_ls);

  return(0);
}
*/

static int IDM_Set(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}

static int IDM_Get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}
/*
static int IDM_Get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double t;
  N_Vector ewt;
  double *this, *next;
  int key, k, which, status, i, nfields;
  mxArray *mx_yy, *mx_yp;

  IDAadjCheckPointRec *ckpnt;
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
  case 1:

    break;

  case 2:

    ewt = N_VClone(yy);

    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    status = IDAGetErrWeights(ida_mem, ewt);
    GetData(ewt, mxGetPr(plhs[0]), N);

    N_VDestroy(ewt);

    break;

  case 3:

    break;

  case 4:

    ckpnt = (IDAadjCheckPointRec *) malloc ( (Nc+1)*sizeof(IDAadjCheckPointRec));
    IDAadjGetCheckPointsInfo(idaadj_mem, ckpnt);
    nfields = sizeof(fnames_ckpnt)/sizeof(*fnames_ckpnt);
    plhs[0] = mxCreateStructMatrix(Nc+1, 1, nfields, fnames_ckpnt);
    for (i=0; i<=Nc; i++) {
      this = (double *)(ckpnt[Nc-i].my_addr);
      next = (double *)(ckpnt[Nc-i].next_addr);
      mxSetField(plhs[0], i, "addr",  mxCreateScalarDouble(*this));
      mxSetField(plhs[0], i, "next",  mxCreateScalarDouble(*next));
      mxSetField(plhs[0], i, "t0",    mxCreateScalarDouble((double)(ckpnt[Nc-i].t0)));
      mxSetField(plhs[0], i, "t1",    mxCreateScalarDouble((double)(ckpnt[Nc-i].t1)));
      mxSetField(plhs[0], i, "nstep", mxCreateScalarDouble((double)(ckpnt[Nc-i].nstep)));
      mxSetField(plhs[0], i, "order", mxCreateScalarDouble((double)(ckpnt[Nc-i].order)));
      mxSetField(plhs[0], i, "step",  mxCreateScalarDouble((double)(ckpnt[Nc-i].step)));
    }

    break;

  case 5:

    status = IDAadjGetCurrentCheckPoint(idaadj_mem, (void **)(&this));

    plhs[0] = mxCreateScalarDouble(*this);

    break;

  case 6:
    
    which = (int) (*mxGetPr(prhs[1]));

    if (interp == IDA_HERMITE) {
    
      status = IDAadjGetDataPointHermite(idaadj_mem, which, &t, y, yd);

      plhs[0] = mxCreateCellMatrix(1, 3);

      mxSetCell(plhs[0],1,mxCreateScalarDouble(t));

      mx_yy  = mxCreateDoubleMatrix(N,1,mxREAL);
      GetData(yy, mxGetPr(mx_yy), N);
      mxSetCell(plhs[0],2,mx_yy);

      mx_yp = mxCreateDoubleMatrix(N,1,mxREAL);
      GetData(yp, mxGetPr(mx_yp), N);
      mxSetCell(plhs[0],3,mx_yp);

    }

    break;

  }

  return(0);
}
*/

static int IDM_Free(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  if (idm_mon) {
    mtlb_IdaMonitor(2, 0.0, NULL, NULL, NULL, NULL, NULL);
  }

  N_VDestroy(yy);
  N_VDestroy(yp);
  /*
  if (idm_quad) N_VDestroy(yQ);
  */
  if (pm == PM_BBDPRE)  IDABBDPrecFree(&bp_data);

  /*
  if (idm_fsa) {
    N_VDestroyVectorArray(yyS, Ns);
    N_VDestroyVectorArray(ypS, Ns);
  }
  */

  IDAFree(&ida_mem);

  /*
  if (idm_asa) {

    if (idm_monB) {
      mtlb_IdaMonitorB(2, 0.0, NULL, NULL, NULL);
    }

    N_VDestroy(yyB);
    N_VDestroy(ypB);

    if (idm_quadB) N_VDestroy(yQB);
   
    IDAadjFree(&idaadj_mem);

  }
  */

  return(0);
}


static void IDM_init()
{
  mxArray *empty;

  /* Allocate space for global IDAS and MATLAB data structures */
  
  idm_Cdata = (idm_IDASdata) mxMalloc(sizeof(struct idm_IDASdataStruct));
  idm_Mdata = (idm_MATLABdata) mxMalloc(sizeof(struct idm_MATLABdataStruct));

  /* Initialize global IDAS data */

  ida_mem = NULL;
  idaadj_mem = NULL;
  bp_data   = NULL;

  yy  = NULL;
  yp  = NULL;
  yQ  = NULL;
  yyS = NULL;
  ypS = NULL;
  yyB = NULL;
  ypB = NULL;
  yQB = NULL;

  N   = 0;
  Nq  = 0;
  Ng  = 0;
  Ns  = 0;
  Nd  = 0;
  Nc  = 0;
  NB  = 0;
  NqB = 0;

  /*
  ism = IDA_STAGGERED;
  interp = IDA_POLYNOMIAL;
  */

  ls  = LS_DENSE;
  lsB = LS_DENSE;
  pm  = PM_NONE;
  pmB = PM_NONE;

  /* Initialize global control variables */

  idm_quad  = FALSE;
  idm_quadB = FALSE;
  idm_fsa   = FALSE;
  idm_asa   = FALSE;
  idm_mon   = FALSE;
  idm_monB  = FALSE;

  /* Initialize global MATLAB data */

  empty = mxCreateDoubleMatrix(0,0,mxREAL);

  mx_data     = mxDuplicateArray(empty);

  mx_RESfct   = mxDuplicateArray(empty);
  mx_Gfct     = mxDuplicateArray(empty);
  mx_QUADfct  = mxDuplicateArray(empty);
  mx_SRESfct  = mxDuplicateArray(empty);
  mx_JACfct   = mxDuplicateArray(empty);
  mx_PSETfct  = mxDuplicateArray(empty);
  mx_PSOLfct  = mxDuplicateArray(empty);
  mx_GLOCfct  = mxDuplicateArray(empty);
  mx_GCOMfct  = mxDuplicateArray(empty);

  mx_RESfctB  = mxDuplicateArray(empty);
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

static void IDM_makePersistent()
{
  /* Make global memory persistent */

  mexMakeArrayPersistent(mx_data);

  mexMakeArrayPersistent(mx_RESfct);
  mexMakeArrayPersistent(mx_Gfct);
  mexMakeArrayPersistent(mx_QUADfct);
  mexMakeArrayPersistent(mx_SRESfct);
  mexMakeArrayPersistent(mx_JACfct);
  mexMakeArrayPersistent(mx_PSETfct);
  mexMakeArrayPersistent(mx_PSOLfct);
  mexMakeArrayPersistent(mx_GLOCfct);
  mexMakeArrayPersistent(mx_GCOMfct);

  mexMakeArrayPersistent(mx_RESfctB);
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

  mexMakeMemoryPersistent(idm_Cdata);
  mexMakeMemoryPersistent(idm_Mdata);

}


static void IDM_final()
{  
  mxDestroyArray(mx_data);
  
  mxDestroyArray(mx_RESfct);
  mxDestroyArray(mx_Gfct);
  mxDestroyArray(mx_QUADfct);
  mxDestroyArray(mx_SRESfct);
  mxDestroyArray(mx_JACfct);
  mxDestroyArray(mx_PSETfct);
  mxDestroyArray(mx_PSOLfct);
  mxDestroyArray(mx_GLOCfct);
  mxDestroyArray(mx_GCOMfct);
  
  mxDestroyArray(mx_RESfctB);
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

  mxFree(idm_Cdata);
  mxFree(idm_Mdata);
}
