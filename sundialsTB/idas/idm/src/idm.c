/*
 * -----------------------------------------------------------------
 * $Revision: 1.14 $
 * $Date: 2007-08-21 17:38:45 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials-x.y.z/src/idas/LICENSE.
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
 * Global interface data variable
 * ---------------------------------------------------------------------------------
 */

idmInterfaceData idmData = NULL;

/*
 * ---------------------------------------------------------------------------------
 * Static function prototypes
 * ---------------------------------------------------------------------------------
 */

static void idmInitIDASdata();
static void idmPersistIDASdata();
static void idmFinalIDASdata();

static void idmInitPbData(idmPbData pb);
static void idmPersistPbData(idmPbData pb);
static void idmFinalPbData(idmPbData pb);


static int IDM_Initialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_QuadInitialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_SensInitialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_AdjInitialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

static int IDM_InitializationB(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_QuadInitializationB(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

static int IDM_SensToggleOff(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

static int IDM_CalcIC(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_CalcICB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

static int IDM_Solve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_SolveB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

static void idmSolveB_one(mxArray *plhs[], int NtoutB, double *toutB, int itaskB);
static void idmSolveB_more(mxArray *plhs[], int NtoutB, double *toutB, int itaskB,
                           booleantype any_quadrB, booleantype any_monB);

static int IDM_Stats(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_StatsB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

static int IDM_Set(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static int IDM_SetB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

static int IDM_Get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

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
     2 - initialize quadratures
     3 - initialize forward sensitivity calculations
     4 - initialize adjoint sensitivity calculations
     5 - initialize backward solver
     6 - initialize backward quadratures

    11 - reinitialize IDAS solver
    12 - reinitialize quadratures
    13 - reinitialize forward sensitivity calculations
    14 - reinitialize adjoint sensitivity calculations
    15 - reinitialize backward solver
    16 - reinitialize backward quadratures

    18 - toggle FSA off

    20 - solve problem
    21 - solve backward problem                TODO

    25 - calculate consistent IC
    26 - calculate backward consistent IC      TODO

    30 - get integrator stats
    31 - get backward integrator stats
    32 - extract data from ida_mem

    33 - set one optional input at a time
    34 - set one optional input at a time for backward problems

    40 - finalize
  */

  mode = (int)mxGetScalar(prhs[0]);

  mexUnlock();

  if ( (mode != 1) && (idmData == NULL) ) {
    idmErrHandler(-999, "IDAS", "-",
                  "Illegal attempt to call before IDAInit.", NULL);
  }


  switch(mode) {

    /* Initialization functions */

  case 1:
    if (idmData != NULL) {
      IDM_Free(nlhs, plhs, nrhs-1, &prhs[1]);
      idmFinalIDASdata();
    }
    idmInitIDASdata();
    IDM_Initialization(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 2:
    IDM_QuadInitialization(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 3:
    IDM_SensInitialization(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 4:
    IDM_AdjInitialization(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 5:
    IDM_InitializationB(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 6:
    IDM_QuadInitializationB(0, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

    /* Re-initialization functions */

  case 11:
    IDM_Initialization(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 12:
    IDM_QuadInitialization(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 13:
    IDM_SensInitialization(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 14:
    IDM_AdjInitialization(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 15:
    IDM_InitializationB(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 16:
    IDM_QuadInitializationB(1, nlhs, plhs, nrhs-1, &prhs[1]);
    break;

    /* Sensitivity toggle function */

  case 18:
    IDM_SensToggleOff(nlhs, plhs, nrhs-1, &prhs[1]);
    break;
    
    /* Solve functions */

  case 20:
    IDM_Solve(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 21:
    IDM_SolveB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

    /* Consistent IC calculation functions */

  case 25:
    IDM_CalcIC(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 26:
    IDM_CalcICB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

    /* Optional output extraction functions */

  case 30:
    IDM_Stats(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 31:
    IDM_StatsB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 32:
    IDM_Get(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 33:
    IDM_Set(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

  case 34:
    IDM_SetB(nlhs, plhs, nrhs-1, &prhs[1]);
    break;

    /* Memory deallocation function */

  case 40:
    IDM_Free(nlhs, plhs, nrhs-1, &prhs[1]);
    idmFinalIDASdata();
    return;

  }

  /* Unless this was the IDAFree call,
   * make data persistent and lock the MEX file */
  if (mode != 40) {
    idmPersistIDASdata();
    mexLock();
  }

  return;
}


/*
 * ---------------------------------------------------------------------------------
 * Private functions
 * ---------------------------------------------------------------------------------
 */


static void idmInitIDASdata()
{
  /* Allocate space for global IDAS data structure */

  idmData = (idmInterfaceData) mxMalloc(sizeof(struct idmInterfaceData_));

  /* Initialize global IDAS data */

  idmData->ida_mem = NULL;

  idmData->fwdPb     = NULL;
  idmData->bckPb     = NULL;

  idmData->NbckPb    = 0;

  idmData->Nd        = 0;
  idmData->Nc        = 0;
  idmData->asa       = FALSE;

  return;
}


static void idmInitPbData(idmPbData pb)
{
  mxArray *empty;

  pb->n  = 0;
  pb->nq = 0;
  pb->ng = 0;
  pb->ns = 0;

  pb->YY = NULL;
  pb->YP = NULL;

  pb->YQ = NULL;

  pb->YYS = NULL;
  pb->YPS = NULL;

  pb->Quadr = FALSE;
  pb->Fsa   = FALSE;
  pb->Mon   = FALSE;

  pb->LS = LS_DENSE;
  pb->PM = PM_NONE;

  empty = mxCreateDoubleMatrix(0,0,mxREAL);

  pb->RESfct   = mxDuplicateArray(empty);
  pb->Gfct     = mxDuplicateArray(empty);
  pb->QUADfct  = mxDuplicateArray(empty);
  pb->SRESfct  = mxDuplicateArray(empty);
  pb->JACfct   = mxDuplicateArray(empty);
  pb->PSETfct  = mxDuplicateArray(empty);
  pb->PSOLfct  = mxDuplicateArray(empty);
  pb->GLOCfct  = mxDuplicateArray(empty);
  pb->GCOMfct  = mxDuplicateArray(empty);

  pb->MONfct   = mxDuplicateArray(empty);
  pb->MONdata  = mxDuplicateArray(empty);

  pb->mtlb_data = mxDuplicateArray(empty);

  pb->fwd = idmData->fwdPb;

  pb->index = 0;
  pb->next  = NULL;

  mxDestroyArray(empty);
}


static void idmPersistIDASdata()
{
  idmPbData tmpPb;

  /* Make global memory persistent */

  if (idmData->fwdPb != NULL) {
    idmPersistPbData(idmData->fwdPb);
    mexMakeMemoryPersistent(idmData->fwdPb);
  }

  tmpPb = idmData->bckPb;
  while(tmpPb != NULL) {
    idmPersistPbData(tmpPb);
    mexMakeMemoryPersistent(tmpPb);
    tmpPb = tmpPb->next;
  }
  
  mexMakeMemoryPersistent(idmData);

  return;
}


static void idmPersistPbData(idmPbData pb)
{
  mexMakeArrayPersistent(pb->mtlb_data);

  mexMakeArrayPersistent(pb->RESfct);
  mexMakeArrayPersistent(pb->Gfct);
  mexMakeArrayPersistent(pb->QUADfct);
  mexMakeArrayPersistent(pb->SRESfct);
  mexMakeArrayPersistent(pb->JACfct);
  mexMakeArrayPersistent(pb->PSETfct);
  mexMakeArrayPersistent(pb->PSOLfct);
  mexMakeArrayPersistent(pb->GLOCfct);
  mexMakeArrayPersistent(pb->GCOMfct);

  mexMakeArrayPersistent(pb->MONfct);
  mexMakeArrayPersistent(pb->MONdata);
}

static void idmFinalIDASdata()
{  
  idmPbData tmpPb;

  if (idmData == NULL) return;

  if (idmData->fwdPb != NULL) {
    idmFinalPbData(idmData->fwdPb);
    mxFree(idmData->fwdPb);
    idmData->fwdPb = NULL;
  }

  while(idmData->bckPb != NULL) {
    tmpPb = idmData->bckPb->next;
    mxFree(idmData->bckPb);
    idmData->bckPb = tmpPb;
  }

  mxFree(idmData);
  idmData = NULL;

  return;
}


static void idmFinalPbData(idmPbData pb)
{

  if (pb->YY != NULL) N_VDestroy(pb->YY);
  if (pb->YP != NULL) N_VDestroy(pb->YP);

  if (pb->YQ != NULL) N_VDestroy(pb->YQ);

  if (pb->YYS != NULL) N_VDestroyVectorArray(pb->YYS, pb->ns);
  if (pb->YPS != NULL) N_VDestroyVectorArray(pb->YPS, pb->ns);

  mxDestroyArray(pb->mtlb_data);

  mxDestroyArray(pb->RESfct);
  mxDestroyArray(pb->Gfct);
  mxDestroyArray(pb->QUADfct);
  mxDestroyArray(pb->SRESfct);
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
 * This function is both passed as the IDAS error handler and used throughout
 * the Matlab interface.
 *
 * If called directly by one of the interface functions, error_code = -999 to
 * indicate an error and err_code = +999 to indicate a warning. Otherwise,
 * err_code is set by the calling IDAS function.
 *
 * NOTE: mexErrMsgTxt will end the execution of the MEX file. Therefore we do
 *       not have to intercept any of the IDAS error return flags.
 *       The only return flags we intercept are those from IDASolve() and IDASolveB()
 *       which are passed back to the user (only positive values will make it). 
 * ---------------------------------------------------------------------------------
 */

void idmErrHandler(int error_code, 
                   const char *module, const char *function, 
                   char *msg, void *f_data)
{
  char err_msg[256];

  if (error_code > 0) {
    sprintf(err_msg,"Warning in ==> %s\n%s",function,msg);
    mexWarnMsgTxt(err_msg);    
  } else if (error_code < 0) {
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

#define ida_mem     (idmData->ida_mem)

#define asa         (idmData->asa)
#define Nd          (idmData->Nd) 
#define Nc          (idmData->Nc) 
#define NbckPb      (idmData->NbckPb)

#define fsa         (fwdPb->Fsa)
#define quadr       (fwdPb->Quadr)
#define mon         (fwdPb->Mon)
#define rootSet     (fwdPb->RootSet)
#define tstopSet    (fwdPb->TstopSet)

#define yy          (fwdPb->YY) 
#define yp          (fwdPb->YP) 
#define yQ          (fwdPb->YQ) 
#define yyS         (fwdPb->YYS) 
#define ypS         (fwdPb->YPS) 
#define N           (fwdPb->n) 
#define Nq          (fwdPb->nq) 
#define Ng          (fwdPb->ng) 
#define Ns          (fwdPb->ns) 
#define ls          (fwdPb->LS) 
#define pm          (fwdPb->PM)

#define mtlb_data     (fwdPb->mtlb_data)

#define mtlb_RESfct   (fwdPb->RESfct)
#define mtlb_QUADfct  (fwdPb->QUADfct)
#define mtlb_JACfct   (fwdPb->JACfct)
#define mtlb_PSETfct  (fwdPb->PSETfct)
#define mtlb_PSOLfct  (fwdPb->PSOLfct)
#define mtlb_GLOCfct  (fwdPb->GLOCfct)
#define mtlb_GCOMfct  (fwdPb->GCOMfct)
#define mtlb_Gfct     (fwdPb->Gfct)
#define mtlb_SRESfct  (fwdPb->SRESfct)

#define mtlb_MONfct   (fwdPb->MONfct)
#define mtlb_MONdata  (fwdPb->MONdata)


#define indexB      (bckPb->index)

#define quadrB      (bckPb->Quadr)
#define monB        (bckPb->Mon)

#define yyB         (bckPb->YY) 
#define ypB         (bckPb->YP) 
#define yQB         (bckPb->YQ) 
#define NB          (bckPb->n) 
#define NqB         (bckPb->nq) 
#define lsB         (bckPb->LS) 
#define pmB         (bckPb->PM) 

#define mtlb_dataB    (bckPb->mtlb_data)

#define mtlb_RESfctB  (bckPb->RESfct)
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

/* IDM_Initialization
 *
 * action = 0   -> IDACreate + IDAMalloc
 * action = 1   -> IDAReInit
 *
 * prhs contains:
 *   res
 *   t0
 *   yy0
 *   yp0
 *   options
 *   data
 *
 * plhs contains:
 *   status
 *
 */

static int IDM_Initialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  idmPbData fwdPb;

  const mxArray *options;

  double t0, *yy0, *yp0;

  int maxord;
  long int mxsteps;

  int itol;
  realtype reltol, Sabstol, *Vabstol;
  N_Vector NV_abstol;

  double hin, hmax;
  double tstop;

  booleantype suppress;

  int mupper, mlower;
  int gstype, maxl;
  int mudq, mldq;
  double dqrely;

  double *id, *cnstr;
  N_Vector NV_id, NV_cnstr;

  booleantype res_s; /* ignored */

  /* 
   * ------------------------------------
   * Process inputs based on action
   * ------------------------------------
   */

  switch (action) {

  case 0:     /* SOLVER INITIALIZATION */

    /* Create and initialize a new problem */

    fwdPb = (idmPbData) mxMalloc(sizeof(struct idmPbData_));
    idmInitPbData(fwdPb);

    idmData->fwdPb = fwdPb;

    /* Initialize appropriate vector module */

    InitVectors();

    /* Extract user-provided RES function */

    mxDestroyArray(mtlb_RESfct);
    mtlb_RESfct = mxDuplicateArray(prhs[0]);

    /* Extract initial time */

    t0 = (double)mxGetScalar(prhs[1]);

    /* Extract initial conditions */

    yy0 = mxGetPr(prhs[2]);
    yp0 = mxGetPr(prhs[3]);
    N = mxGetM(prhs[2]);

    /* Create the solution N_Vectors */

    yy = NewVector(N);
    yp = NewVector(N);

    /* Load initial conditions */

    PutData(yy, yy0, N);
    PutData(yp, yp0, N);

    /* Extract options structure */
    
    options = prhs[4];

    break;

  case 1:    /* SOLVER RE-INITIALIZATION */

    fwdPb = idmData->fwdPb;

    /* If monitoring was enabled, finalize it now. */

    if (mon) mxW_IDAMonitor(2, 0.0, NULL, NULL, NULL, fwdPb);

    /* Extract initial time */

    t0 = (double)mxGetScalar(prhs[0]);

    /* Extract initial conditions */

    yy0 = mxGetPr(prhs[1]);

    if (mxGetM(prhs[1]) != N) idmErrHandler(-999, "IDAS", "IDAReInit",
                                            "Size of yy0 changed from IDAInit call.", NULL);
    yp0 = mxGetPr(prhs[2]);

    if (mxGetM(prhs[2]) != N) idmErrHandler(-999, "IDAS", "IDAReInit",
                                            "Size of yp0 changed from IDAInit call.", NULL);

    /* Load initial conditions */

    PutData(yy, yy0, N);
    PutData(yp, yp0, N);

    /* Extract options structure */
    
    options = prhs[3];

    break;

  }

  /* Process the options structure */

  get_IntgrOptions(options, fwdPb, TRUE,
                   &maxord, &mxsteps,
                   &itol, &reltol, &Sabstol, &Vabstol,
                   &hin, &hmax, &tstop,
                   &suppress,
                   &id, &cnstr,
                   &res_s);


  /* 
   * ----------------------------------------
   * Call appropriate IDAS functions
   *
   * If action = 0
   *    Create IDAS object and allocate memory
   *    Attach error handler function
   *    Redirect output
   * If action = 1
   *    Reinitialize solver
   * ----------------------------------------
   */

  switch (action) {

  case 0:

    /* Create IDAS object */
    ida_mem = IDACreate();
    /* Attach the global IDAS data as 'user-data' */
    IDASetUserData(ida_mem, fwdPb);
    /* Attach error handler function */
    IDASetErrHandlerFn(ida_mem, idmErrHandler);
    /* Call IDAInit */
    IDAInit(ida_mem, mxW_IDARes, t0, yy, yp);
    /* Redirect output */
    IDASetErrFile(ida_mem, stdout);

    break;

  case 1:

    /* Reinitialize solver */
    IDAReInit(ida_mem, t0, yy, yp);

    break;

  }

  /*
   * ----------------------------------------
   * Set tolerances
   * ----------------------------------------
   */

  switch (itol) {
    case IDA_SS:
      IDASStolerances(ida_mem, reltol, Sabstol);
      break;
    case IDA_SV:
      NV_abstol = N_VClone(yy);
      PutData(NV_abstol, Vabstol, N);
      IDASVtolerances(ida_mem, reltol, NV_abstol);
      N_VDestroy(NV_abstol);
      break;
    }


  /*
   * --------------------------------
   * Set various optional inputs
   * --------------------------------
   */

  /* set maxorder (default is 5) */
  IDASetMaxOrd(ida_mem, maxord);

  /* set initial step size (the default value of 0.0 is ignored by IDAS) */
  IDASetInitStep(ida_mem, hin);

  /* set max step (default is infinity) */
  IDASetMaxStep(ida_mem, hmax);

  /* set number of max steps */
  IDASetMaxNumSteps(ida_mem, mxsteps);

  /* set suppressAlg */
  IDASetSuppressAlg(ida_mem, suppress);

  /* set tstop? */
  if (tstopSet)
    IDASetStopTime(ida_mem, tstop);
  
  /* Rootfinding? */
  if ( !mxIsEmpty(mtlb_Gfct) && (Ng > 0) ) {
    IDARootInit(ida_mem, Ng, mxW_IDAGfct);
    rootSet = TRUE;
  } else {
    rootSet = FALSE;
  }

  /* ID vector specified? */
  if (id != NULL) {
    NV_id = N_VClone(yy);
    PutData(NV_id, id, N);
    IDASetId(ida_mem, NV_id);
    N_VDestroy(NV_id);
  }

  /* Constraint vector specified? */
  if (cnstr != NULL) {
    NV_cnstr = N_VClone(yy);
    PutData(NV_cnstr, cnstr, N);
    IDASetConstraints(ida_mem, NV_cnstr);
    N_VDestroy(NV_cnstr);
  }


  /*
   * ----------------------------------------
   * Linear solver
   * ----------------------------------------
   */

  get_LinSolvOptions(options, fwdPb, TRUE,
                     &mupper, &mlower,
                     &mudq, &mldq, &dqrely,
                     &gstype, &maxl);
  
  switch (ls) {

  case LS_DENSE:

    IDADense(ida_mem, N);
    if (!mxIsEmpty(mtlb_JACfct)) IDADlsSetDenseJacFn(ida_mem, mxW_IDADenseJac);

    break;

  case LS_BAND:

    IDABand(ida_mem, N, mupper, mlower);
    if (!mxIsEmpty(mtlb_JACfct)) IDADlsSetBandJacFn(ida_mem, mxW_IDABandJac);

    break;

  case LS_SPGMR:

    IDASpgmr(ida_mem, maxl);
    IDASpilsSetGSType(ida_mem, gstype);

    break;

  case LS_SPBCG:

    IDASpbcg(ida_mem, maxl);

    break;

  case LS_SPTFQMR:

    IDASptfqmr(ida_mem, maxl);

    break;
    
  }

  /* Jacobian * vector and preconditioner for SPILS linear solvers */

  if ( (ls==LS_SPGMR) || (ls==LS_SPBCG) || (ls==LS_SPTFQMR) ) {

    if (!mxIsEmpty(mtlb_JACfct)) IDASpilsSetJacTimesVecFn(ida_mem, mxW_IDASpilsJac);

    switch (pm) {

    case PM_NONE:

      if (!mxIsEmpty(mtlb_PSOLfct)) {
        if (!mxIsEmpty(mtlb_PSETfct)) IDASpilsSetPreconditioner(ida_mem, mxW_IDASpilsPset, mxW_IDASpilsPsol);
        else                          IDASpilsSetPreconditioner(ida_mem, NULL, mxW_IDASpilsPsol);
      }

      break;

    case PM_BBDPRE:

      if (!mxIsEmpty(mtlb_GCOMfct)) IDABBDPrecInit(ida_mem, N, mudq, mldq, mupper, mlower, dqrely, mxW_IDABBDgloc, mxW_IDABBDgcom);
      else                          IDABBDPrecInit(ida_mem, N, mudq, mldq, mupper, mlower, dqrely, mxW_IDABBDgloc, NULL);

      break;

    }

  }

  /* Do we monitor? */
  
  if (mon) mxW_IDAMonitor(0, t0, NULL, NULL, NULL, fwdPb);

  return(0);

}

/* IDM_QuadInitialization
 *
 * action = 0   -> IDAQuadInit
 * prhs contains:
 *   fQ
 *   y0
 *   options
 *
 * action = 1   -> IDAQuadReInit
 * prhs contains:
 *   y0
 *   options
 *
 */

static int IDM_QuadInitialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  idmPbData fwdPb;

  const mxArray *options;

  double *yQ0;

  booleantype rhs_s; /* ignored */

  booleantype errconQ;
  int itolQ;
  realtype reltolQ, SabstolQ, *VabstolQ;
  N_Vector NV_abstolQ;

  fwdPb = idmData->fwdPb;

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

    if (mxGetM(prhs[0]) != Nq) idmErrHandler(-999, "IDAS", "IDAQuadReInit",
                                             "Size of yQ0 changed from IDAQuadInit call.", NULL);

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
   * Call appropriate IDAS functions
   *
   * If action = 0
   *    Initialize quadratures
   * If action = 1
   *    Reinitialize quadratures
   * ----------------------------------------
   */

  switch (action) {
  case 0:
    IDAQuadInit(ida_mem, mxW_IDAQuadFct, yQ);
    break;
  case 1:
    IDAQuadReInit(ida_mem, yQ);
    break;
  }

  /*
   * ----------------------------------------
   * Set tolerances for quadrature variables
   * ----------------------------------------
   */

  IDASetQuadErrCon(ida_mem, errconQ);

  if (errconQ) {
    
    switch (itolQ) {
    case IDA_SS:
      IDAQuadSStolerances(ida_mem, reltolQ, SabstolQ);
      break;
    case IDA_SV:
      NV_abstolQ = N_VClone(yQ);
      PutData(NV_abstolQ, VabstolQ, Nq);
      IDAQuadSVtolerances(ida_mem, reltolQ, NV_abstolQ);
      N_VDestroy(NV_abstolQ);
      break;
    }
    
  }

  /* Quadratures will be integrated */

  quadr = TRUE;

  return(0);

}


/* IDM_SensInitialization
 * action = 0 -> IDASensMalloc
 * action = 1 -> IDASensReInit
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

static int IDM_SensInitialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  idmPbData fwdPb;

  const mxArray *options;

  booleantype fS_DQ;
  IDASensResFn resS;

  double *yyS0, *ypS0;

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

  fwdPb = idmData->fwdPb;

  /* 
   * ------------------------------------
   * Process inputs based on action
   * ------------------------------------
   */

  switch (action) {

  case 0:     /* FSA INITIALIZATION */

    /* Extract number of sensitivities */

    Ns = (int)mxGetScalar(prhs[0]);

    /* Extract user-provided sensitivity residual function */

    if ( mxIsEmpty(prhs[1]) ) {
      resS = NULL;
      fS_DQ = TRUE;
    } else {
      mxDestroyArray(mtlb_SRESfct);
      mtlb_SRESfct = mxDuplicateArray(prhs[1]);
      resS = mxW_IDASensRes;
      fS_DQ = FALSE;
    }

    /* Extract sensitivity initial conditions */

    yyS0 = mxGetPr(prhs[2]);
    ypS0 = mxGetPr(prhs[3]);

    /* Create the sensitivity N_Vectors */

    yyS = N_VCloneVectorArray(Ns, yy);
    ypS = N_VCloneVectorArray(Ns, yy);

    /* Load sensitivity initial conditions */

    for (is=0;is<Ns;is++) {
      PutData(yyS[is], &yyS0[is*N], N);
      PutData(ypS[is], &ypS0[is*N], N);
    }

    /* Extract FSA options structure */

    options = prhs[4];

    break;

  case 1:     /* FSA RE-INITIALIZATION */

    /* Extract sensitivity initial condition */

    yyS0 = mxGetPr(prhs[0]);
    ypS0 = mxGetPr(prhs[1]);

    if ( (mxGetM(prhs[0]) != N) || (mxGetN(prhs[0]) != Ns) )  idmErrHandler(-999, "IDAS", "IDASensReInit",
                                                                            "Size of yyS0 changed from IDASensInit call.", NULL);
    if ( (mxGetM(prhs[1]) != N) || (mxGetN(prhs[1]) != Ns) )  idmErrHandler(-999, "IDAS", "IDASensReInit",
                                                                            "Size of ypS0 changed from IDASensInit call.", NULL);

    /* Load sensitivity initial conditions */

    for (is=0;is<Ns;is++) {
      PutData(yyS[is], &yyS0[is*N], N);
      PutData(ypS[is], &ypS0[is*N], N);
    }

    /* Extract qFSA options structure */

    options = prhs[2];

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
   * Call appropriate IDAS functions
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

    if (fS_DQ) {

      if (pfield_name == NULL) idmErrHandler(-999, "IDAS", "IDASensInit/IDASensReInit",
                                             "pfield required but was not provided.", NULL);

      pfield = mxGetField(mtlb_data,0,pfield_name);

      if (pfield == NULL) idmErrHandler(-999, "IDAS", "IDASensInit/IDASensReInit",
                                        "illegal pfield input.", NULL);

      p = mxGetPr(pfield);

    }

    IDASensInit(ida_mem, Ns, ism, resS, yyS, ypS);

    break;

  case 1:

    IDASensReInit(ida_mem, ism, yyS, ypS);

    break;

  }

  /*
   * ----------------------------------------
   * Set tolerances for sensitivity variables
   * ----------------------------------------
   */

  switch (itolS) {
  case IDA_SS:
    IDASensSStolerances(ida_mem, reltolS, SabstolS);
    break;
  case IDA_SV:
    NV_abstolS = N_VCloneVectorArray(Ns, yy);
    for (is=0;is<Ns;is++)
      PutData(NV_abstolS[is], &VabstolS[is*N], N);
    IDASensSVtolerances(ida_mem, reltolS, NV_abstolS);
    N_VDestroyVectorArray(NV_abstolS, Ns);
    break;
  case IDA_EE:
    IDASensEEtolerances(ida_mem);
    break;
  }

  /*
   * --------------------------------
   * Set various optional inputs
   * --------------------------------
   */
  
  IDASetSensParams(ida_mem, p, pbar, plist);

  IDASetSensDQMethod(ida_mem, dqtype, rho);
  
  IDASetSensErrCon(ida_mem, errconS);

  fsa = TRUE;

  return(0);
}

/*
 * IDM_SensToggleOff
 *
 * deactivates FSA
 */

static int IDM_SensToggleOff(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  idmPbData fwdPb;

  fwdPb = idmData->fwdPb;

  IDASensToggleOff(ida_mem);
  
  fsa = FALSE;

  return(0);
}


/* IDM_AdjInitialization
 *
 * prhs contains:
 *
 * plhs contains:
 *   status
 */

static int IDM_AdjInitialization(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
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
    if(status != 0) idmErrHandler(-999, "IDAS", "IDAAdjInit", 
                                  "Could not parse InterpType.", NULL);

    if(!strcmp(bufval,"Hermite"))         interp = IDA_HERMITE;
    else if(!strcmp(bufval,"Polynomial")) interp = IDA_POLYNOMIAL;
    else idmErrHandler(-999, "IDAS", "IDAAdjInit",
                       "Interp. type has an illegal value.", NULL);

    IDAAdjInit(ida_mem, Nd, interp);

    break;

  case 1:

    IDAAdjReInit(ida_mem);

    break;

  }

  asa = TRUE;

  return(0);
}


/* IDM_InitializationB
 *
 * action = 0   -> IDACreateB + IDAInitB
 * action = 1   -> IDAReInitB
 *
 * prhs contains:
 *   resB
 *   tF
 *   yyB0
 *   ypB0
 *   options
 *   data
 *
 * plhs contains:
 *   status
 *
 */
static int IDM_InitializationB(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  idmPbData bckPb;

  const mxArray *options;

  int idxB;

  double tB0, *yyB0, *ypB0;

  int maxordB;
  long int mxstepsB;

  int itolB;
  realtype reltolB, SabstolB, *VabstolB;
  N_Vector NV_abstolB;

  double hinB, hmaxB;
  double tstopB;            /* ignored */

  booleantype suppressB;

  int mupperB, mlowerB;
  int gstypeB, maxlB;
  int mudqB, mldqB;
  double dqrelyB;

  double *idB, *cnstrB;
  N_Vector NV_idB;

  booleantype res_s;

  booleantype found_bck;

  /* 
   * -----------------------------
   * Finalize Forward monitoring
   * -----------------------------
   */

  if (idmData->fwdPb->Mon) {
    mxW_IDAMonitor(2, 0.0, NULL, NULL, NULL, idmData->fwdPb);
    idmData->fwdPb->Mon = FALSE;
  }

  /* 
   * ------------------------------------
   * Process inputs based on action
   * ------------------------------------
   */

  switch (action) {

  case 0:     /* BACKWARD SOLVER INITIALIZATION */

    /* Create and initialize a new problem */

    bckPb = (idmPbData) mxMalloc(sizeof(struct idmPbData_));
    idmInitPbData(bckPb);

    bckPb->next = idmData->bckPb;
    idmData->bckPb = bckPb;

    /* Extract user-provided RHS function */

    mxDestroyArray(mtlb_RESfctB);
    mtlb_RESfctB = mxDuplicateArray(prhs[0]);

    /* Extract final time */

    tB0 = (double)mxGetScalar(prhs[1]);

    /* Extract final conditions */

    yyB0 = mxGetPr(prhs[2]);
    ypB0 = mxGetPr(prhs[3]);
    NB = mxGetM(prhs[2]);

    /* Create the solution N_Vectors */

    yyB = NewVector(NB);
    ypB = NewVector(NB);

    /* Load final conditions */

    PutData(yyB, yyB0, NB);
    PutData(ypB, ypB0, NB);

    /* Extract options structure */
    
    options = prhs[4];

    break;

  case 1:     /* BACKWARD SOLVER RE-INITIALIZATION */

    /* Extract index of current backward problem */
    
    idxB = (int)mxGetScalar(prhs[0]);

    /* Find current backward problem */

    found_bck = FALSE;
    bckPb = idmData->bckPb;
    while (bckPb != NULL) {
      if (indexB == idxB) {
        found_bck = TRUE;
        break;
      }
      bckPb = bckPb->next;
    }

    if (!found_bck) idmErrHandler(-999, "IDAS", "IDAReInitB",
                                  "idxB has an illegal value.", NULL);

    /* If backward monitoring was enabled, finalize it now. */

    if (monB) mxW_IDAMonitorB(2, indexB, 0.0, NULL, NULL, bckPb);

    /* Extract final time */

    tB0 = (double)mxGetScalar(prhs[1]);

    /* Extract final conditions */

    yyB0 = mxGetPr(prhs[2]);

    if (mxGetM(prhs[2]) != NB) idmErrHandler(-999, "IDAS", "IDAReInitB",
                                             "Size of yyB0 changed from IDAInitB call.", NULL);

    yyB0 = mxGetPr(prhs[3]);

    if (mxGetM(prhs[3]) != NB) idmErrHandler(-999, "IDAS", "IDAReInitB",
                                             "Size of ypB0 changed from IDAInitB call.", NULL);

    /* Load final conditions */

    PutData(yyB, yyB0, NB);
    PutData(ypB, ypB0, NB);

    /* Extract options structure */
    
    options = prhs[4];

    break;

  }

  /* Process the options structure */

  get_IntgrOptions(options, bckPb, FALSE,
                   &maxordB, &mxstepsB,
                   &itolB, &reltolB, &SabstolB, &VabstolB,
                   &hinB, &hmaxB, &tstopB,
                   &suppressB,
                   &idB, &cnstrB,
                   &res_s);

  /* 
   * ----------------------------------------
   * Call appropriate IDAS functions
   *
   * If action = 0
   *    Create IDAS object and allocate memory
   *    Initialize and allocate memory
   * If action = 1
   *    Reinitialize solver
   * ----------------------------------------
   */

  switch (action) {

  case 0:

    IDACreateB(ida_mem, &idxB);
    IDASetUserDataB(ida_mem, idxB, bckPb);

    if (res_s) IDAInitBS(ida_mem, idxB, mxW_IDAResBS, tB0, yyB, ypB);
    else       IDAInitB(ida_mem, idxB, mxW_IDAResB, tB0, yyB, ypB);

    /* Return idxB */

    plhs[0] = mxCreateScalarDouble((double)idxB);

    indexB = idxB;

    NbckPb++;

    break;

  case 1:

    IDAReInitB(ida_mem, idxB, tB0, yyB, ypB);

    break;

  }

  /*
   * ----------------------------------------
   * Set tolerances
   * ----------------------------------------
   */

  switch (itolB) {
  case IDA_SS:
    IDASStolerancesB(ida_mem, idxB, reltolB, SabstolB);
    break;
  case IDA_SV:
    NV_abstolB = N_VClone(yyB);
    PutData(NV_abstolB, VabstolB, NB);
    IDASVtolerancesB(ida_mem, idxB, reltolB, NV_abstolB);
    N_VDestroy(NV_abstolB);
    break;
  }

  /*
   * --------------------------------
   * Set various optional inputs
   * --------------------------------
   */

  /* set maxorder (default is consistent with LMM) */
  IDASetMaxOrdB(ida_mem, idxB, maxordB);

  /* set initial step size (the default value of 0.0 is ignored by IDAS) */
  IDASetInitStepB(ida_mem, idxB, hinB);

  /* set max step (default is infinity) */
  IDASetMaxStepB(ida_mem, idxB, hmaxB);

  /* set number of max steps */
  IDASetMaxNumStepsB(ida_mem, idxB, mxstepsB);

  /* set suppressAlg */
  IDASetSuppressAlgB(ida_mem, idxB, suppressB);

  /* ID vector specified? */
  if (idB != NULL) {
    NV_idB = N_VClone(yyB);
    PutData(NV_idB, idB, NB);
    IDASetIdB(ida_mem, idxB, NV_idB);
    N_VDestroy(NV_idB);
  }

  /*
   * ----------------------------------------
   * Linear solver
   * ----------------------------------------
   */

  get_LinSolvOptions(options, bckPb, FALSE,
                     &mupperB, &mlowerB,
                     &mudqB, &mldqB, &dqrelyB,
                     &gstypeB, &maxlB);

  switch(lsB) {

  case LS_DENSE:

    IDADenseB(ida_mem, idxB, NB);
    if (!mxIsEmpty(mtlb_JACfctB)) IDADlsSetDenseJacFnB(ida_mem, idxB, mxW_IDADenseJacB);

    break;

  case LS_BAND:

    IDABandB(ida_mem, idxB, NB, mupperB, mlowerB);
    if (!mxIsEmpty(mtlb_JACfctB)) IDADlsSetBandJacFnB(ida_mem, idxB, mxW_IDABandJacB);

    break;

  case LS_SPGMR:

    IDASpgmrB(ida_mem, idxB, maxlB);
    IDASpilsSetGSTypeB(ida_mem, idxB, gstypeB);

    break;

  case LS_SPBCG:

    IDASpbcgB(ida_mem, idxB, maxlB);

    break;

  case LS_SPTFQMR:

    IDASptfqmrB(ida_mem, idxB, maxlB);

    break;

  }

  /* Jacobian * vector and preconditioner for SPILS linear solvers */

  if ( (lsB==LS_SPGMR) || (lsB==LS_SPBCG) || (lsB==LS_SPTFQMR) ) {

    if (!mxIsEmpty(mtlb_JACfctB)) IDASpilsSetJacTimesVecFnB(ida_mem, idxB, mxW_IDASpilsJacB);

    switch (pmB) {

    case PM_NONE:

      if (!mxIsEmpty(mtlb_PSOLfctB)) {
        if (!mxIsEmpty(mtlb_PSETfctB)) IDASpilsSetPreconditionerB(ida_mem, idxB, mxW_IDASpilsPsetB, mxW_IDASpilsPsolB);
        else                           IDASpilsSetPreconditionerB(ida_mem, idxB, NULL, mxW_IDASpilsPsolB);
      }

      break;

    case PM_BBDPRE:

      if (!mxIsEmpty(mtlb_GCOMfctB)) IDABBDPrecInitB(ida_mem, idxB, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mxW_IDABBDglocB, mxW_IDABBDgcomB);
      else                           IDABBDPrecInitB(ida_mem, idxB, NB, mudqB, mldqB, mupperB, mlowerB, dqrelyB, mxW_IDABBDglocB, NULL);

      break;

    }

  }

  /* Do we monitor? */

  if (monB) mxW_IDAMonitorB(0, idxB, tB0, NULL, NULL, bckPb);

  return(0);
}

/* IDM_QuadInitializationB
 *
 * action = 0   -> IDAQuadInitB
 * prhs contains:
 *   idxB
 *   fQB
 *   yQB0
 *   options
 *
 * action = 1   -> IDAQuadReInitB
 *   idxB
 *   yQB0
 *   options
 *
 */

static int IDM_QuadInitializationB(int action, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  idmPbData bckPb;

  const mxArray *options;

  int idxB;

  double *yQB0;

  booleantype rhs_s;
  
  booleantype errconQB;
  int itolQB;
  realtype reltolQB, SabstolQB, *VabstolQB;
  N_Vector NV_abstolQB;

  booleantype found_bck;

  /* Extract index of current backward problem */
    
  idxB = (int)mxGetScalar(prhs[0]);

  /* Find current backward problem */

  found_bck = FALSE;
  bckPb = idmData->bckPb;
  while (bckPb != NULL) {
    if (indexB == idxB) {
      found_bck = TRUE;
      break;
    }
    bckPb = bckPb->next;
  }

  if (!found_bck) idmErrHandler(-999, "IDAS", "IDAQuadInitB/IDAQuadReInitB",
                                "idxB has an illegal value.", NULL);

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

    if (mxGetM(prhs[1]) != NqB) idmErrHandler(-999, "IDAS", "IDAQuadReInitB",
                                              "Size of yQB0 changed from IDAQuadInitB call.", NULL);

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
   * Call appropriate IDAS functions
   *
   * If action = 0
   *    Initialize backward quadratures
   * If action = 1
   *    Reinitialize backward quadratures
   * ----------------------------------------
   */

  switch (action) {
  case 0:
    if (rhs_s) IDAQuadInitBS(ida_mem, idxB, mxW_IDAQuadFctBS, yQB);
    else       IDAQuadInitB(ida_mem, idxB, mxW_IDAQuadFctB, yQB);
    break;
  case 1:
    IDAQuadReInitB(ida_mem, idxB, yQB);
    break;
  }

  /*
   * ----------------------------------------
   * Set tolerances for quadrature variables
   * ----------------------------------------
   */
  
  IDASetQuadErrConB(ida_mem, idxB, errconQB);

  if (errconQB) {

    switch (itolQB) {
    case IDA_SS:
      IDAQuadSStolerancesB(ida_mem, idxB, reltolQB, SabstolQB);
      break;
    case IDA_SV:
      NV_abstolQB = N_VClone(yQB);
      PutData(NV_abstolQB, VabstolQB, NqB);
      IDAQuadSVtolerancesB(ida_mem, idxB, reltolQB, NV_abstolQB);
      N_VDestroy(NV_abstolQB);
      break;
    }
    
  }

  quadrB = TRUE;

  return(0);
}

static int IDM_CalcIC(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  idmPbData fwdPb;

  double tout;
  int buflen, status, icopt;
  char *bufval;


  fwdPb = idmData->fwdPb;


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
  return(0);
}

static int IDM_Solve(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  idmPbData fwdPb;

  int buflen;
  char *bufval;

  int nlhs_needed, dims[3];

  int itask, is, Ntout, itout, status, s_idx;
  double *tout, tret, h;
  double *tdata, *yydata, *yQdata, *yySdata;
  long int nst;


  fwdPb = idmData->fwdPb;


  /*
   * ----------------------------------------------------------------
   * Verify if number of output arguments agrees with current options
   * ----------------------------------------------------------------
   */

  nlhs_needed = 3;

  if (quadr) nlhs_needed++;
  if (fsa)   nlhs_needed++;

  if (nlhs < nlhs_needed) idmErrHandler(-999, "IDAS", "IDASolve",
                                        "Too few output arguments.", NULL);
  if (nlhs > nlhs_needed) idmErrHandler(-999, "IDAS", "IDASolve",
                                        "Too many output arguments.", NULL);

  /*
   * ----------------------------------------------------------------
   * Extract input arguments
   * ----------------------------------------------------------------
   */

  /* Extract tout */

  Ntout = mxGetM(prhs[0]) * mxGetN(prhs[0]);
  tout = mxGetPr(prhs[0]);

  /* If rootfinding or tstop are enabled, we do not allow multiple output times */

  if (rootSet && (Ntout>1)) idmErrHandler(-999, "IDAS", "IDASolve",
                                          "More than one tout value prohibited with rootfinding enabled.", NULL);

  if (tstopSet && (Ntout>1)) idmErrHandler(-999, "IDAS", "IDASolve",
                                           "More than one tout value prohibited with tstop enabled.", NULL);

  /* Extract itask */

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Normal")) itask = IDA_NORMAL;
  else if(!strcmp(bufval,"OneStep")) itask = IDA_ONE_STEP;
  else idmErrHandler(-999, "IDAS", "IDASolve",
                     "Illegal value for itask.", NULL); 



  if (itask==IDA_ONE_STEP) {

    /* If itask==IDA_ONE_STEP, we do not allow multiple output times and we do not monitor */

    if (Ntout>1) idmErrHandler(-999, "IDAS", "IDASolve",
                               "More than one tout value prohibited in ONE_STEP mode.", NULL); 

    if (mon) {
      idmErrHandler(+999, "IDAS", "IDASolve",
                    "Monitoring disabled in ONE_STEP mode.", NULL);
      mon = FALSE;
    }

  } else {

    /* Check if tout values are legal */

    IDAGetCurrentTime(ida_mem, &tret);
    IDAGetNumSteps(ida_mem, &nst);

    /* h is used throughout this function as integration direction only */
    if (nst == 0) {
      h = tout[0] - tret;
    } else {
      IDAGetLastStep(ida_mem, &h);
      if ( (tout[0] - tret + h)*h < 0.0 ) idmErrHandler(-999, "IDAS", "IDASolve",
                                                        "Illegal value of tout.", NULL);
    }
    
    for (itout=1; itout<Ntout; itout++) 
      if ( (tout[itout] - tout[itout-1])*h < 0.0 ) idmErrHandler(-999, "IDAS", "IDASolve",
                                                                 "tout values are not monotonic.", NULL);

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
  yydata = mxGetPr(plhs[2]);

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
    yySdata = mxGetPr(plhs[s_idx]);
  }

  /*
   * ----------------------------------------------------------------
   * Call the IDAS main solver function
   * ----------------------------------------------------------------
   */

  if (!mon) {

    /* No monitoring. itask can be either IDA_ONE_STEP or IDA_NORMAL */

    for (itout=0; itout<Ntout; itout++) {

      if (!asa) {
        status = IDASolve(ida_mem, tout[itout], &tret, yy, yp, itask);
      } else {
        status = IDASolveF(ida_mem, tout[itout], &tret, yy, yp, itask, &Nc);
      }

      tdata[itout] = tret;

      GetData(yy, &yydata[itout*N], N);

      if (quadr) {
        IDAGetQuad(ida_mem, &tret, yQ);
        GetData(yQ, &yQdata[itout*Nq], Nq);
      }

      if (fsa) {
        IDAGetSens(ida_mem, &tret, yyS);
        for (is=0; is<Ns; is++)
          GetData(yyS[is], &yySdata[itout*Ns*N+is*N], N);
      }

    }

  } else {

    /* Monitoring. itask = IDA_NORMAL */

    for (itout=0; itout<Ntout; itout++) {
      
      /* In ONE_STEP mode, IDAS reads tout only at the first step.
       * We must therefore check here whether we need to take additional steps,
       * or simply return interpolated solution at tout. */

      IDAGetNumSteps(ida_mem, &nst);
      IDAGetCurrentTime(ida_mem, &tret);

      if ( (nst>0) && ((tret - tout[itout])*h >= 0.0) ) {

        /* No need to take an additional step */
        status = IDA_SUCCESS;

      } else {

        /* Take additional steps */
        while(1) {

          if (!asa) status = IDASolve(ida_mem, tout[itout], &tret, yy, yp, IDA_ONE_STEP);
          else      status = IDASolveF(ida_mem, tout[itout], &tret, yy, yp, IDA_ONE_STEP, &Nc);

          /* Call the monitoring function */
          if (quadr) IDAGetQuad(ida_mem, &tret, yQ);
          if (fsa)   IDAGetSens(ida_mem, &tret, yyS);
          mxW_IDAMonitor(1, tret, yy, yQ, yyS, fwdPb);

          /* If a root was found or tstop was reached, break out of while loop */
          if (status == IDA_TSTOP_RETURN || status == IDA_ROOT_RETURN) break;

          /* If current tout was reached break out of while loop */
          if ( (tret - tout[itout])*h >= 0.0 )  break;

        }

      }
      
      /* On a tstop or root return, return solution at tret.
       * Otherwise (status=IDA_SUCCESS), return solution at tout[itout]. */

      if (status == IDA_TSTOP_RETURN || status == IDA_ROOT_RETURN) {

        if (quadr) IDAGetQuad(ida_mem, &tret, yQ);

        if (fsa)   IDAGetSens(ida_mem, &tret, yyS);

      } else {

        tret = tout[itout];
        
        IDAGetDky(ida_mem, tret, 0, yy);

        if (quadr) IDAGetQuadDky(ida_mem, tret, 0, yQ);

        if (fsa)   IDAGetSensDky(ida_mem, tret, 0, yyS);

      }

      tdata[itout] = tret;

      GetData(yy, &yydata[itout*N], N);

      if (quadr)  GetData(yQ, &yQdata[itout*Nq], Nq);

      if (fsa)
        for (is=0; is<Ns; is++)
          GetData(yyS[is], &yySdata[itout*Ns*N+is*N], N);

    }

  }

  /* IDASolve return flag */

  plhs[0] = mxCreateScalarDouble((double)status);

  return(0);

}


static int IDM_SolveB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  idmPbData bckPb;

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
  bckPb = idmData->bckPb;
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

  if (nlhs_bad < 0) idmErrHandler(-999, "IDAS", "IDASolveB",
                                  "Too few output arguments.", NULL);
  if (nlhs_bad > 0) idmErrHandler(-999, "IDAS", "IDASolveB",
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

  IDAGetLastStep(ida_mem, &h);
  IDAGetCurrentTime(ida_mem, &tret);

  /* The stepsize of the forward problem is used to indicate the integration direction */
  if ( (tret - toutB[0])*h < 0.0 ) idmErrHandler(-999, "IDAS", "IDASolveB",
                                                 "tout value in wrong direction.", NULL);
  
  /* Extract itaskB */

  buflen = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  bufval = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], bufval, buflen);
  if(!strcmp(bufval,"Normal"))       itaskB = IDA_NORMAL;
  else if(!strcmp(bufval,"OneStep")) itaskB = IDA_ONE_STEP;
  else idmErrHandler(-999, "IDAS", "IDASolveB",
                     "Illegal value for itask.", NULL); 

  /* If itask == IDA_ONE_STEP, then
   * - we do not allow multiple output times
   * - we disable monitoring 
   */

  if ( itaskB == IDA_ONE_STEP ) {
    
    if (NtoutB > 1) idmErrHandler(-999, "IDAS", "IDASolveB",
                                  "More than one tout value prohibited in ONE_STEP mode.", NULL); 

    if (any_monB) {
      idmErrHandler(+999, "IDAS", "IDASolveB",
                    "Monitoring disabled in itask=ONE_STEP", NULL);
      bckPb = idmData->bckPb;
      while(bckPb != NULL) {
        monB = FALSE;
        bckPb = bckPb->next;
      }
      any_monB = FALSE;
    }

  }

  if (NbckPb == 1) idmSolveB_one(plhs, NtoutB, toutB, itaskB);
  else             idmSolveB_more(plhs, NtoutB, toutB, itaskB, any_quadrB, any_monB);

  return(0);
}

static void idmSolveB_one(mxArray *plhs[], int NtoutB, double *toutB, int itaskB)
{
  idmPbData bckPb;
  
  void *ida_memB;

  double tretB, hB;
  double *tdata, *ydata, *yQdata;
  int itout, status;
  long int nstB;

  bckPb = idmData->bckPb;

  ida_memB = IDAGetAdjIDABmem(ida_mem, indexB);

  /* Check if tout values are legal */

  IDAGetCurrentTime(ida_memB, &tretB);
  IDAGetNumSteps(ida_memB, &nstB);

  /* hB is used throughout this function as integration direction only */
  if (nstB == 0) {
    hB = toutB[0] - tretB;
  } else {
    IDAGetLastStep(ida_memB, &hB);
    if ( (toutB[0] - tretB + hB)*hB < 0.0 ) idmErrHandler(-999, "IDAS", "IDASolveB",
                                                          "Illegal value of tout.", NULL);
  }

  for (itout=1; itout<NtoutB; itout++) 
    if ( (toutB[itout] - toutB[itout-1])*hB < 0.0 ) idmErrHandler(-999, "IDAS", "IDASolveB",
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
   * Call the IDASolveB main solver function
   * ----------------------------------------------------------------
   */

  if (!monB) {

    /* No monitoring. itaskB can be either IDA_ONE_STEP or IDA_NORMAL */

    for (itout=0; itout<NtoutB; itout++) {
      
      status = IDASolveB(ida_mem, toutB[itout], itaskB);
      
      IDAGetB(ida_mem, indexB, &tretB, yyB, ypB);
      
      tdata[itout] = tretB;
      
      GetData(yyB, &ydata[itout*NB], NB);
      
      if (quadrB) {
        IDAGetQuadB(ida_mem, indexB, &tretB, yQB);
        GetData(yQB, &yQdata[itout*NqB], NqB);
      }
      
    }
    
    
  } else {

    /* Monitoring. itask = IDA_NORMAL */

    for (itout=0; itout<NtoutB; itout++) {

      /* In ONE_STEP mode, IDAS reads tout only at the first step.
       * We must therefore check here whether we need to take additional steps,
       * or simply return interpolated solution at tout. */

      IDAGetNumSteps(ida_memB, &nstB);
      IDAGetCurrentTime(ida_memB, &tretB);

      if ( (nstB>0) && ((tretB - toutB[itout])*hB >= 0.0) ) {

        /* No need to take an additional step */
        status = IDA_SUCCESS;

      } else {

        /* Take additional steps */
        while(1) {
        
          status = IDASolveB(ida_mem, toutB[itout], IDA_ONE_STEP);
          
          /* Call the monitoring function */          
          IDAGetB(ida_mem, indexB, &tretB, yyB, ypB);
          if (quadrB) IDAGetQuadB(ida_mem, indexB, &tretB, yQB);
          mxW_IDAMonitorB(1, indexB, tretB, yyB, yQB, bckPb);
          
          /* If current tout was reached break out of while loop */
          if ( (tretB - toutB[itout])*hB >= 0.0 ) break;
          
        }

      }

      tretB = toutB[itout];

      tdata[itout] = tretB;

      IDAGetDky(ida_memB, tretB, 0, yyB);

      GetData(yyB, &ydata[itout*NB], NB);
      
      if (quadrB) {
        IDAGetQuadDky(ida_memB, tretB, 0, yQB);
        GetData(yQB, &yQdata[itout*NqB], NqB);
      }

    }

  }

  /* IDASolve return flag */

  plhs[0] = mxCreateScalarDouble((double)status);
    
  return;
}


static void idmSolveB_more(mxArray *plhs[], int NtoutB, double *toutB, int itaskB,
                           booleantype any_quadrB, booleantype any_monB)
{
  return;
}


static int IDM_Stats(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  idmPbData fwdPb;

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
    "nrSe",
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
  long int njeB, nfeB;
  long int nli, npe, nps, ncfl, njeSG, nfeSG;

  long int nfQe, netfQ;

  long int nrSe, nfeS, netfS, nsetupsS;
  long int nniS, ncfnS;

  int i, flag;
  mxArray *mxS_root, *mxS_ls, *mxS_quad, *mxS_fsa;
  mxArray *mxS_rootsfound;
  double *tmp;
  int nfields;

  if (idmData == NULL) return(0);

  fwdPb = idmData->fwdPb;

  flag = IDAGetIntegratorStats(ida_mem, &nst, &nfe, &nsetups, 
                               &netf, &qlast, &qcur, &h0used, &hlast, &hcur, &tcur);

  flag = IDAGetNonlinSolvStats(ida_mem, &nni, &ncfn);

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

    flag = IDAGetNumGEvals(ida_mem, &nge);

    nfields = sizeof(fnames_root)/sizeof(*fnames_root);
    mxS_root = mxCreateStructMatrix(1, 1, nfields, fnames_root);

    mxSetField(mxS_root, 0, "nge", mxCreateScalarDouble((double)nge));

    rootsfound = (int *) malloc(nge*sizeof(int));
    flag = IDAGetRootInfo(ida_mem, rootsfound);
    mxS_rootsfound = mxCreateDoubleMatrix(Ng,1,mxREAL);
    tmp = mxGetPr(mxS_rootsfound);
    for (i=0;i<Ng;i++)
      tmp[i] = (double)rootsfound[i];
    mxSetField(mxS_root, 0, "roots", mxS_rootsfound);

  } else {

    mxS_root = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "RootInfo", mxS_root);

  /* Quadrature Statistics */

  if (quadr) {

    flag = IDAGetQuadStats(ida_mem, &nfQe, &netfQ);

    nfields = sizeof(fnames_quad)/sizeof(*fnames_quad);
    mxS_quad = mxCreateStructMatrix(1, 1, nfields, fnames_quad);

    mxSetField(mxS_quad, 0, "nfQe",  mxCreateScalarDouble((double)nfQe));
    mxSetField(mxS_quad, 0, "netfQ", mxCreateScalarDouble((double)netfQ));

  } else {

    mxS_quad = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "QuadInfo", mxS_quad);

  /* Linear Solver Statistics */

  switch(ls){

  case LS_DENSE:
    
    flag = IDADlsGetNumJacEvals(ida_mem, &njeD);
    flag = IDADlsGetNumResEvals(ida_mem, &nfeD);
    
    nfields = sizeof(fnames_dense)/sizeof(*fnames_dense);
    mxS_ls = mxCreateStructMatrix(1, 1, nfields, fnames_dense);
    
    mxSetField(mxS_ls, 0, "name", mxCreateString("Dense"));
    mxSetField(mxS_ls, 0, "njeD", mxCreateScalarDouble((double)njeD));
    mxSetField(mxS_ls, 0, "nfeD", mxCreateScalarDouble((double)nfeD));
    
    break;

  case LS_BAND:
      
    flag = IDADlsGetNumJacEvals(ida_mem, &njeB);
    flag = IDADlsGetNumResEvals(ida_mem, &nfeB);
      
    nfields = sizeof(fnames_band)/sizeof(*fnames_band);
    mxS_ls = mxCreateStructMatrix(1, 1, nfields, fnames_band);
 
    mxSetField(mxS_ls, 0, "name", mxCreateString("Band"));
    mxSetField(mxS_ls, 0, "njeB", mxCreateScalarDouble((double)njeB));
    mxSetField(mxS_ls, 0, "nfeB", mxCreateScalarDouble((double)nfeB));
      
    break;

  case LS_SPGMR:
  case LS_SPBCG:
  case LS_SPTFQMR:

    flag = IDASpilsGetNumLinIters(ida_mem, &nli);
    flag = IDASpilsGetNumPrecEvals(ida_mem, &npe);
    flag = IDASpilsGetNumPrecSolves(ida_mem, &nps);
    flag = IDASpilsGetNumConvFails(ida_mem, &ncfl);
    flag = IDASpilsGetNumJtimesEvals(ida_mem, &njeSG);
    flag = IDASpilsGetNumResEvals(ida_mem, &nfeSG);
    
    nfields = sizeof(fnames_spils)/sizeof(*fnames_spils);
    mxS_ls = mxCreateStructMatrix(1, 1, nfields, fnames_spils);
    
    if (ls == LS_SPGMR)
      mxSetField(mxS_ls, 0, "name",  mxCreateString("GMRES"));
    else if (ls == LS_SPBCG)
      mxSetField(mxS_ls, 0, "name",  mxCreateString("BiCGStab"));
    else
      mxSetField(mxS_ls, 0, "name",  mxCreateString("TFQMR"));

    mxSetField(mxS_ls, 0, "nli",   mxCreateScalarDouble((double)nli));
    mxSetField(mxS_ls, 0, "npe",   mxCreateScalarDouble((double)npe));
    mxSetField(mxS_ls, 0, "nps",   mxCreateScalarDouble((double)nps));
    mxSetField(mxS_ls, 0, "ncfl",  mxCreateScalarDouble((double)ncfl));
    mxSetField(mxS_ls, 0, "njeSG", mxCreateScalarDouble((double)njeSG));
    mxSetField(mxS_ls, 0, "nfeSG", mxCreateScalarDouble((double)nfeSG));
    
    break;
    
  }

  mxSetField(plhs[0], 0, "LSInfo", mxS_ls);

  /* Forward Sensitivity Statistics */

  if (fsa) {

    flag = IDAGetSensStats(ida_mem, &nrSe, &nfeS, &netfS, &nsetupsS); 

    flag = IDAGetSensNonlinSolvStats(ida_mem, &nniS, &ncfnS);

    nfields = sizeof(fnames_sens)/sizeof(*fnames_sens);
    mxS_fsa = mxCreateStructMatrix(1, 1, nfields, fnames_sens);
    
    mxSetField(mxS_fsa, 0, "nrSe",     mxCreateScalarDouble((double)nrSe));
    mxSetField(mxS_fsa, 0, "nfeS",     mxCreateScalarDouble((double)nfeS));
    mxSetField(mxS_fsa, 0, "nsetupsS", mxCreateScalarDouble((double)nsetupsS));
    mxSetField(mxS_fsa, 0, "netfS",    mxCreateScalarDouble((double)netfS));
    mxSetField(mxS_fsa, 0, "nniS",     mxCreateScalarDouble((double)nniS));
    mxSetField(mxS_fsa, 0, "ncfnS",    mxCreateScalarDouble((double)ncfnS));
    
  } else {

    mxS_fsa = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "FSAInfo", mxS_fsa);
  
  return(0);
}

static int IDM_StatsB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  idmPbData bckPb;
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

  void *ida_memB;

  long int nst, nfe, nsetups, nni, ncfn, netf;
  int qlast, qcur;
  double h0used, hlast, hcur, tcur;

  long int njeD, nfeD;
  long int njeB, nfeB;
  long int nli, npe, nps, ncfl, njeSG, nfeSG;

  long int nfQe, netfQ;

  int flag;
  mxArray *mxS_ls, *mxS_quad;
  int nfields;

  booleantype found_bck;
  

  /* Extract index of current backward problem */
    
  idxB = (int)mxGetScalar(prhs[0]);

  /* Find current backward problem */

  found_bck = FALSE;
  bckPb = idmData->bckPb;
  while (bckPb != NULL) {
    if (indexB == idxB) {
      found_bck = TRUE;
      break;
    }
    bckPb = bckPb->next;
  }

  if (!found_bck) idmErrHandler(-999, "IDAS", "IDAGetStatsB",
                                "idxB has an illegal value.", NULL);

  ida_memB = IDAGetAdjIDABmem(ida_mem, indexB);

  flag = IDAGetIntegratorStats(ida_memB, &nst, &nfe, &nsetups, 
                               &netf, &qlast, &qcur, &h0used, &hlast, &hcur, &tcur);

  flag = IDAGetNonlinSolvStats(ida_memB, &nni, &ncfn);

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

    flag = IDAGetQuadStats(ida_memB, &nfQe, &netfQ);

    nfields = sizeof(fnames_quad)/sizeof(*fnames_quad);
    mxS_quad = mxCreateStructMatrix(1, 1, nfields, fnames_quad);

    mxSetField(mxS_quad, 0, "nfQe",  mxCreateScalarDouble((double)nfQe));
    mxSetField(mxS_quad, 0, "netfQ", mxCreateScalarDouble((double)netfQ));

  } else {

    mxS_quad = mxCreateDoubleMatrix(0,0,mxREAL);

  }
  
  mxSetField(plhs[0], 0, "QuadInfo", mxS_quad);

  /* Linear Solver Statistics */

  switch(lsB){

  case LS_DENSE:
    
    flag = IDADlsGetNumJacEvals(ida_memB, &njeD);
    flag = IDADlsGetNumResEvals(ida_memB, &nfeD);
    
    nfields = sizeof(fnames_dense)/sizeof(*fnames_dense);
    mxS_ls = mxCreateStructMatrix(1, 1, nfields, fnames_dense);

    mxSetField(mxS_ls, 0, "name", mxCreateString("Dense"));
    mxSetField(mxS_ls, 0, "njeD", mxCreateScalarDouble((double)njeD));
    mxSetField(mxS_ls, 0, "nfeD", mxCreateScalarDouble((double)nfeD));
    
    break;

  case LS_BAND:
      
    flag = IDADlsGetNumJacEvals(ida_memB, &njeB);
    flag = IDADlsGetNumResEvals(ida_memB, &nfeB);
      
    nfields = sizeof(fnames_band)/sizeof(*fnames_band);
    mxS_ls = mxCreateStructMatrix(1, 1, nfields, fnames_band);
 
    mxSetField(mxS_ls, 0, "name", mxCreateString("Band"));
    mxSetField(mxS_ls, 0, "njeB", mxCreateScalarDouble((double)njeB));
    mxSetField(mxS_ls, 0, "nfeB", mxCreateScalarDouble((double)nfeB));
      
    break;

  case LS_SPGMR:
  case LS_SPBCG:
  case LS_SPTFQMR:

    flag = IDASpilsGetNumLinIters(ida_memB, &nli);
    flag = IDASpilsGetNumPrecEvals(ida_memB, &npe);
    flag = IDASpilsGetNumPrecSolves(ida_memB, &nps);
    flag = IDASpilsGetNumConvFails(ida_memB, &ncfl);
    flag = IDASpilsGetNumJtimesEvals(ida_memB, &njeSG);
    flag = IDASpilsGetNumResEvals(ida_memB, &nfeSG);
    
    nfields = sizeof(fnames_spils)/sizeof(*fnames_spils);
    mxS_ls = mxCreateStructMatrix(1, 1, nfields, fnames_spils);
 
    if (lsB == LS_SPGMR)
      mxSetField(mxS_ls, 0, "name",  mxCreateString("GMRES"));
    else if (lsB == LS_SPBCG)
      mxSetField(mxS_ls, 0, "name",  mxCreateString("BiCGStab"));
    else
      mxSetField(mxS_ls, 0, "name",  mxCreateString("TFQMR"));

    mxSetField(mxS_ls, 0, "nli",   mxCreateScalarDouble((double)nli));
    mxSetField(mxS_ls, 0, "npe",   mxCreateScalarDouble((double)npe));
    mxSetField(mxS_ls, 0, "nps",   mxCreateScalarDouble((double)nps));
    mxSetField(mxS_ls, 0, "ncfl",  mxCreateScalarDouble((double)ncfl));
    mxSetField(mxS_ls, 0, "njeSG", mxCreateScalarDouble((double)njeSG));
    mxSetField(mxS_ls, 0, "nfeSG", mxCreateScalarDouble((double)nfeSG));
    
    break;
  }

  mxSetField(plhs[0], 0, "LSInfo", mxS_ls);

  return(0);
}

static int IDM_Set(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  idmPbData fwdPb;

  const mxArray *options;
  mxArray *opt;

  double tstop;


  fwdPb = idmData->fwdPb;
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
    IDASetStopTime(ida_mem, tstop);
  }

  return(0);
}

static int IDM_SetB(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  return(0);
}

static int IDM_Get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  idmPbData fwdPb;

  double t;
  N_Vector ewt;
  double *this, *next;
  int key, k, i, nfields;

  IDAadjCheckPointRec *ckpnt;
  const char *fnames_ckpnt[]={
    "t0",
    "t1",
    "nstep",
    "order",
    "step"
  };


  fwdPb = idmData->fwdPb;

  key = (int) (*mxGetPr(prhs[0]));

  switch (key) {

  case 1:    /* DerivSolution */

    t = *mxGetPr(prhs[1]);
    k = (int) (*mxGetPr(prhs[2]));

    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    IDAGetDky(ida_mem, t, k, yy);
    GetData(yy, mxGetPr(plhs[0]), N);

    break;

  case 2:    /* ErrorWeights */

    ewt = N_VClone(yy);

    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    IDAGetErrWeights(ida_mem, ewt);
    GetData(ewt, mxGetPr(plhs[0]), N);

    N_VDestroy(ewt);

    break;

  case 3:    /* not used */

    break;

  case 4:    /* CheckPointsInfo */

    ckpnt = (IDAadjCheckPointRec *) malloc ( (Nc+1)*sizeof(IDAadjCheckPointRec));
    IDAGetAdjCheckPointsInfo(ida_mem, ckpnt);
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

  }

  return(0);
}

static int IDM_Free(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  idmPbData fwdPb, bckPb;

  if (idmData == NULL) return(0);
  
  fwdPb = idmData->fwdPb;
  if (mon) mxW_IDAMonitor(2, 0.0, NULL, NULL, NULL, fwdPb);

  bckPb = idmData->bckPb;
  while (bckPb != NULL) {
    if (monB) mxW_IDAMonitorB(2, indexB, 0.0, NULL, NULL, bckPb);   
    bckPb = bckPb->next;
  }

  IDAFree(&ida_mem);

  return(0);
}

