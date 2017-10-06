/*
 * -----------------------------------------------------------------
 * $Revision: 4845 $
 * $Date: 2016-08-03 15:45:09 -0700 (Wed, 03 Aug 2016) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan Hindmarsh, Radu Serban and Aaron Collier @ LLNL
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
 * This is the implementation file for the main IDA solver.
 * It is independent of the linear solver in use.
 * -----------------------------------------------------------------
 *
 * EXPORTED FUNCTIONS
 * ------------------
 *   Creation, allocation and re-initialization functions
 *       IDACreate
 *       IDAInit
 *       IDAReInit
 *       IDARootInit
 *   Main solver function
 *       IDASolve
 *   Interpolated output and extraction functions
 *       IDAGetDky
 *   Deallocation functions
 *       IDAFree
 *
 * PRIVATE FUNCTIONS
 * -----------------
 *       IDACheckNvector
 *   Memory allocation/deallocation
 *       IDAAllocVectors
 *       IDAFreeVectors
 *   Initial setup
 *       IDAInitialSetup
 *       IDAEwtSet
 *       IDAEwtSetSS
 *       IDAEwtSetSV
 *   Stopping tests
 *       IDAStopTest1
 *       IDAStopTest2
 *   Error handler
 *       IDAHandleFailure
 *   Main IDAStep function
 *       IDAStep
 *       IDASetCoeffs
 *   Nonlinear solver functions
 *       IDANls
 *       IDAPredict
 *       IDANewtonIter
 *   Error test
 *       IDATestError
 *       IDARestore
 *   Handler for convergence and/or error test failures
 *       IDAHandleNFlag
 *       IDAReset
 *   Function called after a successful step
 *       IDACompleteStep
 *   Get solution
 *       IDAGetSolution
 *   Norm functions
 *       IDAWrmsNorm
 *   Functions for rootfinding
 *       IDARcheck1
 *       IDARcheck2
 *       IDARcheck3
 *       IDARootfind
 *   IDA Error message handling functions 
 *       IDAProcessError
 *       IDAErrHandler
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "ida_impl.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * MACRO DEFINITIONS
 * =================================================================
 */

/* Macro: loop */
#define loop for(;;)

/* 
 * =================================================================
 * IDAS PRIVATE CONSTANTS
 * =================================================================
 */

#define ZERO      RCONST(0.0)    /* real 0.0    */
#define HALF      RCONST(0.5)    /* real 0.5    */
#define QUARTER   RCONST(0.25)   /* real 0.25   */
#define TWOTHIRDS RCONST(0.667)  /* real 2/3    */
#define ONE       RCONST(1.0)    /* real 1.0    */
#define ONEPT5    RCONST(1.5)    /* real 1.5    */
#define TWO       RCONST(2.0)    /* real 2.0    */
#define FOUR      RCONST(4.0)    /* real 4.0    */
#define FIVE      RCONST(5.0)    /* real 5.0    */
#define TEN       RCONST(10.0)   /* real 10.0   */
#define TWELVE    RCONST(12.0)   /* real 12.0   */
#define TWENTY    RCONST(20.0)   /* real 20.0   */
#define HUNDRED   RCONST(100.0)  /* real 100.0  */
#define PT9       RCONST(0.9)    /* real 0.9    */
#define PT99      RCONST(0.99)   /* real 0.99   */
#define PT1       RCONST(0.1)    /* real 0.1    */
#define PT01      RCONST(0.01)   /* real 0.01   */
#define PT001     RCONST(0.001)  /* real 0.001  */
#define PT0001    RCONST(0.0001) /* real 0.0001 */

/* 
 * =================================================================
 * IDAS ROUTINE-SPECIFIC CONSTANTS
 * =================================================================
 */

/* 
 * Control constants for lower-level functions used by IDASolve 
 * ------------------------------------------------------------
 */

/* IDAStep control constants */

#define PREDICT_AGAIN 20

/* Return values for lower level routines used by IDASolve */

#define IDA_RES_RECVR    +1
#define IDA_LSETUP_RECVR +2
#define IDA_LSOLVE_RECVR +3

#define IDA_NCONV_RECVR  +4
#define IDA_CONSTR_RECVR +5
#define CONTINUE_STEPS   +99

/* IDACompleteStep constants */

#define UNSET            -1
#define LOWER            +1 
#define RAISE            +2 
#define MAINTAIN         +3

/* IDATestError constants */

#define ERROR_TEST_FAIL  +7

/*
 * Control constants for lower-level rootfinding functions
 * -------------------------------------------------------
 */

#define RTFOUND          +1
#define CLOSERT          +3

/*
 * Control constants for tolerances
 * --------------------------------
 */

#define IDA_NN  0
#define IDA_SS  1
#define IDA_SV  2
#define IDA_WF  3

/*
 * Algorithmic constants
 * ---------------------
 */

#define MXNCF           10  /* max number of convergence failures allowed */
#define MXNEF           10  /* max number of error test failures allowed  */
#define MAXNH            5  /* max. number of h tries in IC calc. */
#define MAXNJ            4  /* max. number of J tries in IC calc. */
#define MAXNI           10  /* max. Newton iterations in IC calc. */
#define EPCON RCONST(0.33)  /* Newton convergence test constant */
#define MAXBACKS       100  /* max backtracks per Newton step in IDACalcIC */

/* IDANewtonIter constants */

#define MAXIT   4
#define RATEMAX RCONST(0.9)
#define XRATE   RCONST(0.25)        

/* 
 * =================================================================
 * PRIVATE FUNCTION PROTOTYPES
 * =================================================================
 */

static booleantype IDACheckNvector(N_Vector tmpl);

/* Memory allocation/deallocation */

static booleantype IDAAllocVectors(IDAMem IDA_mem, N_Vector tmpl);
static void IDAFreeVectors(IDAMem IDA_mem);

/* Initial setup */

int IDAInitialSetup(IDAMem IDA_mem);
static int IDAEwtSetSS(IDAMem IDA_mem, N_Vector ycur, N_Vector weight);
static int IDAEwtSetSV(IDAMem IDA_mem, N_Vector ycur, N_Vector weight);

/* Main IDAStep function */

static int IDAStep(IDAMem IDA_mem);

/* Function called at beginning of step */

static void IDASetCoeffs(IDAMem IDA_mem, realtype *ck);

/* Nonlinear solver functions */

static void IDAPredict(IDAMem IDA_mem);
static int IDANls(IDAMem IDA_mem);
static int IDANewtonIter(IDAMem IDA_mem);

/* Error test */

static int IDATestError(IDAMem IDA_mem, realtype ck, 
                        realtype *err_k, realtype *err_km1);

/* Handling of convergence and/or error test failures */

static void IDARestore(IDAMem IDA_mem, realtype saved_t);
static int IDAHandleNFlag(IDAMem IDA_mem, int nflag, realtype err_k, realtype err_km1,
                          long int *ncfnPtr, int *ncfPtr, long int *netfPtr, int *nefPtr);
static void IDAReset(IDAMem IDA_mem);

/* Function called after a successful step */

static void IDACompleteStep(IDAMem IDA_mem, realtype err_k, realtype err_km1);

/* Function called to evaluate the solutions y(t) and y'(t) at t */

int IDAGetSolution(void *ida_mem, realtype t, N_Vector yret, N_Vector ypret);

/* Stopping tests and failure handling */

static int IDAStopTest1(IDAMem IDA_mem, realtype tout,realtype *tret, 
                        N_Vector yret, N_Vector ypret, int itask);
static int IDAStopTest2(IDAMem IDA_mem, realtype tout, realtype *tret, 
                        N_Vector yret, N_Vector ypret, int itask);
static int IDAHandleFailure(IDAMem IDA_mem, int sflag);

/* Functions for rootfinding */

static int IDARcheck1(IDAMem IDA_mem);
static int IDARcheck2(IDAMem IDA_mem);
static int IDARcheck3(IDAMem IDA_mem);
static int IDARootfind(IDAMem IDA_mem);

/* Norm functions */

realtype IDAWrmsNorm(IDAMem IDA_mem, N_Vector x, N_Vector w, booleantype mask);

/* 
 * =================================================================
 * EXPORTED FUNCTIONS IMPLEMENTATION
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * Creation, allocation and re-initialization functions
 * -----------------------------------------------------------------
 */

/* 
 * IDACreate
 *
 * IDACreate creates an internal memory block for a problem to 
 * be solved by IDA.
 * If successful, IDACreate returns a pointer to the problem memory. 
 * This pointer should be passed to IDAInit.  
 * If an initialization error occurs, IDACreate prints an error 
 * message to standard err and returns NULL. 
 */

void *IDACreate(void)
{
  IDAMem IDA_mem;

  IDA_mem = NULL;
  IDA_mem = (IDAMem) malloc(sizeof(struct IDAMemRec));
  if (IDA_mem == NULL) {
    IDAProcessError(NULL, 0, "IDA", "IDACreate", MSG_MEM_FAIL);
    return (NULL);
  }

  /* Zero out ida_mem */
  memset(IDA_mem, 0, sizeof(struct IDAMemRec));

  /* Set unit roundoff in IDA_mem */
  IDA_mem->ida_uround = UNIT_ROUNDOFF;

  /* Set default values for integrator optional inputs */
  IDA_mem->ida_res         = NULL;
  IDA_mem->ida_user_data   = NULL;
  IDA_mem->ida_itol        = IDA_NN;
  IDA_mem->ida_user_efun   = FALSE;
  IDA_mem->ida_efun        = NULL;
  IDA_mem->ida_edata       = NULL;
  IDA_mem->ida_ehfun       = IDAErrHandler;
  IDA_mem->ida_eh_data     = IDA_mem;
  IDA_mem->ida_errfp       = stderr;
  IDA_mem->ida_maxord      = MAXORD_DEFAULT;
  IDA_mem->ida_mxstep      = MXSTEP_DEFAULT;
  IDA_mem->ida_hmax_inv    = HMAX_INV_DEFAULT;
  IDA_mem->ida_hin         = ZERO;
  IDA_mem->ida_epcon       = EPCON;
  IDA_mem->ida_maxnef      = MXNEF;
  IDA_mem->ida_maxncf      = MXNCF;
  IDA_mem->ida_maxcor      = MAXIT;
  IDA_mem->ida_suppressalg = FALSE;
  IDA_mem->ida_id          = NULL;
  IDA_mem->ida_constraints = NULL;
  IDA_mem->ida_constraintsSet = FALSE;
  IDA_mem->ida_tstopset    = FALSE;

  /* set the saved value maxord_alloc */
  IDA_mem->ida_maxord_alloc = MAXORD_DEFAULT;

  /* Set default values for IC optional inputs */
  IDA_mem->ida_epiccon = PT01 * EPCON;
  IDA_mem->ida_maxnh   = MAXNH;
  IDA_mem->ida_maxnj   = MAXNJ;
  IDA_mem->ida_maxnit  = MAXNI;
  IDA_mem->ida_maxbacks  = MAXBACKS;
  IDA_mem->ida_lsoff   = FALSE;
  IDA_mem->ida_steptol = SUNRpowerR(IDA_mem->ida_uround, TWOTHIRDS);

  /* Initialize lrw and liw */
  IDA_mem->ida_lrw = 25 + 5*MXORDP1;
  IDA_mem->ida_liw = 38;

  /* No mallocs have been done yet */
  IDA_mem->ida_VatolMallocDone = FALSE;
  IDA_mem->ida_constraintsMallocDone = FALSE;
  IDA_mem->ida_idMallocDone = FALSE;
  IDA_mem->ida_MallocDone = FALSE;

  /* Return pointer to IDA memory block */
  return((void *)IDA_mem);
}

/*-----------------------------------------------------------------*/

#define lrw   (IDA_mem->ida_lrw)
#define liw   (IDA_mem->ida_liw)

/*-----------------------------------------------------------------*/

/*
 * IDAInit
 *
 * IDAInit allocates and initializes memory for a problem. All
 * problem specification inputs are checked for errors. If any
 * error occurs during initialization, it is reported to the 
 * error handler function.
 */

int IDAInit(void *ida_mem, IDAResFn res,
            realtype t0, N_Vector yy0, N_Vector yp0)
{
  IDAMem IDA_mem;
  booleantype nvectorOK, allocOK;
  long int lrw1, liw1;

  /* Check ida_mem */

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAInit", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;
  
  /* Check for legal input parameters */
  
  if (yy0 == NULL) { 
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAInit", MSG_Y0_NULL);
    return(IDA_ILL_INPUT); 
  }
  
  if (yp0 == NULL) { 
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAInit", MSG_YP0_NULL);
    return(IDA_ILL_INPUT); 
  }

  if (res == NULL) { 
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAInit", MSG_RES_NULL);
    return(IDA_ILL_INPUT); 
  }

  /* Test if all required vector operations are implemented */

  nvectorOK = IDACheckNvector(yy0);
  if (!nvectorOK) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAInit", MSG_BAD_NVECTOR);
    return(IDA_ILL_INPUT);
  }

  /* Set space requirements for one N_Vector */

  if (yy0->ops->nvspace != NULL) {
    N_VSpace(yy0, &lrw1, &liw1);
  } else {
    lrw1 = 0;
    liw1 = 0;
  }
  IDA_mem->ida_lrw1 = lrw1;
  IDA_mem->ida_liw1 = liw1;

  /* Allocate the vectors (using yy0 as a template) */

  allocOK = IDAAllocVectors(IDA_mem, yy0);
  if (!allocOK) {
    IDAProcessError(IDA_mem, IDA_MEM_FAIL, "IDA", "IDAInit", MSG_MEM_FAIL);
    return(IDA_MEM_FAIL);
  }
 
  /* All error checking is complete at this point */

  /* Copy the input parameters into IDA memory block */

  IDA_mem->ida_res = res;
  IDA_mem->ida_tn  = t0;

  /* Set the linear solver addresses to NULL */

  IDA_mem->ida_linit  = NULL;
  IDA_mem->ida_lsetup = NULL;
  IDA_mem->ida_lsolve = NULL;
  IDA_mem->ida_lperf  = NULL;
  IDA_mem->ida_lfree  = NULL;
  IDA_mem->ida_lmem   = NULL;

  /* Initialize the phi array */

  N_VScale(ONE, yy0, IDA_mem->ida_phi[0]);  
  N_VScale(ONE, yp0, IDA_mem->ida_phi[1]);  
 
  /* Initialize all the counters and other optional output values */

  IDA_mem->ida_nst     = 0;
  IDA_mem->ida_nre     = 0;
  IDA_mem->ida_ncfn    = 0;
  IDA_mem->ida_netf    = 0;
  IDA_mem->ida_nni     = 0;
  IDA_mem->ida_nsetups = 0;
  
  IDA_mem->ida_kused = 0;
  IDA_mem->ida_hused = ZERO;
  IDA_mem->ida_tolsf = ONE;

  IDA_mem->ida_nge = 0;

  IDA_mem->ida_irfnd = 0;

  /* Initialize root-finding variables */

  IDA_mem->ida_glo     = NULL;
  IDA_mem->ida_ghi     = NULL;
  IDA_mem->ida_grout   = NULL;
  IDA_mem->ida_iroots  = NULL;
  IDA_mem->ida_rootdir = NULL;
  IDA_mem->ida_gfun    = NULL;
  IDA_mem->ida_nrtfn   = 0;
  IDA_mem->ida_gactive  = NULL;
  IDA_mem->ida_mxgnull  = 1;

  /* Initial setup not done yet */

  IDA_mem->ida_SetupDone = FALSE;

  /* Problem memory has been successfully allocated */

  IDA_mem->ida_MallocDone = TRUE;

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

#define lrw1 (IDA_mem->ida_lrw1)
#define liw1 (IDA_mem->ida_liw1)

/*-----------------------------------------------------------------*/

/*
 * IDAReInit
 *
 * IDAReInit re-initializes IDA's memory for a problem, assuming
 * it has already beeen allocated in a prior IDAInit call.
 * All problem specification inputs are checked for errors.
 * The problem size Neq is assumed to be unchanged since the call
 * to IDAInit, and the maximum order maxord must not be larger.
 * If any error occurs during reinitialization, it is reported to
 * the error handler function.
 * The return value is IDA_SUCCESS = 0 if no errors occurred, or
 * a negative value otherwise.
 */

int IDAReInit(void *ida_mem,
              realtype t0, N_Vector yy0, N_Vector yp0)
{
  IDAMem IDA_mem;

  /* Check for legal input parameters */
  
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAReInit", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Check if problem was malloc'ed */
  
  if (IDA_mem->ida_MallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_MALLOC, "IDA", "IDAReInit", MSG_NO_MALLOC);
    return(IDA_NO_MALLOC);
  }

  /* Check for legal input parameters */
  
  if (yy0 == NULL) { 
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAReInit", MSG_Y0_NULL);
    return(IDA_ILL_INPUT); 
  }
  
  if (yp0 == NULL) { 
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAReInit", MSG_YP0_NULL);
    return(IDA_ILL_INPUT); 
  }

  /* Copy the input parameters into IDA memory block */

  IDA_mem->ida_tn  = t0;

  /* Initialize the phi array */

  N_VScale(ONE, yy0, IDA_mem->ida_phi[0]);  
  N_VScale(ONE, yp0, IDA_mem->ida_phi[1]);  
 
  /* Initialize all the counters and other optional output values */
 
  IDA_mem->ida_nst     = 0;
  IDA_mem->ida_nre     = 0;
  IDA_mem->ida_ncfn    = 0;
  IDA_mem->ida_netf    = 0;
  IDA_mem->ida_nni     = 0;
  IDA_mem->ida_nsetups = 0;
  
  IDA_mem->ida_kused = 0;
  IDA_mem->ida_hused = ZERO;
  IDA_mem->ida_tolsf = ONE;

  IDA_mem->ida_nge = 0;

  IDA_mem->ida_irfnd = 0;

  /* Initial setup not done yet */

  IDA_mem->ida_SetupDone = FALSE;
      
  /* Problem has been successfully re-initialized */

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

/*
 * IDASStolerances
 * IDASVtolerances
 * IDAWFtolerances
 *
 * These functions specify the integration tolerances. One of them
 * MUST be called before the first call to IDA.
 *
 * IDASStolerances specifies scalar relative and absolute tolerances.
 * IDASVtolerances specifies scalar relative tolerance and a vector
 *   absolute tolerance (a potentially different absolute tolerance 
 *   for each vector component).
 * IDAWFtolerances specifies a user-provides function (of type IDAEwtFn)
 *   which will be called to set the error weight vector.
 */

int IDASStolerances(void *ida_mem, realtype reltol, realtype abstol)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASStolerances", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_MallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_MALLOC, "IDA", "IDASStolerances", MSG_NO_MALLOC);
    return(IDA_NO_MALLOC);
  }

  /* Check inputs */

  if (reltol < ZERO) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASStolerances", MSG_BAD_RTOL);
    return(IDA_ILL_INPUT);
  }

  if (abstol < ZERO) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASStolerances", MSG_BAD_ATOL);
    return(IDA_ILL_INPUT);
  }

  /* Copy tolerances into memory */
  
  IDA_mem->ida_rtol  = reltol;
  IDA_mem->ida_Satol = abstol;

  IDA_mem->ida_itol = IDA_SS;

  IDA_mem->ida_user_efun = FALSE;
  IDA_mem->ida_efun = IDAEwtSet;
  IDA_mem->ida_edata = NULL; /* will be set to ida_mem in InitialSetup; */

  return(IDA_SUCCESS);
}


int IDASVtolerances(void *ida_mem, realtype reltol, N_Vector abstol)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASVtolerances", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_MallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_MALLOC, "IDA", "IDASVtolerances", MSG_NO_MALLOC);
    return(IDA_NO_MALLOC);
  }

  /* Check inputs */

  if (reltol < ZERO) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASVtolerances", MSG_BAD_RTOL);
    return(IDA_ILL_INPUT);
  }

  if (N_VMin(abstol) < ZERO) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASVtolerances", MSG_BAD_ATOL);
    return(IDA_ILL_INPUT);
  }

  /* Copy tolerances into memory */
  
  if ( !(IDA_mem->ida_VatolMallocDone) ) {
    IDA_mem->ida_Vatol = N_VClone(IDA_mem->ida_ewt);
    lrw += lrw1;
    liw += liw1;
    IDA_mem->ida_VatolMallocDone = TRUE;
  }

  IDA_mem->ida_rtol = reltol;
  N_VScale(ONE, abstol, IDA_mem->ida_Vatol);

  IDA_mem->ida_itol = IDA_SV;

  IDA_mem->ida_user_efun = FALSE;
  IDA_mem->ida_efun = IDAEwtSet;
  IDA_mem->ida_edata = NULL; /* will be set to ida_mem in InitialSetup; */

  return(IDA_SUCCESS);
}


int IDAWFtolerances(void *ida_mem, IDAEwtFn efun)
{
  IDAMem IDA_mem;

  if (ida_mem==NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAWFtolerances", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  if (IDA_mem->ida_MallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_MALLOC, "IDA", "IDAWFtolerances", MSG_NO_MALLOC);
    return(IDA_NO_MALLOC);
  }

  IDA_mem->ida_itol = IDA_WF;

  IDA_mem->ida_user_efun = TRUE;
  IDA_mem->ida_efun = efun;
  IDA_mem->ida_edata = NULL; /* will be set to user_data in InitialSetup */ 

  return(IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

#define gfun    (IDA_mem->ida_gfun)
#define glo     (IDA_mem->ida_glo)
#define ghi     (IDA_mem->ida_ghi)
#define grout   (IDA_mem->ida_grout)
#define iroots  (IDA_mem->ida_iroots)
#define rootdir (IDA_mem->ida_rootdir)
#define gactive (IDA_mem->ida_gactive)

/*-----------------------------------------------------------------*/

/*
 * IDARootInit
 *
 * IDARootInit initializes a rootfinding problem to be solved
 * during the integration of the DAE system.  It loads the root
 * function pointer and the number of root functions, and allocates
 * workspace memory.  The return value is IDA_SUCCESS = 0 if no
 * errors occurred, or a negative value otherwise.
 */

int IDARootInit(void *ida_mem, int nrtfn, IDARootFn g)
{
  IDAMem IDA_mem;
  int i, nrt;

  /* Check ida_mem pointer */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDARootInit", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  nrt = (nrtfn < 0) ? 0 : nrtfn;

  /* If rerunning IDARootInit() with a different number of root
     functions (changing number of gfun components), then free
     currently held memory resources */
  if ((nrt != IDA_mem->ida_nrtfn) && (IDA_mem->ida_nrtfn > 0)) {

    free(glo); glo = NULL;
    free(ghi); ghi = NULL;
    free(grout); grout = NULL;
    free(iroots); iroots = NULL;
    free(rootdir); iroots = NULL;
    free(gactive); gactive = NULL;

    lrw -= 3 * (IDA_mem->ida_nrtfn);
    liw -= 3 * (IDA_mem->ida_nrtfn);

  }

  /* If IDARootInit() was called with nrtfn == 0, then set ida_nrtfn to
     zero and ida_gfun to NULL before returning */
  if (nrt == 0) {
    IDA_mem->ida_nrtfn = nrt;
    gfun = NULL;
    return(IDA_SUCCESS);
  }

  /* If rerunning IDARootInit() with the same number of root functions
     (not changing number of gfun components), then check if the root
     function argument has changed */
  /* If g != NULL then return as currently reserved memory resources
     will suffice */
  if (nrt == IDA_mem->ida_nrtfn) {
    if (g != gfun) {
      if (g == NULL) {
	free(glo); glo = NULL;
	free(ghi); ghi = NULL;
	free(grout); grout = NULL;
	free(iroots); iroots = NULL;
        free(rootdir); iroots = NULL;
        free(gactive); gactive = NULL;

        lrw -= 3*nrt;
        liw -= 3*nrt;

        IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDARootInit", MSG_ROOT_FUNC_NULL);
        return(IDA_ILL_INPUT);
      }
      else {
        gfun = g;
        return(IDA_SUCCESS);
      }
    }
    else return(IDA_SUCCESS);
  }

  /* Set variable values in IDA memory block */
  IDA_mem->ida_nrtfn = nrt;
  if (g == NULL) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDARootInit", MSG_ROOT_FUNC_NULL);
    return(IDA_ILL_INPUT);
  }
  else gfun = g;

  /* Allocate necessary memory and return */
  glo = NULL;
  glo = (realtype *) malloc(nrt*sizeof(realtype));
  if (glo == NULL) {
    IDAProcessError(IDA_mem, IDA_MEM_FAIL, "IDA", "IDARootInit", MSG_MEM_FAIL);
    return(IDA_MEM_FAIL);
  }

  ghi = NULL;
  ghi = (realtype *) malloc(nrt*sizeof(realtype));
  if (ghi == NULL) {
    free(glo); glo = NULL;
    IDAProcessError(IDA_mem, IDA_MEM_FAIL, "IDA", "IDARootInit", MSG_MEM_FAIL);
    return(IDA_MEM_FAIL);
  }

  grout = NULL;
  grout = (realtype *) malloc(nrt*sizeof(realtype));
  if (grout == NULL) {
    free(glo); glo = NULL;
    free(ghi); ghi = NULL;
    IDAProcessError(IDA_mem, IDA_MEM_FAIL, "IDA", "IDARootInit", MSG_MEM_FAIL);
    return(IDA_MEM_FAIL);
  }

  iroots = NULL;
  iroots = (int *) malloc(nrt*sizeof(int));
  if (iroots == NULL) {
    free(glo); glo = NULL;
    free(ghi); ghi = NULL;
    free(grout); grout = NULL;
    IDAProcessError(IDA_mem, IDA_MEM_FAIL, "IDA", "IDARootInit", MSG_MEM_FAIL);
    return(IDA_MEM_FAIL);
  }

  rootdir = NULL;
  rootdir = (int *) malloc(nrt*sizeof(int));
  if (rootdir == NULL) {
    free(glo); glo = NULL;
    free(ghi); ghi = NULL;
    free(grout); grout = NULL;
    free(iroots); iroots = NULL;
    IDAProcessError(IDA_mem, IDA_MEM_FAIL, "IDA", "IDARootInit", MSG_MEM_FAIL);
    return(IDA_MEM_FAIL);
  }

  gactive = NULL;
  gactive = (booleantype *) malloc(nrt*sizeof(booleantype));
  if (gactive == NULL) {
    free(glo); glo = NULL; 
    free(ghi); ghi = NULL;
    free(grout); grout = NULL;
    free(iroots); iroots = NULL;
    free(rootdir); rootdir = NULL;
    IDAProcessError(IDA_mem, IDA_MEM_FAIL, "IDA", "IDARootInit", MSG_MEM_FAIL);
    return(IDA_MEM_FAIL);
  }

  /* Set default values for rootdir (both directions) */
  for(i=0; i<nrt; i++) rootdir[i] = 0;

  /* Set default values for gactive (all active) */
  for(i=0; i<nrt; i++) gactive[i] = TRUE;

  lrw += 3*nrt;
  liw += 3*nrt;

  return(IDA_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * readability constants
 * -----------------------------------------------------------------
 */

#define res (IDA_mem->ida_res)
#define y0  (IDA_mem->ida_y0)
#define yp0 (IDA_mem->ida_yp0)

#define itol  (IDA_mem->ida_itol)
#define rtol  (IDA_mem->ida_rtol)
#define Satol (IDA_mem->ida_Satol)
#define Vatol (IDA_mem->ida_Vatol)
#define efun  (IDA_mem->ida_efun)
#define edata (IDA_mem->ida_edata)

#define user_data   (IDA_mem->ida_user_data)
#define maxord      (IDA_mem->ida_maxord)
#define mxstep      (IDA_mem->ida_mxstep)
#define hin         (IDA_mem->ida_hin)
#define hmax_inv    (IDA_mem->ida_hmax_inv)
#define tstop       (IDA_mem->ida_tstop)
#define tstopset    (IDA_mem->ida_tstopset)
#define epcon       (IDA_mem->ida_epcon)
#define maxnef      (IDA_mem->ida_maxnef)
#define maxncf      (IDA_mem->ida_maxncf)
#define maxcor      (IDA_mem->ida_maxcor)
#define suppressalg (IDA_mem->ida_suppressalg)
#define id          (IDA_mem->ida_id)
#define constraints (IDA_mem->ida_constraints)

#define epiccon (IDA_mem->ida_epiccon)
#define maxnh   (IDA_mem->ida_maxnh)
#define maxnj   (IDA_mem->ida_maxnj)
#define maxnit  (IDA_mem->ida_maxnit)
#define lsoff   (IDA_mem->ida_lsoff)
#define steptol (IDA_mem->ida_steptol)

#define uround         (IDA_mem->ida_uround)  
#define phi            (IDA_mem->ida_phi) 
#define ewt            (IDA_mem->ida_ewt)  
#define yy             (IDA_mem->ida_yy)
#define yp             (IDA_mem->ida_yp)
#define delta          (IDA_mem->ida_delta)
#define mm             (IDA_mem->ida_mm)
#define ee             (IDA_mem->ida_ee)
#define savres         (IDA_mem->ida_savres)
#define tempv1         (IDA_mem->ida_tempv1)
#define tempv2         (IDA_mem->ida_tempv2) 
#define kk             (IDA_mem->ida_kk)
#define hh             (IDA_mem->ida_hh)
#define h0u            (IDA_mem->ida_h0u)
#define tn             (IDA_mem->ida_tn)
#define tretlast       (IDA_mem->ida_tretlast)
#define cj             (IDA_mem->ida_cj)
#define cjold          (IDA_mem->ida_cjold)
#define cjratio        (IDA_mem->ida_cjratio)
#define cjlast         (IDA_mem->ida_cjlast)
#define nbacktr        (IDA_mem->ida_nbacktr)
#define nst            (IDA_mem->ida_nst)
#define nre            (IDA_mem->ida_nre)
#define ncfn           (IDA_mem->ida_ncfn)
#define netf           (IDA_mem->ida_netf)
#define nni            (IDA_mem->ida_nni)
#define nsetups        (IDA_mem->ida_nsetups)

#define ns             (IDA_mem->ida_ns)
#define linit          (IDA_mem->ida_linit)
#define lsetup         (IDA_mem->ida_lsetup)
#define lsolve         (IDA_mem->ida_lsolve) 
#define lperf          (IDA_mem->ida_lperf)
#define lfree          (IDA_mem->ida_lfree) 
#define lmem           (IDA_mem->ida_lmem) 
#define knew           (IDA_mem->ida_knew)
#define kused          (IDA_mem->ida_kused)          
#define hused          (IDA_mem->ida_hused)         
#define tolsf          (IDA_mem->ida_tolsf)      
#define phase          (IDA_mem->ida_phase)
#define epsNewt        (IDA_mem->ida_epsNewt)
#define toldel         (IDA_mem->ida_toldel)
#define ss             (IDA_mem->ida_ss)
#define rr             (IDA_mem->ida_rr)
#define psi            (IDA_mem->ida_psi)
#define alpha          (IDA_mem->ida_alpha)
#define beta           (IDA_mem->ida_beta)
#define sigma          (IDA_mem->ida_sigma)
#define gamma          (IDA_mem->ida_gamma)
#define setupNonNull   (IDA_mem->ida_setupNonNull) 
#define constraintsSet (IDA_mem->ida_constraintsSet)
#define nrtfn          (IDA_mem->ida_nrtfn)
#define tlo            (IDA_mem->ida_tlo)
#define thi            (IDA_mem->ida_thi)
#define toutc          (IDA_mem->ida_toutc)
#define trout          (IDA_mem->ida_trout)
#define ttol           (IDA_mem->ida_ttol)
#define taskc          (IDA_mem->ida_taskc)
#define irfnd          (IDA_mem->ida_irfnd)
#define nge            (IDA_mem->ida_nge)

/* 
 * -----------------------------------------------------------------
 * Main solver function
 * -----------------------------------------------------------------
 */

/*
 * IDASolve
 *
 * This routine is the main driver of the IDA package. 
 *
 * It integrates over an independent variable interval defined by the user, 
 * by calling IDAStep to take internal independent variable steps.
 *
 * The first time that IDASolve is called for a successfully initialized
 * problem, it computes a tentative initial step size.
 *
 * IDASolve supports two modes, specified by itask:
 * In the IDA_NORMAL mode, the solver steps until it passes tout and then
 * interpolates to obtain y(tout) and yp(tout).
 * In the IDA_ONE_STEP mode, it takes one internal step and returns.
 *
 * IDASolve returns integer values corresponding to success and failure as below:
 *
 * successful returns: 
 *
 * IDA_SUCCESS        
 * IDA_TSTOP_RETURN   
 *
 * failed returns:
 *
 * IDA_ILL_INPUT
 * IDA_TOO_MUCH_WORK
 * IDA_MEM_NULL
 * IDA_TOO_MUCH_ACC
 * IDA_CONV_FAIL
 * IDA_LSETUP_FAIL
 * IDA_LSOLVE_FAIL    
 * IDA_CONSTR_FAIL
 * IDA_ERR_FAIL   
 * IDA_REP_RES_ERR
 * IDA_RES_FAIL
 */

int IDASolve(void *ida_mem, realtype tout, realtype *tret,
             N_Vector yret, N_Vector ypret, int itask)
{
  long int nstloc;
  int sflag, istate, ier, irfndp, ir;
  realtype tdist, troundoff, ypnorm, rh, nrm;
  IDAMem IDA_mem;
  booleantype inactive_roots;

  /* Check for legal inputs in all cases. */

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASolve", MSG_NO_MEM);
    return(IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Check if problem was malloc'ed */
  
  if (IDA_mem->ida_MallocDone == FALSE) {
    IDAProcessError(IDA_mem, IDA_NO_MALLOC, "IDA", "IDASolve", MSG_NO_MALLOC);
    return(IDA_NO_MALLOC);
  }

  /* Check for legal arguments */

  if (yret == NULL) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_YRET_NULL);
    return(IDA_ILL_INPUT);
  }
  yy = yret;  

  if (ypret == NULL) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_YPRET_NULL);
    return(IDA_ILL_INPUT);
  }
  yp = ypret;
  
  if (tret == NULL) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_TRET_NULL);
    return(IDA_ILL_INPUT);
  }

  if ((itask != IDA_NORMAL) && (itask != IDA_ONE_STEP)) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_BAD_ITASK);
    return(IDA_ILL_INPUT);
  }
  
  if (itask == IDA_NORMAL) toutc = tout;
  taskc = itask;

  if (nst == 0) {       /* This is the first call */

    /* Check inputs to IDA for correctness and consistency */

    if (IDA_mem->ida_SetupDone == FALSE) {
      ier = IDAInitialSetup(IDA_mem);
      if (ier != IDA_SUCCESS) return(IDA_ILL_INPUT);
      IDA_mem->ida_SetupDone = TRUE;
    }

    /* On first call, check for tout - tn too small, set initial hh,
       check for approach to tstop, and scale phi[1] by hh.
       Also check for zeros of root function g at and near t0.    */

    tdist = SUNRabs(tout - tn);
    if (tdist == ZERO) {
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_TOO_CLOSE);
      return(IDA_ILL_INPUT);
    }
    troundoff = TWO*uround*(SUNRabs(tn) + SUNRabs(tout));
    if (tdist < troundoff) {
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_TOO_CLOSE);
      return(IDA_ILL_INPUT);
    }

    hh = hin;
    if ( (hh != ZERO) && ((tout-tn)*hh < ZERO) ) {
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_BAD_HINIT);
      return(IDA_ILL_INPUT);
    }

    if (hh == ZERO) {
      hh = PT001*tdist;
      ypnorm = IDAWrmsNorm(IDA_mem, phi[1], ewt, suppressalg);
      if (ypnorm > HALF/hh) hh = HALF/ypnorm;
      if (tout < tn) hh = -hh;
    }

    rh = SUNRabs(hh)*hmax_inv;
    if (rh > ONE) hh /= rh;

    if (tstopset) {
      if ( (tstop - tn)*hh <= ZERO) {
        IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_BAD_TSTOP, tstop, tn);
        return(IDA_ILL_INPUT);
      }
      if ( (tn + hh - tstop)*hh > ZERO) 
        hh = (tstop - tn)*(ONE-FOUR*uround);
    }

    h0u = hh;
    kk = 0; kused = 0;  /* set in case of an error return before a step */

    /* Check for exact zeros of the root functions at or near t0. */
    if (nrtfn > 0) {
      ier = IDARcheck1(IDA_mem);
      if (ier == IDA_RTFUNC_FAIL) {
        IDAProcessError(IDA_mem, IDA_RTFUNC_FAIL, "IDA", "IDARcheck1", MSG_RTFUNC_FAILED, tn);
        return(IDA_RTFUNC_FAIL);
      }
    }

    N_VScale(hh, phi[1], phi[1]);  /* set phi[1] = hh*y' */

    /* Set the convergence test constants epsNewt and toldel */
    epsNewt = epcon;
    toldel = PT0001 * epsNewt;

  } /* end of first-call block. */

  /* Call lperf function and set nstloc for later performance testing. */

  if (lperf != NULL) lperf(IDA_mem, 0);
  nstloc = 0;

  /* If not the first call, perform all stopping tests. */

  if (nst > 0) {

    /* First, check for a root in the last step taken, other than the
       last root found, if any.  If itask = IDA_ONE_STEP and y(tn) was not
       returned because of an intervening root, return y(tn) now.     */

    if (nrtfn > 0) {

      irfndp = irfnd;
      
      ier = IDARcheck2(IDA_mem);

      if (ier == CLOSERT) {
        IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDARcheck2", MSG_CLOSE_ROOTS, tlo);
        return(IDA_ILL_INPUT);
      } else if (ier == IDA_RTFUNC_FAIL) {
        IDAProcessError(IDA_mem, IDA_RTFUNC_FAIL, "IDA", "IDARcheck2", MSG_RTFUNC_FAILED, tlo);
        return(IDA_RTFUNC_FAIL);
      } else if (ier == RTFOUND) {
        tretlast = *tret = tlo;
        return(IDA_ROOT_RETURN);
      }

      /* If tn is distinct from tretlast (within roundoff),
         check remaining interval for roots */
      troundoff = HUNDRED*uround*(SUNRabs(tn) + SUNRabs(hh));
      if ( SUNRabs(tn - tretlast) > troundoff ) {
        ier = IDARcheck3(IDA_mem);
        if (ier == IDA_SUCCESS) {     /* no root found */
          irfnd = 0;
          if ((irfndp == 1) && (itask == IDA_ONE_STEP)) {
            tretlast = *tret = tn;
            ier = IDAGetSolution(IDA_mem, tn, yret, ypret);
            return(IDA_SUCCESS);
          }
        } else if (ier == RTFOUND) {  /* a new root was found */
          irfnd = 1;
          tretlast = *tret = tlo;
          return(IDA_ROOT_RETURN);
        } else if (ier == IDA_RTFUNC_FAIL) {  /* g failed */
          IDAProcessError(IDA_mem, IDA_RTFUNC_FAIL, "IDA", "IDARcheck3", MSG_RTFUNC_FAILED, tlo);
          return(IDA_RTFUNC_FAIL);
        }
      }

    } /* end of root stop check */


  /* Now test for all other stop conditions. */

    istate = IDAStopTest1(IDA_mem, tout, tret, yret, ypret, itask);
    if (istate != CONTINUE_STEPS) return(istate);
  }

  /* Looping point for internal steps. */

  loop {
   
    /* Check for too many steps taken. */
    
    if ( (mxstep>0) && (nstloc >= mxstep) ) {
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_MAX_STEPS, tn);
      istate = IDA_TOO_MUCH_WORK;
      *tret = tretlast = tn;
      break; /* Here yy=yret and yp=ypret already have the current solution. */
    }

    /* Call lperf to generate warnings of poor performance. */

    if (lperf != NULL) lperf(IDA_mem, 1);

    /* Reset and check ewt (if not first call). */

    if (nst > 0) {

      ier = efun(phi[0], ewt, edata);

      if (ier != 0) {

        if (itol == IDA_WF) 
          IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_EWT_NOW_FAIL, tn);
        else
          IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_EWT_NOW_BAD, tn);
	
        istate = IDA_ILL_INPUT;
        ier = IDAGetSolution(IDA_mem, tn, yret, ypret);
        *tret = tretlast = tn;
        break;

      }

    }
    
    /* Check for too much accuracy requested. */
    
    nrm = IDAWrmsNorm(IDA_mem, phi[0], ewt, suppressalg);
    tolsf = uround * nrm;
    if (tolsf > ONE) {
      tolsf *= TEN;
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_TOO_MUCH_ACC, tn);
      istate = IDA_TOO_MUCH_ACC;
      *tret = tretlast = tn;
      if (nst > 0) ier = IDAGetSolution(IDA_mem, tn, yret, ypret);
      break;
    }

    /* Call IDAStep to take a step. */

    sflag = IDAStep(IDA_mem);

    /* Process all failed-step cases, and exit loop. */

    if (sflag != IDA_SUCCESS) {
      istate = IDAHandleFailure(IDA_mem, sflag);
      *tret = tretlast = tn;
      ier = IDAGetSolution(IDA_mem, tn, yret, ypret);
      break;
    }
    
    nstloc++;

    /* After successful step, check for stop conditions; continue or break. */

    /* First check for root in the last step taken. */

    if (nrtfn > 0) {

      ier = IDARcheck3(IDA_mem);

      if (ier == RTFOUND) {  /* A new root was found */
        irfnd = 1;
        istate = IDA_ROOT_RETURN;
        tretlast = *tret = tlo;
        break;
      } else if (ier == IDA_RTFUNC_FAIL) { /* g failed */
        IDAProcessError(IDA_mem, IDA_RTFUNC_FAIL, "IDA", "IDARcheck3", MSG_RTFUNC_FAILED, tlo);
        istate = IDA_RTFUNC_FAIL;
        break;
      }

      /* If we are at the end of the first step and we still have
       * some event functions that are inactive, issue a warning
       * as this may indicate a user error in the implementation
       * of the root function. */

      if (nst==1) {
        inactive_roots = FALSE;
        for (ir=0; ir<nrtfn; ir++) { 
          if (!gactive[ir]) {
            inactive_roots = TRUE;
            break;
          }
        }
        if ((IDA_mem->ida_mxgnull > 0) && inactive_roots) {
          IDAProcessError(IDA_mem, IDA_WARNING, "IDA", "IDASolve", MSG_INACTIVE_ROOTS);
        }
      }

    }

    /* Now check all other stop conditions. */

    istate = IDAStopTest2(IDA_mem, tout, tret, yret, ypret, itask);
    if (istate != CONTINUE_STEPS) break;

  } /* End of step loop */

  return(istate);    
}

/* 
 * -----------------------------------------------------------------
 * Interpolated output
 * -----------------------------------------------------------------
 */

/* 
 * IDAGetDky
 *
 * This routine evaluates the k-th derivative of y(t) as the value of 
 * the k-th derivative of the interpolating polynomial at the independent 
 * variable t, and stores the results in the vector dky.  It uses the current
 * independent variable value, tn, and the method order last used, kused.
 * 
 * The return values are:
 *   IDA_SUCCESS  if t is legal, or
 *   IDA_BAD_T    if t is not within the interval of the last step taken.
 *   IDA_BAD_DKY  if the dky vector is NULL.
 *   IDA_BAD_K    if the requested k is not in the range 0,1,...,order used 
 *
 */

int IDAGetDky(void *ida_mem, realtype t, int k, N_Vector dky)
{
  IDAMem IDA_mem;
  realtype tfuzz, tp, delt, psij_1;
  int i, j;
  realtype cjk  [MXORDP1];
  realtype cjk_1[MXORDP1];

  /* Check ida_mem */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetDky", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  if (dky == NULL) {
    IDAProcessError(IDA_mem, IDA_BAD_DKY, "IDA", "IDAGetDky", MSG_NULL_DKY);
    return(IDA_BAD_DKY);
  }
  
  if ((k < 0) || (k > kused)) {
    IDAProcessError(IDA_mem, IDA_BAD_K, "IDA", "IDAGetDky", MSG_BAD_K);
    return(IDA_BAD_K);
  }

  /* Check t for legality.  Here tn - hused is t_{n-1}. */

  tfuzz = HUNDRED * uround * (SUNRabs(tn) + SUNRabs(hh));
  if (hh < ZERO) tfuzz = - tfuzz;
  tp = tn - hused - tfuzz;
  if ((t - tp)*hh < ZERO) {
    IDAProcessError(IDA_mem, IDA_BAD_T, "IDA", "IDAGetDky", MSG_BAD_T, t, tn-hused, tn);
    return(IDA_BAD_T);
  }

  /* Initialize the c_j^(k) and c_k^(k-1) */
  for(i=0; i<MXORDP1; i++) {
    cjk  [i] = 0;
    cjk_1[i] = 0;
  }

  delt = t-tn;

  for(i=0; i<=k; i++) {

    /* The below reccurence is used to compute the k-th derivative of the solution:
       c_j^(k) = ( k * c_{j-1}^(k-1) + c_{j-1}^{k} (Delta+psi_{j-1}) ) / psi_j
       
       Translated in indexes notation:
       cjk[j] = ( k*cjk_1[j-1] + cjk[j-1]*(delt+psi[j-2]) ) / psi[j-1]

       For k=0, j=1: c_1 = c_0^(-1) + (delt+psi[-1]) / psi[0]

       In order to be able to deal with k=0 in the same way as for k>0, the
       following conventions were adopted:
         - c_0(t) = 1 , c_0^(-1)(t)=0 
         - psij_1 stands for psi[-1]=0 when j=1 
                         for psi[j-2]  when j>1
    */
    if(i==0) {

      cjk[i] = 1;
      psij_1 = 0;
    }else {
      /*                                                i       i-1          1
        c_i^(i) can be always updated since c_i^(i) = -----  --------  ... -----
                                                      psi_j  psi_{j-1}     psi_1
      */
      cjk[i] = cjk[i-1]*i/psi[i-1];
      psij_1 = psi[i-1];
    }

    /* update c_j^(i) */

    /*j does not need to go till kused */
    for(j=i+1; j<=kused-k+i; j++) {

      cjk[j] = ( i* cjk_1[j-1] + cjk[j-1] * (delt + psij_1) ) / psi[j-1];      
      psij_1 = psi[j-1];
    }

    /* save existing c_j^(i)'s */
    for(j=i+1; j<=kused-k+i; j++) cjk_1[j] = cjk[j];
  }

  /* Compute sum (c_j(t) * phi(t)) */

  N_VConst(ZERO, dky);
  for(j=k; j<=kused; j++)
  {
    N_VLinearSum(ONE, dky, cjk[j], phi[j], dky);
  }

  return(IDA_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * Deallocation function
 * -----------------------------------------------------------------
 */

/*
 * IDAFree
 *
 * This routine frees the problem memory allocated by IDAInit
 * Such memory includes all the vectors allocated by IDAAllocVectors,
 * and the memory lmem for the linear solver (deallocated by a call
 * to lfree).
 */

void IDAFree(void **ida_mem)
{
  IDAMem IDA_mem;

  if (*ida_mem == NULL) return;

  IDA_mem = (IDAMem) (*ida_mem);
  
  IDAFreeVectors(IDA_mem);

  if (lfree != NULL) lfree(IDA_mem);

  if (nrtfn > 0) {
    free(glo); glo = NULL; 
    free(ghi);  ghi = NULL;
    free(grout);  grout = NULL;
    free(iroots); iroots = NULL;
    free(rootdir); rootdir = NULL;
    free(gactive); gactive = NULL;
  }

  free(*ida_mem);
  *ida_mem = NULL;
}

/* 
 * =================================================================
 * PRIVATE FUNCTIONS
 * =================================================================
 */

/*
 * IDACheckNvector
 *
 * This routine checks if all required vector operations are present.
 * If any of them is missing it returns FALSE.
 */

static booleantype IDACheckNvector(N_Vector tmpl)
{
  if ((tmpl->ops->nvclone        == NULL) ||
     (tmpl->ops->nvdestroy      == NULL) ||
     (tmpl->ops->nvlinearsum    == NULL) ||
     (tmpl->ops->nvconst        == NULL) ||
     (tmpl->ops->nvprod         == NULL) ||
     (tmpl->ops->nvscale        == NULL) ||
     (tmpl->ops->nvabs          == NULL) ||
     (tmpl->ops->nvinv          == NULL) ||
     (tmpl->ops->nvaddconst     == NULL) ||
     (tmpl->ops->nvwrmsnorm     == NULL) ||
     (tmpl->ops->nvmin          == NULL))
    return(FALSE);
  else
    return(TRUE);
}

/* 
 * -----------------------------------------------------------------
 * Memory allocation/deallocation
 * -----------------------------------------------------------------
 */

/*
 * IDAAllocVectors
 *
 * This routine allocates the IDA vectors ewt, tempv1, tempv2, and
 * phi[0], ..., phi[maxord].
 * If all memory allocations are successful, IDAAllocVectors returns 
 * TRUE. Otherwise all allocated memory is freed and IDAAllocVectors 
 * returns FALSE.
 * This routine also sets the optional outputs lrw and liw, which are
 * (respectively) the lengths of the real and integer work spaces
 * allocated here.
 */

static booleantype IDAAllocVectors(IDAMem IDA_mem, N_Vector tmpl)
{
  int i, j, maxcol;

  /* Allocate ewt, ee, delta, tempv1, tempv2 */
  
  ewt = N_VClone(tmpl);
  if (ewt == NULL) return(FALSE);

  ee = N_VClone(tmpl);
  if (ee == NULL) {
    N_VDestroy(ewt);
    return(FALSE);
  }

  delta = N_VClone(tmpl);
  if (delta == NULL) {
    N_VDestroy(ewt);
    N_VDestroy(ee);
    return(FALSE);
  }

  tempv1 = N_VClone(tmpl);
  if (tempv1 == NULL) {
    N_VDestroy(ewt);
    N_VDestroy(ee);
    N_VDestroy(delta);
    return(FALSE);
  }

  tempv2= N_VClone(tmpl);
  if (tempv2 == NULL) {
    N_VDestroy(ewt);
    N_VDestroy(ee);
    N_VDestroy(delta);
    N_VDestroy(tempv1);
    return(FALSE);
  }

  savres = tempv1;

  /* Allocate phi[0] ... phi[maxord].  Make sure phi[2] and phi[3] are
  allocated (for use as temporary vectors), regardless of maxord.       */

  maxcol = SUNMAX(maxord,3);
  for (j=0; j <= maxcol; j++) {
    phi[j] = N_VClone(tmpl);
    if (phi[j] == NULL) {
      N_VDestroy(ewt);
      N_VDestroy(ee);
      N_VDestroy(delta);
      N_VDestroy(tempv1);
      N_VDestroy(tempv2);
      for (i=0; i < j; i++) N_VDestroy(phi[i]);
      return(FALSE);
    }
  }

  /* Update solver workspace lengths  */
  lrw += (maxcol + 6)*lrw1;
  liw += (maxcol + 6)*liw1;

  /* Store the value of maxord used here */
  IDA_mem->ida_maxord_alloc = maxord;

  return(TRUE);
}

/*
 * IDAfreeVectors
 *
 * This routine frees the IDA vectors allocated for IDA.
 */

static void IDAFreeVectors(IDAMem IDA_mem)
{
  int j, maxcol;
  
  N_VDestroy(ewt);
  N_VDestroy(ee);
  N_VDestroy(delta);
  N_VDestroy(tempv1);
  N_VDestroy(tempv2);
  maxcol = SUNMAX(IDA_mem->ida_maxord_alloc,3);
  for(j=0; j <= maxcol; j++) N_VDestroy(phi[j]);

  lrw -= (maxcol + 6)*lrw1;
  liw -= (maxcol + 6)*liw1;

  if (IDA_mem->ida_VatolMallocDone) {
    N_VDestroy(Vatol);
    lrw -= lrw1;
    liw -= liw1;
  }

  if (IDA_mem->ida_constraintsMallocDone) {
    N_VDestroy(constraints);
    lrw -= lrw1;
    liw -= liw1;
  }

  if (IDA_mem->ida_idMallocDone) {
    N_VDestroy(id);
    lrw -= lrw1;
    liw -= liw1;
  }

}

/* 
 * -----------------------------------------------------------------
 * Initial setup
 * -----------------------------------------------------------------
 */

/*
 * IDAInitialSetup
 *
 * This routine is called by IDASolve once at the first step. 
 * It performs all checks on optional inputs and inputs to 
 * IDAInit/IDAReInit that could not be done before.
 *
 * If no error is encountered, IDAInitialSetup returns IDA_SUCCESS. 
 * Otherwise, it returns an error flag and reported to the error 
 * handler function.
 */

int IDAInitialSetup(IDAMem IDA_mem)
{
  booleantype conOK;
  int ier;

  /* Test for more vector operations, depending on options */
  if (suppressalg)
    if (phi[0]->ops->nvwrmsnormmask == NULL) {
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAInitialSetup", MSG_BAD_NVECTOR);
      return(IDA_ILL_INPUT);
  }

  /* Test id vector for legality */
  if (suppressalg && (id==NULL)){ 
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAInitialSetup", MSG_MISSING_ID);
    return(IDA_ILL_INPUT); 
  }

  /* Did the user specify tolerances? */
  if (itol == IDA_NN) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAInitialSetup", MSG_NO_TOLS);
    return(IDA_ILL_INPUT);
  }

  /* Set data for efun */
  if (IDA_mem->ida_user_efun) edata = user_data;
  else                        edata = IDA_mem;

  /* Initial error weight vector */
  ier = efun(phi[0], ewt, edata);
  if (ier != 0) {
    if (itol == IDA_WF) 
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAInitialSetup", MSG_FAIL_EWT);
    else
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAInitialSetup", MSG_BAD_EWT);
    return(IDA_ILL_INPUT);
  }

  /* Check to see if y0 satisfies constraints. */
  if (constraintsSet) {
    conOK = N_VConstrMask(constraints, phi[0], tempv2);
    if (!conOK) { 
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAInitialSetup", MSG_Y0_FAIL_CONSTR);
      return(IDA_ILL_INPUT); 
    }
  }

  /* Check that lsolve exists and call linit function if it exists. */
  if (lsolve == NULL) {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAInitialSetup", MSG_LSOLVE_NULL);
    return(IDA_ILL_INPUT);
  }

  if (linit != NULL) {
    ier = linit(IDA_mem);
    if (ier != 0) {
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAInitialSetup", MSG_LINIT_FAIL);
      return(IDA_LINIT_FAIL);
    }
  }

  return(IDA_SUCCESS);
}

/*  
 * IDAEwtSet
 *
 * This routine is responsible for loading the error weight vector
 * ewt, according to itol, as follows:
 * (1) ewt[i] = 1 / (rtol * SUNRabs(ycur[i]) + atol), i=0,...,Neq-1
 *     if itol = IDA_SS
 * (2) ewt[i] = 1 / (rtol * SUNRabs(ycur[i]) + atol[i]), i=0,...,Neq-1
 *     if itol = IDA_SV
 *
 *  IDAEwtSet returns 0 if ewt is successfully set as above to a
 *  positive vector and -1 otherwise. In the latter case, ewt is
 *  considered undefined.
 *
 * All the real work is done in the routines IDAEwtSetSS, IDAEwtSetSV.
 */

int IDAEwtSet(N_Vector ycur, N_Vector weight, void *data)
{
  IDAMem IDA_mem;
  int flag = 0;

  /* data points to IDA_mem here */

  IDA_mem = (IDAMem) data;

  switch(itol) {
  case IDA_SS: 
    flag = IDAEwtSetSS(IDA_mem, ycur, weight); 
    break;
  case IDA_SV: 
    flag = IDAEwtSetSV(IDA_mem, ycur, weight); 
    break;
  }
  return(flag);
}

/*
 * IDAEwtSetSS
 *
 * This routine sets ewt as decribed above in the case itol=IDA_SS.
 * It tests for non-positive components before inverting. IDAEwtSetSS
 * returns 0 if ewt is successfully set to a positive vector
 * and -1 otherwise. In the latter case, ewt is considered
 * undefined.
 */

static int IDAEwtSetSS(IDAMem IDA_mem, N_Vector ycur, N_Vector weight)
{
  N_VAbs(ycur, tempv1);
  N_VScale(rtol, tempv1, tempv1);
  N_VAddConst(tempv1, Satol, tempv1);
  if (N_VMin(tempv1) <= ZERO) return(-1);
  N_VInv(tempv1, weight);
  return(0);
}

/*
 * IDAEwtSetSV
 *
 * This routine sets ewt as decribed above in the case itol=IDA_SV.
 * It tests for non-positive components before inverting. IDAEwtSetSV
 * returns 0 if ewt is successfully set to a positive vector
 * and -1 otherwise. In the latter case, ewt is considered
 * undefined.
 */

static int IDAEwtSetSV(IDAMem IDA_mem, N_Vector ycur, N_Vector weight)
{
  N_VAbs(ycur, tempv1);
  N_VLinearSum(rtol, tempv1, ONE, Vatol, tempv1);
  if (N_VMin(tempv1) <= ZERO) return(-1);
  N_VInv(tempv1, weight);
  return(0);
}

/* 
 * -----------------------------------------------------------------
 * Stopping tests
 * -----------------------------------------------------------------
 */

/*
 * IDAStopTest1
 *
 * This routine tests for stop conditions before taking a step.
 * The tests depend on the value of itask.
 * The variable tretlast is the previously returned value of tret.
 *
 * The return values are:
 * CONTINUE_STEPS       if no stop conditions were found
 * IDA_SUCCESS          for a normal return to the user
 * IDA_TSTOP_RETURN     for a tstop-reached return to the user
 * IDA_ILL_INPUT        for an illegal-input return to the user 
 *
 * In the tstop cases, this routine may adjust the stepsize hh to cause
 * the next step to reach tstop exactly.
 */

static int IDAStopTest1(IDAMem IDA_mem, realtype tout, realtype *tret, 
                        N_Vector yret, N_Vector ypret, int itask)
{
  int ier;
  realtype troundoff;

  switch (itask) {
    
  case IDA_NORMAL:

    if (tstopset) {
      /* Test for tn past tstop, tn = tretlast, tn past tout, tn near tstop. */
      if ( (tn - tstop)*hh > ZERO) {
        IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_BAD_TSTOP, tstop, tn);
        return(IDA_ILL_INPUT);
      }
    }

    /* Test for tout = tretlast, and for tn past tout. */
    if (tout == tretlast) {
      *tret = tretlast = tout;
      return(IDA_SUCCESS);
    }
    if ((tn - tout)*hh >= ZERO) {
      ier = IDAGetSolution(IDA_mem, tout, yret, ypret);
      if (ier != IDA_SUCCESS) {
        IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_BAD_TOUT, tout);
        return(IDA_ILL_INPUT);
      }
      *tret = tretlast = tout;
      return(IDA_SUCCESS);
    }

    if (tstopset) {
      troundoff = HUNDRED*uround*(SUNRabs(tn) + SUNRabs(hh));
      if (SUNRabs(tn - tstop) <= troundoff) {
        ier = IDAGetSolution(IDA_mem, tstop, yret, ypret);
        if (ier != IDA_SUCCESS) {
          IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_BAD_TSTOP, tstop, tn);
          return(IDA_ILL_INPUT);
        }
        *tret = tretlast = tstop;
        tstopset = FALSE;
        return(IDA_TSTOP_RETURN);
      }
      if ((tn + hh - tstop)*hh > ZERO) 
        hh = (tstop - tn)*(ONE-FOUR*uround);
    }

    return(CONTINUE_STEPS);
    
  case IDA_ONE_STEP:

    if (tstopset) {
      /* Test for tn past tstop, tn past tretlast, and tn near tstop. */
      if ((tn - tstop)*hh > ZERO) {
        IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_BAD_TSTOP, tstop, tn);
        return(IDA_ILL_INPUT);
      }
    }

    /* Test for tn past tretlast. */
    if ((tn - tretlast)*hh > ZERO) {
      ier = IDAGetSolution(IDA_mem, tn, yret, ypret);
      *tret = tretlast = tn;
      return(IDA_SUCCESS);
    }

    if (tstopset) {
      troundoff = HUNDRED*uround*(SUNRabs(tn) + SUNRabs(hh));
      if (SUNRabs(tn - tstop) <= troundoff) {
        ier = IDAGetSolution(IDA_mem, tstop, yret, ypret);
        if (ier != IDA_SUCCESS) {
          IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASolve", MSG_BAD_TSTOP, tstop, tn);
          return(IDA_ILL_INPUT);
        }
        *tret = tretlast = tstop;
        tstopset = FALSE;
        return(IDA_TSTOP_RETURN);
      }
      if ((tn + hh - tstop)*hh > ZERO) 
        hh = (tstop - tn)*(ONE-FOUR*uround);
    }

    return(CONTINUE_STEPS);
        
  }
  return(IDA_ILL_INPUT);
}

/*
 * IDAStopTest2
 *
 * This routine tests for stop conditions after taking a step.
 * The tests depend on the value of itask.
 *
 * The return values are:
 *  CONTINUE_STEPS     if no stop conditions were found
 *  IDA_SUCCESS        for a normal return to the user
 *  IDA_TSTOP_RETURN   for a tstop-reached return to the user
 *
 * In the two cases with tstop, this routine may reset the stepsize hh
 * to cause the next step to reach tstop exactly.
 *
 * In the two cases with ONE_STEP mode, no interpolation to tn is needed
 * because yret and ypret already contain the current y and y' values.
 *
 * Note: No test is made for an error return from IDAGetSolution here,
 * because the same test was made prior to the step.
 */

static int IDAStopTest2(IDAMem IDA_mem, realtype tout, realtype *tret, 
                        N_Vector yret, N_Vector ypret, int itask)
{
  /* int ier; */
  realtype troundoff;

  switch (itask) {

    case IDA_NORMAL:  

      /* Test for tn past tout. */
      if ((tn - tout)*hh >= ZERO) {
        /* ier = */ IDAGetSolution(IDA_mem, tout, yret, ypret);
        *tret = tretlast = tout;
        return(IDA_SUCCESS);
      }

      if (tstopset) {
        /* Test for tn at tstop and for tn near tstop */
        troundoff = HUNDRED*uround*(SUNRabs(tn) + SUNRabs(hh));
        if (SUNRabs(tn - tstop) <= troundoff) {
          /* ier = */ IDAGetSolution(IDA_mem, tstop, yret, ypret);
          *tret = tretlast = tstop;
          tstopset = FALSE;
          return(IDA_TSTOP_RETURN);
        }
        if ((tn + hh - tstop)*hh > ZERO) 
          hh = (tstop - tn)*(ONE-FOUR*uround);
      }

      return(CONTINUE_STEPS);

    case IDA_ONE_STEP:

      if (tstopset) {
        /* Test for tn at tstop and for tn near tstop */
        troundoff = HUNDRED*uround*(SUNRabs(tn) + SUNRabs(hh));
        if (SUNRabs(tn - tstop) <= troundoff) {
          /* ier = */ IDAGetSolution(IDA_mem, tstop, yret, ypret);
          *tret = tretlast = tstop;
          tstopset = FALSE;
          return(IDA_TSTOP_RETURN);
        }
        if ((tn + hh - tstop)*hh > ZERO) 
          hh = (tstop - tn)*(ONE-FOUR*uround);
      }

      *tret = tretlast = tn;
      return(IDA_SUCCESS);

  }
  return IDA_ILL_INPUT;
}

/* 
 * -----------------------------------------------------------------
 * Error handler
 * -----------------------------------------------------------------
 */

/*
 * IDAHandleFailure
 *
 * This routine prints error messages for all cases of failure by
 * IDAStep.  It returns to IDASolve the value that it is to return to
 * the user.
 */

static int IDAHandleFailure(IDAMem IDA_mem, int sflag)
{
  /* Depending on sflag, print error message and return error flag */
  switch (sflag) {

    case IDA_ERR_FAIL:
      IDAProcessError(IDA_mem, IDA_ERR_FAIL, "IDA", "IDASolve", MSG_ERR_FAILS, tn, hh);
      return(IDA_ERR_FAIL);

    case IDA_CONV_FAIL:
      IDAProcessError(IDA_mem, IDA_CONV_FAIL, "IDA", "IDASolve", MSG_CONV_FAILS, tn, hh);
      return(IDA_CONV_FAIL);

    case IDA_LSETUP_FAIL:  
      IDAProcessError(IDA_mem, IDA_LSETUP_FAIL, "IDA", "IDASolve", MSG_SETUP_FAILED, tn);
      return(IDA_LSETUP_FAIL);

    case IDA_LSOLVE_FAIL: 
      IDAProcessError(IDA_mem, IDA_LSOLVE_FAIL, "IDA", "IDASolve", MSG_SOLVE_FAILED, tn);
      return(IDA_LSOLVE_FAIL);

    case IDA_REP_RES_ERR:
      IDAProcessError(IDA_mem, IDA_REP_RES_ERR, "IDA", "IDASolve", MSG_REP_RES_ERR, tn);
      return(IDA_REP_RES_ERR);

    case IDA_RES_FAIL: 
      IDAProcessError(IDA_mem, IDA_RES_FAIL, "IDA", "IDASolve", MSG_RES_NONRECOV, tn);
      return(IDA_RES_FAIL);

    case IDA_CONSTR_FAIL: 
      IDAProcessError(IDA_mem, IDA_CONSTR_FAIL, "IDA", "IDASolve", MSG_FAILED_CONSTR, tn);
      return(IDA_CONSTR_FAIL);

  }

  return (IDA_UNRECOGNISED_ERROR);
}

/* 
 * -----------------------------------------------------------------
 * Main IDAStep function
 * -----------------------------------------------------------------
 */

/*
 * IDAStep
 *
 * This routine performs one internal IDA step, from tn to tn + hh.
 * It calls other routines to do all the work.
 *
 * It solves a system of differential/algebraic equations of the form
 *       F(t,y,y') = 0, for one step. In IDA, tt is used for t,
 * yy is used for y, and yp is used for y'. The function F is supplied as 'res'
 * by the user.
 *
 * The methods used are modified divided difference, fixed leading 
 * coefficient forms of backward differentiation formulas.
 * The code adjusts the stepsize and order to control the local error per step.
 *
 * The main operations done here are as follows:
 *  * initialize various quantities;
 *  * setting of multistep method coefficients;
 *  * solution of the nonlinear system for yy at t = tn + hh;
 *  * deciding on order reduction and testing the local error;
 *  * attempting to recover from failure in nonlinear solver or error test;
 *  * resetting stepsize and order for the next step.
 *  * updating phi and other state data if successful;
 *
 * On a failure in the nonlinear system solution or error test, the
 * step may be reattempted, depending on the nature of the failure.
 *
 * Variables or arrays (all in the IDAMem structure) used in IDAStep are:
 *
 * tt -- Independent variable.
 * yy -- Solution vector at tt.
 * yp -- Derivative of solution vector after successful stelp.
 * res -- User-supplied function to evaluate the residual. See the 
 *        description given in file ida.h .
 * lsetup -- Routine to prepare for the linear solver call. It may either
 *        save or recalculate quantities used by lsolve. (Optional)
 * lsolve -- Routine to solve a linear system. A prior call to lsetup
 *        may be required. 
 * hh  -- Appropriate step size for next step.
 * ewt -- Vector of weights used in all convergence tests.
 * phi -- Array of divided differences used by IDAStep. This array is composed 
 *       of  (maxord+1) nvectors (each of size Neq). (maxord+1) is the maximum 
 *       order for the problem, maxord, plus 1.
 *
 *       Return values are:
 *       IDA_SUCCESS   IDA_RES_FAIL        LSETUP_ERROR_NONRECVR       
 *                     IDA_LSOLVE_FAIL   IDA_ERR_FAIL            
 *                     IDA_CONSTR_FAIL               IDA_CONV_FAIL          
 *                     IDA_REP_RES_ERR            
 */

static int IDAStep(IDAMem IDA_mem)
{
  realtype saved_t, ck;
  realtype err_k, err_km1;
  int ncf, nef;
  int nflag, kflag;

  saved_t = tn;
  ncf = nef = 0;

  if (nst == ZERO){
    kk = 1;
    kused = 0;
    hused = ZERO;
    psi[0] = hh;
    cj = ONE/hh;
    phase = 0;
    ns = 0;
  }

  /* To prevent 'unintialized variable' warnings */
  err_k = ZERO;
  err_km1 = ZERO;

  /* Looping point for attempts to take a step */

  loop {  

    /*-----------------------
      Set method coefficients
      -----------------------*/

    IDASetCoeffs(IDA_mem, &ck);

    kflag = IDA_SUCCESS;

    /*----------------------------------------------------
      If tn is past tstop (by roundoff), reset it to tstop.
      -----------------------------------------------------*/
    
    tn = tn + hh;
    if (tstopset) {
      if ((tn - tstop)*hh > ZERO) tn = tstop;
    }

    /*-----------------------
      Advance state variables
      -----------------------*/

    /* Nonlinear system solution */
    nflag = IDANls(IDA_mem);

    /* If NLS was successful, perform error test */
    if (nflag == IDA_SUCCESS)
      nflag = IDATestError(IDA_mem, ck, &err_k, &err_km1);

    /* Test for convergence or error test failures */
    if (nflag != IDA_SUCCESS) {

      /* restore and decide what to do */
      IDARestore(IDA_mem, saved_t);
      kflag = IDAHandleNFlag(IDA_mem, nflag, err_k, err_km1, 
                             &ncfn, &ncf, &netf, &nef);

      /* exit on nonrecoverable failure */ 
      if (kflag != PREDICT_AGAIN) return(kflag);

      /* recoverable error; predict again */
      if(nst==0) IDAReset(IDA_mem);
      continue;

    }

    /* kflag == IDA_SUCCESS */
    break;

  }

  /* Nonlinear system solve and error test were both successful;
     update data, and consider change of step and/or order */

  IDACompleteStep(IDA_mem, err_k, err_km1);

  /* 
     Rescale ee vector to be the estimated local error
     Notes:
       (1) altering the value of ee is permissible since
           it will be re-initialized to the zero vector by
           IDASolve()->IDAStep()->IDANls()->IDANewtonIter()
           before it is needed again
       (2) the value of ee is only valid if IDAHandleNFlag()
           returns either PREDICT_AGAIN or IDA_SUCCESS
   */

  N_VScale(ck, ee, ee);

  return(IDA_SUCCESS);
}

/*
 * IDASetCoeffs
 *
 *  This routine computes the coefficients relevant to the current step.
 *  The counter ns counts the number of consecutive steps taken at
 *  constant stepsize h and order k, up to a maximum of k + 2.
 *  Then the first ns components of beta will be one, and on a step  
 *  with ns = k + 2, the coefficients alpha, etc. need not be reset here.
 *  Also, IDACompleteStep prohibits an order increase until ns = k + 2.
 */

static void IDASetCoeffs(IDAMem IDA_mem, realtype *ck)
{
  int i;
  realtype temp1, temp2, alpha0, alphas;

  /* Set coefficients for the current stepsize h */

  if (hh != hused || kk != kused) ns = 0;
  ns = SUNMIN(ns+1,kused+2);
  if (kk+1 >= ns){
    beta[0] = ONE;
    alpha[0] = ONE;
    temp1 = hh;
    gamma[0] = ZERO;
    sigma[0] = ONE;
    for(i=1;i<=kk;i++){
      temp2 = psi[i-1];
      psi[i-1] = temp1;
      beta[i] = beta[i-1] * psi[i-1] / temp2;
      temp1 = temp2 + hh;
      alpha[i] = hh / temp1;
      sigma[i] = i * sigma[i-1] * alpha[i]; 
      gamma[i] = gamma[i-1] + alpha[i-1] / hh;
   }
    psi[kk] = temp1;
  }
  /* compute alphas, alpha0 */
  alphas = ZERO;
  alpha0 = ZERO;
  for(i=0;i<kk;i++){
    alphas = alphas - ONE/(i+1);
    alpha0 = alpha0 - alpha[i];
  }

  /* compute leading coefficient cj  */
  cjlast = cj;
  cj = -alphas/hh;
  
  /* compute variable stepsize error coefficient ck */

  *ck = SUNRabs(alpha[kk] + alphas - alpha0);
  *ck = SUNMAX(*ck, alpha[kk]);

 /* change phi to phi-star  */

  for(i=ns;i<=kk;i++) N_VScale(beta[i], phi[i], phi[i]);

}

/* 
 * -----------------------------------------------------------------
 * Nonlinear solver functions
 * -----------------------------------------------------------------
 */

/*
 * IDANls
 *
 * This routine attempts to solve the nonlinear system using the linear
 * solver specified. NOTE: this routine uses N_Vector ee as the scratch
 * vector tempv3 passed to lsetup.
 *
 *  Possible return values:
 *
 *  IDA_SUCCESS
 *
 *  IDA_RES_RECVR       IDA_RES_FAIL
 *  IDA_LSETUP_RECVR    IDA_LSETUP_FAIL
 *  IDA_LSOLVE_RECVR    IDA_LSOLVE_FAIL
 *
 *  IDA_CONSTR_RECVR
 *  IDA_NCONV_RECVR
 */

static int IDANls(IDAMem IDA_mem)
{
  int retval;
  booleantype constraintsPassed, callSetup, tryAgain;
  realtype temp1, temp2, vnorm;
  N_Vector tempv3;

  callSetup = FALSE;

  /* Initialize if the first time called */

  if (nst == 0){
    cjold = cj;
    ss = TWENTY;
    if (setupNonNull) callSetup = TRUE;
  }

  mm = tempv2;
  tempv3 = ee;

  /* Decide if lsetup is to be called */

  if (setupNonNull){
    cjratio = cj / cjold;
    temp1 = (ONE - XRATE) / (ONE + XRATE);
    temp2 = ONE/temp1;
    {if (cjratio < temp1 || cjratio > temp2) callSetup = TRUE;}
    {if (cj != cjlast) ss=HUNDRED;}
  }

  /* Begin the main loop. This loop is traversed at most twice. 
     The second pass only occurs when the first pass had a recoverable
     failure with old Jacobian data */
  loop{

    /* Compute predicted values for yy and yp, and compute residual there. */
    IDAPredict(IDA_mem);

    retval = res(tn, yy, yp, delta, user_data);
    nre++;
    if (retval < 0) return(IDA_RES_FAIL);
    if (retval > 0) return(IDA_RES_RECVR);

    /* If indicated, call linear solver setup function and reset parameters. */
    if (callSetup){
      nsetups++;
      retval = lsetup(IDA_mem, yy, yp, delta, tempv1, tempv2, tempv3);
      cjold = cj;
      cjratio = ONE;
      ss = TWENTY;
      if (retval < 0) return(IDA_LSETUP_FAIL);
      if (retval > 0) return(IDA_LSETUP_RECVR);
    }

    /* Call the Newton iteration routine.  */

    retval = IDANewtonIter(IDA_mem);

    /* Retry the current step on recoverable failure with old Jacobian data. */

    tryAgain = (retval>0)&&(setupNonNull) &&(!callSetup);

    if (tryAgain){
      callSetup = TRUE;
      continue;
    }
    else break;

  }  /* end of loop */

  if (retval != IDA_SUCCESS) return(retval);

  /* If otherwise successful, check and enforce inequality constraints. */

  if (constraintsSet){  /* Check constraints and get mask vector mm, 
                          set where constraints failed */
    constraintsPassed = N_VConstrMask(constraints,yy,mm);
    if (constraintsPassed) return(IDA_SUCCESS);
    else {
      N_VCompare(ONEPT5, constraints, tempv1);  
      /* a , where a[i] =1. when |c[i]| = 2 ,  c the vector of constraints */
      N_VProd(tempv1, constraints, tempv1);       /* a * c */
      N_VDiv(tempv1, ewt, tempv1);                /* a * c * wt */
      N_VLinearSum(ONE, yy, -PT1, tempv1, tempv1);/* y - 0.1 * a * c * wt */
      N_VProd(tempv1, mm, tempv1);               /*  v = mm*(y-.1*a*c*wt) */
      vnorm = IDAWrmsNorm(IDA_mem, tempv1, ewt, FALSE); /*  ||v|| */
      
      /* If vector v of constraint corrections is small
         in norm, correct and accept this step */      
      if (vnorm <= epsNewt){  
        N_VLinearSum(ONE, ee, -ONE, tempv1, ee);  /* ee <- ee - v */
        return(IDA_SUCCESS);
      }
      else {
        /* Constraints not met -- reduce h by computing rr = h'/h */
        N_VLinearSum(ONE, phi[0], -ONE, yy, tempv1);
        N_VProd(mm, tempv1, tempv1);
        rr = PT9*N_VMinQuotient(phi[0], tempv1);
        rr = SUNMAX(rr,PT1);
        return(IDA_CONSTR_RECVR);
      }
    }
  }
  return(IDA_SUCCESS);
}


/*
 * IDAPredict
 *
 * This routine predicts the new values for vectors yy and yp.
 */

static void IDAPredict(IDAMem IDA_mem)
{
  int j;

  N_VScale(ONE, phi[0], yy);
  N_VConst(ZERO, yp);
  
  for(j=1; j<=kk; j++) {
    N_VLinearSum(ONE,      phi[j], ONE, yy, yy);
    N_VLinearSum(gamma[j], phi[j], ONE, yp, yp);
  }
}

/*
 * IDANewtonIter
 *
 * This routine performs the Newton iteration.  
 * It assumes that delta contains the initial residual vector on entry.
 * If the iteration succeeds, it returns the value IDA_SUCCESS = 0.
 * If not, it returns either:
 *   a positive value (for a recoverable failure), namely one of:
 *     IDA_RES_RECVR
 *     IDA_LSOLVE_RECVR
 *     IDA_NCONV_RECVR
 * or
 *   a negative value (for a nonrecoverable failure), namely one of:
 *     IDA_RES_FAIL
 *     IDA_LSOLVE_FAIL
 *
 * NOTE: This routine uses N_Vector savres, which is preset to tempv1.
 */

static int IDANewtonIter(IDAMem IDA_mem)
{
  int mnewt, retval;
  realtype delnrm, oldnrm, rate;

  /* Initialize counter mnewt and cumulative correction vector ee. */
  mnewt = 0;
  N_VConst(ZERO, ee);

  /* Initialize oldnrm to avoid compiler warning message */
  oldnrm = ZERO;

  /* Looping point for Newton iteration.  Break out on any error. */
  loop {

    nni++;

    /* Save a copy of the residual vector in savres. */
    N_VScale(ONE, delta, savres);

    /* Call the lsolve function to get correction vector delta. */
    retval = lsolve(IDA_mem, delta, ewt, yy, yp, savres); 
    if (retval < 0) return(IDA_LSOLVE_FAIL);
    if (retval > 0) return(IDA_LSOLVE_RECVR);

    /* Apply delta to yy, yp, and ee, and get norm(delta). */
    N_VLinearSum(ONE, yy, -ONE, delta, yy);
    N_VLinearSum(ONE, ee, -ONE, delta, ee);
    N_VLinearSum(ONE, yp, -cj,  delta, yp);
    delnrm = IDAWrmsNorm(IDA_mem, delta, ewt, FALSE);

    /* Test for convergence, first directly, then with rate estimate. */

    if (mnewt == 0){ 
       oldnrm = delnrm;
       if (delnrm <= toldel) return(IDA_SUCCESS);
    }
    else {
      rate = SUNRpowerR( delnrm/oldnrm, ONE/mnewt );
      if (rate > RATEMAX) return(IDA_NCONV_RECVR); 
      ss = rate/(ONE - rate);
    }

    if (ss*delnrm <= epsNewt) return(IDA_SUCCESS);

    /* Not yet converged.  Increment mnewt and test for max allowed. */
    mnewt++;
    if (mnewt >= maxcor) {retval = IDA_NCONV_RECVR; break;}

    /* Call res for new residual and check error flag from res. */
    retval = res(tn, yy, yp, delta, user_data);
    nre++;
    if (retval < 0) return(IDA_RES_FAIL);
    if (retval > 0) return(IDA_RES_RECVR);

    /* Loop for next iteration. */

  } /* end of Newton iteration loop */

  /* All error returns exit here. */
  return(retval);

}

/* 
 * -----------------------------------------------------------------
 * Error test
 * -----------------------------------------------------------------
 */

/*
 * IDATestError
 *
 * This routine estimates errors at orders k, k-1, k-2, decides 
 * whether or not to suggest an order decrease, and performs 
 * the local error test. 
 *
 * IDATestError returns either IDA_SUCCESS or ERROR_TEST_FAIL.
 */

static int IDATestError(IDAMem IDA_mem, realtype ck, 
                        realtype *err_k, realtype *err_km1)
{
  realtype err_km2;                         /* estimated error at k-2 */
  realtype enorm_k, enorm_km1, enorm_km2;   /* error norms */
  realtype terr_k, terr_km1, terr_km2;      /* local truncation error norms */

  /* Compute error for order k. */
  enorm_k = IDAWrmsNorm(IDA_mem, ee, ewt, suppressalg);
  *err_k = sigma[kk] * enorm_k;
  terr_k = (kk+1) * (*err_k);

  knew = kk;

  if ( kk > 1 ) {

    /* Compute error at order k-1 */
    N_VLinearSum(ONE, phi[kk], ONE, ee, delta);
    enorm_km1 = IDAWrmsNorm(IDA_mem, delta, ewt, suppressalg);
    *err_km1 = sigma[kk-1] * enorm_km1;
    terr_km1 = kk * (*err_km1);

    if ( kk > 2 ) {

      /* Compute error at order k-2 */
      N_VLinearSum(ONE, phi[kk-1], ONE, delta, delta);
      enorm_km2 = IDAWrmsNorm(IDA_mem, delta, ewt, suppressalg);
      err_km2 = sigma[kk-2] * enorm_km2;
      terr_km2 = (kk-1) * err_km2;

      /* Decrease order if errors are reduced */
      if (SUNMAX(terr_km1, terr_km2) <= terr_k)  knew = kk - 1;

    } else {

      /* Decrease order to 1 if errors are reduced by at least 1/2 */
      if (terr_km1 <= (HALF * terr_k) )  knew = kk - 1; 

    }

  }

  /* Perform error test */
  if (ck * enorm_k > ONE) return(ERROR_TEST_FAIL);
  else                    return(IDA_SUCCESS);
}

/*
 * IDARestore
 *
 * This routine restores tn, psi, and phi in the event of a failure.
 * It changes back phi-star to phi (changed in IDASetCoeffs)
 */

static void IDARestore(IDAMem IDA_mem, realtype saved_t)
{
  int j;

  tn = saved_t;
  
  for (j = 1; j <= kk; j++) 
    psi[j-1] = psi[j] - hh;

  for (j = ns; j <= kk; j++) 
    N_VScale(ONE/beta[j], phi[j], phi[j]);

}

/* 
 * -----------------------------------------------------------------
 * Handler for convergence and/or error test failures
 * -----------------------------------------------------------------
 */

/*
 * IDAHandleNFlag
 *
 * This routine handles failures indicated by the input variable nflag. 
 * Positive values indicate various recoverable failures while negative
 * values indicate nonrecoverable failures. This routine adjusts the
 * step size for recoverable failures. 
 *
 *  Possible nflag values (input):
 *
 *   --convergence failures--
 *   IDA_RES_RECVR              > 0
 *   IDA_LSOLVE_RECVR           > 0
 *   IDA_CONSTR_RECVR           > 0
 *   IDA_NCONV_RECVR            > 0
 *   IDA_RES_FAIL               < 0
 *   IDA_LSOLVE_FAIL            < 0
 *   IDA_LSETUP_FAIL            < 0
 *
 *   --error test failure--
 *   ERROR_TEST_FAIL            > 0
 *
 *  Possible kflag values (output):
 *
 *   --recoverable--
 *   PREDICT_AGAIN
 *
 *   --nonrecoverable--
 *   IDA_CONSTR_FAIL   
 *   IDA_REP_RES_ERR    
 *   IDA_ERR_FAIL  
 *   IDA_CONV_FAIL 
 *   IDA_RES_FAIL
 *   IDA_LSETUP_FAIL
 *   IDA_LSOLVE_FAIL
 */

static int IDAHandleNFlag(IDAMem IDA_mem, int nflag, realtype err_k, realtype err_km1,
                          long int *ncfnPtr, int *ncfPtr, long int *netfPtr, int *nefPtr)
{
  realtype err_knew;

  phase = 1;
    
  if (nflag != ERROR_TEST_FAIL) {

    /*-----------------------
      Nonlinear solver failed 
      -----------------------*/

    (*ncfPtr)++;      /* local counter for convergence failures */
    (*ncfnPtr)++;     /* global counter for convergence failures */
    
    if (nflag < 0) {  /* nonrecoverable failure */

      return(nflag);

    } else {          /* recoverable failure    */
      
      /* Reduce step size for a new prediction
         Note that if nflag=IDA_CONSTR_RECVR then rr was already set in IDANls */
      if (nflag != IDA_CONSTR_RECVR) rr = QUARTER;
      hh *= rr;

      /* Test if there were too many convergence failures */
      if (*ncfPtr < maxncf)               return(PREDICT_AGAIN);
      else if (nflag == IDA_RES_RECVR)    return(IDA_REP_RES_ERR);
      else if (nflag == IDA_CONSTR_RECVR) return(IDA_CONSTR_FAIL);
      else                                return(IDA_CONV_FAIL);
    }
    
  } else { 

    /*-----------------
      Error Test failed 
      -----------------*/

    (*nefPtr)++;      /* local counter for error test failures */
    (*netfPtr)++;     /* global counter for error test failures */
    
    if (*nefPtr == 1) {
      
      /* On first error test failure, keep current order or lower order by one. 
         Compute new stepsize based on differences of the solution. */

      err_knew = (kk==knew)? err_k : err_km1;

      kk = knew;      
      rr = PT9 * SUNRpowerR( TWO * err_knew + PT0001,(-ONE/(kk+1)) );
      rr = SUNMAX(QUARTER, SUNMIN(PT9,rr));
      hh *=rr;
      return(PREDICT_AGAIN);
      
    } else if (*nefPtr == 2) {
      
      /* On second error test failure, use current order or decrease order by one. 
         Reduce stepsize by factor of 1/4. */

      kk = knew;
      rr = QUARTER;
      hh *= rr;
      return(PREDICT_AGAIN);
      
    } else if (*nefPtr < maxnef) {
      
      /* On third and subsequent error test failures, set order to 1.
         Reduce stepsize by factor of 1/4. */
      kk = 1;
      rr = QUARTER;
      hh *= rr;
      return(PREDICT_AGAIN);

    } else {

      /* Too many error test failures */
      return(IDA_ERR_FAIL);

    }    

  }
      
}

/*
 * IDAReset
 *
 * This routine is called only if we need to predict again at the 
 * very first step. In such a case, reset phi[1] and psi[0].
 */

static void IDAReset(IDAMem IDA_mem)
{
  psi[0] = hh;

  N_VScale(rr, phi[1], phi[1]);
}

/* 
 * -----------------------------------------------------------------
 * Function called after a successful step
 * -----------------------------------------------------------------
 */

/*
 * IDACompleteStep
 *
 * This routine completes a successful step.  It increments nst,
 * saves the stepsize and order used, makes the final selection of
 * stepsize and order for the next step, and updates the phi array.
 */

static void IDACompleteStep(IDAMem IDA_mem, realtype err_k, realtype err_km1)
{
  int j, kdiff, action;
  realtype terr_k, terr_km1, terr_kp1;
  realtype err_knew, err_kp1;
  realtype enorm, tmp, hnew;

  nst++;
  kdiff = kk - kused;
  kused = kk;
  hused = hh;

  if ( (knew == kk-1) || (kk == maxord) ) phase = 1;

  /* For the first few steps, until either a step fails, or the order is 
     reduced, or the order reaches its maximum, we raise the order and double 
     the stepsize. During these steps, phase = 0. Thereafter, phase = 1, and
     stepsize and order are set by the usual local error algorithm.         
     
     Note that, after the first step, the order is not increased, as not all 
     of the neccessary information is available yet. */
  
  if (phase == 0) {

    if(nst > 1) {
      kk++;
      hnew = TWO * hh;
      if( (tmp = SUNRabs(hnew)*hmax_inv) > ONE ) hnew /= tmp;
      hh = hnew;
    }

  } else {

    action = UNSET;
    
    /* Set action = LOWER/MAINTAIN/RAISE to specify order decision */
    
    if (knew == kk-1)                   {action = LOWER;    goto takeaction;}
    if (kk == maxord)                   {action = MAINTAIN; goto takeaction;}
    if ( (kk+1 >= ns ) || (kdiff == 1)) {action = MAINTAIN; goto takeaction;}
    
    /* Estimate the error at order k+1, unless already decided to
       reduce order, or already using maximum order, or stepsize has not
       been constant, or order was just raised. */
    
    N_VLinearSum (ONE, ee, -ONE, phi[kk+1], tempv1);
    enorm = IDAWrmsNorm(IDA_mem, tempv1, ewt, suppressalg);
    err_kp1= enorm/(kk+2);

    /* Choose among orders k-1, k, k+1 using local truncation error norms. */

    terr_k   = (kk+1) * err_k;
    terr_kp1 = (kk+2) * err_kp1;

    if (kk == 1) {
      if (terr_kp1 >= HALF * terr_k)         {action = MAINTAIN; goto takeaction;}
      else                                   {action = RAISE;    goto takeaction;}
    } else {
      terr_km1 = kk * err_km1;
      if (terr_km1 <= SUNMIN(terr_k, terr_kp1)) {action = LOWER;    goto takeaction;}
      else if (terr_kp1 >= terr_k)           {action = MAINTAIN; goto takeaction;}
      else                                   {action = RAISE;    goto takeaction;}
    }
    
  takeaction:
    
    /* Set the estimated error norm and, on change of order, reset kk. */
    if      (action == RAISE) { kk++; err_knew = err_kp1; }
    else if (action == LOWER) { kk--; err_knew = err_km1; }
    else                      {       err_knew = err_k;   }  

    /* Compute rr = tentative ratio hnew/hh from error norm estimate.
       Reduce hh if rr <= 1, double hh if rr >= 2, else leave hh as is.
       If hh is reduced, hnew/hh is restricted to be between .5 and .9. */
    
    hnew = hh;
    rr = SUNRpowerR( (TWO * err_knew + PT0001) , (-ONE/(kk+1) ) );
    
    if (rr >= TWO) {
      hnew = TWO * hh;
      if( (tmp = SUNRabs(hnew)*hmax_inv) > ONE ) hnew /= tmp;
    } else if (rr <= ONE ) { 
      rr = SUNMAX(HALF, SUNMIN(PT9,rr));
      hnew = hh * rr;
    }
    
    hh = hnew;
    
  } /* end of phase if block */
  
  /* Save ee for possible order increase on next step */
  if (kused < maxord) {
    N_VScale(ONE, ee, phi[kused+1]);
  }

  /* Update phi arrays */
  N_VLinearSum(ONE, ee, ONE, phi[kused], phi[kused]);
  for (j= kused-1; j>=0; j--)
    N_VLinearSum(ONE, phi[j], ONE, phi[j+1], phi[j]);

}

/* 
 * -----------------------------------------------------------------
 * Interpolated output
 * -----------------------------------------------------------------
 */

/* 
 * IDAGetSolution
 *
 * This routine evaluates y(t) and y'(t) as the value and derivative of 
 * the interpolating polynomial at the independent variable t, and stores
 * the results in the vectors yret and ypret.  It uses the current
 * independent variable value, tn, and the method order last used, kused.
 * This function is called by IDASolve with t = tout, t = tn, or t = tstop.
 * 
 * If kused = 0 (no step has been taken), or if t = tn, then the order used
 * here is taken to be 1, giving yret = phi[0], ypret = phi[1]/psi[0].
 * 
 * The return values are:
 *   IDA_SUCCESS  if t is legal, or
 *   IDA_BAD_T    if t is not within the interval of the last step taken.
 */

int IDAGetSolution(void *ida_mem, realtype t, N_Vector yret, N_Vector ypret)
{
  IDAMem IDA_mem;
  realtype tfuzz, tp, delt, c, d, gam;
  int j, kord;

  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetSolution", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem; 

  /* Check t for legality.  Here tn - hused is t_{n-1}. */
 
  tfuzz = HUNDRED * uround * (SUNRabs(tn) + SUNRabs(hh));
  if (hh < ZERO) tfuzz = - tfuzz;
  tp = tn - hused - tfuzz;
  if ((t - tp)*hh < ZERO) {
    IDAProcessError(IDA_mem, IDA_BAD_T, "IDA", "IDAGetSolution", MSG_BAD_T, t, tn-hused, tn);
    return(IDA_BAD_T);
  }

  /* Initialize yret = phi[0], ypret = 0, and kord = (kused or 1). */

  N_VScale (ONE, phi[0], yret);
  N_VConst (ZERO, ypret);
  kord = kused; 
  if (kused == 0) kord = 1;

  /* Accumulate multiples of columns phi[j] into yret and ypret. */

  delt = t - tn;
  c = ONE; d = ZERO;
  gam = delt/psi[0];
  for (j=1; j <= kord; j++) {
    d = d*gam + c/psi[j-1];
    c = c*gam;
    gam = (delt + psi[j-1])/psi[j];
    N_VLinearSum(ONE,  yret, c, phi[j],  yret);
    N_VLinearSum(ONE, ypret, d, phi[j], ypret);
  }
  return(IDA_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * Norm function
 * -----------------------------------------------------------------
 */

/*
 * IDAWrmsNorm
 *
 *  Returns the WRMS norm of vector x with weights w.
 *  If mask = TRUE, the weight vector w is masked by id, i.e.,
 *      nrm = N_VWrmsNormMask(x,w,id);
 *  Otherwise,
 *      nrm = N_VWrmsNorm(x,w);
 * 
 * mask = FALSE       when the call is made from the nonlinear solver.
 * mask = suppressalg otherwise.
 */

realtype IDAWrmsNorm(IDAMem IDA_mem, N_Vector x, N_Vector w, 
                     booleantype mask)
{
  realtype nrm;

  if (mask) nrm = N_VWrmsNormMask(x, w, id);
  else      nrm = N_VWrmsNorm(x, w);

  return(nrm);
}

/* 
 * -----------------------------------------------------------------
 * Functions for rootfinding
 * -----------------------------------------------------------------
 */

/*
 * IDARcheck1
 *
 * This routine completes the initialization of rootfinding memory
 * information, and checks whether g has a zero both at and very near
 * the initial point of the IVP.
 *
 * This routine returns an int equal to:
 *  IDA_RTFUNC_FAIL < 0 if the g function failed, or
 *  IDA_SUCCESS     = 0 otherwise.
 */

static int IDARcheck1(IDAMem IDA_mem)
{
  int i, retval;
  realtype smallh, hratio, tplus;
  booleantype zroot;

  for (i = 0; i < nrtfn; i++) iroots[i] = 0;
  tlo = tn;
  ttol = (SUNRabs(tn) + SUNRabs(hh))*uround*HUNDRED;

  /* Evaluate g at initial t and check for zero values. */
  retval = gfun (tlo, phi[0], phi[1], glo, user_data);
  nge = 1;
  if (retval != 0) return(IDA_RTFUNC_FAIL);

  zroot = FALSE;
  for (i = 0; i < nrtfn; i++) {
    if (SUNRabs(glo[i]) == ZERO) {
      zroot = TRUE;
      gactive[i] = FALSE;
    }
  }
  if (!zroot) return(IDA_SUCCESS);

  /* Some g_i is zero at t0; look at g at t0+(small increment). */
  hratio = SUNMAX(ttol/SUNRabs(hh), PT1);
  smallh = hratio*hh;
  tplus = tlo + smallh;
  N_VLinearSum(ONE, phi[0], smallh, phi[1], yy);
  retval = gfun (tplus, yy, phi[1], ghi, user_data);  
  nge++;
  if (retval != 0) return(IDA_RTFUNC_FAIL);

  /* We check now only the components of g which were exactly 0.0 at t0
   * to see if we can 'activate' them. */
  for (i = 0; i < nrtfn; i++) {
    if (!gactive[i] && SUNRabs(ghi[i]) != ZERO) {
      gactive[i] = TRUE;
      glo[i] = ghi[i];
    }
  }
  return(IDA_SUCCESS);
}

/*
 * IDARcheck2
 *
 * This routine checks for exact zeros of g at the last root found,
 * if the last return was a root.  It then checks for a close pair of
 * zeros (an error condition), and for a new root at a nearby point.
 * The array glo = g(tlo) at the left endpoint of the search interval
 * is adjusted if necessary to assure that all g_i are nonzero
 * there, before returning to do a root search in the interval.
 *
 * On entry, tlo = tretlast is the last value of tret returned by
 * IDASolve.  This may be the previous tn, the previous tout value,
 * or the last root location.
 *
 * This routine returns an int equal to:
 *     IDA_RTFUNC_FAIL < 0 if the g function failed, or
 *     CLOSERT         = 3 if a close pair of zeros was found, or
 *     RTFOUND         = 1 if a new zero of g was found near tlo, or
 *     IDA_SUCCESS     = 0 otherwise.
 */

static int IDARcheck2(IDAMem IDA_mem)
{
  int i, retval;
  realtype smallh, hratio, tplus;
  booleantype zroot;

  if (irfnd == 0) return(IDA_SUCCESS);

  (void) IDAGetSolution(IDA_mem, tlo, yy, yp);
  retval = gfun (tlo, yy, yp, glo, user_data);  
  nge++;
  if (retval != 0) return(IDA_RTFUNC_FAIL);

  zroot = FALSE;
  for (i = 0; i < nrtfn; i++) iroots[i] = 0;
  for (i = 0; i < nrtfn; i++) {
    if (!gactive[i]) continue;
    if (SUNRabs(glo[i]) == ZERO) {
      zroot = TRUE;
      iroots[i] = 1;
    }
  }
  if (!zroot) return(IDA_SUCCESS);

  /* One or more g_i has a zero at tlo.  Check g at tlo+smallh. */
  ttol = (SUNRabs(tn) + SUNRabs(hh))*uround*HUNDRED;
  smallh = (hh > ZERO) ? ttol : -ttol;
  tplus = tlo + smallh;
  if ( (tplus - tn)*hh >= ZERO) {
    hratio = smallh/hh;
    N_VLinearSum(ONE, yy, hratio, phi[1], yy);
  } else {
    (void) IDAGetSolution(IDA_mem, tplus, yy, yp);
  }
  retval = gfun (tplus, yy, yp, ghi, user_data);  
  nge++;
  if (retval != 0) return(IDA_RTFUNC_FAIL);

  /* Check for close roots (error return), for a new zero at tlo+smallh,
  and for a g_i that changed from zero to nonzero. */
  zroot = FALSE;
  for (i = 0; i < nrtfn; i++) {
    if (!gactive[i]) continue;
    if (SUNRabs(ghi[i]) == ZERO) {
      if (iroots[i] == 1) return(CLOSERT);
      zroot = TRUE;
      iroots[i] = 1;
    } else {
      if (iroots[i] == 1) glo[i] = ghi[i];
    }
  }
  if (zroot) return(RTFOUND);
  return(IDA_SUCCESS);
}

/*
 * IDARcheck3
 *
 * This routine interfaces to IDARootfind to look for a root of g
 * between tlo and either tn or tout, whichever comes first.
 * Only roots beyond tlo in the direction of integration are sought.
 *
 * This routine returns an int equal to:
 *     IDA_RTFUNC_FAIL < 0 if the g function failed, or
 *     RTFOUND         = 1 if a root of g was found, or
 *     IDA_SUCCESS     = 0 otherwise.
 */

static int IDARcheck3(IDAMem IDA_mem)
{
  int i, ier, retval;

  /* Set thi = tn or tout, whichever comes first. */
  if (taskc == IDA_ONE_STEP) thi = tn;
  if (taskc == IDA_NORMAL) {
    thi = ( (toutc - tn)*hh >= ZERO) ? tn : toutc;
  }

  /* Get y and y' at thi. */
  (void) IDAGetSolution(IDA_mem, thi, yy, yp);


  /* Set ghi = g(thi) and call IDARootfind to search (tlo,thi) for roots. */
  retval = gfun (thi, yy, yp, ghi, user_data);  
  nge++;
  if (retval != 0) return(IDA_RTFUNC_FAIL);

  ttol = (SUNRabs(tn) + SUNRabs(hh))*uround*HUNDRED;
  ier = IDARootfind(IDA_mem);
  if (ier == IDA_RTFUNC_FAIL) return(IDA_RTFUNC_FAIL);
  for(i=0; i<nrtfn; i++) {
    if(!gactive[i] && grout[i] != ZERO) gactive[i] = TRUE;
  }
  tlo = trout;
  for (i = 0; i < nrtfn; i++) glo[i] = grout[i];

  /* If no root found, return IDA_SUCCESS. */  
  if (ier == IDA_SUCCESS) return(IDA_SUCCESS);

  /* If a root was found, interpolate to get y(trout) and return.  */
  (void) IDAGetSolution(IDA_mem, trout, yy, yp);
  return(RTFOUND);
}

/*
 * IDARootfind
 *
 * This routine solves for a root of g(t) between tlo and thi, if
 * one exists.  Only roots of odd multiplicity (i.e. with a change
 * of sign in one of the g_i), or exact zeros, are found.
 * Here the sign of tlo - thi is arbitrary, but if multiple roots
 * are found, the one closest to tlo is returned.
 *
 * The method used is the Illinois algorithm, a modified secant method.
 * Reference: Kathie L. Hiebert and Lawrence F. Shampine, Implicitly
 * Defined Output Points for Solutions of ODEs, Sandia National
 * Laboratory Report SAND80-0180, February 1980.
 *
 * This routine uses the following parameters for communication:
 *
 * nrtfn    = number of functions g_i, or number of components of
 *            the vector-valued function g(t).  Input only.
 *
 * gfun     = user-defined function for g(t).  Its form is
 *            (void) gfun(t, y, yp, gt, user_data)
 *
 * rootdir  = in array specifying the direction of zero-crossings.
 *            If rootdir[i] > 0, search for roots of g_i only if
 *            g_i is increasing; if rootdir[i] < 0, search for
 *            roots of g_i only if g_i is decreasing; otherwise
 *            always search for roots of g_i.
 *
 * gactive  = array specifying whether a component of g should
 *            or should not be monitored. gactive[i] is initially
 *            set to TRUE for all i=0,...,nrtfn-1, but it may be
 *            reset to FALSE if at the first step g[i] is 0.0
 *            both at the I.C. and at a small perturbation of them.
 *            gactive[i] is then set back on TRUE only after the 
 *            corresponding g function moves away from 0.0.
 *
 * nge      = cumulative counter for gfun calls.
 *
 * ttol     = a convergence tolerance for trout.  Input only.
 *            When a root at trout is found, it is located only to
 *            within a tolerance of ttol.  Typically, ttol should
 *            be set to a value on the order of
 *               100 * UROUND * max (SUNRabs(tlo), SUNRabs(thi))
 *            where UROUND is the unit roundoff of the machine.
 *
 * tlo, thi = endpoints of the interval in which roots are sought.
 *            On input, these must be distinct, but tlo - thi may
 *            be of either sign.  The direction of integration is
 *            assumed to be from tlo to thi.  On return, tlo and thi
 *            are the endpoints of the final relevant interval.
 *
 * glo, ghi = arrays of length nrtfn containing the vectors g(tlo)
 *            and g(thi) respectively.  Input and output.  On input,
 *            none of the glo[i] should be zero.
 *
 * trout    = root location, if a root was found, or thi if not.
 *            Output only.  If a root was found other than an exact
 *            zero of g, trout is the endpoint thi of the final
 *            interval bracketing the root, with size at most ttol.
 *
 * grout    = array of length nrtfn containing g(trout) on return.
 *
 * iroots   = int array of length nrtfn with root information.
 *            Output only.  If a root was found, iroots indicates
 *            which components g_i have a root at trout.  For
 *            i = 0, ..., nrtfn-1, iroots[i] = 1 if g_i has a root
 *            and g_i is increasing, iroots[i] = -1 if g_i has a
 *            root and g_i is decreasing, and iroots[i] = 0 if g_i
 *            has no roots or g_i varies in the direction opposite
 *            to that indicated by rootdir[i].
 *
 * This routine returns an int equal to:
 *      IDA_RTFUNC_FAIL < 0 if the g function failed, or
 *      RTFOUND         = 1 if a root of g was found, or
 *      IDA_SUCCESS     = 0 otherwise.
 *
 */

static int IDARootfind(IDAMem IDA_mem)
{
  realtype alph, tmid, gfrac, maxfrac, fracint, fracsub;
  int i, retval, imax, side, sideprev;
  booleantype zroot, sgnchg;

  imax = 0;

  /* First check for change in sign in ghi or for a zero in ghi. */
  maxfrac = ZERO;
  zroot = FALSE;
  sgnchg = FALSE;
  for (i = 0;  i < nrtfn; i++) {
    if(!gactive[i]) continue;
    if (SUNRabs(ghi[i]) == ZERO) {
      if(rootdir[i]*glo[i] <= ZERO) {
        zroot = TRUE;
      }
    } else {
      if ( (glo[i]*ghi[i] < ZERO) && (rootdir[i]*glo[i] <= ZERO) ) {
        gfrac = SUNRabs(ghi[i]/(ghi[i] - glo[i]));
        if (gfrac > maxfrac) {
          sgnchg = TRUE;
          maxfrac = gfrac;
          imax = i;
        }
      }
    }
  }

  /* If no sign change was found, reset trout and grout.  Then return
     IDA_SUCCESS if no zero was found, or set iroots and return RTFOUND.  */ 
  if (!sgnchg) {
    trout = thi;
    for (i = 0; i < nrtfn; i++) grout[i] = ghi[i];
    if (!zroot) return(IDA_SUCCESS);
    for (i = 0; i < nrtfn; i++) {
      iroots[i] = 0;
      if(!gactive[i]) continue;
      if ( (SUNRabs(ghi[i]) == ZERO) && (rootdir[i]*glo[i] <= ZERO) )
        iroots[i] = glo[i] > 0 ? -1:1;
    }
    return(RTFOUND);
  }

  /* Initialize alph to avoid compiler warning */
  alph = ONE;

  /* A sign change was found.  Loop to locate nearest root. */

  side = 0;  sideprev = -1;
  loop {                                    /* Looping point */

    /* If interval size is already less than tolerance ttol, break. */
      if (SUNRabs(thi - tlo) <= ttol) break;

    /* Set weight alph.
       On the first two passes, set alph = 1.  Thereafter, reset alph
       according to the side (low vs high) of the subinterval in which
       the sign change was found in the previous two passes.
       If the sides were opposite, set alph = 1.
       If the sides were the same, then double alph (if high side),
       or halve alph (if low side).
       The next guess tmid is the secant method value if alph = 1, but
       is closer to tlo if alph < 1, and closer to thi if alph > 1.    */

    if (sideprev == side) {
      alph = (side == 2) ? alph*TWO : alph*HALF;
    } else {
      alph = ONE;
    }

    /* Set next root approximation tmid and get g(tmid).
       If tmid is too close to tlo or thi, adjust it inward,
       by a fractional distance that is between 0.1 and 0.5.  */
    tmid = thi - (thi - tlo)*ghi[imax]/(ghi[imax] - alph*glo[imax]);
    if (SUNRabs(tmid - tlo) < HALF*ttol) {
      fracint = SUNRabs(thi - tlo)/ttol;
      fracsub = (fracint > FIVE) ? PT1 : HALF/fracint;
      tmid = tlo + fracsub*(thi - tlo);
    }
    if (SUNRabs(thi - tmid) < HALF*ttol) {
      fracint = SUNRabs(thi - tlo)/ttol;
      fracsub = (fracint > FIVE) ? PT1 : HALF/fracint;
      tmid = thi - fracsub*(thi - tlo);
    }

    (void) IDAGetSolution(IDA_mem, tmid, yy, yp);
    retval = gfun (tmid, yy, yp, grout, user_data);  
    nge++;
    if (retval != 0) return(IDA_RTFUNC_FAIL);

    /* Check to see in which subinterval g changes sign, and reset imax.
       Set side = 1 if sign change is on low side, or 2 if on high side.  */  
    maxfrac = ZERO;
    zroot = FALSE;
    sgnchg = FALSE;
    sideprev = side;
    for (i = 0;  i < nrtfn; i++) {
      if(!gactive[i]) continue;
      if (SUNRabs(grout[i]) == ZERO) {
        if(rootdir[i]*glo[i] <= ZERO) zroot = TRUE;
      } else {
        if ( (glo[i]*grout[i] < ZERO) && (rootdir[i]*glo[i] <= ZERO) ) {
          gfrac = SUNRabs(grout[i]/(grout[i] - glo[i]));
          if (gfrac > maxfrac) {
            sgnchg = TRUE;
            maxfrac = gfrac;
            imax = i;
          }
        }
      }
    }
    if (sgnchg) {
      /* Sign change found in (tlo,tmid); replace thi with tmid. */
      thi = tmid;
      for (i = 0; i < nrtfn; i++) ghi[i] = grout[i];
      side = 1;
      /* Stop at root thi if converged; otherwise loop. */
      if (SUNRabs(thi - tlo) <= ttol) break;
      continue;  /* Return to looping point. */
    }

    if (zroot) {
      /* No sign change in (tlo,tmid), but g = 0 at tmid; return root tmid. */
      thi = tmid;
      for (i = 0; i < nrtfn; i++) ghi[i] = grout[i];
      break;
    }

    /* No sign change in (tlo,tmid), and no zero at tmid.
       Sign change must be in (tmid,thi).  Replace tlo with tmid. */
    tlo = tmid;
    for (i = 0; i < nrtfn; i++) glo[i] = grout[i];
    side = 2;
    /* Stop at root thi if converged; otherwise loop back. */
    if (SUNRabs(thi - tlo) <= ttol) break;

  } /* End of root-search loop */

  /* Reset trout and grout, set iroots, and return RTFOUND. */
  trout = thi;
  for (i = 0; i < nrtfn; i++) {
    grout[i] = ghi[i];
    iroots[i] = 0;
    if(!gactive[i]) continue;
    if ( (SUNRabs(ghi[i]) == ZERO) && (rootdir[i]*glo[i] <= ZERO) )
      iroots[i] = glo[i] > 0 ? -1:1;
    if ( (glo[i]*ghi[i] < ZERO) && (rootdir[i]*glo[i] <= ZERO) ) 
      iroots[i] = glo[i] > 0 ? -1:1;
  }
  return(RTFOUND);
}

/* 
 * =================================================================
 * IDA error message handling functions   
 * =================================================================
 */

/*
 * IDAProcessError is a high level error handling function.
 * - If ida_mem==NULL it prints the error message to stderr.
 * - Otherwise, it sets up and calls the error handling function 
 *   pointed to by ida_ehfun.
 */

#define ehfun   (IDA_mem->ida_ehfun)
#define eh_data (IDA_mem->ida_eh_data)

void IDAProcessError(IDAMem IDA_mem, 
                    int error_code, const char *module, const char *fname, 
                    const char *msgfmt, ...)
{
  va_list ap;
  char msg[256];

  /* Initialize the argument pointer variable 
     (msgfmt is the last required argument to IDAProcessError) */

  va_start(ap, msgfmt);

  /* Compose the message */

  vsprintf(msg, msgfmt, ap);

  if (IDA_mem == NULL) {    /* We write to stderr */
#ifndef NO_FPRINTF_OUTPUT
    fprintf(stderr, "\n[%s ERROR]  %s\n  ", module, fname);
    fprintf(stderr, "%s\n\n", msg);
#endif

  } else {                 /* We can call ehfun */
    ehfun(error_code, module, fname, msg, eh_data);
  }

  /* Finalize argument processing */
  va_end(ap);

  return;
}

/* IDAErrHandler is the default error handling function.
   It sends the error message to the stream pointed to by ida_errfp */

#define errfp (IDA_mem->ida_errfp)

void IDAErrHandler(int error_code, const char *module,
                   const char *function, char *msg, void *data)
{
  IDAMem IDA_mem;
  char err_type[10];

  /* data points to IDA_mem here */

  IDA_mem = (IDAMem) data;

  if (error_code == IDA_WARNING)
    sprintf(err_type,"WARNING");
  else
    sprintf(err_type,"ERROR");

#ifndef NO_FPRINTF_OUTPUT
  if (errfp!=NULL) {
    fprintf(errfp,"\n[%s %s]  %s\n",module,err_type,function);
    fprintf(errfp,"  %s\n\n",msg);
  }
#endif

  return;
}
