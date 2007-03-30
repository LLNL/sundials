/*
 * -----------------------------------------------------------------
 * $Revision: 1.10 $
 * $Date: 2007-03-30 15:05:02 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the CVODEA adjoint integrator.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvodes_impl.h"

#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

/* 
 * =================================================================
 * MACRO DEFINITIONS
 * =================================================================
 */

#define loop for(;;)

/* 
 * =================================================================
 * CVODEA PRIVATE CONSTANTS
 * =================================================================
 */

#define ZERO        RCONST(0.0)        /* real 0.0   */
#define ONE         RCONST(1.0)        /* real 1.0   */
#define TWO         RCONST(2.0)        /* real 2.0   */
#define HUNDRED     RCONST(100.0)      /* real 100.0 */
#define FUZZ_FACTOR RCONST(1000000.0)  /* fuzz factor for getY */

/* 
 * =================================================================
 * PRIVATE FUNCTION PROTOTYPES
 * =================================================================
 */

static CkpntMem CVAckpntInit(CVodeMem cv_mem);
static CkpntMem CVAckpntNew(CVodeMem cv_mem);
static void CVAckpntDelete(CkpntMem *ck_memPtr);

static void CVAbckpbDelete(CVodeBMem *cvB_memPtr);

static int  CVAdataStore(CVodeMem cv_mem, CkpntMem ck_mem);
static int  CVAckpntGet(CVodeMem cv_mem, CkpntMem ck_mem); 

static int CVAfindIndex(CVodeMem cv_mem, realtype t, 
                        long int *indx, booleantype *newpoint);

static booleantype CVAhermiteMalloc(CVodeMem cv_mem, long int steps);
static void CVAhermiteFree(DtpntMem *dt_mem, long int steps);
static int CVAhermiteGetY(CVodeMem cv_mem, realtype t, N_Vector y);
static int CVAhermiteStorePnt(CVodeMem cv_mem, DtpntMem d);

static booleantype CVApolynomialMalloc(CVodeMem cv_mem, long int steps);
static void CVApolynomialFree(DtpntMem *dt_mem, long int steps);
static int CVApolynomialGetY(CVodeMem cv_mem, realtype t, N_Vector y);
static int CVApolynomialStorePnt(CVodeMem cv_mem, DtpntMem d);

/* Wrappers */

static int CVArhs(realtype t, N_Vector yB, 
                  N_Vector yBdot, void *cvode_mem);

static int CVArhsQ(realtype t, N_Vector yB, 
                   N_Vector qBdot, void *cvode_mem);

/* 
 * =================================================================
 * EXPORTED FUNCTIONS IMPLEMENTATION
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * Readibility Constants
 * -----------------------------------------------------------------
 */

#define tinitial    (ca_mem->ca_tinitial)
#define tfinal      (ca_mem->ca_tfinal)
#define nckpnts     (ca_mem->ca_nckpnts)
#define nsteps      (ca_mem->ca_nsteps)
#define nbckpbs     (ca_mem->ca_nbckpbs)
#define ckpntData   (ca_mem->ca_ckpntData)
#define newData     (ca_mem->ca_newData)
#define np          (ca_mem->ca_np)
#define ytmp        (ca_mem->ca_ytmp)
#define Y0          (ca_mem->ca_Y0)
#define Y1          (ca_mem->ca_Y1)
#define Y           (ca_mem->ca_Y)
#define T           (ca_mem->ca_T)
#define interpType  (ca_mem->ca_interpType)
#define getY        (ca_mem->ca_getY)
#define storePnt    (ca_mem->ca_storePnt)

#define uround     (cv_mem->cv_uround)
#define zn         (cv_mem->cv_zn)
#define nst        (cv_mem->cv_nst)
#define q          (cv_mem->cv_q)
#define qu         (cv_mem->cv_qu)
#define qprime     (cv_mem->cv_qprime)
#define qwait      (cv_mem->cv_qwait)
#define L          (cv_mem->cv_L)
#define gammap     (cv_mem->cv_gammap)
#define h          (cv_mem->cv_h)
#define hprime     (cv_mem->cv_hprime)
#define hscale     (cv_mem->cv_hscale)
#define eta        (cv_mem->cv_eta)
#define etamax     (cv_mem->cv_etamax)
#define tn         (cv_mem->cv_tn)
#define tretlast   (cv_mem->cv_tretlast)
#define tau        (cv_mem->cv_tau)
#define tq         (cv_mem->cv_tq)
#define l          (cv_mem->cv_l)
#define saved_tq5  (cv_mem->cv_saved_tq5)
#define forceSetup (cv_mem->cv_forceSetup)
#define f          (cv_mem->cv_f)
#define lmm        (cv_mem->cv_lmm)
#define iter       (cv_mem->cv_iter)
#define itol       (cv_mem->cv_itol)
#define reltol     (cv_mem->cv_reltol)
#define Sabstol    (cv_mem->cv_Sabstol)
#define Vabstol    (cv_mem->cv_Vabstol)
#define efun       (cv_mem->cv_efun)
#define f_data     (cv_mem->cv_f_data)
#define errfp      (cv_mem->cv_errfp)
#define h0u        (cv_mem->cv_h0u)
#define quadr      (cv_mem->cv_quadr)
#define errconQ    (cv_mem->cv_errconQ)
#define znQ        (cv_mem->cv_znQ)
#define itolQ      (cv_mem->cv_itolQ)
#define reltolQ    (cv_mem->cv_reltolQ)
#define SabstolQ   (cv_mem->cv_SabstolQ)
#define VabstolQ   (cv_mem->cv_VabstolQ)
#define fQ         (cv_mem->cv_fQ)
#define tempv      (cv_mem->cv_tempv)
#define tempvQ     (cv_mem->cv_tempvQ)

#define t0_        (ck_mem->ck_t0)
#define t1_        (ck_mem->ck_t1)
#define zn_        (ck_mem->ck_zn)
#define znQ_       (ck_mem->ck_znQ)
#define quadr_     (ck_mem->ck_quadr)
#define zqm_       (ck_mem->ck_zqm)
#define nst_       (ck_mem->ck_nst)
#define tretlast_  (ck_mem->ck_tretlast)
#define q_         (ck_mem->ck_q)
#define qprime_    (ck_mem->ck_qprime)
#define qwait_     (ck_mem->ck_qwait)
#define L_         (ck_mem->ck_L)
#define gammap_    (ck_mem->ck_gammap)
#define h_         (ck_mem->ck_h)
#define hprime_    (ck_mem->ck_hprime)
#define hscale_    (ck_mem->ck_hscale)
#define eta_       (ck_mem->ck_eta)
#define etamax_    (ck_mem->ck_etamax)
#define tau_       (ck_mem->ck_tau)
#define tq_        (ck_mem->ck_tq)
#define l_         (ck_mem->ck_l)
#define saved_tq5_ (ck_mem->ck_saved_tq5)
#define next_      (ck_mem->ck_next)

/*
 * CVodeAdjMalloc
 *
 * This routine initializes ASA and allocates space for the adjoint 
 * memory structure.
 */

int CVodeAdjMalloc(void *cvode_mem, long int steps, int interp)
{
  CVadjMem ca_mem;
  CVodeMem cv_mem;
  booleantype allocOK;
  long int i, ii;

  /* ---------------
   * Check arguments
   * --------------- */

  if (cvode_mem == NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeAdjMalloc", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem)cvode_mem;

  if (steps <= 0) {
    CVProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeAdjMalloc", MSGCV_BAD_STEPS);
    return(CV_ILL_INPUT);
  }

  if ( (interp != CV_HERMITE) && (interp != CV_POLYNOMIAL) ) {
    CVProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeAdjMalloc", MSGCV_BAD_INTERP);
    return(CV_ILL_INPUT);
  } 

  /* ----------------------------
   * Allocate CVODEA memory block
   * ---------------------------- */

  ca_mem = NULL;
  ca_mem = (CVadjMem) malloc(sizeof(struct CVadjMemRec));
  if (ca_mem == NULL) {
    CVProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeAdjMalloc", MSGCV_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  /* Attach ca_mem to CVodeMem structure */

  cv_mem->cv_adj_mem = ca_mem;

  /* Allocate memory for the workspace vector ytmp */

  ytmp = NULL;
  ytmp = N_VClone(tempv);
  if (ytmp == NULL) {
    free(ca_mem); ca_mem = NULL;
    CVProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeAdjMalloc", MSGCV_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  /* ------------------------------
   * Initialization of check points
   * ------------------------------ */

  /* Initialize Check Points linked list */

  ca_mem->ck_mem = NULL;
  ca_mem->ck_mem = CVAckpntInit(cv_mem);
  if (ca_mem->ck_mem == NULL) {
    N_VDestroy(ytmp);
    free(ca_mem); ca_mem = NULL;
    CVProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeAdjMalloc", MSGCV_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  /* Initialize nckpnts to ZERO */

  nckpnts = 0;
  ca_mem->ca_ckpntData = NULL;

  /* ------------------------------------
   * Initialization of interpolation data
   * ------------------------------------ */

  /* Interpolation type */

  interpType = interp;

  /* Number of steps between check points */

  nsteps = steps;

  /* Allocate space for the array of Data Point structures */

  ca_mem->dt_mem = NULL;
  ca_mem->dt_mem = (DtpntMem *) malloc((steps+1)*sizeof(struct DtpntMemRec *));
  if (ca_mem->dt_mem == NULL) {
    while (ca_mem->ck_mem != NULL) CVAckpntDelete(&(ca_mem->ck_mem));
    N_VDestroy(ytmp);
    free(ca_mem); ca_mem = NULL;
    CVProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeAdjMalloc", MSGCV_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  for (i=0; i<=steps; i++) { 
    ca_mem->dt_mem[i] = NULL;
    ca_mem->dt_mem[i] = (DtpntMem) malloc(sizeof(struct DtpntMemRec));
    if (ca_mem->dt_mem[i] == NULL) {
      for(ii=0; ii<i; ii++) {free(ca_mem->dt_mem[ii]); ca_mem->dt_mem[ii] = NULL;}
      free(ca_mem->dt_mem); ca_mem->dt_mem = NULL;
      while (ca_mem->ck_mem != NULL) CVAckpntDelete(&(ca_mem->ck_mem));
      N_VDestroy(ytmp);
      free(ca_mem); ca_mem = NULL;
      CVProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeAdjMalloc", MSGCV_MEM_FAIL);
      return(CV_MEM_FAIL);
    }
  }

  /* Workspace for interpolation functions */

  switch (interpType) {

  case CV_HERMITE:
    /* Allocate Data Points memory */
    allocOK = CVAhermiteMalloc(cv_mem, steps);
    if (!allocOK) {
      for(i=0; i<=steps; i++) {free(ca_mem->dt_mem[i]); ca_mem->dt_mem[i] = NULL;}
      free(ca_mem->dt_mem); ca_mem->dt_mem = NULL;
      while (ca_mem->ck_mem != NULL) CVAckpntDelete(&(ca_mem->ck_mem));
      N_VDestroy(ytmp);
      free(ca_mem); ca_mem = NULL;
      CVProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeAdjMalloc", MSGCV_MEM_FAIL);
      return(CV_MEM_FAIL);
    }
    /* Attach interpolation functions getY and storePnt */
    getY = CVAhermiteGetY;
    storePnt = CVAhermiteStorePnt;
    /* Rename zn[0] and zn[1] for use in interpolation */
    Y0 = zn[0];
    Y1 = zn[1];
    break;
  case CV_POLYNOMIAL:
    /* Allocate Data Points memory */
    allocOK = CVApolynomialMalloc(cv_mem, steps);
    if (!allocOK) {
      for(i=0; i<=steps; i++) {free(ca_mem->dt_mem[i]); ca_mem->dt_mem[i] = NULL;}
      free(ca_mem->dt_mem); ca_mem->dt_mem = NULL;
      while (ca_mem->ck_mem != NULL) CVAckpntDelete(&(ca_mem->ck_mem));
      N_VDestroy(ytmp);
      free(ca_mem); ca_mem = NULL;
      CVProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeAdjMalloc", MSGCV_MEM_FAIL);
      return(CV_MEM_FAIL);
    }
    /* Attach interpolation functions getY and storePnt */
    getY = CVApolynomialGetY;
    storePnt = CVApolynomialStorePnt;
    /* Rename zn for use in interpolation */
    for (i=0;i<L_MAX;i++) Y[i] = zn[i];
    break;
  }

  /* ------------------------------------
   * Initialize list of backward problems
   * ------------------------------------ */

  ca_mem->cvB_mem = NULL;
  ca_mem->ca_bckpbCrt = NULL;
  nbckpbs = 0;

  /* --------------------------------
   * CVodeF and CVodeB not called yet
   * -------------------------------- */

  ca_mem->ca_firstCVodeFcall = TRUE;
  ca_mem->ca_tstopCVodeFcall = FALSE;

  ca_mem->ca_firstCVodeBcall = TRUE;

  /* ---------------------------------------------
   * ASA initialized and allocated
   * --------------------------------------------- */

  cv_mem->cv_adj = TRUE;
  cv_mem->cv_adjMallocDone = TRUE;

  return(CV_SUCCESS);
} 

/* CVodeAdjReInit
 *
 * This routine reinitializes the CVODEA memory structure assuming that the
 * the number of steps between check points and the type of interpolation
 * remain unchanged.
 * The list of check points (and associated memory) is deleted.
 * The list of backward problems is kept (however, new backward problems can 
 * be added to this list by calling CVodeCreateB).
 * The CVODES memory for the forward and backward problems can be reinitialized
 * separately by calling CVodeReInit and CVodeReInitB, respectively.
 * NOTE: if a completely new list of backward problems is also needed, then
 *       simply free the adjoint memory (by calling CVodeAdjFree) and reinitialize
 *       ASA with CVodeAdjMalloc.
 */

int CVodeAdjReInit(void *cvode_mem)
{
  CVadjMem ca_mem;
  CVodeMem cv_mem;

  /* Check cvode_mem */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeAdjReInit", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    CVProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeAdjReInit", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 

  ca_mem = cv_mem->cv_adj_mem;

  /* Free current list of Check Points */

  while (ca_mem->ck_mem != NULL) CVAckpntDelete(&(ca_mem->ck_mem));

  /* Re-initialize Check Points linked list */

  ca_mem->ck_mem = NULL;
  ca_mem->ck_mem = CVAckpntInit(cv_mem);
  if (ca_mem->ck_mem == NULL) {
    CVProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeAdjReInit", MSGCV_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  nckpnts = 0;
  ca_mem->ca_ckpntData = NULL;

  /* CVodeF and CVodeB not called yet */
 
  ca_mem->ca_firstCVodeFcall = TRUE;
  ca_mem->ca_tstopCVodeFcall = FALSE;
  ca_mem->ca_firstCVodeBcall = TRUE;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetAdjInterpType
 *
 * Changes the interpolation type.
 */

int CVodeSetAdjInterpType(void *cvode_mem, int interp)
{
  CVadjMem ca_mem;
  CVodeMem cv_mem;
  booleantype allocOK;
  long int i;
  
  /* Check cvode_mem */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeSetAdjInterpType", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    CVProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeSetAdjInterpType", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 

  ca_mem = cv_mem->cv_adj_mem;

  if ( (interp != CV_HERMITE) && (interp != CV_POLYNOMIAL) ) {
    CVProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeSetAdjInterpType", MSGCV_BAD_INTERP);
    return(CV_ILL_INPUT);
  } 

  if (interp == interpType) return(CV_SUCCESS);

  interpType = interp;

  switch (interpType) {

  case CV_HERMITE:
    /* Delete Data Points memory */
    CVApolynomialFree(ca_mem->dt_mem, nsteps);
    /* Allocate Data Points memory */
    allocOK = CVAhermiteMalloc(cv_mem, nsteps);
    if (!allocOK) {
      CVProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeSetAdjInterpType", MSGCV_MEM_FAIL);
      return(CV_MEM_FAIL);
    }
    /* Attach interpolation functions getY and storePnt */
    getY = CVAhermiteGetY;
    storePnt = CVAhermiteStorePnt;
    /* Rename zn[0] and zn[1] for use in interpolation */
    Y0 = zn[0];
    Y1 = zn[1];
    break;
  case CV_POLYNOMIAL:
    /* Delete Data Points memory */
    CVAhermiteFree(ca_mem->dt_mem, nsteps);
    /* Allocate Data Points memory */
    allocOK = CVApolynomialMalloc(cv_mem, nsteps);
    if (!allocOK) {
      CVProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeSetAdjInterpType", MSGCV_MEM_FAIL);
      return(CV_MEM_FAIL);
    }
    /* Attach interpolation functions getY and storePnt */
    getY = CVApolynomialGetY;
    storePnt = CVApolynomialStorePnt;
    /* Rename zn for use in interpolation */
    for (i=0;i<L_MAX;i++) Y[i] = zn[i];
    break;
  }

  return(CV_SUCCESS);

}

/*
 * CVodeAdjFree
 *
 * This routine frees the memory allocated by CVodeAdjMalloc.
 */

void CVodeAdjFree(void *cvode_mem)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  void *cvode_bmem;
  long int i;
  
  if (cvode_mem == NULL) return;
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_adjMallocDone) {

    ca_mem = cv_mem->cv_adj_mem;

    /* Delete check points one by one */
    while (ca_mem->ck_mem != NULL) CVAckpntDelete(&(ca_mem->ck_mem));

    /* Free vectors at each data point */
    switch (interpType) {
    case CV_HERMITE:
      CVAhermiteFree(ca_mem->dt_mem, nsteps);
      break;
    case CV_POLYNOMIAL:
      CVApolynomialFree(ca_mem->dt_mem, nsteps);
      break;
    }
    for(i=0; i<=nsteps; i++) {free(ca_mem->dt_mem[i]); ca_mem->dt_mem[i] = NULL;}
    free(ca_mem->dt_mem); ca_mem->dt_mem = NULL;

    /* Delete backward problems one by one */
    while (ca_mem->cvB_mem != NULL) CVAbckpbDelete(&(ca_mem->cvB_mem));

    /* Free workspace vector ytmp in ca_mem */
    N_VDestroy(ytmp);

    /* Free CVODEA memory */
    free(ca_mem);
    cv_mem->cv_adj_mem = NULL;

  }

}

/*
 * CVodeF
 *
 * This routine integrates to tout and returns solution into yout.
 * In the same time, it stores check point data every 'steps' steps. 
 * 
 * CVodeF can be called repeatedly by the user.
 *
 * ncheckPtr points to the number of check points stored so far.
 */

int CVodeF(void *cvode_mem, realtype tout, N_Vector yout, 
           realtype *tret, int itask, int *ncheckPtr)
{
  CVadjMem ca_mem;
  CVodeMem cv_mem;
  CkpntMem tmp;
  DtpntMem *dt_mem;
  int cv_itask, flag;
  booleantype iret;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeF", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    CVProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeF", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 

  ca_mem = cv_mem->cv_adj_mem;

  /* Check for yout != NULL */
  if (yout == NULL) {
    CVProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeF", MSGCV_YOUT_NULL);
    return(CV_ILL_INPUT);
  }
  
  /* Check for tret != NULL */
  if (tret == NULL) {
    CVProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeF", MSGCV_TRET_NULL);
    return(CV_ILL_INPUT);
  }

  /* Check for valid itask */
  if ((itask != CV_NORMAL)       && 
      (itask != CV_ONE_STEP)     &&
      (itask != CV_NORMAL_TSTOP) &&
      (itask != CV_ONE_STEP_TSTOP) ) {
    CVProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeF", MSGCV_BAD_ITASK);
    return(CV_ILL_INPUT);
  }

  /* All error checking done */

  dt_mem = ca_mem->dt_mem;

  iret = TRUE;
  cv_itask = CV_ONE_STEP;

  /* Interpret itask */
  switch (itask) {
  case CV_NORMAL:
    iret = FALSE;
    cv_itask = CV_ONE_STEP;
    break;
  case CV_ONE_STEP:
    iret = TRUE;
    cv_itask = CV_ONE_STEP;
    break;
  case CV_NORMAL_TSTOP:
    iret = FALSE;
    cv_itask = CV_ONE_STEP_TSTOP;
    ca_mem->ca_tstopCVodeFcall = TRUE;
    ca_mem->ca_tstopCVodeF = cv_mem->cv_tstop;
    break;
  case CV_ONE_STEP_TSTOP:
    iret = TRUE;
    cv_itask = CV_ONE_STEP_TSTOP;
    ca_mem->ca_tstopCVodeFcall = TRUE;
    ca_mem->ca_tstopCVodeF = cv_mem->cv_tstop;
    break;
  }

  /* On the first step, set tinitial and load dt_mem[0].
     On subsequent steps, test if taking a new step is necessary. */
  if ( ca_mem->ca_firstCVodeFcall ) {

    tinitial = tn;

    dt_mem[0]->t = ca_mem->ck_mem->ck_t0;
    storePnt(cv_mem, dt_mem[0]);

    ca_mem->ca_firstCVodeFcall = FALSE;

  } else if ( (tn - tout)*h >= ZERO ) {

    /* If tout was passed, return interpolated solution. 
       No changes to ck_mem or dt_mem are needed. */
    *tret = tout;
    flag = CVodeGetDky(cv_mem, tout, 0, yout);
    *ncheckPtr = nckpnts;
    newData = TRUE;
    ckpntData = ca_mem->ck_mem;
    np = nst % nsteps + 1;
    return(flag);

  }

  /* Integrate to tout (in CV_ONE_STEP mode) while loading check points */
  loop {

    /* Perform one step of the integration */

    flag = CVode(cv_mem, tout, yout, tret, cv_itask);
    if (flag < 0) break;

    /* Test if a new check point is needed */

    if ( nst % nsteps == 0 ) {

      ca_mem->ck_mem->ck_t1 = *tret;

      /* Create a new check point, load it, and append it to the list */
      tmp = CVAckpntNew(cv_mem);
      if (tmp == NULL) {
        CVProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeF", MSGCV_MEM_FAIL);
        flag = CV_MEM_FAIL;
        break;
      }
      tmp->ck_next = ca_mem->ck_mem;
      ca_mem->ck_mem = tmp;
      nckpnts++;
      forceSetup = TRUE;
      
      /* Reset i=0 and load dt_mem[0] */
      dt_mem[0]->t = ca_mem->ck_mem->ck_t0;
      storePnt(cv_mem, dt_mem[0]);

    } else {
      
      /* Load next point in dt_mem */
      dt_mem[nst%nsteps]->t = *tret;
      storePnt(cv_mem, dt_mem[nst%nsteps]);

    }

    /* Set t1 field of the current ckeck point structure
       for the case in which there will be no future
       check points */
    ca_mem->ck_mem->ck_t1 = *tret;

    /* tfinal is now set to *tret */
    tfinal = *tret;

    /* Return if in CV_ONE_STEP mode */
    if (iret) break;

    /* Return if tout reached */
    if ( (*tret - tout)*h >= ZERO ) {
      *tret = tout;
      CVodeGetDky(cv_mem, tout, 0, yout);
      /* Reset tretlast in cv_mem so that CVodeGetQuad 
       * evaluates quadratures at the proper time */
      cv_mem->cv_tretlast = tout;
      break;
    }

  } /* end of loop() */

  /* Get ncheck from ca_mem */ 
  *ncheckPtr = nckpnts;

  /* Data is available for the last interval */
  newData = TRUE;
  ckpntData = ca_mem->ck_mem;
  np = nst % nsteps + 1;

  return(flag);
}



/* 
 * =================================================================
 * FUNCTIONS FOR BACKWARD PROBLEMS
 * =================================================================
 */


int CVodeCreateB(void *cvode_mem, int lmmB, int iterB, int *which)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem new_cvB_mem;
  void *cvodeB_mem;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeCreateB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    CVProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeCreateB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  }
  ca_mem = cv_mem->cv_adj_mem;

  /* Allocate space for new CVodeBMem object */

  new_cvB_mem = NULL;
  new_cvB_mem = (CVodeBMem) malloc(sizeof(struct CVodeBMemRec));
  if (new_cvB_mem == NULL) {
    CVProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeCreateB", MSGCV_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  /* Create and set a new CVODES object for the backward problem */

  cvodeB_mem = CVodeCreate(lmmB, iterB);
  if (cvodeB_mem == NULL) {
    CVProcessError(cv_mem, CV_MEM_FAIL, "CVODEA", "CVodeCreateB", MSGCV_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  CVodeSetFdata(cvodeB_mem, cvode_mem);

  CVodeSetMaxHnilWarns(cvodeB_mem, -1);

  /* Set/initialize fields in the new CVodeBMem object, new_cvB_mem */

  new_cvB_mem->cv_index   = nbckpbs;

  new_cvB_mem->cv_mem     = (CVodeMem) cvodeB_mem;

  new_cvB_mem->cv_f       = NULL;

  new_cvB_mem->cv_fQ      = NULL;
  new_cvB_mem->cv_f_data  = NULL;
  new_cvB_mem->cv_fQ_data = NULL;
  new_cvB_mem->cv_lmem    = NULL;
  new_cvB_mem->cv_lfree   = NULL;
  new_cvB_mem->cv_pmem    = NULL;
  new_cvB_mem->cv_pfree   = NULL;

  new_cvB_mem->cv_y       = NULL;

  /* Attach the new object to the linked list cvB_mem */

  new_cvB_mem->cv_next = ca_mem->cvB_mem;
  ca_mem->cvB_mem = new_cvB_mem;
  
  /* Return the index of the newly created CVodeBMem object.
   * This must be passed to CVodeMallocB and to other ***B 
   * functions to set optional inputs for this backward problem */

  *which = nbckpbs;

  nbckpbs++;

  return(CV_SUCCESS);
}


int CVodeMallocB(void *cvode_mem, int which, 
                 CVRhsFnB fB,
                 realtype tB0, N_Vector yB0,
                 int itolB, realtype reltolB, void *abstolB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeMallocB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    CVProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeMallocB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check the value of which which */
  if ( which >= nbckpbs ) {
    CVProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeMallocB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);
  
  /* Allocate and set the CVODES object */

  flag = CVodeMalloc(cvodeB_mem, CVArhs, tB0, yB0,
                     itolB, reltolB, abstolB);
  if (flag != CV_SUCCESS) return(flag);

  /* Copy fB function in cvB_mem */

  cvB_mem->cv_f = fB;

  /* Allocate space and initialize for the y Nvector in cvB_mem */
  cvB_mem->cv_t0 = tB0;
  cvB_mem->cv_y = N_VClone(yB0);
  N_VScale(ONE, yB0, cvB_mem->cv_y);

  return(CV_SUCCESS);
}


int CVodeReInitB(void *cvode_mem, int which,
                 CVRhsFnB fB, 
                 realtype tB0, N_Vector yB0, 
                 int itolB, realtype reltolB, void *abstolB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeReInitB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    CVProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeReInitB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  }
  ca_mem = cv_mem->cv_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    CVProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeReInitB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  /* Reinitialize CVODES object */

  flag = CVodeReInit(cvodeB_mem, CVArhs, tB0, yB0,
                     itolB, reltolB, abstolB);
  if (flag != CV_SUCCESS)  return(flag);

  /* Copy fB function in cvB_mem */

  cvB_mem->cv_f = fB;

  return(CV_SUCCESS);
}

/*
 * CVodeQuadMallocB, and CVodeQuadReInitB
 *
 * Wrappers for the backward phase around the corresponding 
 * CVODES functions
 */

int CVodeQuadMallocB(void *cvode_mem, int which,
                     CVQuadRhsFnB fQB, N_Vector yQB0)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeQuadMallocB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    CVProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeQuadMallocB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= nbckpbs ) {
    CVProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeQuadMallocB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  flag = CVodeQuadMalloc(cvodeB_mem, CVArhsQ, yQB0);
  if (flag != CV_SUCCESS) return(flag);

  flag = CVodeSetQuadFdata(cvodeB_mem, cvode_mem);

  cvB_mem->cv_fQ = fQB;

  return(CV_SUCCESS);
}

int CVodeQuadReInitB(void *cvode_mem, int which,
                     CVQuadRhsFnB fQB, N_Vector yQB0)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeQuadReInitB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    CVProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeQuadReInitB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    CVProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeQuadReInitB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  flag = CVodeQuadReInit(cvodeB_mem, CVArhsQ, yQB0);
  if (flag != CV_SUCCESS) return(flag);

  cvB_mem->cv_fQ = fQB;

  return(CV_SUCCESS);
}


/*
 * CVodeB
 *
 * This routine performs the backward integration towards tBout
 * of all backward problems that were defined.
 * When necessary, it performs a forward integration between two 
 * consecutive check points to update interpolation data.
 * itask can be CV_NORMAL or CV_ONE_STEP only.
 *
 * On a successful return, CVodeB returns either CV_SUCCESS
 * (in ONE_STEP mode or if tBout was reached in NORMAL mode)
 * unless the return time happens to be a checkpoint, in which
 * case it returns CV_TSTOP_RETURN)
 *
 * NOTE that CVodeB DOES NOT return the solution for the backward
 * problem(s). Use CVodeGetB to extract the solution at tBret
 * for any given backward problem.
 */

int CVodeB(void *cvode_mem, realtype tBout, int itaskB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem, tmp_cvB_mem;
  CkpntMem ck_mem;
  int sign, flag, cv_itask;
  realtype tBret, tBn, hB, troundoff, tmp_tBn;
  booleantype gotCheckpoint, isActive, reachedTBout;
  
  /* Check if cvode_mem exists */

  if (cvode_mem == NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */

  if (cv_mem->cv_adjMallocDone == FALSE) {
    CVProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  }
  ca_mem = cv_mem->cv_adj_mem;

  /* Check if any backward problem has been defined */

  if ( nbckpbs == 0 ) {
    CVProcessError(cv_mem, CV_NO_BCK, "CVODEA", "CVodeB", MSGCV_NO_BCK);
    return(CV_NO_BCK);
  }
  cvB_mem = ca_mem->cvB_mem;

  /* Check whether CVodeF has been called */

  if ( ca_mem->ca_firstCVodeFcall ) {
    CVProcessError(cv_mem, CV_NO_FWD, "CVODEA", "CVodeB", MSGCV_NO_FWD);
    return(CV_NO_FWD);
  }
  sign = (tfinal - tinitial > ZERO) ? 1 : -1;

  /* If this is the first call, check the tB0 time for each backward problem */

  if ( ca_mem->ca_firstCVodeBcall ) {

    tmp_cvB_mem = cvB_mem;
    while(tmp_cvB_mem != NULL) {
      tBn = tmp_cvB_mem->cv_mem->cv_tn;
      if ( (sign*(tBn-tinitial) < ZERO) || (sign*(tfinal-tBn) < ZERO) ) {
        CVProcessError(cv_mem, CV_BAD_TB0, "CVODEA", "CVodeB", MSGCV_BAD_TB0, tmp_cvB_mem->cv_index);
        return(CV_BAD_TB0);
      }
      tmp_cvB_mem = tmp_cvB_mem->cv_next;
    }

    ca_mem->ca_firstCVodeBcall = FALSE;
  }

  /* Check if itaskB is legal */

  if (itaskB == CV_NORMAL)
    cv_itask = CV_NORMAL_TSTOP;
  else if (itaskB == CV_ONE_STEP)
    cv_itask = CV_ONE_STEP_TSTOP;
  else {
    CVProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeB", MSGCV_BAD_ITASKB);
    return(CV_ILL_INPUT);
  }

  /* Check if tBout is legal */

  if ( (sign*(tBout-tinitial) < ZERO) || (sign*(tfinal-tBout) < ZERO) ) {
    if ( ABS(tBout-tinitial) < HUNDRED*uround ) {
      tBout = tinitial;
    } else {
      CVProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeB", MSGCV_BAD_TBOUT);
      return(CV_ILL_INPUT);
    }
  }

  /* Loop through the check points and stop as soon as a backward
   * problem has its tn value larger than the current check point's
   * t0_ value (taking into account the direction of integration) */

  ck_mem = ca_mem->ck_mem;

  gotCheckpoint = FALSE;

  loop {
    
    tmp_cvB_mem = cvB_mem;
    while(tmp_cvB_mem != NULL) {
      tBn = tmp_cvB_mem->cv_mem->cv_tn;
      hB  = tmp_cvB_mem->cv_mem->cv_hu;
      troundoff = HUNDRED*uround*(ABS(tBn) + ABS(hB));
      if ( sign * (tBn-t0_) > troundoff ) {
        gotCheckpoint = TRUE;
        break;
      }
      tmp_cvB_mem = tmp_cvB_mem->cv_next;
    }

    if (gotCheckpoint) break;

    if (next_ == NULL) break;

    ck_mem = next_;
  }

  /* Loop while propagating backward problems */

  loop {

    /* Store interpolation data if not available.
       This is the 2nd forward integration pass */

    if (ck_mem != ckpntData) {
      flag = CVAdataStore(cv_mem, ck_mem);
      if (flag != CV_SUCCESS) break;
    }

    /* Loop through all backward problems and, if needed,
     * propagate their solution towards tBout */

    tmp_cvB_mem = cvB_mem;
    while (tmp_cvB_mem != NULL) {

      /* Decide if current backward problem is "active" */

      isActive = TRUE;

      tBn = tmp_cvB_mem->cv_mem->cv_tn;
      hB  = tmp_cvB_mem->cv_mem->cv_hu;
      troundoff = HUNDRED*uround*(ABS(tBn) + ABS(hB));

      if ( sign * (tBn - t0_)   < troundoff ) isActive = FALSE;
      if ( sign * (tBn - tBout) < troundoff ) isActive = FALSE;

      if ( isActive ) {

        /* Store the address of current backward problem memory 
         * in ca_mem to be used in the wrapper functions */
        ca_mem->ca_bckpbCrt = tmp_cvB_mem;

        /* Integrate current backward problem */
        CVodeSetStopTime(tmp_cvB_mem->cv_mem, t0_);
        flag = CVode(tmp_cvB_mem->cv_mem, tBout, tmp_cvB_mem->cv_y, &tBret, cv_itask);

        /* If an error occured, exit while loop */
        if (flag < 0) break;

        /* Set the time at which we will report solution and/or quadratures */
        tmp_cvB_mem->cv_tout = tBret;

      } else {

        flag = CV_SUCCESS;
        tmp_cvB_mem->cv_tout = tBn;

      }

      /* Move to next backward problem */

      tmp_cvB_mem = tmp_cvB_mem->cv_next;
    }

    /* If an error occured, return now */

    if (flag <0) {
      CVProcessError(cv_mem, flag, "CVODEA", "CVodeB", MSGCV_BACK_ERROR, tmp_cvB_mem->cv_index);
      return(flag);
    }

    /* If in CV_ONE_STEP mode, return now (flag=CV_SUCCESS or flag=CV_TSTOP_RETURN) */

    if (itaskB == CV_ONE_STEP) break;

    /* If all backward problems have succesfully reached tBout, return now */

    reachedTBout = TRUE;

    tmp_cvB_mem = cvB_mem;
    while(tmp_cvB_mem != NULL) {
      if ( sign*(tmp_cvB_mem->cv_tout - tBout) > 0 ) {
        reachedTBout = FALSE;
        break;
      }
      tmp_cvB_mem = tmp_cvB_mem->cv_next;
    }

    if ( reachedTBout ) break;

    /* Move check point in linked list to next one */

    ck_mem = next_;

  } 

  return(flag);
}



int CVodeGetB(void *cvode_mem, int which, realtype *tret, N_Vector yB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeGetB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    CVProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeGetB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 

  ca_mem = cv_mem->cv_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    CVProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeGetB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  } 

  N_VScale(ONE, cvB_mem->cv_y, yB);
  *tret = cvB_mem->cv_tout;

  return(CV_SUCCESS);
}


/*
 * CVodeGetQuadB
 */

int CVodeGetQuadB(void *cvode_mem, int which, realtype *tret, N_Vector qB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeGetQuadB", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    CVProcessError(cv_mem, CV_NO_ADJ, "CVODEA", "CVodeGetQuadB", MSGCV_NO_ADJ);
    return(CV_NO_ADJ);
  } 

  ca_mem = cv_mem->cv_adj_mem;

  /* Check the value of which */
  if ( which >= nbckpbs ) {
    CVProcessError(cv_mem, CV_ILL_INPUT, "CVODEA", "CVodeGetQuadB", MSGCV_BAD_WHICH);
    return(CV_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  } 

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  flag = CVodeGetQuad(cvodeB_mem, tret, qB);
  /*  *tret = cvB_mem->cv_tout; */

  return(flag);
}


/* 
 * =================================================================
 * PRIVATE FUNCTIONS FOR CHECK POINTS
 * =================================================================
 */

/*
 * CVAckpntInit
 *
 * This routine initializes the check point linked list with 
 * information from the initial time.
 */

static CkpntMem CVAckpntInit(CVodeMem cv_mem)
{
  CkpntMem ck_mem;

  /* Allocate space for ckdata */
  ck_mem = NULL;
  ck_mem = (CkpntMem) malloc(sizeof(struct CkpntMemRec));
  if (ck_mem == NULL) return(NULL);

  zn_[0] = NULL;
  zn_[0] = N_VClone(tempv);
  if (zn_[0] == NULL) {
    free(ck_mem); ck_mem = NULL;
    return(NULL);
  }
  
  zn_[1] = NULL;
  zn_[1] = N_VClone(tempv);
  if (zn_[1] == NULL) {
    N_VDestroy(zn_[0]);
    free(ck_mem); ck_mem = NULL;
    return(NULL);
  }

  /* zn_[qmax] was not allocated */
  zqm_ = 0;

  /* Load ckdata from cv_mem */
  N_VScale(ONE, zn[0], zn_[0]);
  t0_    = tn;
  nst_   = 0;
  q_     = 1;
  h_     = 0.0;
  
  /* Do we need to carry quadratures */
  quadr_ = quadr && errconQ;

  if (quadr_) {
    znQ_[0] = NULL;
    znQ_[0] = N_VClone(tempvQ);
    if (znQ_[0] == NULL) {
      N_VDestroy(zn_[0]);
      N_VDestroy(zn_[1]);
      free(ck_mem); ck_mem = NULL;
      return(NULL);
    }
    N_VScale(ONE, znQ[0], znQ_[0]);
  }

  /* Next in list */
  next_  = NULL;

  return(ck_mem);
}

/*
 * CVAckpntNew
 *
 * This routine allocates space for a new check point and sets 
 * its data from current values in cv_mem.
 */

static CkpntMem CVAckpntNew(CVodeMem cv_mem)
{
  CkpntMem ck_mem;
  int j, jj;
  int qmax; 

  /* Allocate space for ckdata */
  ck_mem = NULL;
  ck_mem = (CkpntMem) malloc(sizeof(struct CkpntMemRec));
  if (ck_mem == NULL) return(NULL);

  /* Test if we need to allocate space for the last zn.
     NOTE: zn(qmax) may be needed for a hot restart, if an order
     increase is deemed necessary at the first step after a check 
     point */
  qmax = cv_mem->cv_qmax;
  zqm_ = (q < qmax) ? qmax : 0;

  for (j=0; j<=q; j++) {
    zn_[j] = NULL;
    zn_[j] = N_VClone(tempv);
    if (zn_[j] == NULL) {
      for (jj=0; jj<j; jj++) N_VDestroy(zn_[jj]);
      free(ck_mem); ck_mem = NULL;
      return(NULL);
    }
  }

  if (q < qmax) {
    zn_[qmax] = NULL;
    zn_[qmax] = N_VClone(tempv);
    if (zn_[qmax] == NULL) {
      for (j=0; j<=q; j++) N_VDestroy(zn_[j]);
      free(ck_mem); ck_mem = NULL;
      return(NULL);
    }
  }

  /* Test if we need to carry quadratures */
  quadr_ = quadr && errconQ;

  if (quadr_) {

    for (j=0; j<=q; j++) {
      znQ_[j] = NULL;
      znQ_[j] = N_VClone(tempvQ);
      if(znQ_[j] == NULL) {
        for (jj=0; jj<j; jj++) N_VDestroy(znQ_[jj]);
        if (q < qmax) N_VDestroy(zn_[qmax]);
        for (j=0; j<=q; j++) N_VDestroy(zn_[j]);
        free(ck_mem); ck_mem = NULL;
        return(NULL);
      }
    }

    if ( q < qmax) {
      znQ_[qmax] = NULL;
      znQ_[qmax] = N_VClone(tempvQ);
      if (znQ_[qmax] == NULL) {
        for (j=0; j<=q; j++) N_VDestroy(znQ_[j]);
        N_VDestroy(zn_[qmax]);
        for (j=0; j<=q; j++) N_VDestroy(zn_[j]);
        free(ck_mem); ck_mem = NULL;
        return(NULL);
      }
    }
  }

  /* Load check point data from cv_mem */

  for (j=0; j<=q; j++) N_VScale(ONE, zn[j], zn_[j]);
  if ( q < qmax ) N_VScale(ONE, zn[qmax], zn_[qmax]);

  if(quadr_) {
    for (j=0; j<=q; j++) N_VScale(ONE, znQ[j], znQ_[j]);
    if ( q < qmax ) N_VScale(ONE, znQ[qmax], znQ_[qmax]);
  }

  for (j=0; j<=L_MAX; j++)     tau_[j] = tau[j];
  for (j=0; j<=NUM_TESTS; j++) tq_[j] = tq[j];
  for (j=0; j<=q; j++)         l_[j] = l[j];
  nst_       = nst;
  tretlast_  = tretlast;
  q_         = q;
  qprime_    = qprime;
  qwait_     = qwait;
  L_         = L;
  gammap_    = gammap;
  h_         = h;
  hprime_    = hprime;
  hscale_    = hscale;
  eta_       = eta;
  etamax_    = etamax;
  t0_        = tn;
  saved_tq5_ = saved_tq5;

  return(ck_mem);
}

/*
 * CVAckpntDelete
 *
 * This routine deletes the first check point in list.
 */

static void CVAckpntDelete(CkpntMem *ck_memPtr)
{
  CkpntMem tmp;
  int j;

  if (*ck_memPtr != NULL) {

    /* store head of list */
    tmp = *ck_memPtr;

    /* move head of list */
    *ck_memPtr = (*ck_memPtr)->ck_next;

    /* free N_Vectors in tmp */
    for (j=0;j<=tmp->ck_q;j++) N_VDestroy(tmp->ck_zn[j]);
    if (tmp->ck_zqm != 0) N_VDestroy(tmp->ck_zn[tmp->ck_zqm]);

    /* free N_Vectors for quadratures in tmp 
       Note that at the check point at t_initial, only znQ_[0] 
       was allocated*/
    if(tmp->ck_quadr) {
      if(tmp->ck_next != NULL) {
        for (j=0;j<=tmp->ck_q;j++) N_VDestroy(tmp->ck_znQ[j]);
        if (tmp->ck_zqm != 0) N_VDestroy(tmp->ck_znQ[tmp->ck_zqm]);
      } else {
        N_VDestroy(tmp->ck_znQ[0]);
      }
    }

    free(tmp); tmp = NULL;

  }
}

/* 
 * =================================================================
 * PRIVATE FUNCTIONS FOR BACKWARD PROBLEMS
 * =================================================================
 */

static void CVAbckpbDelete(CVodeBMem *cvB_memPtr)
{
  CVodeBMem tmp;
  void *cvode_mem;

  if (*cvB_memPtr != NULL) {

    /* Save head of the list */
    tmp = *cvB_memPtr;

    /* Move head of the list */
    *cvB_memPtr = (*cvB_memPtr)->cv_next;

    /* Free CVODES memory in tmp */
    cvode_mem = (void *)(tmp->cv_mem);
    CVodeFree(&cvode_mem);

    /* Free linear solver memory */
    if (tmp->cv_lfree != NULL) tmp->cv_lfree(tmp);

    /* Free preconditioner memory */
    if (tmp->cv_pfree != NULL) tmp->cv_pfree(tmp);

    /* Free workspace Nvector */
    N_VDestroy(tmp->cv_y);

    free(tmp); tmp = NULL;

  }

}

/* 
 * =================================================================
 * PRIVATE FUNCTIONS FOR INTERPOLATION
 * =================================================================
 */

/*
 * CVAdataStore
 *
 * This routine integrates the forward model starting at the check
 * point ck_mem and stores y and yprime at all intermediate steps.
 *
 * Return values:
 * CV_SUCCESS
 * CV_REIFWD_FAIL
 * CV_FWD_FAIL
 */

static int CVAdataStore(CVodeMem cv_mem, CkpntMem ck_mem)
{
  CVadjMem ca_mem;
  DtpntMem *dt_mem;
  realtype t;
  long int i;
  int cv_itask, flag;

  ca_mem = cv_mem->cv_adj_mem;
  dt_mem = ca_mem->dt_mem;

  /* Initialize cv_mem with data from ck_mem */
  flag = CVAckpntGet(cv_mem, ck_mem);
  if (flag != CV_SUCCESS)
    return(CV_REIFWD_FAIL);

  /* Set first structure in dt_mem[0] */
  dt_mem[0]->t = t0_;
  storePnt(cv_mem, dt_mem[0]);

  /* Decide whether TSTOP must be activated */
  if (ca_mem->ca_tstopCVodeFcall) {
    CVodeSetStopTime(cv_mem, ca_mem->ca_tstopCVodeF);
    cv_itask = CV_ONE_STEP_TSTOP;
  } else {
    cv_itask = CV_ONE_STEP;
  }

  /* Run CVode to set following structures in dt_mem[i] */
  i = 1;
  do {
    flag = CVode(cv_mem, t1_, ytmp, &t, cv_itask);
    if (flag < 0) return(CV_FWD_FAIL);
    dt_mem[i]->t = t;
    storePnt(cv_mem, dt_mem[i]);
    i++;
  } while (t<t1_);


  newData = TRUE;       /* New data is now available    */
  ckpntData = ck_mem;   /* starting at this check point */
  np = i;               /* and we have this many points */

  return(CV_SUCCESS);
}

/*
 * CVAckpntGet
 *
 * This routine prepares CVODES for a hot restart from
 * the check point ck_mem
 */

static int CVAckpntGet(CVodeMem cv_mem, CkpntMem ck_mem) 
{
  int j;
  int flag;
  int qmax;
  void *abstol;

  abstol = NULL;

  if (next_ == NULL) {

    /* In this case, we just call the reinitialization routine,
       but make sure we use the same initial stepsize as on 
       the first run. */

    CVodeSetInitStep(cv_mem, h0u);

    switch (itol) {
    case CV_SS:
      abstol = (void *) &Sabstol;
      break;
    case CV_SV:
      abstol = (void *)Vabstol;
      break;
    case CV_EE:
      abstol = (void *)efun;
      break;
    }
    flag = CVodeReInit(cv_mem, f, t0_, zn_[0], itol, reltol, abstol);
    if (flag != CV_SUCCESS) return(flag);

    if(quadr_) {
      flag = CVodeQuadReInit(cv_mem, fQ, znQ_[0]);
      if (flag != CV_SUCCESS) return(flag);
    }

  } else {
    
    qmax = cv_mem->cv_qmax;

    /* Copy parameters from check point data structure */
    nst       = nst_;
    tretlast  = tretlast_;
    q         = q_;
    qprime    = qprime_;
    qwait     = qwait_;
    L         = L_;
    gammap    = gammap_;
    h         = h_;
    hprime    = hprime_;
    hscale    = hscale_;
    eta       = eta_;
    etamax    = etamax_;
    tn        = t0_;
    saved_tq5 = saved_tq5_;
    
    /* Copy the arrays from check point data structure */
    for (j=0; j<=q; j++) N_VScale(ONE, zn_[j], zn[j]);
    if ( q < qmax ) N_VScale(ONE, zn_[qmax], zn[qmax]);
    if(quadr_) {
      for (j=0; j<=q; j++) N_VScale(ONE, znQ_[j], znQ[j]);
      if ( q < qmax ) N_VScale(ONE, znQ_[qmax], znQ[qmax]);
    }
    for (j=0; j<=L_MAX; j++)     tau[j] = tau_[j];
    for (j=0; j<=NUM_TESTS; j++) tq[j] = tq_[j];
    for (j=0; j<=q; j++)         l[j] = l_[j];
    
    /* Force a call to setup */
    forceSetup = TRUE;

  }

  return(CV_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * Functions for interpolation
 * -----------------------------------------------------------------
 */

/*
 * CVAfindIndex
 *
 * Finds the index in the array of data point strctures such that
 *     dt_mem[indx-1].t <= t < dt_mem[indx].t
 * If indx is changed from the previous invocation, then newpoint = TRUE
 *
 * If t is beyond the leftmost limit, but close enough, indx=0.
 *
 * Returns CV_SUCCESS if successful and CV_GETY_BADT if unable to
 * find indx (t is too far beyond limits).
 */

static int CVAfindIndex(CVodeMem cv_mem, realtype t, 
                        long int *indx, booleantype *newpoint)
{
  CVadjMem ca_mem;
  static long int ilast;
  DtpntMem *dt_mem;
  int sign;
  booleantype to_left, to_right;

  ca_mem = cv_mem->cv_adj_mem;
  dt_mem = ca_mem->dt_mem;

  *newpoint = FALSE;

  /* Find the direction of integration */
  sign = (tfinal - tinitial > ZERO) ? 1 : -1;

  /* If this is the first time we use new data */
  if (newData) {
    ilast     = np-1;
    *newpoint = TRUE;
    newData   = FALSE;
  }

  /* Search for indx starting from ilast */
  to_left  = ( sign*(t - dt_mem[ilast-1]->t) < ZERO);
  to_right = ( sign*(t - dt_mem[ilast]->t)   > ZERO);

  if ( to_left ) {
    /* look for a new indx to the left */

    *newpoint = TRUE;
    
    *indx = ilast;
    loop {
      if ( *indx == 0 ) break;
      if ( sign*(t - dt_mem[*indx-1]->t) <= ZERO ) (*indx)--;
      else                                         break;
    }

    if ( *indx == 0 )
      ilast = 1;
    else
      ilast = *indx;

    if ( *indx == 0 ) {
      /* t is beyond leftmost limit. Is it too far? */  
      if ( ABS(t - dt_mem[0]->t) > FUZZ_FACTOR * uround ) {
        return(CV_GETY_BADT);
      }
    }

  } else if ( to_right ) {
    /* look for a new indx to the right */

    *newpoint = TRUE;

    *indx = ilast;
    loop {
      if ( sign*(t - dt_mem[*indx]->t) > ZERO) (*indx)++;
      else                                     break;
    }

    ilast = *indx;


  } else {
    /* ilast is still OK */

    *indx = ilast;

  }

  return(CV_SUCCESS);


}

/*
 * CVodeGetAdjY
 *
 * This routine returns the interpolated forward solution at time t.
 * The user must allocate space for y.
 */

int CVodeGetAdjY(void *cvode_mem, realtype t, N_Vector y)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  int flag;

  if (cvode_mem == NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODEA", "CVodeGetAdjY", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  flag = getY(cv_mem, t, y);

  return(flag);
}

/* 
 * -----------------------------------------------------------------
 * Functions specific to cubic Hermite interpolation
 * -----------------------------------------------------------------
 */

/*
 * CVAhermiteMalloc
 *
 * This routine allocates memory for storing information at all
 * intermediate points between two consecutive check points. 
 * This data is then used to interpolate the forward solution 
 * at any other time.
 */

static booleantype CVAhermiteMalloc(CVodeMem cv_mem, long int steps)
{
  CVadjMem ca_mem;
  DtpntMem *dt_mem;
  HermiteDataMem content;
  long int i, ii=0;
  booleantype allocOK;

  allocOK = TRUE;

  ca_mem = cv_mem->cv_adj_mem;

  dt_mem = ca_mem->dt_mem;

  for (i=0; i<=steps; i++) {

    content = NULL;
    content = (HermiteDataMem) malloc(sizeof(struct HermiteDataMemRec));
    if (content == NULL) {
      ii = i;
      allocOK = FALSE;
      break;
    }

    content->y = NULL;
    content->y = N_VClone(tempv);
    if (content->y == NULL) {
      free(content); content = NULL;
      ii = i;
      allocOK = FALSE;
      break;
    }

    content->yd = NULL;
    content->yd = N_VClone(tempv);
    if (content->yd == NULL) {
      N_VDestroy(content->y);
      free(content); content = NULL;
      ii = i;
      allocOK = FALSE;
      break;
    }
    
    dt_mem[i]->content = content;

  } 

  if (!allocOK) {
    for (i=0; i<ii; i++) {
      content = (HermiteDataMem) (dt_mem[i]->content);
      N_VDestroy(content->y);
      N_VDestroy(content->yd);
      free(dt_mem[i]->content); dt_mem[i]->content = NULL;
    }
  }

  return(allocOK);
}

/*
 * CVAhermiteFree
 *
 * This routine frees the memory allocated for data storage.
 */

static void CVAhermiteFree(DtpntMem *dt_mem, long int steps)
{
  long int i;
  HermiteDataMem content;

  for (i=0; i<=steps; i++) {
    content = (HermiteDataMem) (dt_mem[i]->content);
    N_VDestroy(content->y);
    N_VDestroy(content->yd);
    free(dt_mem[i]->content); dt_mem[i]->content = NULL;
  }
}

/*
 * CVAhermiteStorePnt ( -> storePnt )
 *
 * This routine stores a new point (y,yd) in the structure d for use
 * in the cubic Hermite interpolation.
 * Note that the time is already stored.
 */

static int CVAhermiteStorePnt(CVodeMem cv_mem, DtpntMem d)
{
  int retval;
  HermiteDataMem content;

  content = (HermiteDataMem) d->content;

  N_VScale(ONE, zn[0], content->y);
  
  if (nst == 0)
    retval = f(tn, content->y, content->yd, f_data);
  else
    N_VScale(ONE/h, zn[1], content->yd);

  return(0);
}

/*
 * CVAhermiteGetY ( -> getY )
 *
 * This routine uses cubic piece-wise Hermite interpolation for 
 * the forward solution vector. 
 * It is typically called by the wrapper routines before calling
 * user provided routines (fB, djacB, bjacB, jtimesB, psolB) but
 * can be directly called by the user through CVodeGetAdjY
 */

static int CVAhermiteGetY(CVodeMem cv_mem, realtype t, N_Vector y)
{
  CVadjMem ca_mem;
  DtpntMem *dt_mem;
  HermiteDataMem content0, content1;
  realtype t0, t1, delta, factor;
  N_Vector y0, yd0, y1, yd1;
  int flag;
  long int indx;
  booleantype newpoint;

  ca_mem = cv_mem->cv_adj_mem;
  dt_mem = ca_mem->dt_mem;
  
  /* Get the index in dt_mem */

  flag = CVAfindIndex(cv_mem, t, &indx, &newpoint);
  if (flag != CV_SUCCESS) return(flag);

  /* If we are beyond the left limit but close enough,
     then return y at the left limit. */

  if (indx == 0) {
    content0 = (HermiteDataMem) (dt_mem[0]->content);
    N_VScale(ONE, content0->y, y);
    return(CV_SUCCESS);
  }

  /* Extract stuff from the appropriate data points */

  t0 = dt_mem[indx-1]->t;
  t1 = dt_mem[indx]->t;
  delta = t1 - t0;

  content0 = (HermiteDataMem) (dt_mem[indx-1]->content);
  y0  = content0->y;
  yd0 = content0->yd;

  if (newpoint) {
    
    /* Recompute Y0 and Y1 */

    content1 = (HermiteDataMem) (dt_mem[indx]->content);
    y1  = content1->y;
    yd1 = content1->yd;

    N_VLinearSum(ONE, y1, -ONE, y0, Y0);
    N_VLinearSum(ONE, yd1,  ONE, yd0, Y1);
    N_VLinearSum(delta, Y1, -TWO, Y0, Y1);
    N_VLinearSum(ONE, Y0, -delta, yd0, Y0);

  }

  /* Perform the actual interpolation. */

  factor = t - t0;
  N_VLinearSum(ONE, y0, factor, yd0, y);

  factor = factor/delta;
  factor = factor*factor;
  N_VLinearSum(ONE, y, factor, Y0, y);

  factor = factor*(t-t1)/delta;
  N_VLinearSum(ONE, y, factor, Y1, y);

  return(CV_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * Functions specific to Polynomial interpolation
 * -----------------------------------------------------------------
 */

/*
 * CVApolynomialMalloc
 *
 * This routine allocates memory for storing information at all
 * intermediate points between two consecutive check points. 
 * This data is then used to interpolate the forward solution 
 * at any other time.
 */

static booleantype CVApolynomialMalloc(CVodeMem cv_mem, long int steps)
{
  CVadjMem ca_mem;
  DtpntMem *dt_mem;
  PolynomialDataMem content;
  long int i, ii=0;
  booleantype allocOK;

  allocOK = TRUE;

  ca_mem = cv_mem->cv_adj_mem;

  dt_mem = ca_mem->dt_mem;

  for (i=0; i<=steps; i++) {
    content = NULL;
    content = (PolynomialDataMem) malloc(sizeof(struct PolynomialDataMemRec));
    if (content == NULL) {
      ii = i;
      allocOK = FALSE;
      break;
    }

    content->y = NULL;
    content->y = N_VClone(tempv);
    if (content->y == NULL) {
      free(content); content = NULL;
      ii = i;
      allocOK = FALSE;
      break;
    }

    dt_mem[i]->content = content;
  } 

  if (!allocOK) {
    for (i=0; i<ii; i++) {
      content = (PolynomialDataMem) (dt_mem[i]->content);
      N_VDestroy(content->y);
      free(dt_mem[i]->content); dt_mem[i]->content = NULL;
    }
  }

  return(allocOK);

}

/*
 * CVApolynomialFree
 *
 * This routine frees the memeory allocated for data storage.
 */

static void CVApolynomialFree(DtpntMem *dt_mem, long int steps)
{
  long int i;
  PolynomialDataMem content;

  for (i=0; i<=steps; i++) {
    content = (PolynomialDataMem) (dt_mem[i]->content);
    N_VDestroy(content->y);
    free(dt_mem[i]->content); dt_mem[i]->content = NULL;
  }
}

/*
 * CVApolynomialStorePnt ( -> storePnt )
 *
 * This routine stores a new point y in the structure d for use
 * in the Polynomial interpolation.
 * Note that the time is already stored.
 */

static int CVApolynomialStorePnt(CVodeMem cv_mem, DtpntMem d)
{
  PolynomialDataMem content;

  content = (PolynomialDataMem) d->content;

  N_VScale(ONE, zn[0], content->y);
  content->order = qu;

  return(0);
}

/*
 * CVApolynomialGetY ( -> getY )
 *
 * This routine uses polynomial interpolation for the forward solution vector. 
 * It is typically called by the wrapper routines before calling
 * user provided routines (fB, djacB, bjacB, jtimesB, psolB)) but
 * can be directly called by the user through CVodeGetAdjY.
 */

static int CVApolynomialGetY(CVodeMem cv_mem, realtype t, N_Vector y)
{
  CVadjMem ca_mem;
  DtpntMem *dt_mem;
  PolynomialDataMem content;
  int flag, dir, order, i, j;
  long int indx, base;
  booleantype newpoint;
  realtype dt, factor;

  ca_mem = cv_mem->cv_adj_mem;
  dt_mem = ca_mem->dt_mem;
  
  /* Get the index in dt_mem */

  flag = CVAfindIndex(cv_mem, t, &indx, &newpoint);
  if (flag != CV_SUCCESS) return(flag);

  /* If we are beyond the left limit but close enough,
     then return y at the left limit. */

  if (indx == 0) {
    content = (PolynomialDataMem) (dt_mem[0]->content);
    N_VScale(ONE, content->y, y);
    return(CV_SUCCESS);
  }

  /* Scaling factor */

  dt = ABS(dt_mem[indx]->t - dt_mem[indx-1]->t);

  /* Find the direction of the forward integration */

  dir = (tfinal - tinitial > ZERO) ? 1 : -1;

  /* Establish the base point depending on the integration direction.
     Modify the base if there are not enough points for the current order */

  if (dir == 1) {
    base = indx;
    content = (PolynomialDataMem) (dt_mem[base]->content);
    order = content->order;
    if(indx < order) base += order-indx;
  } else {
    base = indx-1;
    content = (PolynomialDataMem) (dt_mem[base]->content);
    order = content->order;
    if (np-indx > order) base -= indx+order-np;
  }

  /* Recompute Y (divided differences for Newton polynomial) if needed */

  if (newpoint) {

    /* Store 0-th order DD */
    if (dir == 1) {
      for(j=0;j<=order;j++) {
        T[j] = dt_mem[base-j]->t;
        content = (PolynomialDataMem) (dt_mem[base-j]->content);
        N_VScale(ONE, content->y, Y[j]);
      }
    } else {
      for(j=0;j<=order;j++) {
        T[j] = dt_mem[base-1+j]->t;
        content = (PolynomialDataMem) (dt_mem[base-1+j]->content);
        N_VScale(ONE, content->y, Y[j]);
      }
    }

    /* Compute higher-order DD */
    for(i=1;i<=order;i++) {
      for(j=order;j>=i;j--) {
        factor = dt/(T[j]-T[j-i]);
        N_VLinearSum(factor, Y[j], -factor, Y[j-1], Y[j]);
      }
    }
  }

  /* Perform the actual interpolation using nested multiplications */

  N_VScale(ONE, Y[order], y);
  for (i=order-1; i>=0; i--) {
    factor = (t-T[i])/dt;
    N_VLinearSum(factor, y, ONE, Y[i], y);
  }

  return(CV_SUCCESS);

}

/* 
 * =================================================================
 * WRAPPERS FOR ADJOINT SYSTEM
 * =================================================================
 */
/*
 * CVArhs
 *
 * This routine interfaces to the CVRhsFnB routine provided by
 * the user.
 * NOTE: f_data actually contains cvadj_mem
 */

static int CVArhs(realtype t, N_Vector yB, 
                  N_Vector yBdot, void *cvode_mem)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  int flag, retval;

  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  cvB_mem = ca_mem->ca_bckpbCrt;

  /* Forward solution from interpolation */
  flag = getY(cv_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    CVProcessError(cv_mem, -1, "CVODEA", "CVArhs", MSGCV_BAD_TINTERP, t);
    return(-1);
  }

  /* Call user's adjoint RHS routine */
  retval = (cvB_mem->cv_f)(t, ytmp, yB, yBdot, cvB_mem->cv_f_data);

  return(retval);
}

/*
 * CVArhsQ
 *
 * This routine interfaces to the CVQuadRhsFnB routine provided by
 * the user.
 * NOTE: fQ_data actually contains cvadj_mem
 */

static int CVArhsQ(realtype t, N_Vector yB, 
                   N_Vector qBdot, void *cvode_mem)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  int flag, retval;

  cv_mem = (CVodeMem) cvode_mem;

  ca_mem = cv_mem->cv_adj_mem;

  cvB_mem = ca_mem->ca_bckpbCrt;

  /* Forward solution from interpolation */
  flag = getY(cv_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    CVProcessError(cv_mem, -1, "CVODEA", "CVArhsQ", MSGCV_BAD_TINTERP, t);
    return(-1);
  }

  /* Call user's adjoint RHS routine */
  retval = (cvB_mem->cv_fQ)(t, ytmp, yB, qBdot, cvB_mem->cv_fQ_data);

  return(retval);
}



