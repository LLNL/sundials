/**************************************************************
 *                                                            *
 * File          : cvodea.c                                   *
 * Programmers   : Radu Serban @ LLNL                         *
 * Version of    : 9 January 2002                             *
 *------------------------------------------------------------*
 * This is the implementation file for the CVODEA adjoint     *
 * integrator.                                                *
 *                                                            *
 **************************************************************/

/******************* BEGIN Imports ****************************/

#include <stdio.h>
#include <stdlib.h>
#include "cvodea.h"
#include "llnlmath.h"

/******************** END Imports *****************************/


/*********************** BEGIN Macros *************************/

/* Macro: loop */

#define loop for(;;)

/************************ END Macros **************************/


/**************** BEGIN CVODEA Private Constants **************/

#define ZERO        RCONST(0.0)     /* real 0.0 */
#define ONE         RCONST(1.0)     /* real 1.0 */
#define TWO         RCONST(2.0)     /* real 2.0 */
#define FUZZ_FACTOR RCONST(10000.0) /* fuzz factor for CVadjGetY */

/***************** END CVODEA Private Constants ***************/


/********* BEGIN Private Helper Functions Prototypes **********/

static CkpntMem CVAckpntInit(CVodeMem cv_mem);
static CkpntMem CVAckpntNew(CVodeMem cv_mem);
static int  CVAckpntAdd(CVadjMem ca_mem, real tout, N_Vector yout, 
                        real *t, int itask);
static void CVAckpntDelete(CkpntMem *ck_memPtr);

static DtpntMem *CVAdataMalloc(CVodeMem cv_mem, long int steps);
static void CVAdataFree(DtpntMem *dt_mem, long int steps);

static int  CVAdataStore(CVadjMem ca_mem, int icheck, 
                         real *t0Ptr, real *t1Ptr, long int *npPtr);
static int  CVAckpntGet(CVodeMem cv_mem, CkpntMem ck_mem, 
                        DtpntMem dt_mem, int ncheck, int icheck, 
                        real *t0Ptr, real *t1Ptr);
static void CVAhermitePrepare(CVadjMem ca_mem, DtpntMem *dt_mem, 
                              long int i);
static void CVAhermiteInterpolate(CVadjMem ca_mem, DtpntMem *dt_mem,
                                  long int i, real t, N_Vector y);

static void CVArhs(integer NB, real t, N_Vector yB, 
                   N_Vector yBdot, void *passed_data);
static void CVAdenseJac(integer NB, DenseMat JB, RhsFn fB, 
                        void *cvadj_mem, real t, 
                        N_Vector yB, N_Vector fyB, N_Vector ewtB,
                        real hB, real uroundB, void *cvadj_mem_bis,
                        long int *nfePtrB, N_Vector vtemp1B,
                        N_Vector vtemp2B, N_Vector vtemp3B);
static void CVAbandJac(integer NB, integer mupperB, integer mlowerB,
                       BandMat JB, RhsFn fB, void *cvadj_mem, real t,
                       N_Vector yB, N_Vector fyB, N_Vector ewtB, 
                       real hB, real uroundB, void *cvadj_mem_bis, 
                       long int *nfePtrB, N_Vector vtemp1B, 
                       N_Vector vtemp2B, N_Vector vtemp3B);
static int CVAspgmrPrecond(integer NB, real t, N_Vector yB, 
                           N_Vector fyB, boole jokB, 
                           boole *jcurPtrB, real gammaB,
                           N_Vector ewtB, real hB, real uroundB,
                           long int *nfePtrB, void *cvadj_mem,
                           N_Vector vtemp1B, N_Vector vtemp2B,
                           N_Vector vtemp3B);
static int CVAspgmrPsolve(integer NB, real t, N_Vector yB, 
                          N_Vector fyB, N_Vector vtempB, 
                          real gammaB, N_Vector ewtB,
                          real deltaB, long int *nfePtrB, 
                          N_Vector rB, int lrB, void *cvadj_mem, 
                          N_Vector zB);
static int CVAspgmrJtimes(integer NB, N_Vector vB, N_Vector JvB, 
                          RhsFn fB, void *cvadj_mem, real t, 
                          N_Vector yB, N_Vector fyB,
                          real vnrmB, N_Vector ewtB, real hB, 
                          real uroundB, void *cvadj_mem_bis, 
                          long int *nfePtrB, N_Vector workB);

/********** END Private Helper Functions Prototypes ***********/


/**************** BEGIN Readability Constants *****************/

#define uround     (ca_mem->ca_uround)
#define tinitial   (ca_mem->ca_tinitial)
#define tfinal     (ca_mem->ca_tfinal)
#define nckpnts    (ca_mem->ca_nckpnts)
#define nsteps     (ca_mem->ca_nsteps)
#define newData    (ca_mem->ca_newData)
#define np         (ca_mem->ca_np)
#define delta      (ca_mem->ca_delta)
#define Y0         (ca_mem->ca_Y0)
#define Y1         (ca_mem->ca_Y1)
#define ytmp       (ca_mem->ca_ytmp)
#define f_B        (ca_mem->ca_fB)
#define f_data_B   (ca_mem->ca_f_dataB)
#define djac_B     (ca_mem->ca_djacB)
#define bjac_B     (ca_mem->ca_bjacB)
#define precond_B  (ca_mem->ca_precondB)
#define psolve_B   (ca_mem->ca_psolveB)
#define jtimes_B   (ca_mem->ca_jtimesB)
#define jac_data_B (ca_mem->ca_jac_dataB)
#define P_data_B   (ca_mem->ca_P_dataB)
#define ioptBalloc (ca_mem->ca_ioptBalloc)
#define roptBalloc (ca_mem->ca_roptBalloc)

#define N          (cv_mem->cv_N)
#define zn         (cv_mem->cv_zn)
#define nst        (cv_mem->cv_nst)
#define q          (cv_mem->cv_q)
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
#define tau        (cv_mem->cv_tau)
#define tq         (cv_mem->cv_tq)
#define l          (cv_mem->cv_l)
#define saved_tq5  (cv_mem->cv_saved_tq5)
#define forceSetup (cv_mem->cv_forceSetup)
#define machenv    (cv_mem->cv_machenv)
#define f          (cv_mem->cv_f)
#define lmm        (cv_mem->cv_lmm)
#define iter       (cv_mem->cv_iter)
#define itol       (cv_mem->cv_itol)
#define reltol     (cv_mem->cv_reltol)
#define abstol     (cv_mem->cv_abstol)
#define f_data     (cv_mem->cv_f_data)
#define errfp      (cv_mem->cv_errfp)
#define optIn      (cv_mem->cv_optIn)
#define iopt       (cv_mem->cv_iopt) 
#define ropt       (cv_mem->cv_ropt) 
#define h0u        (cv_mem->cv_h0u)

#define N_B        (cvb_mem->cv_N)
#define machenv_B  (cvb_mem->cv_machenv)
#define optIn_B    (cvb_mem->cv_iopt)
#define iopt_B     (cvb_mem->cv_iopt)
#define ropt_B     (cvb_mem->cv_ropt)

#define t0_        (ck_mem->ck_t0)
#define t1_        (ck_mem->ck_t1)
#define zn_        (ck_mem->ck_zn)
#define nst_       (ck_mem->ck_nst)
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

/***************** END Readability Constants ******************/


/************* BEGIN CVODE Implementation *********************/


/********* BEGIN Exported Functions Implementation ************/

/*--------------------- CVadjMalloc ---------------------------
  This routine allocates space for the global CVODEA memory
  structure.
--------------------------------------------------------------*/

void *CVadjMalloc(void *cvode_mem, long int steps)
{
  CVadjMem ca_mem;
  CVodeMem cv_mem;

  ca_mem = (CVadjMem) malloc(sizeof(struct CVadjMemRec));

  /* Attach CVODE memory for forward runs */
  cv_mem = (CVodeMem)cvode_mem;
  ca_mem->cv_mem = cv_mem;

  /* Initialize Check Points linked list */
  ca_mem->ck_mem = CVAckpntInit(cv_mem);

  /* Allocate Data Points memory */
  ca_mem->dt_mem = CVAdataMalloc(cv_mem, steps);

  /* Workspace memory */
  Y0   = N_VNew(N,machenv);
  Y1   = N_VNew(N,machenv);
  ytmp = N_VNew(N,machenv);

  /* Other entries in ca_mem */
  uround   = cv_mem->cv_uround;
  nsteps   = steps;
  tinitial = tn; 

  /* Initialize nckpnts to ZERO */
  nckpnts = 0;

  return((void *)ca_mem);
} 

/*--------------------- CVodeF -------------------------------
 This routine integrates to tout and returns solution into yout.
 In the same time, it stores check point data every 'steps' steps. 

 CVodeF can be called repeatedly by the user. The last tout
 will be used as the starting time for the backward integration.

 ncheckPtr points to the number of check points stored so far.
--------------------------------------------------------------*/

int CVodeF(void *cvadj_mem, real tout, N_Vector yout, real *t, 
           int itask, int *ncheckPtr)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Integrate to tout while loading check points */
  flag = CVAckpntAdd(ca_mem, tout, yout, t, itask);

  /* tfinal is now tout */
  tfinal = tout;

  /* Get ncheck from ca_mem */ 
  *ncheckPtr = nckpnts;

  return(flag);

}

/*--------------------- CVodeMallocB --------------------------
  CVodeMallocB allocates memory for the backward run.            
  It is essentailly a call to CVodeMalloc but with some          
  particularizations for backward integration:                   
    - it takes as an argument a pointer to the adjoint memory    
      structure and sets the cvode memory for the backward run   
      as a field of this structure. CVodeMallocB returns a 
      pointer to the CVodeMem structure to be used in setting
      the linear solver;
    - no 'initial time' is required as the backward integration  
      starts at the last tout with which CVodeF was called;      
    - the routine that provides the ODE right hand side is of a  
      different type (it also gets the solution of the forward   
      integration at the current time).    
--------------------------------------------------------------*/

int CVodeMallocB(void *cvadj_mem, integer NB, RhsFnB fB, 
                 N_Vector yB0, int lmmB, int iterB, int itolB, 
                 real *reltolB, void *abstolB, void *f_dataB, 
                 FILE *errfpB, boole optInB, 
                 long int ioptB[], real roptB[], M_Env machEnv)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  real tout;
  int j;

  ca_mem = (CVadjMem) cvadj_mem;

  if (ioptB == NULL) {
    ioptB = (long int *)malloc(OPT_SIZE*sizeof(long int));
    ioptBalloc = TRUE;
  } else {
    ioptBalloc = FALSE;
  }
  if (roptB == NULL) {
    roptB = (real *)malloc(OPT_SIZE*sizeof(real));
    roptBalloc = TRUE;
  } else {
    roptBalloc = FALSE;
  }
  if (optInB == FALSE)
    for (j=0; j<OPT_SIZE; j++) {
      ioptB[j] = 0;
      roptB[j] = 0.0;
    }
  optInB = TRUE;
  ioptB[MXHNIL] = -1;
  ioptB[ISTOP]  =  1;

  tout   = tfinal;

  cvode_mem = CVodeMalloc(NB, CVArhs, tout, yB0, lmmB, iterB, 
                          itolB, reltolB, abstolB, cvadj_mem,
                          errfpB, optInB, ioptB, roptB, machEnv);

  ca_mem->cvb_mem = (CVodeMem) cvode_mem;

  f_B      = fB;
  f_data_B = f_dataB;

  return(0);

}

/*------------------------- CVDenseB ---------------------------
 CVDenseB links the main CVODE integrator with the CVDENSE
 linear solver for the backward integration.
--------------------------------------------------------------*/

int CVDenseB(void *cvadj_mem, CVDenseJacFnB djacB, void *jac_dataB)
{
  CVadjMem ca_mem;
  CVodeMem cvb_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvb_mem = ca_mem->cvb_mem;
  
  if (djacB == NULL) {
    flag = CVDense(cvb_mem, NULL, NULL);
  } else {
    flag = CVDense(cvb_mem, CVAdenseJac, cvadj_mem);
    djac_B     = djacB;
    jac_data_B = jac_dataB;
  }

  return(flag);

}


/*------------------------- CVBandB ---------------------------
 CVBandB links the main CVODE integrator with the CVBAND
 linear solver for the backward integration.
--------------------------------------------------------------*/

int CVBandB(void *cvadj_mem, integer mupperB, integer mlowerB,
            CVBandJacFnB bjacB, void *jac_dataB)
{
  CVadjMem ca_mem;
  CVodeMem cvb_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvb_mem = ca_mem->cvb_mem;
  
  if (bjacB == NULL) {
    flag = CVBand(cvb_mem, mupperB, mlowerB, NULL, NULL);
  } else {
    flag = CVBand(cvb_mem, mupperB, mlowerB, CVAbandJac, cvadj_mem);
    bjac_B     = bjacB;
    jac_data_B = jac_dataB;
  }

  return(flag);

}

/*------------------------- CVBandPreAllocB -------------------
  CVBandPreAllocB interfaces to the CVBandPre preconditioner 
  for the backward integration. The pointer to the CVBandPreData
  structure returned by this routine can then be used together
  with the functions CVBandPrecond and CVBandPSolve in a call 
  to CVSpgmrB. 
--------------------------------------------------------------*/

CVBandPreData CVBandPreAllocB(void *cvadj_mem, integer NB, 
                              integer muB, integer mlB)
{
  CVadjMem ca_mem;
  CVodeMem cvb_mem;
  CVBandPreData bpdata;

  ca_mem = (CVadjMem) cvadj_mem;
  cvb_mem = ca_mem->cvb_mem;

  bpdata = CVBandPreAlloc(NB, CVArhs, cvadj_mem, muB, mlB, cvb_mem);
  return(bpdata);
}

/*------------------------- CVBandPrecondB --------------------
  CVBandPrecondB 
--------------------------------------------------------------*/

int CVBandPrecondB(integer NB, real t, N_Vector y, 
                   N_Vector yB, N_Vector fyB, boole jokB, 
                   boole *jcurPtrB, real gammaB,
                   N_Vector ewtB, real hB, real uroundB,
                   long int *nfePtrB, void *P_dataB,
                   N_Vector vtemp1B, N_Vector vtemp2B,
                   N_Vector vtemp3B)
{
  int flag;
  flag = CVBandPrecond(NB, t, yB, fyB, jokB, jcurPtrB, gammaB,
                       ewtB, hB, uroundB, nfePtrB, P_dataB,
                       vtemp1B, vtemp2B, vtemp3B);
  return(flag);
}

/*------------------------- CVBandPSolveB ---------------------
  CVBandPSolveB 
--------------------------------------------------------------*/

int CVBandPSolveB(integer NB, real t, N_Vector y, 
                  N_Vector yB, N_Vector fyB, 
                  N_Vector vtempB,  real gammaB, 
                  N_Vector ewtB, real deltaB, 
                  long int *nfePtrB, N_Vector rB, 
                  int lrB, void *P_dataB, N_Vector zB)
{
  int flag;
  flag = CVBandPSolve(NB, t, yB, fyB, vtempB, gammaB, ewtB, 
                      deltaB, nfePtrB, rB, lrB, P_dataB, zB);
  return(flag);
}

/*------------------------- CVSpgmrB --------------------------
 CVSpgmrB links the main CVODE integrator with the CVSPGMR
 linear solver for the backward integration.
--------------------------------------------------------------*/

int CVSpgmrB(void *cvadj_mem, int pretypeB, int gstypeB, 
             int maxlB, real deltB, CVSpgmrPrecondFnB precondB, 
             CVSpgmrPSolveFnB psolveB, void *P_dataB,
             CVSpgmrJtimesFnB jtimesB, void *jac_dataB)
{
  CVadjMem ca_mem;
  CVodeMem cvb_mem;
  int flag;
 
  CVSpgmrPrecondFn my_precond;
  CVSpgmrPSolveFn  my_psolve;
  CVSpgmrJtimesFn  my_jtimes;

  ca_mem = (CVadjMem) cvadj_mem;
  cvb_mem = ca_mem->cvb_mem;

  if (precondB == NULL) {
    my_precond = NULL;
  } else {
    my_precond = CVAspgmrPrecond;
    precond_B  = precondB;
  }

  if (psolveB == NULL) {
    my_psolve = NULL;
  } else {
    my_psolve = CVAspgmrPsolve;
    psolve_B  = psolveB;
  }

  if (jtimesB == NULL) {
    my_jtimes = NULL;
  } else {
    my_jtimes = CVAspgmrJtimes;
    jtimes_B  = jtimesB;
  }

  flag = CVSpgmr(cvb_mem, pretypeB, gstypeB, maxlB, deltB,
                 my_precond, my_psolve, cvadj_mem, my_jtimes, cvadj_mem);

  P_data_B   = P_dataB;
  jac_data_B = jac_dataB;

  return(flag);

}

/*------------------------- CVodeB ----------------------------
 This routine performs the backward integration from tfinal 
 to tinitial through a sequence of forward-backward runs in
 between consecutive check points. It returns the values of
 the adjoint variables and any existing quadrature variables
 at tinitial.
--------------------------------------------------------------*/

int CVodeB(void *cvadj_mem, N_Vector yB)
{
  CVadjMem ca_mem;
  CVodeMem cvb_mem;
  int start, flag;
  real t, t0, t1;
  long int np_actual;
  
  ca_mem  = (CVadjMem) cvadj_mem;
  cvb_mem = ca_mem->cvb_mem;

  for(start = nckpnts; start >= 0; start--) {

    /* Rerun forward from check point 'start' */
    flag = CVAdataStore(ca_mem, start, &t0, &t1, &np_actual);

    /* Run backwards */
    ropt_B[TSTOP] = t0;
    flag = CVode(cvb_mem, t0, yB, &t, NORMAL);

  } 

  return(flag);

}

/*------------------------- CVadjFree -------------------------
 This routine frees the memory allocated by CVadjMalloc.
--------------------------------------------------------------*/

void CVadjFree(void *cvadj_mem)
{
  CVadjMem ca_mem;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Delete check points one by one */
  while (ca_mem->ck_mem != NULL) {
    CVAckpntDelete(&(ca_mem->ck_mem));
  }

  /* Free vectors at each data point */
  CVAdataFree(ca_mem->dt_mem, nsteps);
  free(ca_mem->dt_mem);

  /* Free vectors in ca_mem */
  N_VFree(Y0);
  N_VFree(Y1);
  N_VFree(ytmp);

  /* Free CVODE memory for backward run */
  if(ioptBalloc) 
    free(ca_mem->cvb_mem->cv_iopt);
  if(roptBalloc)
    free(ca_mem->cvb_mem->cv_ropt);
  CVodeFree(ca_mem->cvb_mem);

  /* Free CVODEA memory */
  free(ca_mem);

}

/*---------------------- CVadjGetY -----------------------------
 This routine uses cubic piece-wise Hermite interpolation for 
 the forward solution vector. 
 It is typically called by the wrapper routines before calling
 user provided routines (fB, djacB, bjacB, jtimesB, psolB) but
 can be directly called by the user if memory for the bacward 
 run is allocated through CVODE calls and not through CVODEA
 calls.
--------------------------------------------------------------*/

int CVadjGetY(void *cvadj_mem, real t, N_Vector y)
{
  CVadjMem ca_mem;
  DtpntMem *dt_mem;
  static long int i;
  long int inew;
  boole to_left, to_right;
  real troundoff;

  ca_mem = (CVadjMem) cvadj_mem;
  dt_mem = ca_mem->dt_mem;

  if ( newData ) {
    i = np-1;
    CVAhermitePrepare(ca_mem, dt_mem, i); 
    newData = FALSE;
  }

  /* Search for inew starting from last i */ 
  to_left  = ( t < dt_mem[i-1]->t );
  to_right = ( t > dt_mem[i]->t );
  
  /* Test if t is beyond left limit */
  if ( (to_left) && (i==1) ) {
    /*troundoff = FUZZ_FACTOR*uround*(ABS(dt_mem[0]->t)+ABS(dt_mem[1]->t));*/
    troundoff = FUZZ_FACTOR*uround;
    if ( ABS(t-dt_mem[0]->t) <= troundoff ) {
      N_VScale(ONE, dt_mem[0]->y, y);
      return(GETY_OK);
    }
    else {
      printf("\n TROUBLE IN GETY \n ");
      printf("%g = ABS(t-dt_mem[0]->t) > troundoff = %g  uround = %g\n",
             ABS(t-dt_mem[0]->t), troundoff, uround);
      return(GETY_BADT);
    }
  }

  inew = i;
  if ( to_left ) {
    /* Search to the left */
    inew--;
    loop {
      if ( inew == 1 ) break;
      if ( t <= dt_mem[inew-1]->t ) inew--;
      else                          break;
    }
  } else if ( to_right ) {
    /* Search to the right */
    inew++;
    loop {
      if ( t > dt_mem[inew]->t ) inew++;
      else                       break;
    }
  }
  
  if ( inew != i )
    CVAhermitePrepare(ca_mem, dt_mem, inew);

  CVAhermiteInterpolate(ca_mem, dt_mem, inew, t, y);

  i = inew;

  return(GETY_OK);

}

/*--------------------- CVadjCheckPointsList ------------------
 This routine lists the linked list of check point structures.
 For debugging....
--------------------------------------------------------------*/

void CVadjCheckPointsList(void *cvadj_mem)
{
  CVadjMem ca_mem;
  CkpntMem ck_mem;
  int i;

  ca_mem = (CVadjMem) cvadj_mem;
  ck_mem = ca_mem->ck_mem;
  i = 0;

  while (ck_mem != NULL) {
    printf("Check point %d\n", nckpnts-i);
    printf("  address   %x\n", (int)ck_mem);
    printf("  t0        %f\n", t0_);
    printf("  t1        %f\n", t1_);
    printf("  next      %x\n", (int)next_);
    ck_mem = next_;
    i++;
  }

}

/*--------------------- CVadjDataExtract ----------------------
 This routine returns the solution stored in the data structure
 at the 'which' data point.
 For debugging....
--------------------------------------------------------------*/

void CVadjDataExtract(void *cvadj_mem, long int which, 
                      real *t, N_Vector yout, N_Vector ydout)
{
  CVadjMem ca_mem;
  DtpntMem *dt_mem;

  ca_mem = (CVadjMem) cvadj_mem;
  dt_mem = ca_mem->dt_mem;

  *t = dt_mem[which]->t;

  if (yout != NULL)
    N_VScale(ONE, dt_mem[which]->y, yout);

  if (ydout != NULL)
    N_VScale(ONE, dt_mem[which]->yd, ydout);

}

/********** END Exported Functions Implementation *************/


/******** BEGIN Private Helper Functions Implementation *******/

/*--------------------- CVAckpntInit --------------------------
 This routine initializes the check point linked list with 
 information from the initial time.
--------------------------------------------------------------*/

static CkpntMem CVAckpntInit(CVodeMem cv_mem)
{
  CkpntMem ck_mem;

  /* Allocate space for ckdata */
  ck_mem = (CkpntMem) malloc(sizeof(struct CkpntMemRec));
  zn_[0] = N_VNew(N,machenv);

  /* Load ckdata from cv_mem */
  N_VScale(ONE, zn[0], zn_[0]);
  t0_    = tn;
  q_     = 0;
  /* Next in list */
  next_  = NULL;

  return(ck_mem);

}

/*--------------------- CVAckpntAdd --------------------------
 This routine integrates the forward model from intial time to
 tout and stores check point information every 'steps'
 integration steps. It builds a linked list of check point 
 data structures with the head of the list being the last check
 point.
--------------------------------------------------------------*/

static int CVAckpntAdd(CVadjMem ca_mem, real tout, N_Vector yout, 
                       real *t, int itask) 
{
  CVodeMem cv_mem;
  CkpntMem tmp;
  int flag;

  cv_mem = ca_mem->cv_mem;

  loop {

    flag = CVode(cv_mem, tout, yout, t, ONE_STEP);
    if (flag < 0) break;

    if ( nst % nsteps == 0 ) {
      ca_mem->ck_mem->ck_t1 = *t;
      tmp = CVAckpntNew(cv_mem);
      tmp->ck_next = ca_mem->ck_mem;
      ca_mem->ck_mem = tmp;
      nckpnts++;
      forceSetup = TRUE;
    }

    if (itask == ONE_STEP) {
      ca_mem->ck_mem->ck_t1 = *t;
      break;
    }
    if ((itask==NORMAL) && (*t >= tout)) {
      CVodeDky(cv_mem, tout, 0, yout);
      ca_mem->ck_mem->ck_t1 = tout;
      break;
    }
  }

  return(flag);
 
}

/*--------------------- CVAckpntNew ---------------------------

 This routine allocates space for a new check point and sets 
 its data from current values in cv_mem.

--------------------------------------------------------------*/

static CkpntMem CVAckpntNew(CVodeMem cv_mem)
{
  CkpntMem ck_mem;
  int j;

  /* Allocate space for ckdata */
  ck_mem = (CkpntMem) malloc(sizeof(struct CkpntMemRec));
  for (j=0; j<=q; j++)
    zn_[j] = N_VNew(N,machenv);

  /* Load ckdata from cv_mem */
  for (j=0; j<=q; j++)         N_VScale(ONE, zn[j], zn_[j]);
  for (j=0; j<=L_MAX; j++)     tau_[j] = tau[j];
  for (j=0; j<=NUM_TESTS; j++) tq_[j] = tq[j];
  for (j=0; j<=q; j++)         l_[j] = l[j];
  nst_       = nst;
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

/*--------------------- CVAckpntDelete --------------------------
 This routine deletes the first check point in list.
--------------------------------------------------------------*/

static void CVAckpntDelete(CkpntMem *ck_memPtr)
{
  CkpntMem tmp;
  int j;

  if (*ck_memPtr != NULL) {
    /* store head of list */
    tmp = *ck_memPtr;
    /* move head of list */
    *ck_memPtr = (*ck_memPtr)->ck_next;
    /* free tmp */
    for (j=0;j<=tmp->ck_q;j++) N_VFree(tmp->ck_zn[j]);
    free(tmp);
  }
}

/*--------------------- CVAdataMalloc -------------------------
 This routine allocates memory for storing information at all
 intermediate points between two consecutive check points. 
 This data is then used to interpolate the forward solution 
 at any other time.
--------------------------------------------------------------*/

static DtpntMem *CVAdataMalloc(CVodeMem cv_mem, long int steps)
{
  DtpntMem *dt_mem;
  long int i;

  dt_mem = (DtpntMem *)malloc((steps+1)*sizeof(struct DtpntMemRec *));

  for (i=0; i<=steps; i++) {
    dt_mem[i] = (DtpntMem)malloc(sizeof(struct DtpntMemRec));
    dt_mem[i]->y  = N_VNew(N,machenv);
    dt_mem[i]->yd = N_VNew(N,machenv);
  } 

  return(dt_mem);

}

/*--------------------- CVAdataFree ---------------------------
 This routine frees the memeory allocated for data storage.
--------------------------------------------------------------*/

static void CVAdataFree(DtpntMem *dt_mem, long int steps)
{
  long int i;

  for (i=0; i<=steps; i++) {
    N_VFree(dt_mem[i]->y);
    N_VFree(dt_mem[i]->yd);
    free(dt_mem[i]);
  }

}

/*--------------------- CVAdataStore --------------------------
 This routine integrates the forward model starting at check
 point icheck and stores y and yprime at all intermediate 
 steps. It returns an error flag and the integration interval 
 in t0Ptr and t1Ptr.
--------------------------------------------------------------*/

int CVAdataStore(CVadjMem ca_mem, int icheck, 
                 real *t0Ptr, real *t1Ptr, long int *npPtr)
{
  CVodeMem cv_mem;
  CkpntMem ck_mem;
  DtpntMem *dt_mem;
  N_Vector y, yd;
  real t0, t1, t;
  long int i;
  int flag;

  cv_mem = ca_mem->cv_mem;
  ck_mem = ca_mem->ck_mem;
  dt_mem = ca_mem->dt_mem;

  /* Initialize cv_mem with check point (icheck) data from ck_mem
     and set first structure in the dt_mem array */
  flag = CVAckpntGet(cv_mem, ck_mem, dt_mem[0], nckpnts, icheck, &t0, &t1);

  /* Run CVode to set following structures in dt_mem */
  i = 1;
  do {
    y  = dt_mem[i]->y;
    yd = dt_mem[i]->yd;
    flag = CVode(cv_mem, t1, y, &t, ONE_STEP);
    if (flag < 0) return(flag);
    dt_mem[i]->t = t;
    CVodeDky(cv_mem, t, 1, yd);
    i++;
  } while (t<t1);

  /* Because CVode is called in ONE_STEP mode, we
     might overshoot at TOUT. Here we take care of
     that possibility... */
  if ( t != t1 ) {
    dt_mem[i-1]->t = t1;
    CVodeDky(cv_mem, t1, 0, y);
    CVodeDky(cv_mem, t1, 1, yd);
  }

  *t0Ptr = t0;
  *t1Ptr = t1;
  *npPtr = i;

  newData = TRUE;
  np  = i;

  return(flag);

}

/*--------------------- CVAckpntGet ---------------------------
 This routine extracts data from the 'icheck'-th check point
 structure in the linked list of check points and sets cv_mem
 for a hot start. In addition, it sets the first element of the
 array of data structures dt_mem.
--------------------------------------------------------------*/

static int CVAckpntGet(CVodeMem cv_mem, CkpntMem ck_mem, 
                       DtpntMem dt_mem, int ncheck, int icheck, 
                       real *t0Ptr, real *t1Ptr)
{
  int j;
  real t0;
  N_Vector y0;

  /* First, find the structure for desired check point */

  for (j=ncheck; j>icheck; j--)
    ck_mem = next_;

  *t1Ptr = t1_;
  *t0Ptr = t0_;

  if (icheck == 0) {

    /* In this case, we just call the reinitialization routine,
       but make sure we use the same initial stepsize as on 
       the first run. */

    t0 = t0_;
    y0 = zn_[0];
    if (iopt == NULL)
      iopt = (long int *)malloc(OPT_SIZE*sizeof(long int));
    if (ropt == NULL)
      ropt = (real *)malloc(OPT_SIZE*sizeof(real));
    if (optIn == FALSE)
      for (j=0; j<OPT_SIZE; j++) {
        iopt[j] = 0;
        ropt[j] = 0.0;
      }
    optIn = TRUE;
    ropt[H0] = h0u;
    CVReInit(cv_mem, f, t0, y0,
             lmm, iter, itol, reltol, abstol,
             f_data, errfp, optIn, iopt, ropt, machenv);

    /* Load t, y, and yd members of dt_mem */

    dt_mem->t = t0;
    N_VScale(ONE, y0, dt_mem->y);
    f(N, t0, y0, dt_mem->yd, f_data);

  } else {

    /* Copy parameters from check point data structure */

    nst       = nst_;
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
    
    /* Copy the history array from check point data structure */
    
    for (j=0; j<=q; j++)
      N_VScale(ONE, zn_[j], zn[j]);
    
    /* Copy other arrays from check point data structure */
    
    for (j=0; j<=L_MAX; j++)
      tau[j] = tau_[j];
    
    for (j=0; j<=NUM_TESTS; j++)
      tq[j] = tq_[j];
    
    for (j=0; j<=q; j++)
      l[j] = l_[j];
    
    /* Force a call to setup */
    
    forceSetup = TRUE;

    /* Load t, y, and yd members of dt_mem */

    dt_mem->t = tn;
    CVodeDky(cv_mem, tn, 0, dt_mem->y);
    CVodeDky(cv_mem, tn, 1, dt_mem->yd);
    
  }

  return(0);
    
}

/*--------------------- CVAhermitePrepare ----------------------
 This routine computes quantities required by the Hermite
 interpolation that are independent of the interpolation point.
--------------------------------------------------------------*/

static void CVAhermitePrepare(CVadjMem ca_mem, DtpntMem *dt_mem, 
                              long int i)
{
  real t0, t1; 
  N_Vector y0, y1, yd0, yd1;

  t0  = dt_mem[i-1]->t;
  y0  = dt_mem[i-1]->y;
  yd0 = dt_mem[i-1]->yd;

  t1  = dt_mem[i]->t;
  y1  = dt_mem[i]->y;
  yd1 = dt_mem[i]->yd;

  delta = t1 - t0;

  N_VLinearSum(ONE, y1, -ONE, y0, Y0);
  N_VLinearSum(ONE, yd1,  ONE, yd0, Y1);
  N_VLinearSum(delta, Y1, -TWO, Y0, Y1);
  N_VLinearSum(ONE, Y0, -delta, yd0, Y0);
}

/*---------- CVAhermiteInterpolate ----------------------------
 This routine performs the Hermite interpolation.
--------------------------------------------------------------*/

static void CVAhermiteInterpolate(CVadjMem ca_mem, DtpntMem *dt_mem,
                                  long int i, real t, N_Vector y)
{
  real t0, t1;
  N_Vector y0, yd0;
  real factor;

  t0  = dt_mem[i-1]->t;
  t1  = dt_mem[i]->t;
  y0  = dt_mem[i-1]->y;
  yd0 = dt_mem[i-1]->yd;

  factor = t - t0;
  N_VLinearSum(ONE, y0, factor, yd0, y);

  factor = factor/delta;
  factor = factor*factor;
  N_VLinearSum(ONE, y, factor, Y0, y);

  factor = factor*(t-t1)/delta;
  N_VLinearSum(ONE, y, factor, Y1, y);
}

/********* BEGIN wrappers for adjoint equations ***************/

/*--------------------- CVArhs ---------------------------------
 This routine interfaces to the RhsFnB routine provided by
 the user.
--------------------------------------------------------------*/

static void CVArhs(integer NB, real t, N_Vector yB, 
                   N_Vector yBdot, void *cvadj_mem)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != GETY_OK) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  }

  /* Call user's adjoint RHS routine */
  f_B(NB, t, ytmp, yB, yBdot, f_data_B);

}

/*--------------------- CVAdenseJac ----------------------------
 This routine interfaces to the CVDenseJacFnB routine provided 
 by the user.
--------------------------------------------------------------*/

static void CVAdenseJac(integer NB, DenseMat JB, RhsFn fB, 
                        void *cvadj_mem, real t, 
                        N_Vector yB, N_Vector fyB, N_Vector ewtB,
                        real hB, real uroundB, void *cvadj_mem_bis,
                        long int *nfePtrB, N_Vector vtemp1B,
                        N_Vector vtemp2B, N_Vector vtemp3B)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != GETY_OK) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  }

  /* Call user's adjoint dense djacB routine */
  djac_B(NB, JB, f_B, f_data_B, t, ytmp, yB, fyB, ewtB, hB,
         uroundB, jac_data_B, nfePtrB, vtemp1B, vtemp2B, vtemp3B);

}

/*--------------------- CVAbandJac ----------------------------
 This routine interfaces to the CVBandJacFnB routine provided 
 by the user.
--------------------------------------------------------------*/

static void CVAbandJac(integer NB, integer mupperB, integer mlowerB,
                       BandMat JB, RhsFn fB, void *cvadj_mem, real t,
                       N_Vector yB, N_Vector fyB, N_Vector ewtB, 
                       real hB, real uroundB, void *cvadj_mem_bis, 
                       long int *nfePtrB, N_Vector vtemp1B, 
                       N_Vector vtemp2B, N_Vector vtemp3B)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != GETY_OK) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  }

  /* Call user's adjoint band bjacB routine */
  bjac_B(NB, mupperB, mlowerB, JB, f_B, f_data_B, t, ytmp, 
         yB, fyB, ewtB, hB, uroundB, jac_data_B, nfePtrB, 
         vtemp1B, vtemp2B, vtemp3B);

}

/*--------------------- CVAspgmrPrecond -----------------------
 This routine interfaces to the CVSpgmrPrecondFnB routine 
 provided by the user.
--------------------------------------------------------------*/

static int CVAspgmrPrecond(integer NB, real t, N_Vector yB, 
                           N_Vector fyB, boole jokB, 
                           boole *jcurPtrB, real gammaB,
                           N_Vector ewtB, real hB, real uroundB,
                           long int *nfePtrB, void *cvadj_mem,
                           N_Vector vtemp1B, N_Vector vtemp2B,
                           N_Vector vtemp3B)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != GETY_OK) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  } 

  /* Call user's adjoint precondB routine */
  flag = precond_B(NB, t, ytmp, yB, fyB, jokB, jcurPtrB, gammaB, ewtB, hB, 
                   uroundB, nfePtrB, P_data_B, vtemp1B, vtemp2B, vtemp3B);

  return(flag);
}

/*--------------------- CVAspgmrPsolve ------------------------
 This routine interfaces to the CVSpgmrPsolveFnB routine 
 provided by the user.
--------------------------------------------------------------*/

static int CVAspgmrPsolve(integer NB, real t, N_Vector yB, 
                          N_Vector fyB, N_Vector vtempB, 
                          real gammaB, N_Vector ewtB,
                          real deltaB, long int *nfePtrB, 
                          N_Vector rB, int lrB, void *cvadj_mem, 
                          N_Vector zB)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != GETY_OK) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  } 

  /* Call user's adjoint psolveB routine */
  flag = psolve_B(NB, t, ytmp, yB, fyB, vtempB, gammaB, ewtB, deltaB, 
                  nfePtrB, rB, lrB, P_data_B, zB);

  return(flag);
}

/*--------------------- CVAspgmrJtimes ------------------------
 This routine interfaces to the CVSpgmrJtimesFnB routine 
 provided by the user.
--------------------------------------------------------------*/

static int CVAspgmrJtimes(integer NB, N_Vector vB, N_Vector JvB, 
                          RhsFn fB, void *cvadj_mem, real t, 
                          N_Vector yB, N_Vector fyB,
                          real vnrmB, N_Vector ewtB, real hB, 
                          real uroundB, void *cvadj_mem_bis, 
                          long int *nfePtrB, N_Vector workB)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != GETY_OK) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  } 

  /* Call user's adjoint jtimesB routine */
  flag = jtimes_B(NB, vB, JvB, f_B, f_data_B, t, ytmp, yB, fyB, vnrmB,
                  ewtB, hB, uroundB, jac_data_B, nfePtrB, workB);

  return(flag);
}

/********* END wrappers adjoint equations *********************/

/********* END Private Helper Functions Implementation ********/

/********* END CVODE Implementation ***************************/
