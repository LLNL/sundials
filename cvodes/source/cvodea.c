/*
 * -----------------------------------------------------------------
 * $Revision: 1.37.2.2 $
 * $Date: 2005-01-27 17:39:30 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the CVODEA adjoint integrator.
 * -----------------------------------------------------------------
 */

/*=================================================================*/
/*BEGIN             Import Header Files                            */
/*=================================================================*/

#include <stdio.h>
#include <stdlib.h>

#include "cvdiag.h"
#include "cvodea_impl.h"
#include "sundialsmath.h"
#include "sundialstypes.h"

/*=================================================================*/
/*BEGIN             Macros                                         */
/*=================================================================*/

/* Macro: loop */

#define loop for(;;)

/*=================================================================*/
/*BEGIN             CVODEA Private Constants                       */
/*=================================================================*/

#define ZERO        RCONST(0.0)     /* real 0.0 */
#define ONE         RCONST(1.0)     /* real 1.0 */
#define TWO         RCONST(2.0)     /* real 2.0 */
#define FUZZ_FACTOR RCONST(1000000.0)  /* fuzz factor for CVadjGetY */

/*=================================================================*/
/*BEGIN             Private Functions Prototypes                   */
/*=================================================================*/

static CkpntMem CVAckpntInit(CVodeMem cv_mem);
static CkpntMem CVAckpntNew(CVodeMem cv_mem);
static void CVAckpntDelete(CkpntMem *ck_memPtr);

static DtpntMem *CVAdataMalloc(CVodeMem cv_mem, long int steps);
static void CVAdataFree(DtpntMem *dt_mem, long int steps);

static int  CVAdataStore(CVadjMem ca_mem, CkpntMem ck_mem);
static int  CVAckpntGet(CVodeMem cv_mem, CkpntMem ck_mem); 
static void CVAhermitePrepare(CVadjMem ca_mem, DtpntMem *dt_mem, long int i);
static void CVAhermiteInterpolate(CVadjMem ca_mem, DtpntMem *dt_mem,
                                  long int i, realtype t, N_Vector y);

/* Wrappers */

static void CVArhs(realtype t, N_Vector yB, 
                   N_Vector yBdot, void *cvadj_mem);
static void CVAdenseJac(long int nB, DenseMat JB, realtype t, 
                        N_Vector yB, N_Vector fyB, void *cvadj_mem,
                        N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
static void CVAbandJac(long int nB, long int mupperB, 
                       long int mlowerB, BandMat JB, realtype t, 
                       N_Vector yB, N_Vector fyB, void *cvadj_mem, 
                       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
static int CVAspgmrPrecSetup(realtype t, N_Vector yB, 
                             N_Vector fyB, booleantype jokB, 
                             booleantype *jcurPtrB, realtype gammaB,
                             void *cvadj_mem,
                             N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
static int CVAspgmrPrecSolve(realtype t, N_Vector yB, N_Vector fyB,
                             N_Vector rB, N_Vector zB,
                             realtype gammaB, realtype deltaB,
                             int lrB, void *cvadj_mem, N_Vector tmpB);
static int CVAspgmrJacTimesVec(N_Vector vB, N_Vector JvB, realtype t, 
                               N_Vector yB, N_Vector fyB, 
                               void *cvadj_mem, N_Vector tmpB);
static void CVArhsQ(realtype t, N_Vector yB, 
                    N_Vector qBdot, void *cvadj_mem);

static void CVAgloc(long int NlocalB, realtype t, N_Vector yB, N_Vector gB, 
                    void *cvadj_mem);

static void CVAcfn(long int NlocalB, realtype t, N_Vector yB,
                   void *cvadj_mem);

/*=================================================================*/
/*END               Private Functions Prototypes                   */
/*=================================================================*/

/*=================================================================*/
/*BEGIN             Readibility Constants                          */
/*=================================================================*/

#define uround     (ca_mem->ca_uround)
#define tinitial   (ca_mem->ca_tinitial)
#define tfinal     (ca_mem->ca_tfinal)
#define nckpnts    (ca_mem->ca_nckpnts)
#define nsteps     (ca_mem->ca_nsteps)
#define ckpntData  (ca_mem->ca_ckpntData)
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
#define jtimes_B   (ca_mem->ca_jtimesB)
#define jac_data_B (ca_mem->ca_jac_dataB)
#define pset_B     (ca_mem->ca_psetB)
#define psolve_B   (ca_mem->ca_psolveB)
#define P_data_B   (ca_mem->ca_P_dataB)
#define fQ_B       (ca_mem->ca_fQB)
#define fQ_data_B  (ca_mem->ca_fQ_dataB)
#define gloc_B     (ca_mem->ca_glocB)
#define cfn_B      (ca_mem->ca_cfnB)
#define bbd_data_B (ca_mem->ca_bbd_dataB)
#define bp_data_B  (ca_mem->ca_bp_dataB)
#define t_for_quad (ca_mem->ca_t_for_quad)

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
#define f          (cv_mem->cv_f)
#define lmm        (cv_mem->cv_lmm)
#define iter       (cv_mem->cv_iter)
#define itol       (cv_mem->cv_itol)
#define reltol     (cv_mem->cv_reltol)
#define abstol     (cv_mem->cv_abstol)
#define f_data     (cv_mem->cv_f_data)
#define errfp      (cv_mem->cv_errfp)
#define h0u        (cv_mem->cv_h0u)
#define quadr      (cv_mem->cv_quadr)
#define errconQ    (cv_mem->cv_errconQ)
#define znQ        (cv_mem->cv_znQ)
#define itolQ      (cv_mem->cv_itolQ)
#define reltolQ    (cv_mem->cv_reltolQ)
#define abstolQ    (cv_mem->cv_abstolQ)
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

/*=================================================================*/
/*END               Readibility Constants                          */
/*=================================================================*/

/*=================================================================*/
/*BEGIN             Exported Functions                             */
/*=================================================================*/

/*------------------    CVadjMalloc      --------------------------*/
/*
  This routine allocates space for the global CVODEA memory
  structure.
*/
/*-----------------------------------------------------------------*/

void *CVadjMalloc(void *cvode_mem, long int steps)
{
  CVadjMem ca_mem;
  CVodeMem cv_mem;

  /* Check arguments */
  if (cvode_mem == NULL) {
    fprintf(stderr, MSGAM_NO_MEM);
    return(NULL);
  }
  if (steps <= 0) {
    fprintf(stderr, MSGAM_BAD_STEPS);
    return(NULL);
  }

  /* Allocate memory block */
  ca_mem = (CVadjMem) malloc(sizeof(struct CVadjMemRec));
  if (ca_mem == NULL) {
    fprintf(stderr, MSGAM_MEM_FAIL);
    return(NULL);
  }

  /* Attach CVODE memory for forward runs */
  cv_mem = (CVodeMem)cvode_mem;
  ca_mem->cv_mem = cv_mem;

  /* Initialize Check Points linked list */
  ca_mem->ck_mem = CVAckpntInit(cv_mem);
  if (ca_mem->ck_mem == NULL) {
    free(ca_mem);
    fprintf(stderr, MSGAM_MEM_FAIL);
    return(NULL);
  }

  /* Allocate Data Points memory */
  ca_mem->dt_mem = CVAdataMalloc(cv_mem, steps);
  if (ca_mem->dt_mem == NULL) {
    CVAckpntDelete(&(ca_mem->ck_mem));
    free(ca_mem);
    fprintf(stderr, MSGAM_MEM_FAIL);
    return(NULL);
  }

  /* Workspace memory */
  Y0 = N_VClone(tempv);
  if (Y0 == NULL) {
    CVAdataFree(ca_mem->dt_mem, steps);
    CVAckpntDelete(&(ca_mem->ck_mem));
    free(ca_mem);
    fprintf(stderr, MSGAM_MEM_FAIL);
    return(NULL);
  }

  Y1 = N_VClone(tempv);
  if (Y1 == NULL) {
    N_VDestroy(Y0);
    CVAdataFree(ca_mem->dt_mem, steps);
    CVAckpntDelete(&(ca_mem->ck_mem));
    free(ca_mem);
    fprintf(stderr, MSGAM_MEM_FAIL);
    return(NULL);
  }

  ytmp = N_VClone(tempv);
  if (ytmp == NULL) {
    N_VDestroy(Y0);
    N_VDestroy(Y1);
    CVAdataFree(ca_mem->dt_mem, steps);
    CVAckpntDelete(&(ca_mem->ck_mem));
    free(ca_mem);
    fprintf(stderr, MSGAM_MEM_FAIL);
    return(NULL);
  }

  /* Other entries in ca_mem */
  uround   = cv_mem->cv_uround;
  nsteps   = steps;
  tinitial = tn; 

  /* Initialize nckpnts to ZERO */
  nckpnts = 0;

  /* Initialize backward cvode memory to NULL */
  ca_mem->cvb_mem = NULL;

  ca_mem->ca_f_dataB = NULL;
  ca_mem->ca_fQ_dataB = NULL;
  ca_mem->ca_jac_dataB = NULL;
  ca_mem->ca_P_dataB = NULL;
  ca_mem->ca_bp_dataB = NULL;
  ca_mem->ca_bbd_dataB = NULL;
  

  return((void *)ca_mem);
} 

/*=================================================================*/
/*BEGIN             Wrappers for CVODEA                            */
/*=================================================================*/

/*------------------     CVodeF          --------------------------*/
/*
  This routine integrates to tout and returns solution into yout.
  In the same time, it stores check point data every 'steps' steps. 
  
  CVodeF can be called repeatedly by the user.
  
  ncheckPtr points to the number of check points stored so far.
*/
/*-----------------------------------------------------------------*/

int CVodeF(void *cvadj_mem, realtype tout, N_Vector yout, 
           realtype *tret, int itask, int *ncheckPtr)
{
  CVadjMem ca_mem;
  CVodeMem cv_mem;
  CkpntMem tmp;
  DtpntMem *dt_mem;
  int cv_itask, flag;
  booleantype iret, istop;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  cv_mem = ca_mem->cv_mem;
  dt_mem = ca_mem->dt_mem;

  iret = TRUE;
  cv_itask = CV_ONE_STEP;

  /* Interpret itask */
  switch (itask) {
  case CV_NORMAL:
    iret = FALSE;
    istop = FALSE;
    cv_itask = CV_ONE_STEP;
    break;
  case CV_ONE_STEP:
    iret = TRUE;
    istop = FALSE;
    cv_itask = CV_ONE_STEP;
    break;
  case CV_NORMAL_TSTOP:
    iret = FALSE;
    istop = TRUE;
    cv_itask = CV_ONE_STEP_TSTOP;
    break;
  case CV_ONE_STEP_TSTOP:
    iret = TRUE;
    istop = TRUE;
    cv_itask = CV_ONE_STEP_TSTOP;
    break;
  }

  /* On the first step, load dt_mem[0] */
  if ( nst == 0) {
    dt_mem[0]->t = ca_mem->ck_mem->ck_t0;
    N_VScale(ONE, ca_mem->ck_mem->ck_zn[0], dt_mem[0]->y);
    N_VScale(ONE, ca_mem->ck_mem->ck_zn[1], dt_mem[0]->yd);
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
        flag = CV_MEM_FAIL;
        break;
      }
      tmp->ck_next = ca_mem->ck_mem;
      ca_mem->ck_mem = tmp;
      nckpnts++;
      forceSetup = TRUE;
      
      /* Reset i=0 and load dt_mem[0] */
      dt_mem[0]->t = ca_mem->ck_mem->ck_t0;
      N_VScale(ONE, ca_mem->ck_mem->ck_zn[0], dt_mem[0]->y);
      N_VScale(ONE, ca_mem->ck_mem->ck_zn[1], dt_mem[0]->yd);

    } else {
      
      /* Load next point in dt_mem */
      dt_mem[nst%nsteps]->t = *tret;
      N_VScale(ONE, yout, dt_mem[nst%nsteps]->y);
      CVodeGetDky(cv_mem, *tret, 1, dt_mem[nst%nsteps]->yd);

    }

    /* Set t1 field of the current ckeck point structure
       for the case in which there will be no future
       check points */
    ca_mem->ck_mem->ck_t1 = *tret;

    /* tfinal is now set to *tret */
    tfinal = *tret;

    /* Return if in CV_ONE_STEP mode */
    if (iret) 
      break;

    /* Return if tout reached */
    if ( (*tret - tout)*h >= ZERO ) {
      *tret = tout;
      CVodeGetDky(cv_mem, tout, 0, yout);
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

/*-- CVodeCreateB, CVodeSet*B, CVodeMallocB, and CVodeReInitB -----*/
/*-----------------------------------------------------------------*/

int CVodeCreateB(void *cvadj_mem, int lmmB, int iterB)
{
  CVadjMem ca_mem;
  void *cvode_mem;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);

  ca_mem = (CVadjMem) cvadj_mem;

  cvode_mem = CVodeCreate(lmmB, iterB);

  if (cvode_mem == NULL) return(CV_MEM_FAIL);

  ca_mem->cvb_mem = (CVodeMem) cvode_mem;

  return(CV_SUCCESS);

}

/*-----------------------------------------------------------------*/

int CVodeSetIterTypeB(void *cvadj_mem, int iterB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvode_mem = (void *)ca_mem->cvb_mem;

  flag = CVodeSetIterType(cvode_mem, iterB);
  
  return(flag);
}

int CVodeSetFdataB(void *cvadj_mem, void *f_dataB)
{
  CVadjMem ca_mem;

  ca_mem = (CVadjMem) cvadj_mem;

  f_data_B = f_dataB;

  return(CV_SUCCESS);
}

int CVodeSetErrFileB(void *cvadj_mem, FILE *errfpB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvode_mem = (void *)ca_mem->cvb_mem;

  flag = CVodeSetErrFile(cvode_mem, errfpB);

  return(flag);
}

int CVodeSetMaxOrdB(void *cvadj_mem, int maxordB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvode_mem = (void *)ca_mem->cvb_mem;

  flag = CVodeSetMaxOrd(cvode_mem, maxordB);

  return(flag);
}


int CVodeSetMaxNumStepsB(void *cvadj_mem, long int mxstepsB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvode_mem = (void *)ca_mem->cvb_mem;

  flag = CVodeSetMaxNumSteps(cvode_mem, mxstepsB);

  return(flag);
}

int CVodeSetStabLimDetB(void *cvadj_mem, booleantype stldetB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvode_mem = (void *)ca_mem->cvb_mem;

  flag = CVodeSetStabLimDet(cvode_mem, stldetB);

  return(flag);
}

int CVodeSetInitStepB(void *cvadj_mem, realtype hinB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvode_mem = (void *)ca_mem->cvb_mem;

  flag = CVodeSetInitStep(cvode_mem, hinB);

  return(flag);
}

int CVodeSetMinStepB(void *cvadj_mem, realtype hminB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvode_mem = (void *)ca_mem->cvb_mem;

  flag = CVodeSetMinStep(cvode_mem, hminB);

  return(flag);
}

int CVodeSetMaxStepB(void *cvadj_mem, realtype hmaxB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvode_mem = (void *)ca_mem->cvb_mem;

  flag = CVodeSetMaxStep(cvode_mem, hmaxB);

  return(flag);
}


/*-----------------------------------------------------------------*/

int CVodeMallocB(void *cvadj_mem, CVRhsFnB fB, 
                 realtype tB0, N_Vector yB0,
                 int itolB, realtype *reltolB, void *abstolB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int sign, flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);

  ca_mem = (CVadjMem) cvadj_mem;

  sign = (tfinal - tinitial > ZERO) ? 1 : -1;
  if ( (sign*(tB0-tinitial) < ZERO) || (sign*(tfinal-tB0) < ZERO) )
    return(CV_BAD_TB0);

  f_B = fB;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVodeMalloc(cvode_mem, CVArhs, tB0, yB0,
                     itolB, reltolB, abstolB);

  if (flag != CV_SUCCESS) return(flag);

  CVodeSetMaxHnilWarns(cvode_mem, -1);
  CVodeSetFdata(cvode_mem, cvadj_mem);

  return(CV_SUCCESS);

}

/*-----------------------------------------------------------------*/

int CVodeReInitB(void *cvadj_mem, CVRhsFnB fB, 
                 realtype tB0, N_Vector yB0, 
                 int itolB, realtype *reltolB, void *abstolB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int sign, flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);

  ca_mem = (CVadjMem) cvadj_mem;

  sign = (tfinal - tinitial > ZERO) ? 1 : -1;
  if ( (sign*(tB0-tinitial) < ZERO) || (sign*(tfinal-tB0) < ZERO) )
    return(CV_BAD_TB0);
  
  f_B  = fB;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVodeReInit(cvode_mem, CVArhs, tB0, yB0,
                     itolB, reltolB, abstolB);

  if (flag != CV_SUCCESS) return(flag);

  CVodeSetMaxHnilWarns(cvode_mem, -1);
  CVodeSetFdata(cvode_mem, cvadj_mem);

  return(CV_SUCCESS);

}

/*-- CVodeSetQuad*B, CVodeQuadMallocB, and CVodeQuadReInitB -------*/
/*-----------------------------------------------------------------*/

int CVodeSetQuadErrConB(void *cvadj_mem, booleantype errconQB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvode_mem = (void *)ca_mem->cvb_mem;

  flag = CVodeSetQuadErrCon(cvode_mem, errconQB);

  return(flag);
}

int CVodeSetQuadFdataB(void *cvadj_mem, void *fQ_dataB)
{
  CVadjMem ca_mem;

  ca_mem = (CVadjMem) cvadj_mem;

  fQ_data_B = fQ_dataB;

  return(CV_SUCCESS);
}

int CVodeSetQuadTolerancesB(void *cvadj_mem, int itolQB, 
                            realtype *reltolQB, void *abstolQB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvode_mem = (void *)ca_mem->cvb_mem;

  flag = CVodeSetQuadTolerances(cvode_mem, itolQB, reltolQB, abstolQB);

  return(flag);
}

/*-----------------------------------------------------------------*/

int CVodeQuadMallocB(void *cvadj_mem, CVQuadRhsFnB fQB, N_Vector yQB0)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);

  ca_mem = (CVadjMem) cvadj_mem;

  fQ_B = fQB;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVodeQuadMalloc(cvode_mem, CVArhsQ, yQB0);
  if (flag != CV_SUCCESS) return(flag);

  flag = CVodeSetQuadFdata(cvode_mem, cvadj_mem);

  return(flag);

}

/*-----------------------------------------------------------------*/

int CVodeQuadReInitB(void *cvadj_mem, CVQuadRhsFnB fQB, N_Vector yQB0)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);

  ca_mem = (CVadjMem) cvadj_mem;

  fQ_B = fQB;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVodeQuadReInit(cvode_mem, CVArhsQ, yQB0);

  return(flag);

}

/*---------  CVDenseB and CVdenseSet*B    -------------------------*/
/*-----------------------------------------------------------------*/

int CVDenseB(void *cvadj_mem, long int nB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  cvode_mem = (void *) ca_mem->cvb_mem;
  
  flag = CVDense(cvode_mem, nB);

  return(flag);
}

int CVDenseSetJacFnB(void *cvadj_mem, CVDenseJacFnB djacB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  djac_B     = djacB;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVDenseSetJacData(cvode_mem, cvadj_mem);
  if (flag != CVDENSE_SUCCESS) return(flag);

  CVDenseSetJacFn(cvode_mem, CVAdenseJac);

  return(CVDENSE_SUCCESS);
}

int CVDenseSetJacDataB(void *cvadj_mem, void *jac_dataB)
{
  CVadjMem ca_mem;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  jac_data_B = jac_dataB;

  return(CVDENSE_SUCCESS);
}

/*-----------------  CVDiagB    -----------------------------------*/
/*-----------------------------------------------------------------*/

int CVDiagB(void *cvadj_mem)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  cvode_mem = (void *) ca_mem->cvb_mem;
  
  flag = CVDiag(cvode_mem);

  return(flag);
}

/*-----------  CVBandB and CVBandSet*B      -----------------------*/
/*-----------------------------------------------------------------*/

int CVBandB(void *cvadj_mem, long int nB, 
            long int mupperB, long int mlowerB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  cvode_mem = (void *) ca_mem->cvb_mem;
  
  flag = CVBand(cvode_mem, nB, mupperB, mlowerB);

  return(flag);
}

int CVBandSetJacFnB(void *cvadj_mem, CVBandJacFnB bjacB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  bjac_B     = bjacB;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVBandSetJacData(cvode_mem, cvadj_mem);
  if (flag != CVBAND_SUCCESS) return(flag);

  CVBandSetJacFn(cvode_mem, CVAbandJac);

  return(CVBAND_SUCCESS);
}

int CVBandSetJacDataB(void *cvadj_mem, void *jac_dataB)
{ 
  CVadjMem ca_mem;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  jac_data_B = jac_dataB;

  return(CVBAND_SUCCESS);
}

/*------------   CVSpgmrB and CVSpgmrSet*B    ---------------------*/
/*-----------------------------------------------------------------*/
int CVSpgmrB(void *cvadj_mem, int pretypeB, int maxlB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  cvode_mem = (void *) ca_mem->cvb_mem;
  
  flag = CVSpgmr(cvode_mem, pretypeB, maxlB);

  return(flag);
}

int CVSpgmrSetPrecTypeB(void *cvadj_mem, int pretypeB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVSpgmrSetPrecType(cvode_mem, pretypeB);

  return(flag);
}

int CVSpgmrSetGSTypeB(void *cvadj_mem, int gstypeB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVSpgmrSetGSType(cvode_mem,gstypeB);

  return(flag);
}

int CVSpgmrSetDeltB(void *cvadj_mem, realtype deltB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVSpgmrSetDelt(cvode_mem,deltB);

  return(flag);
}

int CVSpgmrSetPrecSetupFnB(void *cvadj_mem, CVSpgmrPrecSetupFnB psetB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  pset_B = psetB;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVSpgmrSetPrecData(cvode_mem, cvadj_mem);
  if (flag != CVSPGMR_SUCCESS) return(flag);

  CVSpgmrSetPrecSetupFn(cvode_mem, CVAspgmrPrecSetup);

  return(CVSPGMR_SUCCESS);
}

int CVSpgmrSetPrecSolveFnB(void *cvadj_mem, CVSpgmrPrecSolveFnB psolveB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  psolve_B = psolveB;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVSpgmrSetPrecData(cvode_mem, cvadj_mem);
  if (flag != CVSPGMR_SUCCESS) return(flag);

  CVSpgmrSetPrecSolveFn(cvode_mem, CVAspgmrPrecSolve);

  return(CVSPGMR_SUCCESS);
}

int CVSpgmrSetJacTimesVecFnB(void *cvadj_mem, CVSpgmrJacTimesVecFnB jtimesB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  jtimes_B = jtimesB;

  cvode_mem = (void *) ca_mem->cvb_mem;

  flag = CVSpgmrSetJacData(cvode_mem, cvadj_mem);
  if (flag != CVSPGMR_SUCCESS) return(flag);

  CVSpgmrSetJacTimesVecFn(cvode_mem, CVAspgmrJacTimesVec);

  return(CVSPGMR_SUCCESS);
}

int CVSpgmrSetPrecDataB(void *cvadj_mem, void *P_dataB)
{
  CVadjMem ca_mem;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  P_data_B = P_dataB;

  return(CVSPGMR_SUCCESS);
}

int CVSpgmrSetJacDataB(void *cvadj_mem, void *jac_dataB)
{
  CVadjMem ca_mem;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  jac_data_B = jac_dataB;

  return(CVSPGMR_SUCCESS);
}

/*- CVBandPrecAllocB, CVBPSpgmrB, CVBandPrecFreeB                      -*/
/*----------------------------------------------------------------------*/

int CVBandPrecAllocB(void *cvadj_mem, long int nB, 
                     long int muB, long int mlB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  void *bp_dataB;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  cvode_mem = (void *) ca_mem->cvb_mem;

  bp_dataB = CVBandPrecAlloc(cvode_mem, nB, muB, mlB);

  if (bp_dataB == NULL) return(CV_PDATA_NULL);

  bp_data_B = bp_dataB;

  return(CV_SUCCESS);

}

/*-----------------------------------------------------------------*/

int CVBPSpgmrB(void *cvadj_mem, int pretypeB, int maxlB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  cvode_mem = (void *) ca_mem->cvb_mem;
  
  flag = CVBPSpgmr(cvode_mem, pretypeB, maxlB, bp_data_B);

  return(flag);
}

/*- CVBBDPrecAllocB, CVBPSpgmrB, CVBandPrecFreeB                       -*/
/*----------------------------------------------------------------------*/

int CVBBDPrecAllocB(void *cvadj_mem, long int NlocalB, 
                    long int mudqB, long int mldqB, 
                    long int mukeepB, long int mlkeepB, 
                    realtype dqrelyB,
                    CVLocalFnB glocB, CVCommFnB cfnB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  void *bbd_dataB;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;

  cvode_mem = (void *) ca_mem->cvb_mem;

  gloc_B = glocB;
  cfn_B  = cfnB;

  bbd_dataB = CVBBDPrecAlloc(cvode_mem, NlocalB, 
                             mudqB, mldqB,
                             mukeepB, mlkeepB, 
                             dqrelyB, 
                             CVAgloc, CVAcfn);

  if (bbd_dataB == NULL) return(CV_PDATA_NULL);

  bbd_data_B = bbd_dataB;

  return(CV_SUCCESS);

}

int CVBBDSpgmrB(void *cvadj_mem, int pretypeB, int maxlB)
{

  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;
  
  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;
  
  cvode_mem = (void *) ca_mem->cvb_mem;
  
  flag = CVBBDSpgmr(cvode_mem, pretypeB, maxlB, bbd_data_B);

  return(flag);

}

int CVBBDPrecReInitB(void *cvadj_mem, long int mudqB, long int mldqB,
                     realtype dqrelyB, CVLocalFnB glocB, CVCommFnB cfnB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;

  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem = (CVadjMem) cvadj_mem;
  
  cvode_mem = (void *) ca_mem->cvb_mem;
  
  gloc_B = glocB;
  cfn_B  = cfnB;

  flag = CVBBDPrecReInit(bbd_data_B, mudqB, mldqB,
                         dqrelyB, CVAgloc, CVAcfn);

  return(flag);
}

/*------------------     CVodeB          --------------------------*/
/*
  This routine performs the backward integration towards tBout. 
  When necessary, it performs a forward integration between two 
  consecutive check points to update interpolation data.
  itask can be CV_NORMAL or CV_ONE_STEP only.
*/
/*-----------------------------------------------------------------*/

int CVodeB(void *cvadj_mem, realtype tBout, N_Vector yBout, 
           realtype *tBret, int itaskB)
{
  CVadjMem ca_mem;
  CkpntMem ck_mem;
  CVodeMem cvb_mem;
  int sign, flag, cv_itask;
  realtype tBn;
  
  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem  = (CVadjMem) cvadj_mem;

  cvb_mem = ca_mem->cvb_mem;
  if (cvb_mem == NULL) return(CV_BCKMEM_NULL);

  if (itaskB == CV_NORMAL)
    cv_itask = CV_NORMAL_TSTOP;
  else if (itaskB == CV_ONE_STEP)
    cv_itask = CV_ONE_STEP_TSTOP;
  else
    return(CV_BAD_ITASK);

  ck_mem = ca_mem->ck_mem;

  sign = (tfinal - tinitial > ZERO) ? 1 : -1;

  if ( (sign*(tBout-tinitial) < ZERO) || (sign*(tfinal-tBout) < ZERO) )
    return(CV_BAD_TBOUT);

  tBn = cvb_mem->cv_tn;
  while ( sign*(tBn - t0_) <= ZERO ) ck_mem = next_;
  
  loop {

    /* Store interpolation data if not available */
    if (ck_mem != ckpntData) {
      flag = CVAdataStore(ca_mem, ck_mem);
      if (flag != CV_SUCCESS) return(flag);
    }

    /* Backward integration */
    CVodeSetStopTime((void *)cvb_mem, t0_);
    flag = CVode(cvb_mem, tBout, yBout, tBret, cv_itask);

    /* If an error occured, return now */
    if (flag < 0) return(flag);

    /* Set the time at which CVodeGetQuadB will evaluate any quadratures */
    t_for_quad = *tBret;

    /* If in CV_ONE_STEP mode, return now (flag=CV_SUCCESS or flag=CV_TSTOP_RETURN) */
    if (itaskB == CV_ONE_STEP) return(flag);

    /* If succesfully reached tBout, return now */
    if (*tBret == tBout) return(flag);

    /* Move check point in linked list to next one */
    ck_mem = next_;

  } 


  return(CV_SUCCESS);

}

/*------------------  CVodeGetQuadB  ------------------------------*/
/*-----------------------------------------------------------------*/

int CVodeGetQuadB(void *cvadj_mem, N_Vector qB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int flag;
  
  if (cvadj_mem == NULL) return(CV_ADJMEM_NULL);
  ca_mem  = (CVadjMem) cvadj_mem;
  cvode_mem = (void *) ca_mem->cvb_mem;
  
  flag = CVodeGetQuad(cvode_mem, t_for_quad, qB);

  return(flag);
}

/*=================================================================*/
/*END               Wrappers for CVODEA                            */
/*=================================================================*/

/*------------------     CVAdjFree       --------------------------*/
/*
  This routine frees the memory allocated by CVadjMalloc.
*/
/*-----------------------------------------------------------------*/

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
  N_VDestroy(Y0);
  N_VDestroy(Y1);
  N_VDestroy(ytmp);

  /* Free CVODES memory for backward run */
  CVodeFree(ca_mem->cvb_mem);

  /* Free preconditioner data (the routines below check for non-NULL data) */
  CVBandPrecFree(bp_data_B);
  CVBBDPrecFree(bbd_data_B);

  /* Free CVODEA memory */
  free(ca_mem);

}

/*------------------  CVadjGetCVodeBmem  --------------------------*/
/*
  CVadjGetCVodeBmem returns a (void *) pointer to the CVODES     
  memory allocated for the backward problem. This pointer can    
  then be used to call any of the CVodeGet* CVODES routines to  
  extract optional output for the backward integration phase.
*/
/*-----------------------------------------------------------------*/

void *CVadjGetCVodeBmem(void *cvadj_mem)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  
  if (cvadj_mem == NULL) return(NULL);
  ca_mem  = (CVadjMem) cvadj_mem;
  cvode_mem = (void *) ca_mem->cvb_mem;

  return(cvode_mem);
}

/*------------------     CVAdjGetY       --------------------------*/
/*
  This routine uses cubic piece-wise Hermite interpolation for 
  the forward solution vector. 
  It is typically called by the wrapper routines before calling
  user provided routines (fB, djacB, bjacB, jtimesB, psolB) but
  can be directly called by the user if memory for the bacward 
  run is allocated through CVODE calls and not through CVODEA
  calls.
*/
/*-----------------------------------------------------------------*/

int CVadjGetY(void *cvadj_mem, realtype t, N_Vector y)
{
  CVadjMem ca_mem;
  DtpntMem *dt_mem;
  static long int i;
  long int inew;
  int sign;
  booleantype to_left, to_right;
  realtype troundoff;

  ca_mem = (CVadjMem) cvadj_mem;
  dt_mem = ca_mem->dt_mem;

  sign = (tfinal - tinitial > ZERO) ? 1 : -1;

  if ( newData ) {
    i = np-1;
    CVAhermitePrepare(ca_mem, dt_mem, i); 
    newData = FALSE;
  }

  /* Search for inew starting from last i */ 
  to_left  = ( sign*(t - dt_mem[i-1]->t) < ZERO);
  to_right = ( sign*(t - dt_mem[i]->t)   > ZERO);
  
  /* Test if t is beyond left limit */
  if ( (to_left) && (i==1) ) {
    /*troundoff = FUZZ_FACTOR*uround*(ABS(dt_mem[0]->t)+ABS(dt_mem[1]->t));*/
    troundoff = FUZZ_FACTOR*uround;
    if ( ABS(t - dt_mem[0]->t) <= troundoff ) {
      N_VScale(ONE, dt_mem[0]->y, y);
      return(CV_SUCCESS);
    }
    else {
      printf("\n TROUBLE IN GETY\n ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
      printf("%Lg = ABS(t-dt_mem[0]->t) > troundoff = %Lg  uround = %Lg\n",
             ABS(t - dt_mem[0]->t), troundoff, uround);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      printf("%lg = ABS(t-dt_mem[0]->t) > troundoff = %lg  uround = %lg\n",
             ABS(t - dt_mem[0]->t), troundoff, uround);
#else
      printf("%g = ABS(t-dt_mem[0]->t) > troundoff = %g  uround = %g\n",
             ABS(t - dt_mem[0]->t), troundoff, uround);
#endif
      return(CV_GETY_BADT);
    }
  }

  inew = i;
  if ( to_left ) {
    /* Search to the left */
    inew--;
    loop {
      if ( inew == 1 ) break;
      if ( sign*(t - dt_mem[inew-1]->t) <= ZERO) inew--;
      else                          break;
    }
  } else if ( to_right ) {
    /* Search to the right */
    inew++;
    loop {
      if ( sign*(t - dt_mem[inew]->t) > ZERO) inew++;
      else                       break;
    }
  }
  
  if ( inew != i )
    CVAhermitePrepare(ca_mem, dt_mem, inew);

  CVAhermiteInterpolate(ca_mem, dt_mem, inew, t, y);

  i = inew;

  return(CV_SUCCESS);

}

/*------------------  CVAdjGetCheckPointsList ---------------------*/
/*
  This routine lists the linked list of check point structures.
  For debugging....
*/
/*-----------------------------------------------------------------*/

void CVadjGetCheckPointsList(void *cvadj_mem)
{
  CVadjMem ca_mem;
  CkpntMem ck_mem;
  int i;

  ca_mem = (CVadjMem) cvadj_mem;
  ck_mem = ca_mem->ck_mem;
  i = 0;

  while (ck_mem != NULL) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%2d  addr: %p  time = [ %9.3Le %9.3Le ]  next: %p\n", 
           nckpnts-i, (void *)ck_mem, t0_, t1_, (void *)next_ );
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%2d  addr: %p  time = [ %9.3le %9.3le ]  next: %p\n", 
           nckpnts-i, (void *)ck_mem, t0_, t1_, (void *)next_ );
#else
    printf("%2d  addr: %p  time = [ %9.3e %9.3e ]  next: %p\n", 
           nckpnts-i, (void *)ck_mem, t0_, t1_, (void *)next_ );
#endif
    ck_mem = next_;
    i++;
  }

}

/*------------------   CVAdjGetStoredData  ------------------------*/
/*
  This routine returns the solution stored in the data structure
  at the 'which' data point.
  For debugging....
*/
/*-----------------------------------------------------------------*/

void CVadjGetStoredData(void *cvadj_mem, long int which, 
                        realtype *t, N_Vector yout, N_Vector ydout)
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

/*=================================================================*/
/*BEGIN             Exported Functions                             */
/*=================================================================*/

/*=================================================================*/
/*BEGIN             Private Functions Implementation               */
/*=================================================================*/

/*------------------     CVAckpntInit    --------------------------*/
/*
  This routine initializes the check point linked list with 
  information from the initial time.
*/
/*-----------------------------------------------------------------*/

static CkpntMem CVAckpntInit(CVodeMem cv_mem)
{
  CkpntMem ck_mem;

  /* Allocate space for ckdata */
  ck_mem = (CkpntMem) malloc(sizeof(struct CkpntMemRec));

  zn_[0] = N_VClone(tempv);
  zn_[1] = N_VClone(tempv);

  /* zn_[qmax] was not allocated */
  zqm_ = 0;

  /* Load ckdata from cv_mem */
  N_VScale(ONE, zn[0], zn_[0]);
  t0_    = tn;
  q_     = 1;
  /* Compute zn_[1] by calling the user f routine */
  f(t0_, zn_[0], zn_[1], f_data);
  
  /* Do we need to carry quadratures */
  quadr_ = quadr && errconQ;

  if (quadr_) {
    znQ_[0] = N_VClone(tempvQ);
    N_VScale(ONE, znQ[0], znQ_[0]);
  }

  /* Next in list */
  next_  = NULL;

  return(ck_mem);
}

/*------------------    CVAckpntNew      --------------------------*/
/*
  This routine allocates space for a new check point and sets 
  its data from current values in cv_mem.
*/
/*-----------------------------------------------------------------*/

static CkpntMem CVAckpntNew(CVodeMem cv_mem)
{
  CkpntMem ck_mem;
  int j;
  int qmax; 

  /* Allocate space for ckdata */
  ck_mem = (CkpntMem) malloc(sizeof(struct CkpntMemRec));
  if (ck_mem == NULL) return(NULL);

  /* Test if we need to allocate space for the last zn.
     NOTE: zn(qmax) may be needed for a hot restart, if an order
     increase is deemed necessary at the first step after a check 
     point */
  qmax = cv_mem->cv_qmax;
  zqm_ = (q < qmax) ? qmax : 0;

  for (j=0; j<=q; j++) {
    zn_[j] = N_VClone(tempv);
    if(zn_[j] == NULL) return(NULL);
  }

  if ( q < qmax) {
    zn_[qmax] = N_VClone(tempv);
    if ( zn_[qmax] == NULL ) return(NULL);
  }

  /* Test if we need to carry quadratures */
  quadr_ = quadr && errconQ;

  if (quadr_) {
    for (j=0; j<=q; j++) {
      znQ_[j] = N_VClone(tempvQ);
      if(znQ_[j] == NULL) return(NULL);
    }

    if ( q < qmax) {
      znQ_[qmax] = N_VClone(tempvQ);
      if ( znQ_[qmax] == NULL ) return(NULL);
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

/*------------------    CVAckpntDelete   --------------------------*/
/*
  This routine deletes the first check point in list.
*/
/*-----------------------------------------------------------------*/

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

    free(tmp);

  }

}

/*------------------    CVAdataMalloc    --------------------------*/
/*
  This routine allocates memory for storing information at all
  intermediate points between two consecutive check points. 
  This data is then used to interpolate the forward solution 
  at any other time.
*/
/*-----------------------------------------------------------------*/

static DtpntMem *CVAdataMalloc(CVodeMem cv_mem, long int steps)
{
  DtpntMem *dt_mem;
  long int i;

  dt_mem = (DtpntMem *)malloc((steps+1)*sizeof(struct DtpntMemRec *));

  for (i=0; i<=steps; i++) {
    dt_mem[i] = (DtpntMem)malloc(sizeof(struct DtpntMemRec));
    dt_mem[i]->y  = N_VClone(tempv);
    dt_mem[i]->yd = N_VClone(tempv);
  } 

  return(dt_mem);

}

/*------------------    CVAdataFree      --------------------------*/
/*
  This routine frees the memeory allocated for data storage.
*/
/*-----------------------------------------------------------------*/

static void CVAdataFree(DtpntMem *dt_mem, long int steps)
{
  long int i;

  for (i=0; i<=steps; i++) {
    N_VDestroy(dt_mem[i]->y);
    N_VDestroy(dt_mem[i]->yd);
    free(dt_mem[i]);
  }

}

/*------------------    CVAdataStore     --------------------------*/
/*
  This routine integrates the forward model starting at the check
  point ck_mem and stores y and yprime at all intermediate steps.

  Return values:
  CV_SUCCESS
  CV_REIFWD_FAIL
  CV_FWD_FAIL
*/
/*-----------------------------------------------------------------*/

int CVAdataStore(CVadjMem ca_mem, CkpntMem ck_mem)
{
  CVodeMem cv_mem;
  DtpntMem *dt_mem;
  realtype t;
  long int i;
  int flag;

  cv_mem = ca_mem->cv_mem;
  dt_mem = ca_mem->dt_mem;

  /* Initialize cv_mem with data from ck_mem */
  flag = CVAckpntGet(cv_mem, ck_mem);
  if (flag != CV_SUCCESS) return(CV_REIFWD_FAIL);

  /* Set first structure in dt_mem[0] */
  dt_mem[0]->t = t0_;
  N_VScale(ONE, zn_[0], dt_mem[0]->y);
  N_VScale(ONE, zn_[1], dt_mem[0]->yd);

  /* Run CVode to set following structures in dt_mem[i] */

  i = 1;
  do {
    flag = CVode(cv_mem, t1_, dt_mem[i]->y, &t, CV_ONE_STEP);
    if (flag < 0) return(CV_FWD_FAIL);
    dt_mem[i]->t = t;
    flag = CVodeGetDky(cv_mem, t, 1, dt_mem[i]->yd);
    if (flag != CV_SUCCESS) return(CV_FWD_FAIL);
    i++;
  } while (t<t1_);

  /* New data is now available */
  ckpntData = ck_mem;
  newData = TRUE;
  np  = i;

  return(CV_SUCCESS);

}

/*------------------   CVAckpntGet       --------------------------*/
/*
  This routine extracts data from the check point structure at 
  ck_mem and sets cv_mem for a hot start.
*/
/*-----------------------------------------------------------------*/

static int CVAckpntGet(CVodeMem cv_mem, CkpntMem ck_mem) 
{
  int j;
  int flag;
  int qmax;

  if (next_ == NULL) {

    /* In this case, we just call the reinitialization routine,
       but make sure we use the same initial stepsize as on 
       the first run. */

    CVodeSetInitStep(cv_mem, h0u);

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

/*------------------   CVAhermitePrepare --------------------------*/
/*
  This routine computes quantities required by the Hermite
  interpolation that are independent of the interpolation point.
*/
/*-----------------------------------------------------------------*/

static void CVAhermitePrepare(CVadjMem ca_mem, DtpntMem *dt_mem, long int i)
{
  realtype t0, t1; 
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

/*------------------   CVAhermiteInterpolate ----------------------*/
/*
  This routine performs the Hermite interpolation.
*/
/*-----------------------------------------------------------------*/

static void CVAhermiteInterpolate(CVadjMem ca_mem, DtpntMem *dt_mem,
                                  long int i, realtype t, N_Vector y)
{
  realtype t0, t1;
  N_Vector y0, yd0;
  realtype factor;

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

/*=================================================================*/
/*BEGIN        Wrappers for adjoint system                         */
/*=================================================================*/

/*------------------       CVArhs        --------------------------*/
/*
  This routine interfaces to the CVRhsFnB routine provided by
  the user.
  NOTE: f_data actually contains cvadj_mem
*/
/*-----------------------------------------------------------------*/

static void CVArhs(realtype t, N_Vector yB, 
                   N_Vector yBdot, void *cvadj_mem)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  }

  /* Call user's adjoint RHS routine */
  f_B(t, ytmp, yB, yBdot, f_data_B);

}

/*------------------       CVArhsQ       --------------------------*/
/*
  This routine interfaces to the CVQuadRhsFnB routine provided by
  the user.
  NOTE: fQ_data actually contains cvadj_mem
*/
/*-----------------------------------------------------------------*/

static void CVArhsQ(realtype t, N_Vector yB, 
                    N_Vector qBdot, void *cvadj_mem)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  }

  /* Call user's adjoint RHS routine */
  fQ_B(t, ytmp, yB, qBdot, fQ_data_B);

}

/*------------------    CVAdenseJac      --------------------------*/
/*
  This routine interfaces to the CVDenseJacFnB routine provided 
  by the user.
  NOTE: jac_data actually contains cvadj_mem
*/
/*-----------------------------------------------------------------*/

static void CVAdenseJac(long int nB, DenseMat JB, realtype t, 
                        N_Vector yB, N_Vector fyB, void *cvadj_mem,
                        N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  }

  /* Call user's adjoint dense djacB routine */
  djac_B(nB, JB, t, ytmp, yB, fyB, jac_data_B, 
         tmp1B, tmp2B, tmp3B);

}

/*------------------    CVAbandJac       --------------------------*/
/*
  This routine interfaces to the CVBandJacFnB routine provided 
  by the user.
  NOTE: jac_data actually contains cvadj_mem
*/
/*-----------------------------------------------------------------*/

static void CVAbandJac(long int nB, long int mupperB, 
                       long int mlowerB, BandMat JB, realtype t, 
                       N_Vector yB, N_Vector fyB, void *cvadj_mem, 
                       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  }

  /* Call user's adjoint band bjacB routine */
  bjac_B(nB, mupperB, mlowerB, JB, t, ytmp, yB, fyB, jac_data_B,
         tmp1B, tmp2B, tmp3B);

}

/*------------------   CVAspgmrPrecSetup   ------------------------*/
/*
  This routine interfaces to the CVSpgmrPrecSetupFnB routine 
  provided by the user.
  NOTE: p_data actually contains cvadj_mem
*/
/*-----------------------------------------------------------------*/

static int CVAspgmrPrecSetup(realtype t, N_Vector yB, 
                             N_Vector fyB, booleantype jokB, 
                             booleantype *jcurPtrB, realtype gammaB,
                             void *cvadj_mem,
                             N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  } 

  /* Call user's adjoint precondB routine */
  flag = pset_B(t, ytmp, yB, fyB, jokB, jcurPtrB, gammaB,
                P_data_B, tmp1B, tmp2B, tmp3B);

  return(flag);
}

/*----------------   CVAspgmrPrecSolve    -------------------------*/
/*
  This routine interfaces to the CVSpgmrPrecSolveFnB routine 
  provided by the user.
  NOTE: p_data actually contains cvadj_mem
*/
/*-----------------------------------------------------------------*/

static int CVAspgmrPrecSolve(realtype t, N_Vector yB, N_Vector fyB,
                             N_Vector rB, N_Vector zB,
                             realtype gammaB, realtype deltaB,
                             int lrB, void *cvadj_mem, N_Vector tmpB)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  } 

  /* Call user's adjoint psolveB routine */
  flag = psolve_B(t, ytmp, yB, fyB, rB, zB, gammaB, deltaB, 
                  lrB, P_data_B, tmpB);

  return(flag);
}

/*------------------   CVAspgmrJacTimesVec    ---------------------*/
/*
  This routine interfaces to the CVSpgmrJacTimesVecFnB routine 
  provided by the user.
  NOTE: jac_data actually contains cvadj_mem
*/
/*-----------------------------------------------------------------*/

static int CVAspgmrJacTimesVec(N_Vector vB, N_Vector JvB, realtype t, 
                               N_Vector yB, N_Vector fyB, 
                               void *cvadj_mem, N_Vector tmpB)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  } 

  /* Call user's adjoint jtimesB routine */
  flag = jtimes_B(vB, JvB, t, ytmp, yB, fyB, jac_data_B, tmpB);

  return(flag);
}

/*-------------------   CVAgloc    --------------------------------*/
/*
  This routine interfaces to the CVLocalFnB routine 
  provided by the user.
  NOTE: f_data actually contains cvadj_mem
*/
/*-----------------------------------------------------------------*/

static void CVAgloc(long int NlocalB, realtype t, N_Vector yB, N_Vector gB, 
                    void *cvadj_mem)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  } 

  /* Call user's adjoint glocB routine */
  gloc_B(NlocalB, t, ytmp, yB, gB, f_data_B);

}

/*-------------------   CVAcfn    ---------------------------------*/
/*
  This routine interfaces to the CVCommFnB routine 
  provided by the user.
  NOTE: f_data actually contains cvadj_mem
*/
/*-----------------------------------------------------------------*/

static void CVAcfn(long int NlocalB, realtype t, N_Vector yB,
                   void *cvadj_mem)
{
  CVadjMem ca_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;

  /* Forward solution from Hermite interpolation */
  flag = CVadjGetY(ca_mem, t, ytmp);
  if (flag != CV_SUCCESS) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  } 

  /* Call user's adjoint cfnB routine */
  if(cfn_B != NULL)
    cfn_B(NlocalB, t, ytmp, yB, f_data_B);

}


/*=================================================================*/
/*END               Wrappers for adjoint system                    */
/*=================================================================*/

/*=================================================================*/
/*END               Private Functions Implementation               */
/*=================================================================*/
