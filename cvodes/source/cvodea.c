/*******************************************************************
 *                                                                 *
 * File          : cvodea.c                                        *
 * Programmers   : Radu Serban @ LLNL                              *
 * Version of    : 31 March 2003                                   *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/cvodes/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the implementation file for the CVODEA adjoint          *
 * integrator.                                                     *
 *                                                                 *
 *******************************************************************/

/*=================================================================*/
/*BEGIN             Import Header Files                            */
/*=================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "cvodea.h"
#include "sundialsmath.h"

/*=================================================================*/
/*END               Import Header Files                            */
/*=================================================================*/

/*=================================================================*/
/*BEGIN             Macros                                         */
/*=================================================================*/

/* Macro: loop */

#define loop for(;;)

/*=================================================================*/
/*END               Macros                                         */
/*=================================================================*/

/*=================================================================*/
/*BEGIN             CVODEA Private Constants                       */
/*=================================================================*/

#define ZERO        RCONST(0.0)     /* real 0.0 */
#define ONE         RCONST(1.0)     /* real 1.0 */
#define TWO         RCONST(2.0)     /* real 2.0 */
#define FUZZ_FACTOR RCONST(1000000.0)  /* fuzz factor for CVadjGetY */

/*=================================================================*/
/*END               CVODEA Private Constants                       */
/*=================================================================*/

/*=================================================================*/
/*BEGIN             Error Messages                                 */
/*=================================================================*/

#define CVAM                "CVadjMalloc-- "
#define MSG_CVAM_NO_MEM     CVAM "cvode_mem=NULL illegal.\n\n"
#define MSG_CVAM_BAD_STEPS  CVAM "steps non-positive illegal.\n\n"
#define MSG_CVAM_MEM_FAIL   CVAM "a memory request failed.\n\n"

#define CVF                 "CVodeF-- "
#define MSG_CVODEF_MEM_FAIL CVF "a memory request failed.\n\n"

#define CVBM                "CVodeMallocB/CVodeReInitB-- "
#define MSG_CVBM_NO_MEM     CVBM "cvadj_mem=NULL illegal.\n\n"
#define MSG_CVBM_BAD_TB0    CVBM "tB0 out of range.\n\n"
#define MSG_CVBM_MEM_FAIL   CVBM "a memory request failed.\n\n"

#define CVBQM               "CVodeQuadMallocB-- "
#define MSG_CVBQM_NO_MEM    CVBQM "cvadj_mem=NULL illegal.\n\n"
#define MSG_BAD_ECONQB_1    CVBQM "errconQB=%d illegal.\n"
#define MSG_BAD_ECONQB_2    "The legal values are FULL=%d and PARTIAL=%d.\n\n"
#define MSG_BAD_ECONQB      MSG_BAD_ECONQB_1 MSG_BAD_ECONQB_2

#define CVB                 "CVodeB-- "
#define MSG_CVODEB_FWD      CVB "an error occured during the forward phase.\n\n"
#define MSG_CVODEB_BCK      CVB "an error occured during the backward phase.\n\n"

/*=================================================================*/
/*END               Error Messages                                 */
/*=================================================================*/

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
static void CVAhermitePrepare(CVadjMem ca_mem, DtpntMem *dt_mem, 
                              long int i);
static void CVAhermiteInterpolate(CVadjMem ca_mem, DtpntMem *dt_mem,
                                  long int i, realtype t, N_Vector y);

static void CVArhs(realtype t, N_Vector yB, 
                   N_Vector yBdot, void *passed_data);
static void CVAdenseJac(integertype nB, DenseMat JB, realtype t, 
                        N_Vector yB, N_Vector fyB, void *cvadj_mem,
                        N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
static void CVAbandJac(integertype nB, integertype mupperB, 
                       integertype mlowerB, BandMat JB, realtype t, 
                       N_Vector yB, N_Vector fyB, void *cvadj_mem, 
                       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
static int CVAspgmrPrecond(realtype t, N_Vector yB, 
                           N_Vector fyB, booleantype jokB, 
                           booleantype *jcurPtrB, realtype gammaB,
                           N_Vector ewtB, realtype hB, realtype uroundB,
                           long int *nfePtrB, void *cvadj_mem,
                           N_Vector vtemp1B, N_Vector vtemp2B,
                           N_Vector vtemp3B);
static int CVAspgmrPsolve(realtype t, N_Vector yB, 
                          N_Vector fyB, N_Vector vtempB, 
                          realtype gammaB, N_Vector ewtB,
                          realtype deltaB, long int *nfePtrB, 
                          N_Vector rB, int lrB, void *cvadj_mem, 
                          N_Vector zB);
static int CVAspgmrJtimes(N_Vector vB, N_Vector JvB, realtype t, 
                          N_Vector yB, N_Vector fyB, 
                          void *cvadj_mem, N_Vector workB);
static void CVArhsQ(realtype t, N_Vector yB, 
                    N_Vector qBdot, void *passed_data);
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
#define precond_B  (ca_mem->ca_precondB)
#define psolve_B   (ca_mem->ca_psolveB)
#define jtimes_B   (ca_mem->ca_jtimesB)
#define jac_data_B (ca_mem->ca_jac_dataB)
#define P_data_B   (ca_mem->ca_P_dataB)
#define ioptBalloc (ca_mem->ca_ioptBalloc)
#define roptBalloc (ca_mem->ca_roptBalloc)
#define fQ_B       (ca_mem->ca_fQB)
#define fQ_data_B  (ca_mem->ca_fQ_dataB)
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

#define machenv_B  (cvb_mem->cv_machenv)
#define optIn_B    (cvb_mem->cv_iopt)
#define iopt_B     (cvb_mem->cv_iopt)
#define ropt_B     (cvb_mem->cv_ropt)
#define errconQ_B  (cvb_mem->cv_errconQ)
#define machenvQ_B (cvb_mem->cv_machenvQ)

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
    fprintf(stdout, MSG_CVAM_NO_MEM);
    return (NULL);
  }
  if (steps <= 0) {
    fprintf(stdout, MSG_CVAM_BAD_STEPS);
    return (NULL);
  }

  /* Allocate memory block */
  ca_mem = (CVadjMem) malloc(sizeof(struct CVadjMemRec));
  if (ca_mem == NULL) {
    fprintf(stdout, MSG_CVAM_MEM_FAIL);
    return(NULL);
  }

  /* Attach CVODE memory for forward runs */
  cv_mem = (CVodeMem)cvode_mem;
  ca_mem->cv_mem = cv_mem;

  /* Initialize Check Points linked list */
  ca_mem->ck_mem = CVAckpntInit(cv_mem);
  if (ca_mem->ck_mem == NULL) {
    free(ca_mem);
    fprintf(stdout, MSG_CVAM_MEM_FAIL);
    return(NULL);
  }

  /* Allocate Data Points memory */
  ca_mem->dt_mem = CVAdataMalloc(cv_mem, steps);
  if (ca_mem->dt_mem == NULL) {
    CVAckpntDelete(&(ca_mem->ck_mem));
    free(ca_mem);
    fprintf(stdout, MSG_CVAM_MEM_FAIL);
    return(NULL);
  }

  /* Workspace memory */
  Y0 = N_VNew(machenv);
  if (Y0 == NULL) {
    CVAdataFree(ca_mem->dt_mem, nsteps);
    CVAckpntDelete(&(ca_mem->ck_mem));
    free(ca_mem);
    fprintf(stdout, MSG_CVAM_MEM_FAIL);
    return(NULL);
  }

  Y1 = N_VNew(machenv);
  if (Y1 == NULL) {
    N_VFree(Y0);
    CVAdataFree(ca_mem->dt_mem, nsteps);
    CVAckpntDelete(&(ca_mem->ck_mem));
    free(ca_mem);
    fprintf(stdout, MSG_CVAM_MEM_FAIL);
    return(NULL);
  }

  ytmp = N_VNew(machenv);
  if (ytmp == NULL) {
    N_VFree(Y0);
    N_VFree(Y1);
    CVAdataFree(ca_mem->dt_mem, nsteps);
    CVAckpntDelete(&(ca_mem->ck_mem));
    free(ca_mem);
    fprintf(stdout, MSG_CVAM_MEM_FAIL);
    return(NULL);
  }

  /* Other entries in ca_mem */
  uround   = cv_mem->cv_uround;
  nsteps   = steps;
  tinitial = tn; 

  /* Initialize nckpnts to ZERO */
  nckpnts = 0;

  return((void *)ca_mem);
} 

/*=================================================================*/
/*BEGIN             Wrappers for CVODEA                            */
/*=================================================================*/

/*------------------     CVodeF          --------------------------*/
/*
  This routine integrates to tout and returns solution into yout.
  In the same time, it stores check point data every 'steps' steps. 
  
  CVodeF can be called repeatedly by the user. The last tout
  will be used as the starting time for the backward integration.
  
  ncheckPtr points to the number of check points stored so far.
*/
/*-----------------------------------------------------------------*/

int CVodeF(void *cvadj_mem, realtype tout, N_Vector yout, realtype *t, 
           int itask, int *ncheckPtr)
{
  CVadjMem ca_mem;
  CVodeMem cv_mem;
  CkpntMem tmp;
  DtpntMem *dt_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cv_mem = ca_mem->cv_mem;
  dt_mem = ca_mem->dt_mem;

  /* On the first step, load dt_mem[0] */
  if ( nst == 0) {
    dt_mem[0]->t = ca_mem->ck_mem->ck_t0;
    N_VScale(ONE, ca_mem->ck_mem->ck_zn[0], dt_mem[0]->y);
    N_VScale(ONE, ca_mem->ck_mem->ck_zn[1], dt_mem[0]->yd);
  }

  /* Integrate to tout while loading check points */

  loop {

    /* Perform one step of the integration */

    flag = CVode(cv_mem, tout, yout, t, ONE_STEP);
    if (flag < 0) break;

    /* Test if a new check point is needed */

    if ( nst % nsteps == 0 ) {

      ca_mem->ck_mem->ck_t1 = *t;

      /* Create a new check point, load it, and append it to the list */
      tmp = CVAckpntNew(cv_mem);
      if (tmp == NULL) {
        flag = CVODEF_MEM_FAIL;
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
      dt_mem[nst%nsteps]->t = *t;
      N_VScale(ONE, yout, dt_mem[nst%nsteps]->y);
      CVodeDky(cv_mem, *t, 1, dt_mem[nst%nsteps]->yd);

    }

    /* Set t1 field of the current ckeck point structure
       for the case in which there will be no future
       check points */
    ca_mem->ck_mem->ck_t1 = *t;

    /* tfinal is now set to *t */
    tfinal = *t;

    /* In ONE_STEP mode break from loop */
    if (itask == ONE_STEP) break;

    /* In NORMAL mode and if we passed tout,
       evaluate yout at tout, set t=tout,
       then break from the loop */
    if ( (itask == NORMAL) && (*t >= tout) ) {
      *t = tout;
      CVodeDky(cv_mem, tout, 0, yout);
      break;
    }

  }

  if (flag == CVODEF_MEM_FAIL) 
    fprintf(stdout, MSG_CVODEF_MEM_FAIL);

  /* Get ncheck from ca_mem */ 
  *ncheckPtr = nckpnts;

  /* Data is available for the last interval */
  newData = TRUE;
  ckpntData = ca_mem->ck_mem;
  np = nst % nsteps + 1;

  return(flag);

}

/*------------------    CVodeMallocB     --------------------------*/
/*
  CVodeMallocB allocates memory for the backward run.            
  It is essentailly a call to CVodeMalloc but with some          
  particularizations for backward integration:                   
    - it takes as an argument a pointer to the adjoint memory    
      structure and sets the cvode memory for the backward run   
      as a field of this structure.
    - the routine that provides the ODE right hand side is of a  
      different type (it also gets the solution of the forward   
      integration at the current time).    
*/
/*-----------------------------------------------------------------*/

int CVodeMallocB(void *cvadj_mem, RhsFnB fB, 
                 realtype tB0, N_Vector yB0, int lmmB, int iterB, int itolB, 
                 realtype *reltolB, void *abstolB, void *f_dataB, 
                 FILE *errfpB, booleantype optInB, 
                 long int ioptB[], realtype roptB[], M_Env machEnvB)
{
  CVadjMem ca_mem;
  void *cvode_mem;
  int j;

  if (cvadj_mem == NULL) {
    fprintf(stdout, MSG_CVBM_NO_MEM);
    return(CVBM_NO_MEM);
  }

  ca_mem = (CVadjMem) cvadj_mem;

  if ( (tB0 < tinitial) || (tB0 > tfinal) ) {
    fprintf(stdout, MSG_CVBM_BAD_TB0);
    return(CVBM_BAD_TB0);
  }

  if (ioptB == NULL) {
    ioptB = (long int *)malloc(OPT_SIZE*sizeof(long int));
    if (ioptB == NULL) {
      fprintf(stdout, MSG_CVBM_MEM_FAIL);
      return(CVBM_MEM_FAIL);
    }
    ioptBalloc = TRUE;
  } else {
    ioptBalloc = FALSE;
  }
  if (roptB == NULL) {
    roptB = (realtype *)malloc(OPT_SIZE*sizeof(realtype));
    if (roptB == NULL) {
      fprintf(stdout, MSG_CVBM_MEM_FAIL);
      return(CVBM_MEM_FAIL);
    }
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

  cvode_mem = CVodeMalloc(CVArhs, tB0, yB0, lmmB, iterB, 
                          itolB, reltolB, abstolB, cvadj_mem,
                          errfpB, optInB, ioptB, roptB, machEnvB);

  if (cvode_mem == NULL) {
    fprintf(stdout, MSG_CVBM_MEM_FAIL);
    return(CVBM_MEM_FAIL);
  }

  ca_mem->cvb_mem = (CVodeMem) cvode_mem;

  f_B      = fB;
  f_data_B = f_dataB;

  return(SUCCESS);

}

/*------------------    CVodeReInitB     --------------------------*/
/*
  CVodeReInitB reinitializes the backward problem for a new 
  integration using the same forward solution and hence check 
  points. It sets a new final time and final conditions for the 
  backward problem.
*/
/*-----------------------------------------------------------------*/

int CVodeReInitB(void *cvadj_mem, RhsFnB fB, 
                 realtype tB0, N_Vector yB0, int lmmB, int iterB, int itolB, 
                 realtype *reltolB, void *abstolB, void *f_dataB, 
                 FILE *errfpB, booleantype optInB, 
                 long int ioptB[], realtype roptB[], M_Env machEnvB)
{
  CVadjMem ca_mem;
  CVodeMem cvb_mem;
  int flag, j;

  if (cvadj_mem == NULL) {
    fprintf(stdout, MSG_CVBM_NO_MEM);
    return(CVBM_NO_MEM);
  }

  ca_mem = (CVadjMem) cvadj_mem;

  if ( (tB0 < tinitial) || (tB0 > tfinal) ) {
    fprintf(stdout, MSG_CVBM_BAD_TB0);
    return(CVBM_BAD_TB0);
  }

  if (ioptB == NULL) {
    ioptB = (long int *)malloc(OPT_SIZE*sizeof(long int));
    if (ioptB == NULL) {
      fprintf(stdout, MSG_CVBM_MEM_FAIL);
      return(CVBM_MEM_FAIL);
    }
    ioptBalloc = TRUE;
  } else {
    ioptBalloc = FALSE;
  }
  if (roptB == NULL) {
    roptB = (realtype *)malloc(OPT_SIZE*sizeof(realtype));
    if (roptB == NULL) {
      fprintf(stdout, MSG_CVBM_MEM_FAIL);
      return(CVBM_MEM_FAIL);
    }
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

  cvb_mem = ca_mem->cvb_mem;

  flag = CVodeReInit(cvb_mem, CVArhs, tB0, yB0, lmmB, iterB,
                     itolB, reltolB, abstolB, cvadj_mem,
                     errfpB, optInB, ioptB, roptB, machEnvB);
  if (flag != SUCCESS)
    return(flag);

  f_B      = fB;
  f_data_B = f_dataB;

  return(SUCCESS);

}

/*------------------    CVodeQuadMallocB --------------------------*/
/*
  CVodeQuadMallocB initializes quadrature computation during
  hte backward integration phase
*/
/*-----------------------------------------------------------------*/

int CVodeQuadMallocB(void *cvadj_mem, QuadRhsFnB fQB,
                     int errconQB, realtype *reltolQB, void *abstolQB,
                     void *fQ_dataB, M_Env machEnvQB)
{
  CVadjMem ca_mem;
  CVodeMem cvb_mem;
  FILE *fp;
  int flag;

  if (cvadj_mem == NULL) {
    fprintf(stdout, MSG_CVBQM_NO_MEM);
    return(CVBM_NO_MEM);
  }

  ca_mem = (CVadjMem) cvadj_mem;
  cvb_mem = ca_mem->cvb_mem;
  fp = cvb_mem->cv_errfp;

  /* Check if errcon is legal */
  if ((errconQB!=FULL) && (errconQB!=PARTIAL)) {
    fprintf(fp, MSG_BAD_ECONQB,errconQB,FULL,PARTIAL);
    return (CVBM_ILL_INPUT);
  }

  flag = CVodeQuadMalloc(cvb_mem, CVArhsQ, errconQB, 
                         reltolQB, abstolQB, cvadj_mem, machEnvQB);

  if (flag != SUCCESS)
    return(flag);

  fQ_B = fQB;
  fQ_data_B = fQ_dataB;

  return(SUCCESS);
}

/*------------------    CVodeQuadReInitB --------------------------*/
/*
  CVodeQuadReInitB re-initializes quadrature computation during
  the backward integration phase
*/
/*-----------------------------------------------------------------*/

int CVodeQuadReInitB(void *cvadj_mem, QuadRhsFnB fQB,
                     int errconQB, realtype *reltolQB, void *abstolQB,
                     void *fQ_dataB, M_Env machEnvQB)
{
  CVadjMem ca_mem;
  CVodeMem cvb_mem;
  FILE *fp;
  int flag;

  if (cvadj_mem == NULL) {
    fprintf(stdout, MSG_CVBQM_NO_MEM);
    return(CVBM_NO_MEM);
  }

  ca_mem = (CVadjMem) cvadj_mem;
  cvb_mem = ca_mem->cvb_mem;
  fp = cvb_mem->cv_errfp;

  /* Check if errcon is legal */
  if ((errconQB!=FULL) && (errconQB!=PARTIAL)) {
    fprintf(fp, MSG_BAD_ECONQB,errconQB,FULL,PARTIAL);
    return (CVBM_ILL_INPUT);
  }

  flag = CVodeQuadReInit(cvb_mem, CVArhsQ, errconQB, 
                         reltolQB, abstolQB, 
                         cvadj_mem, machEnvQB);

  if (flag != SUCCESS)
    return(flag);

  fQ_B = fQB;
  fQ_data_B = fQ_dataB;

  return(SUCCESS);
}

/*------------------    CVDenseB         --------------------------*/
/*
  CVDenseB links the main CVODE integrator with the CVDENSE
  linear solver for the backward integration.
*/
/*-----------------------------------------------------------------*/

int CVDenseB(void *cvadj_mem, integertype nB, CVDenseJacFnB djacB, void *jac_dataB)
{
  CVadjMem ca_mem;
  CVodeMem cvb_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvb_mem = ca_mem->cvb_mem;
  
  if (djacB == NULL) {
    flag = CVDense(cvb_mem, nB, NULL, NULL);
  } else {
    flag = CVDense(cvb_mem, nB, CVAdenseJac, cvadj_mem);
    djac_B     = djacB;
    jac_data_B = jac_dataB;
  }

  return(flag);
}

/*------------------     CVBandB         --------------------------*/
/*
  CVBandB links the main CVODE integrator with the CVBAND
  linear solver for the backward integration.
*/
/*-----------------------------------------------------------------*/

int CVBandB(void *cvadj_mem, integertype nB, 
            integertype mupperB, integertype mlowerB,
            CVBandJacFnB bjacB, void *jac_dataB)
{
  CVadjMem ca_mem;
  CVodeMem cvb_mem;
  int flag;

  ca_mem = (CVadjMem) cvadj_mem;
  cvb_mem = ca_mem->cvb_mem;
  
  if (bjacB == NULL) {
    flag = CVBand(cvb_mem, nB, mupperB, mlowerB, NULL, NULL);
  } else {
    flag = CVBand(cvb_mem, nB, mupperB, mlowerB, CVAbandJac, cvadj_mem);
    bjac_B     = bjacB;
    jac_data_B = jac_dataB;
  }

  return(flag);
}

/*------------------   CVBandPreAllocB   --------------------------*/
/*
  CVBandPreAllocB interfaces to the CVBandPre preconditioner 
  for the backward integration. The pointer to the CVBandPreData
  structure returned by this routine can then be used together
  with the functions CVBandPrecond and CVBandPSolve in a call 
  to CVSpgmrB. 
*/
/*-----------------------------------------------------------------*/

CVBandPreData CVBandPreAllocB(void *cvadj_mem, integertype nB, 
                              integertype muB, integertype mlB)
{
  CVadjMem ca_mem;
  CVodeMem cvb_mem;
  CVBandPreData bpdata;

  ca_mem = (CVadjMem) cvadj_mem;
  cvb_mem = ca_mem->cvb_mem;

  bpdata = CVBandPreAlloc(nB, CVArhs, cvadj_mem, muB, mlB, cvb_mem);
  return(bpdata);
}

/*------------------    CVBandPrecondB   --------------------------*/
/*-----------------------------------------------------------------*/

int CVBandPrecondB(realtype t, N_Vector y, 
                   N_Vector yB, N_Vector fyB, booleantype jokB, 
                   booleantype *jcurPtrB, realtype gammaB,
                   N_Vector ewtB, realtype hB, realtype uroundB,
                   long int *nfePtrB, void *P_dataB,
                   N_Vector vtemp1B, N_Vector vtemp2B,
                   N_Vector vtemp3B)
{
  int flag;
  flag = CVBandPrecond(t, yB, fyB, jokB, jcurPtrB, gammaB,
                       ewtB, hB, uroundB, nfePtrB, P_dataB,
                       vtemp1B, vtemp2B, vtemp3B);
  return(flag);
}

/*------------------     CVBandPSolveB   --------------------------*/
/*-----------------------------------------------------------------*/

int CVBandPSolveB(realtype t, N_Vector y, 
                  N_Vector yB, N_Vector fyB, 
                  N_Vector vtempB,  realtype gammaB, 
                  N_Vector ewtB, realtype deltaB, 
                  long int *nfePtrB, N_Vector rB, 
                  int lrB, void *P_dataB, N_Vector zB)
{
  int flag;
  flag = CVBandPSolve(t, yB, fyB, vtempB, gammaB, ewtB, 
                      deltaB, nfePtrB, rB, lrB, P_dataB, zB);
  return(flag);
}

/*------------------   CVSpgmrB          --------------------------*/
/*
  CVSpgmrB links the main CVODE integrator with the CVSPGMR
  linear solver for the backward integration.
*/
/*-----------------------------------------------------------------*/

int CVSpgmrB(void *cvadj_mem, int pretypeB, int gstypeB, 
             int maxlB, realtype deltB, CVSpgmrPrecondFnB precondB, 
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

/*------------------     CVodeB          --------------------------*/
/*
  This routine performs the backward integration from tB0 
  to tinitial through a sequence of forward-backward runs in
  between consecutive check points. It returns the values of
  the adjoint variables and any existing quadrature variables
  at tinitial.
*/
/*-----------------------------------------------------------------*/

int CVodeB(void *cvadj_mem, N_Vector yB)
{
  CVadjMem ca_mem;
  CkpntMem ck_mem;
  CVodeMem cvb_mem;
  int flag;
  realtype tB0, t;
  
  ca_mem  = (CVadjMem) cvadj_mem;
  ck_mem = ca_mem->ck_mem;
  cvb_mem = ca_mem->cvb_mem;

  /* First decide which check points tB0 falls in between */
  tB0 = cvb_mem->cv_tn;
  while ( tB0 <= t0_ ) ck_mem = next_;
  
  do {

    /* If data in dt_mem is not available for the 
       current check point, compute it */
    if (ck_mem != ckpntData) {
      flag = CVAdataStore(ca_mem, ck_mem);
      if (flag < 0) {
        fprintf(stdout, MSG_CVODEB_FWD);
        return(flag);
      }
    }

    /* Propagate backward integration to next check point */
    ropt_B[TSTOP] = t0_;
    flag = CVode(cvb_mem, t0_, yB, &t, NORMAL);
    if (flag < 0) {
      fprintf(stdout, MSG_CVODEB_BCK);
      return(flag);
    }

    /* Move check point in linked list to next one */
    ck_mem = next_;

  } while (ck_mem != NULL);

  t_for_quad = t;

  return(flag);

}

/*------------------  CVodeQuadExtractB  --------------------------*/
/*
  This routine extracts values of the quadrature variables
*/
/*-----------------------------------------------------------------*/

int CVodeQuadExtractB(void *cvadj_mem, N_Vector qB)
{
  CVadjMem ca_mem;
  CVodeMem cvb_mem;
  int flag;
  
  ca_mem  = (CVadjMem) cvadj_mem;
  cvb_mem = ca_mem->cvb_mem;
  
  flag = CVodeQuadExtract(cvb_mem, t_for_quad, qB);

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
  booleantype to_left, to_right;
  realtype troundoff;

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
      printf("\n TROUBLE IN GETY\n ");
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

/*------------------  CVAdjCheckPointsList ------------------------*/
/*
  This routine lists the linked list of check point structures.
  For debugging....
*/
/*-----------------------------------------------------------------*/

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

/*------------------   CVAdjDataExtract  --------------------------*/
/*
  This routine returns the solution stored in the data structure
  at the 'which' data point.
  For debugging....
*/
/*-----------------------------------------------------------------*/

void CVadjDataExtract(void *cvadj_mem, long int which, 
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
  zn_[0] = N_VNew(machenv);
  zn_[1] = N_VNew(machenv);

  /* Load ckdata from cv_mem */
  N_VScale(ONE, zn[0], zn_[0]);
  t0_    = tn;
  q_     = 1;
  /* Compute zn_[1] by calling the user f routine */
  f(t0_, zn_[0], zn_[1], f_data);
  
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
  int j, jj;

  /* Allocate space for ckdata */
  ck_mem = (CkpntMem) malloc(sizeof(struct CkpntMemRec));
  if (ck_mem == NULL) return(NULL);

  for (j=0; j<=q; j++) {
    zn_[j] = N_VNew(machenv);
    if(zn_[j] == NULL) {
      for(jj=0; jj<j; jj++) N_VFree(zn_[jj]);
      return(NULL);
    }
  }

  /* Load check point data from cv_mem */
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
    /* free tmp */
    for (j=0;j<=tmp->ck_q;j++) N_VFree(tmp->ck_zn[j]);
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
    dt_mem[i]->y  = N_VNew(machenv);
    dt_mem[i]->yd = N_VNew(machenv);
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
    N_VFree(dt_mem[i]->y);
    N_VFree(dt_mem[i]->yd);
    free(dt_mem[i]);
  }

}

/*------------------    CVAdataStore     --------------------------*/
/*
  This routine integrates the forward model starting at the check
  point ck_mem and stores y and yprime at all intermediate 
  steps. It returns the error flag from CVode.
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

  /* Set first structure in dt_mem[0] */
  dt_mem[0]->t = t0_;
  N_VScale(ONE, zn_[0], dt_mem[0]->y);
  N_VScale(ONE, zn_[1], dt_mem[0]->yd);

  /* Run CVode to set following structures in dt_mem[i] */
  i = 1;
  do {
    flag = CVode(cv_mem, t1_, dt_mem[i]->y, &t, ONE_STEP);
    if (flag < 0) return(flag);
    dt_mem[i]->t = t;
    CVodeDky(cv_mem, t, 1, dt_mem[i]->yd);
    i++;
  } while (t<t1_);

  /* New data is now available */
  ckpntData = ck_mem;
  newData = TRUE;
  np  = i;

  return(flag);

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

  if (next_ == NULL) {

    /* In this case, we just call the reinitialization routine,
       but make sure we use the same initial stepsize as on 
       the first run. */

    if (iopt == NULL)
      iopt = (long int *)malloc(OPT_SIZE*sizeof(long int));
    if (ropt == NULL)
      ropt = (realtype *)malloc(OPT_SIZE*sizeof(realtype));
    if (optIn == FALSE)
      for (j=0; j<OPT_SIZE; j++) {
        iopt[j] = 0;
        ropt[j] = 0.0;
      }
    optIn = TRUE;
    ropt[H0] = h0u;
    CVodeReInit(cv_mem, f, t0_, zn_[0],
                lmm, iter, itol, reltol, abstol,
                f_data, errfp, optIn, iopt, ropt, machenv);
    
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
    
    /* Copy the arrays from check point data structure */
    for (j=0; j<=q; j++)         N_VScale(ONE, zn_[j], zn[j]);
    for (j=0; j<=L_MAX; j++)     tau[j] = tau_[j];
    for (j=0; j<=NUM_TESTS; j++) tq[j] = tq_[j];
    for (j=0; j<=q; j++)         l[j] = l_[j];
    
    /* Force a call to setup */
    forceSetup = TRUE;

  }

  return(0);
}

/*------------------   CVAhermitePrepare --------------------------*/
/*
  This routine computes quantities required by the Hermite
  interpolation that are independent of the interpolation point.
*/
/*-----------------------------------------------------------------*/

static void CVAhermitePrepare(CVadjMem ca_mem, DtpntMem *dt_mem, 
                              long int i)
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
  This routine interfaces to the RhsFnB routine provided by
  the user.
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
  if (flag != GETY_OK) {
    printf("\n\nBad t in interpolation\n\n");
    exit(1);
  }

  /* Call user's adjoint RHS routine */
  f_B(t, ytmp, yB, yBdot, f_data_B);

}

/*------------------       CVArhsQ       --------------------------*/
/*
  This routine interfaces to the QuadRhsFnB routine provided by
  the user.
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
  if (flag != GETY_OK) {
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
*/
/*-----------------------------------------------------------------*/

static void CVAdenseJac(integertype nB, DenseMat JB, realtype t, 
                        N_Vector yB, N_Vector fyB, void *cvadj_mem,
                        N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
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
  djac_B(nB, JB, t, ytmp, yB, fyB, jac_data_B, 
         tmp1B, tmp2B, tmp3B);

}

/*------------------    CVAbandJac       --------------------------*/
/*
  This routine interfaces to the CVBandJacFnB routine provided 
  by the user.
*/
/*-----------------------------------------------------------------*/

static void CVAbandJac(integertype nB, integertype mupperB, 
                       integertype mlowerB, BandMat JB, realtype t, 
                       N_Vector yB, N_Vector fyB, void *cvadj_mem, 
                       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
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
  bjac_B(nB, mupperB, mlowerB, JB, t, ytmp, yB, fyB, jac_data_B,
         tmp1B, tmp2B, tmp3B);

}

/*------------------   CVAspgmrPrecond   --------------------------*/
/*
  This routine interfaces to the CVSpgmrPrecondFnB routine 
  provided by the user.
*/
/*-----------------------------------------------------------------*/

static int CVAspgmrPrecond(realtype t, N_Vector yB, 
                           N_Vector fyB, booleantype jokB, 
                           booleantype *jcurPtrB, realtype gammaB,
                           N_Vector ewtB, realtype hB, realtype uroundB,
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
  flag = precond_B(t, ytmp, yB, fyB, jokB, jcurPtrB, gammaB, ewtB, hB, 
                   uroundB, nfePtrB, P_data_B, vtemp1B, vtemp2B, vtemp3B);

  return(flag);
}

/*------------------   CVAspgmrPSolve    --------------------------*/
/*
  This routine interfaces to the CVSpgmrPsolveFnB routine 
  provided by the user.
*/
/*-----------------------------------------------------------------*/

static int CVAspgmrPsolve(realtype t, N_Vector yB, 
                          N_Vector fyB, N_Vector vtempB, 
                          realtype gammaB, N_Vector ewtB,
                          realtype deltaB, long int *nfePtrB, 
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
  flag = psolve_B(t, ytmp, yB, fyB, vtempB, gammaB, ewtB, deltaB, 
                  nfePtrB, rB, lrB, P_data_B, zB);

  return(flag);
}

/*------------------   CVAspgmrJtimes    --------------------------*/
/*
  This routine interfaces to the CVSpgmrJtimesFnB routine 
  provided by the user.
*/
/*-----------------------------------------------------------------*/

static int CVAspgmrJtimes(N_Vector vB, N_Vector JvB, realtype t, 
                          N_Vector yB, N_Vector fyB, 
                          void *cvadj_mem, N_Vector workB)
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
  flag = jtimes_B(vB, JvB, t, ytmp, yB, fyB, jac_data_B, workB);

  return(flag);
}

/*=================================================================*/
/*END               Wrappers for adjoint system                    */
/*=================================================================*/

/*=================================================================*/
/*END               Private Functions Implementation               */
/*=================================================================*/
