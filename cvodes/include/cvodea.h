/*******************************************************************
 *                                                                 *
 * File          : cvodea.h                                        *
 * Programmers   : Radu Serban @ LLNL                              *
 * Version of    : 27 June 2002                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/cvodes/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the interface file for the CVODEA adjoint integrator.   *
 *                                                                 *
 * Function types:                                                 *
 *    RhsFnB                                                       *
 *    CVDenseJacFnB                                                *
 *    CVBandJacFnB                                                 *
 *    CVSpgmrPrecondFnB                                            *
 *    CVSpgmrPsolveB                                               *
 *    CVSpgmrJtimesB                                               *
 * Exported functions prototypes:                                  *
 *    CVadjMalloc                                                  *
 *    CVodeF                                                       *
 *    CVodeMallocB                                                 *
 *    CVDenseB                                                     *
 *    CVBandB                                                      *
 *    CVBandPreAllocB                                              *
 *    CVSpgmrB                                                     *
 *    CVodeB                                                       *
 *    CVadjFree                                                    *
 *    CVadjGetY                                                    *
 *    CVadjCheckPointsList                                         *
 *    CVadjDataExtract                                             *
 * Type definitions:                                               *
 *    struct CkpntMemRec, CkpntMem                                 *
 *    struct DtpntMemRec, DtpntMem                                 *
 *    struct CVadjMemRec, CVadjMem                                 *
 *                                                                 *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvodea_h
#define _cvodea_h

#include <stdio.h>
#include "nvector.h"
#include "cvodes.h"
#include "cvsdense.h"
#include "cvsband.h"
#include "cvsspgmr.h"
#include "cvsbandpre.h"

/******************************************************************
 *                                                                *
 * Type : RhsFnB                                                  *
 *----------------------------------------------------------------*
 * The fB function which defines the right hand side of the       *
 * ODE systems to be integrated backwards must have type RhsFnB.  *
 *                                                                *
 ******************************************************************/

typedef void (*RhsFnB)(integertype NB, realtype t, N_Vector y, 
                       N_Vector yB, N_Vector yBdot, void *f_dataB);
/******************************************************************
 *                                                                *
 * Type : CVDenseJacFnB                                           *
 *----------------------------------------------------------------*
 * A dense Jacobian approximation function djacB for the backward *
 * integration must have the prototype given below.               *
 *                                                                *
 ******************************************************************/

typedef void (*CVDenseJacFnB)(integertype NB, DenseMat JB, RhsFnB fB, 
                              void *f_dataB, realtype t, N_Vector y, 
                              N_Vector yB, N_Vector fyB, N_Vector ewtB,
                              realtype hB, realtype uroundB, 
                              void *jac_dataB,
                              long int *nfePtrB, N_Vector vtemp1B,
                              N_Vector vtemp2B, N_Vector vtemp3B);

/******************************************************************
 *                                                                *
 * Type : CVBandJacFnB                                            *
 *----------------------------------------------------------------*
 * A band Jacobian approximation function bjacB for the backward  *
 * integration must have the prototype given below.               *
 *                                                                *
 ******************************************************************/

typedef void (*CVBandJacFnB)(integertype NB, integertype mupperB, 
                             integertype mlowerB, BandMat JB, RhsFnB fB, 
                             void *f_dataB, realtype t, N_Vector y, 
                             N_Vector yB, N_Vector fyB, N_Vector ewtB, 
                             realtype hB, realtype uroundB, 
                             void *jac_dataB, 
                             long int *nfePtrB, N_Vector vtemp1B, 
                             N_Vector vtemp2B, N_Vector vtemp3B);

/******************************************************************
 *                                                                *
 * Type : CVSpgmrPrecondFnB                                       *
 *----------------------------------------------------------------*
 * A preconditioner setup function precondB for the backward      *
 * integration must have the prototype given below.               *
 *                                                                *
 ******************************************************************/

typedef int (*CVSpgmrPrecondFnB)(integertype NB, realtype t, N_Vector y, 
                                 N_Vector yB, N_Vector fyB, booleantype jokB, 
                                 booleantype *jcurPtrB, realtype gammaB,
                                 N_Vector ewtB, realtype hB, realtype uroundB,
                                 long int *nfePtrB, void *P_dataB,
                                 N_Vector vtemp1B, N_Vector vtemp2B,
                                 N_Vector vtemp3B);

/******************************************************************
 *                                                                *
 * Type : CVSpgmrPSolveFnB                                        *
 *----------------------------------------------------------------*
 * A preconditioner solve function psolveB for the backward       *
 * integration must have the prototype given below.               *
 *                                                                *
 ******************************************************************/

typedef int (*CVSpgmrPSolveFnB)(integertype NB, realtype t, N_Vector y, 
                                N_Vector yB, N_Vector fyB, 
                                N_Vector vtempB,  realtype gammaB, 
                                N_Vector ewtB, realtype deltaB, 
                                long int *nfePtrB, N_Vector rB, 
                                int lrB, void *P_dataB, N_Vector zB);

/******************************************************************
 *                                                                *
 * Type : CVSpgmrJtimesFnB                                        *
 *----------------------------------------------------------------*
 * A Jacobian times vector function jtimesB for the backward      *
 * integration must have the prototype given below.               *
 *                                                                *
 ******************************************************************/

typedef int (*CVSpgmrJtimesFnB)(integertype NB, N_Vector vB, N_Vector JvB, 
                                RhsFnB fB, void *f_dataB, realtype t, 
                                N_Vector y, N_Vector yB, N_Vector fyB,
                                realtype vnrmB, N_Vector ewtB, realtype hB, 
                                realtype uroundB, void *jac_dataB, 
                                long int *nfePtrB, N_Vector workB);

/******************************************************************
 *                                                                *
 * Function : CVadjMalloc                                         *
 *----------------------------------------------------------------*
 * CVadjMalloc space for the global CVODEA memory structure.      *
 *                                                                *
 ******************************************************************/

void *CVadjMalloc(void *cvode_mem, long int steps);

/******************************************************************
 *                                                                *
 * Function : CVodeF                                              *
 *----------------------------------------------------------------*
 * CVodeF integrates towards tout and returns solution into yout. *
 * In the same time, it stores check point data every 'steps'.    *
 *                                                                *
 * CVodeF can be called repeatedly by the user. The last tout     *
 * will be used as the starting time for the backward integration.*
 *                                                                *
 * ncheckPtr points to the number of check points stored so far.  *
 *                                                                *
 ******************************************************************/

int CVodeF(void *cvadj_mem, realtype tout, N_Vector yout, realtype *t,
           int itask, int *ncheckPtr);

/* CVodeF return values */
#define CVODEF_MEM_FAIL -10
/* or any CVode return value */ 

/******************************************************************
 *                                                                *
 * Function : CVodeMallocB                                        *
 *----------------------------------------------------------------*
 * CVodeMallocB allocates memory for the backward run.            *
 * It is essentailly a call to CVodeMalloc but with some          *
 * particularizations for backward integration.                   *
 *                                                                *
 ******************************************************************/

int CVodeMallocB(void *cvadj_mem, integertype NB, RhsFnB fB, 
                 N_Vector yB0, int lmmB, int iterB, int itolB, 
                 realtype *reltolB, void *abstolB, void *f_dataB, 
                 FILE *errfpB, booleantype optInB, 
                 long int ioptB[], realtype roptB[], M_Env machEnv);

/* CVodeMallocB return values */
/* SUCCESS=0, defined under CVode return values */
#define CVBM_NO_MEM   -1
#define CVBM_MEM_FAIL -2

/******************************************************************
 *                                                                *
 * Function : CVDenseB                                            *
 *----------------------------------------------------------------*
 * CVDenseB links the main CVODE integrator with the CVDENSE      *
 * linear solver for the backward integration.                    *
 *                                                                *
 ******************************************************************/

int CVDenseB(void *cvadj_mem, CVDenseJacFnB djacB, void *jac_dataB);

/******************************************************************
 *                                                                *
 * Function : CVBandB                                             *
 *----------------------------------------------------------------*
 * CVBandB links the main CVODE integrator with the CVBAND        *
 * linear solver for the backward integration.                    *
 *                                                                *
 ******************************************************************/

int CVBandB(void *cvadj_mem, integertype mupperB, integertype mlowerB,
            CVBandJacFnB bjacB, void *jac_dataB);

/******************************************************************
 *                                                                *
 * Function : CVBandPreAllocB, CVBandPrecondB, CVBandPSolveB      *
 *----------------------------------------------------------------*
 * CVBandPreAllocB interfaces to the CVBandPre preconditioner for *
 * the backward integration. The pointer to the CVBandPreData     *
 * structure returned by this routine can then be used together   *
 * with the functions CVBandPrecondB and CVBandPSolveB in a call  *
 * to CVSpgmrB.                                                   *
 *                                                                *
 * CVBandPrecondB is of type CVSpgmrPrecondFnB to agree with the  *
 * expected argument in CVSpgmrB.                                 *
 *                                                                *
 * CVBandPSolveB is of type CVSpgmrPSolveFnB to agree with the    *
 * expected argument in CVSpgmrB.                                 *
 *                                                                *
 ******************************************************************/

CVBandPreData CVBandPreAllocB(void *cvadj_mem, integertype NB, 
                              integertype muB, integertype mlB);
int CVBandPrecondB(integertype NB, realtype t, N_Vector y, 
                   N_Vector yB, N_Vector fyB, booleantype jokB, 
                   booleantype *jcurPtrB, realtype gammaB,
                   N_Vector ewtB, realtype hB, realtype uroundB,
                   long int *nfePtrB, void *P_dataB,
                   N_Vector vtemp1B, N_Vector vtemp2B,
                   N_Vector vtemp3B);
int CVBandPSolveB(integertype NB, realtype t, N_Vector y, 
                  N_Vector yB, N_Vector fyB, 
                  N_Vector vtempB,  realtype gammaB, 
                  N_Vector ewtB, realtype deltaB, 
                  long int *nfePtrB, N_Vector rB, 
                  int lrB, void *P_dataB, N_Vector zB);

/******************************************************************
 *                                                                *
 * Function : CVSpgmrB                                            *
 *----------------------------------------------------------------*
 * CVSpgmrB links the main CVODE integrator with the CVSPGMR      *
 * linear solver for the backward integration.                    *
 *                                                                *
 ******************************************************************/

int CVSpgmrB(void *cvadj_mem, int pretypeB, int gstypeB, 
             int maxlB, realtype deltB, CVSpgmrPrecondFnB precondB, 
             CVSpgmrPSolveFnB psolveB, void *P_dataB,
             CVSpgmrJtimesFnB jtimesB, void *jac_dataB);

/******************************************************************
 *                                                                *
 * Function : CVodeB                                              *
 *----------------------------------------------------------------*
 * CVodeB performs the backward integration from tfinal to        *
 * tinitial through a sequence of forward-backward runs in        *
 * between consecutive check points. It returns the values of     *
 * the adjoint variables and any existing quadrature variables    *
 * at tinitial.                                                   *
 *                                                                *
 ******************************************************************/

int CVodeB(void *cvadj_mem, N_Vector yB);

/******************************************************************
 *                                                                *
 * Function : CVadjFree                                           *
 *----------------------------------------------------------------*
 * CVadjFree frees the memory allocated by CVadjMalloc.           *
 *                                                                *
 ******************************************************************/

void CVadjFree(void *cvadj_mem);

/******************************************************************
 *                                                                *
 * Function : CVadjGetY                                           *
 *----------------------------------------------------------------*
 * This routine uses cubic piece-wise Hermite interpolation for   *
 * the forward solution vector.                                   *
 *                                                                *
 ******************************************************************/

int   CVadjGetY(void *cvadj_mem, realtype t, N_Vector y);

/* CVadjGetY return values */

enum { GETY_OK=0, GETY_BADT=-1 };

/******************************************************************
 *                                                                *
 * Function : CVadjCheckPointsList                                *
 *----------------------------------------------------------------*
 *                                                                *
 ******************************************************************/

void  CVadjCheckPointsList(void *cvadj_mem);

/******************************************************************
 *                                                                *
 * Function : CVadjDataExtract                                    *
 *----------------------------------------------------------------*
 *                                                                *
 ******************************************************************/

void CVadjDataExtract(void *cvadj_mem, long int which, 
                      realtype *t, N_Vector yout, N_Vector ydout);


/******************************************************************
 *                                                                *
 * Types : struct CkpntMemRec, CkpntMem                           *
 *----------------------------------------------------------------*
 * The type CkpntMem is type pointer to struct CkpntMemRec.       *
 * This structure contains fields to store all information at a   *
 * check point that is needed to 'hot' start cvodes.              *
 *                                                                *
 ******************************************************************/

typedef struct CkpntMemRec {

  /* Integration limits */
  realtype     ck_t0;
  realtype     ck_t1;
   
  /* Nordsieck History Array */
  N_Vector ck_zn[L_MAX];
  
  /* Step data */
  int          ck_nst;
  int          ck_q;
  int          ck_qprime;
  int          ck_qwait;
  int          ck_L;
  realtype     ck_gammap;
  realtype     ck_h;
  realtype     ck_hprime;
  realtype     ck_hscale;
  realtype     ck_eta;
  realtype     ck_etamax;
  realtype     ck_tau[L_MAX+1];
  realtype     ck_tq[NUM_TESTS+1];
  realtype     ck_l[L_MAX];
  
  /* Saved values */
  realtype     ck_saved_tq5;

  /* Pointer to next structure in list */
  struct CkpntMemRec *ck_next;

} *CkpntMem;

/******************************************************************
 *                                                                *
 * Types : struct DtpntMemRec, DtpntMem                           *
 *----------------------------------------------------------------*
 * The type DtpntMem is type pointer to struct DtpntMemRec.       *
 * This structure contains fields to store all information at a   *
 * data point that is needed to interpolate solution of forward   *
 * simulations.                                                   *
 *                                                                *
 ******************************************************************/

typedef struct DtpntMemRec {
  
  /* time */
  realtype t;
  
  /* solution */
  N_Vector y;

  /* solution derivative */
  N_Vector yd;

} *DtpntMem;


/******************************************************************
 *                                                                *
 * Types : struct CVadjMemRec, CVadjMem                           *
 *----------------------------------------------------------------*
 * The type CVadjMem is type pointer to struct CVadjMemRec.       *
 * This structure contins fields to store all information         *
 * necessary for adjoint sensitivity analysis.                    *
 *                                                                *
 ******************************************************************/

typedef struct CVadjMemRec {

  /* CVODE memory for forward runs */
  struct CVodeMemRec *cv_mem;
  
  /* CVODE memory for backward run */
  struct CVodeMemRec *cvb_mem;
  
  /* Storage for check point information */
  struct CkpntMemRec *ck_mem;
  
  /* Storage for data from forward runs */
  struct DtpntMemRec **dt_mem;

  /* Does CVodeMalocB allocate memory for ioptB and roptB? */
  booleantype ca_ioptBalloc;
  booleantype ca_roptBalloc;

  /* Right hand side function (fB) for backward run */
  RhsFnB ca_fB;

  /* Dense Jacobian function (djacB) for backward run */
  CVDenseJacFnB ca_djacB;

  /* Banded Jacobian function (bjacB) for backward run */
  CVBandJacFnB ca_bjacB;

  /* Preconditioner routines (precondB and psolveB) for backward run */
  CVSpgmrPrecondFnB ca_precondB;
  CVSpgmrPSolveFnB  ca_psolveB;

  /* Jac times vec routine (jtimesB) for backward run */
  CVSpgmrJtimesFnB ca_jtimesB;

  /* User f_dataB */
  void *ca_f_dataB;
  
  /* User jac_dataB */
  void *ca_jac_dataB;

  /* User P_dataB */
  void *ca_P_dataB;
  
  /* Unit roundoff */
  realtype ca_uround;
  
  /* Integration interval */
  realtype ca_tinitial, ca_tfinal;
  
  /* Number of check points */
  int ca_nckpnts;
  
  /* Number of steps between 2 check points */
  long int ca_nsteps;

  /* Flag to indicate that data in dt_mem is new */
  booleantype ca_newData;

  /* Actual number of data points saved in current dt_mem */
  /* Commonly, np = nsteps+1                              */
  long int ca_np;
  
  /* Temporary space used by the Hermite interpolation */
  realtype ca_delta;
  N_Vector ca_Y0, ca_Y1;
  N_Vector ca_ytmp;

} *CVadjMem;

#endif

#ifdef __cplusplus
}
#endif
