/*******************************************************************
 *                                                                 *
 * File          : cvodea.h                                        *
 * Programmers   : Radu Serban @ LLNL                              *
 * Version of    : 07 January 2004                                 *
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
 *    CVSpgmrPrecSetupFnB                                          *
 *    CVSpgmrPrecSolveB                                            *
 *    CVSpgmrJacTimesVecFnB                                        *
 * Exported functions prototypes:                                  *
 *    CVadjMalloc                                                  *
 *    CVodeF                                                       *
 *    CVodeCreateB                                                 *
 *    CVodeMallocB                                                 *
 *    CVDenseB                                                     *
 *    CVBandB                                                      *
 *    CVBandPrecAllocB                                             *
 *    CVBandPrecFreeB                                              *
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
#include "cvdense.h"
#include "cvband.h"
#include "cvspgmr.h"
#include "cvbandpre.h"

/******************************************************************
 *                                                                *
 * Type : RhsFnB                                                  *
 *----------------------------------------------------------------*
 * The fB function which defines the right hand side of the       *
 * ODE systems to be integrated backwards must have type RhsFnB.  *
 *                                                                *
 ******************************************************************/

typedef void (*RhsFnB)(realtype t, N_Vector y, 
                       N_Vector yB, N_Vector yBdot, 
                       void *f_dataB);

/******************************************************************
 *                                                                *
 * Type : QuadRhsFnB                                              *
 *----------------------------------------------------------------*
 * The fQB function which defines the quadratures to be integrated*
 * backwards must have type QuadRhsFnB.                           *
 *                                                                *
 ******************************************************************/

typedef void (*QuadRhsFnB)(realtype t, N_Vector y, 
                           N_Vector yB, N_Vector qBdot, 
                           void *fQ_dataB);

/******************************************************************
 *                                                                *
 * Type : CVDenseJacFnB                                           *
 *----------------------------------------------------------------*
 * A dense Jacobian approximation function djacB for the backward *
 * integration must have the prototype given below.               *
 *                                                                *
 ******************************************************************/

typedef void (*CVDenseJacFnB)(long int nB, DenseMat JB, realtype t, 
                              N_Vector y, N_Vector yB, N_Vector fyB,
                              void *jac_dataB, N_Vector tmp1B,
                              N_Vector tmp2B, N_Vector tmp3B);

/******************************************************************
 *                                                                *
 * Type : CVBandJacFnB                                            *
 *----------------------------------------------------------------*
 * A band Jacobian approximation function bjacB for the backward  *
 * integration must have the prototype given below.               *
 *                                                                *
 ******************************************************************/

typedef void (*CVBandJacFnB)(long int nB, long int mupperB, 
                             long int mlowerB, BandMat JB,
                             realtype t, N_Vector y, 
                             N_Vector yB, N_Vector fyB,
                             void *jac_dataB, N_Vector tmp1B, 
                             N_Vector tmp2B, N_Vector tmp3B);

/******************************************************************
 *                                                                *
 * Type : CVSpgmrPrecSetupFnB                                     *
 *----------------------------------------------------------------*
 * A preconditioner setup function precondB for the backward      *
 * integration must have the prototype given below.               *
 *                                                                *
 ******************************************************************/

typedef int (*CVSpgmrPrecSetupFnB)(realtype t, N_Vector y, 
                                   N_Vector yB, N_Vector fyB, booleantype jokB, 
                                   booleantype *jcurPtrB, realtype gammaB,
                                   void *P_dataB,
                                   N_Vector tmp1B, N_Vector tmp2B,
                                   N_Vector tmp3B);

/******************************************************************
 *                                                                *
 * Type : CVSpgmrPrecSolveFnB                                     *
 *----------------------------------------------------------------*
 * A preconditioner solve function psolveB for the backward       *
 * integration must have the prototype given below.               *
 *                                                                *
 ******************************************************************/

typedef int (*CVSpgmrPrecSolveFnB)(realtype t, N_Vector y, 
                                   N_Vector yB, N_Vector fyB, 
                                   N_Vector rB, N_Vector zB,
                                   realtype gammaB, realtype deltaB, 
                                   int lrB, void *P_dataB, N_Vector tmpB);
  
/******************************************************************
 *                                                                *
 * Type : CVSpgmrJacTimesVecFnB                                   *
 *----------------------------------------------------------------*
 * A Jacobian times vector function jtimesB for the backward      *
 * integration must have the prototype given below.               *
 *                                                                *
 ******************************************************************/

typedef int (*CVSpgmrJacTimesVecFnB)(N_Vector vB, N_Vector JvB, realtype t, 
                                     N_Vector y, N_Vector yB, N_Vector fyB,
                                     void *jac_dataB, N_Vector tmpB);

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
 * Function : CVodeCreateB, CVodeMallocB, CVodeSet*B              *
 *----------------------------------------------------------------*
 * These functions are just wrappers around the corresponding     *
 * functions in cvodes.h, with some particularizations for the    *
 * backward integration.                                          *
 *                                                                *
 ******************************************************************/

int CVodeCreateB(void *cvadj_mem, int lmmB, int iterB);

int CVodeResetIterTypeB(void *cvadj_mem, int iterB);

int CVodeSetFdataB(void *cvadj_mem, void *f_dataB);
int CVodeSetErrFileB(void *cvadj_mem, FILE *errfpB);
int CVodeSetMaxOrdB(void *cvadj_mem, int maxordB);
int CVodeSetMaxNumStepsB(void *cvadj_mem, long int mxstepsB);
int CVodeSetStabLimDetB(void *cvadj_mem, booleantype stldetB);
int CVodeSetInitStepB(void *cvadj_mem, realtype hinB);
int CVodeSetMinStepB(void *cvadj_mem, realtype hminB);
int CVodeSetMaxStepB(void *cvadj_mem, realtype hmaxB);

int CVodeMallocB(void *cvadj_mem, RhsFnB fB, 
                 realtype tB0, N_Vector yB0,
                 int itolB, realtype *reltolB, void *abstolB, 
                 NV_Spec nvspecB);

/* CVodeMallocB return values */
/* SUCCESS=0, defined under CVode return values */
#define CVBM_NO_MEM   -101
#define CVBM_MEM_FAIL -102
#define CVBM_BAD_TB0  -103

/******************************************************************
 *                                                                *
 * Function : CVodeReInitB                                        *
 *----------------------------------------------------------------*
 * CVodeReInitB resets the final time and final condition for the *
 * backward system, assuming prior calls to CVodeCreateB and      *
 * CVodeMallocB have been made.                                   *
 *                                                                *
 ******************************************************************/

int CVodeReInitB(void *cvadj_mem, RhsFnB fB, 
                 realtype tB0, N_Vector yB0,
                 int itolB, realtype *reltolB, void *abstolB);

/******************************************************************
 *                                                                *
 * Function : CVodeSetQuad*B, CVodeQuadMallocB                    *
 *----------------------------------------------------------------*
 *                                                                *
 ******************************************************************/

int CVodeSetQuadErrConB(void *cvadj_mem, int errconQB);
int CVodeSetQuadFdataB(void *cvadj_mem, void *fQ_dataB);

int CVodeQuadMallocB(void *cvadj_mem, QuadRhsFnB fQB,
                     int itolQB, realtype *reltolQB, void *abstolQB,
                     NV_Spec nvspecQB);

/* CVodeQuadMallocB return values               */
/* SUCCESS=0, defined under CVode return values */
/* CVBM_NO_MEM, defined under CVodeMallocB      */
#define CVBM_ILL_INPUT   -104

/******************************************************************
 *                                                                *
 * Function : CVodeQuadReInitB                                    *
 *----------------------------------------------------------------*
 * CVodeQuadReInitB re-initializaes memory for quadrature         *
 * integration during the backward phase                          *
 *                                                                *
 ******************************************************************/

int CVodeQuadReInitB(void *cvadj_mem, QuadRhsFnB fQB,
                     int itolQB, realtype *reltolQB, void *abstolQB);

/******************************************************************
 *                                                                *
 * Function : CVDenseB, CVDenseSet*B                              *
 *----------------------------------------------------------------*
 * CVDenseB links the main CVODE integrator with the CVDENSE      *
 * linear solver for the backward integration.                    *
 *                                                                *
 ******************************************************************/

int CVDenseB(void *cvadj_mem, long int nB);

int CVDenseSetJacFnB(void *cvadj_mem, CVDenseJacFnB djacB);
int CVDenseSetJacDataB(void *cvadj_mem, void *jac_dataB);

/******************************************************************
 *                                                                *
 * Function : CVBandB, CVBandSet*B                                *
 *----------------------------------------------------------------*
 * CVBandB links the main CVODE integrator with the CVBAND        *
 * linear solver for the backward integration.                    *
 *                                                                *
 ******************************************************************/

int CVBandB(void *cvadj_mem, long int nB,
            long int mupperB, long int mlowerB);

int CVBandSetJacFnB(void *cvadj_mem, CVBandJacFnB bjacB);
int CVBandSetJacDataB(void *cvadj_mem, void *jac_dataB);

/******************************************************************
 *                                                                *
 * Function : CVSpgmrB, CVSpgmrSet*B                              *
 *----------------------------------------------------------------*
 * CVSpgmrB links the main CVODE integrator with the CVSPGMR      *
 * linear solver for the backward integration.                    *
 *                                                                *
 ******************************************************************/

int CVSpgmrB(void *cvadj_mem, int pretypeB, int maxlB);

int CVSpgmrResetPrecTypeB(void *cvadj_mem, int pretypeB);

int CVSpgmrSetGSTypeB(void *cvadj_mem, int gstypeB);
int CVSpgmrSetDeltB(void *cvadj_mem, realtype deltB);
int CVSpgmrSetPrecSetupFnB(void *cvadj_mem, CVSpgmrPrecSetupFnB psetB);
int CVSpgmrSetPrecSolveFnB(void *cvadj_mem, CVSpgmrPrecSolveFnB psolveB);
int CVSpgmrSetJacTimesVecFnB(void *cvadj_mem, CVSpgmrJacTimesVecFnB jtimesB);
int CVSpgmrSetPrecDataB(void *cvadj_mem, void *P_dataB);
int CVSpgmrSetJacDataB(void *cvadj_mem, void *jac_dataB);

/******************************************************************
 *                                                                *
 * Function: CVBandPrecAllocB, CVBandPrecSetupB, CVBandPrecSolveB *
 *----------------------------------------------------------------*
 * CVBandPrecAllocB interfaces to the CVBANDPRE preconditioner for*
 * the backward integration. The pointer to the structure         *
 * returned by this routine can then be used together with the    *
 * functions CVBandPrecSetupB and CVBandPrecSolveB in a call      *
 * to CVSpgmrB.                                                   *
 *                                                                *
 * CVBandPrecSetupB is of type CVSpgmrPrecSetupFnB to agree with  *
 * the expected argument in CVSpgmrB.                             *
 *                                                                *
 * CVBandPrecSolveB is of type CVSpgmrPrecSolveFnB to agree with  *
 * the expected argument in CVSpgmrB.                             *
 *                                                                *
 ******************************************************************/

void *CVBandPrecAllocB(void *cvadj_mem, long int nB, 
                       long int muB, long int mlB);

void CVBandPrecFreeB(void *bp_dataB);

int CVBandPrecSetupB(realtype t, N_Vector y, 
                     N_Vector yB, N_Vector fyB, booleantype jokB, 
                     booleantype *jcurPtrB, realtype gammaB,
                     void *bp_dataB,
                     N_Vector tmp1B, N_Vector tmp2B,
                     N_Vector tmp3B);
int CVBandPrecSolveB(realtype t, N_Vector y, 
                     N_Vector yB, N_Vector fyB, 
                     N_Vector rB, N_Vector zB,
                     realtype gammaB, realtype deltaB, 
                     int lrB, void *bp_dataB, N_Vector tmpB);

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
 * Function : CVodeGetQuadB                                       *
 *----------------------------------------------------------------*
 * CVodeGetQuadB extracts values for quadrature variables in      *
 * the N_Vector qB.                                               *
 *                                                                *
 ******************************************************************/

int CVodeGetQuadB(void *cvadj_mem, N_Vector qB);

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
 * Function : CVadjGetCVodeBmem                                   *
 *----------------------------------------------------------------*
 * CVadjGetCVodeBmem returns a (void *) pointer to the CVODES     *
 * memory allocated for the backward problem. This pointer can    *
 * then be used to call any of the CVodeGet* CVODES routines to   *
 * extract optional output for the backward integration phase.    *
 *                                                                *
 ******************************************************************/

void *CVadjGetCVodeBmem(void *cvadj_mem);

/******************************************************************
 *                                                                *
 * Function : CVadjGetY                                           *
 *----------------------------------------------------------------*
 * This routine uses cubic piece-wise Hermite interpolation for   *
 * the forward solution vector.                                   *
 *                                                                *
 ******************************************************************/

int CVadjGetY(void *cvadj_mem, realtype t, N_Vector y);

/* CVadjGetY return values */

enum { GETY_OK=0, GETY_BADT=-1 };

/******************************************************************
 *                                                                *
 * Function : CVadjGetCheckPointsList                             *
 *----------------------------------------------------------------*
 *                                                                *
 ******************************************************************/

void  CVadjGetCheckPointsList(void *cvadj_mem);

/******************************************************************
 *                                                                *
 * Function : CVadjGetStoredData                                  *
 *----------------------------------------------------------------*
 *                                                                *
 ******************************************************************/

void CVadjGetStoredData(void *cvadj_mem, long int which, 
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

  /* Was ck_zn[qmax] allocated?
     ck_zqm = 0    - no
     ck_zqm = qmax - yes      */
  int ck_zqm;
  
  /* Step data */
  long int     ck_nst;
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

  /* Right hand side function (fB) for backward run */
  RhsFnB ca_fB;

  /* Right hand side quadrature function (fQB) for backward run */
  QuadRhsFnB ca_fQB;

  /* Dense Jacobian function (djacB) for backward run */
  CVDenseJacFnB ca_djacB;

  /* Banded Jacobian function (bjacB) for backward run */
  CVBandJacFnB ca_bjacB;

  /* Preconditioner routines (precondB and psolveB) for backward run */
  CVSpgmrPrecSetupFnB ca_psetB;
  CVSpgmrPrecSolveFnB ca_psolveB;

  /* Jac times vec routine (jtimesB) for backward run */
  CVSpgmrJacTimesVecFnB ca_jtimesB;

  /* User f_dataB */
  void *ca_f_dataB;
  
  /* User fQ_dataB */
  void *ca_fQ_dataB;

  /* User jac_dataB */
  void *ca_jac_dataB;

  /* User P_dataB */
  void *ca_P_dataB;
  
  /* Unit roundoff */
  realtype ca_uround;
  
  /* Integration interval */
  realtype ca_tinitial, ca_tfinal;

  /* Time at which to extract quadratures */
  realtype ca_t_for_quad;
  
  /* Number of check points */
  int ca_nckpnts;
  
  /* Number of steps between 2 check points */
  long int ca_nsteps;

  /* Flag to indicate that data in dt_mem is new */
  booleantype ca_newData;

  /* address of the check point structure for which data is available */
  struct CkpntMemRec *ca_ckpntData;

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
