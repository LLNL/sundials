/*
 * -----------------------------------------------------------------
 * $Revision: 1.33 $
 * $Date: 2005-04-07 23:28:43 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * This is the interface file for the CVODEA adjoint integrator.
 *
 * Function types:
 *    CVRhsFnB
 *    CVQuadRhsFnB
 *    CVDenseJacFnB
 *    CVBandJacFnB
 *    CVSpbcgPrecSetupFnB
 *    CVSpbcgPrecSolveB
 *    CVSpbcgJacTimesVecFnB
 *    CVSpgmrPrecSetupFnB
 *    CVSpgmrPrecSolveB
 *    CVSpgmrJacTimesVecFnB
 * Exported functions prototypes:
 *    CVadjMalloc
 *    CVodeF
 *    CVodeCreateB
 *    CVodeMallocB
 *    CVDenseB
 *    CVBandB
 *    CVSpbcgB
 *    CVSpgmrB
 *    CVBandPrecAllocB
 *    CVBPSpbcgB
 *    CVBPSpgmrB
 *    CVBBDPrecAllocB
 *    CVBBDPrecReInit
 *    CVBBDSpbcgB
 *    CVBBDSpgmrB
 *    CVodeB
 *    CVadjFree
 *    CVadjGetY
 *    CVadjCheckPointsList
 *    CVadjDataExtract
 * Type definitions:
 *    struct CkpntMemRec, CkpntMem
 *    struct DtpntMemRec, DtpntMem
 *    struct CVadjMemRec, CVadjMem
 * -----------------------------------------------------------------
 */

#ifndef _CVODEA_H
#define _CVODEA_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <stdio.h>

#include "dense.h"
#include "band.h"
#include "spbcg.h"
#include "spgmr.h"
#include "sundialstypes.h"
#include "nvector.h"

/*
 * -----------------------------------------------------------------
 * Type : CVRhsFnB
 * -----------------------------------------------------------------
 * The fB function which defines the right hand side of the
 * ODE systems to be integrated backwards must have type CVRhsFnB.
 * -----------------------------------------------------------------
 */

typedef void (*CVRhsFnB)(realtype t, N_Vector y,
                         N_Vector yB, N_Vector yBdot,
                         void *f_dataB);

/*
 * -----------------------------------------------------------------
 * Type : CVQuadRhsFnB
 * -----------------------------------------------------------------
 * The fQB function which defines the quadratures to be integrated
 * backwards must have type CVQuadRhsFnB.
 * -----------------------------------------------------------------
 */

typedef void (*CVQuadRhsFnB)(realtype t, N_Vector y,
                             N_Vector yB, N_Vector qBdot,
                             void *fQ_dataB);

/*
 * -----------------------------------------------------------------
 * Type : CVDenseJacFnB
 * -----------------------------------------------------------------
 * A dense Jacobian approximation function djacB for the backward
 * integration must have the prototype given below.
 * -----------------------------------------------------------------
 */

typedef void (*CVDenseJacFnB)(long int nB, DenseMat JB, realtype t,
                              N_Vector y, N_Vector yB, N_Vector fyB,
                              void *jac_dataB, N_Vector tmp1B,
                              N_Vector tmp2B, N_Vector tmp3B);

/*
 * -----------------------------------------------------------------
 * Type : CVBandJacFnB
 * -----------------------------------------------------------------
 * A band Jacobian approximation function bjacB for the backward
 * integration must have the prototype given below.
 * -----------------------------------------------------------------
 */

typedef void (*CVBandJacFnB)(long int nB, long int mupperB,
                             long int mlowerB, BandMat JB,
                             realtype t, N_Vector y,
                             N_Vector yB, N_Vector fyB,
                             void *jac_dataB, N_Vector tmp1B,
                             N_Vector tmp2B, N_Vector tmp3B);

/*
 * -----------------------------------------------------------------
 * Type : CVSpbcgPrecSetupFnB
 * -----------------------------------------------------------------
 * A preconditioner setup function precondB for the backward
 * integration must have the prototype given below.
 * -----------------------------------------------------------------
 */

typedef int (*CVSpbcgPrecSetupFnB)(realtype t, N_Vector y,
                                   N_Vector yB, N_Vector fyB,
                                   booleantype jokB,
				   booleantype *jcurPtrB, realtype gammaB,
                                   void *P_dataB,
                                   N_Vector tmp1B, N_Vector tmp2B,
                                   N_Vector tmp3B);

/*
 * -----------------------------------------------------------------
 * Type : CVSpgmrPrecSetupFnB
 * -----------------------------------------------------------------
 * A preconditioner setup function precondB for the backward
 * integration must have the prototype given below.
 * -----------------------------------------------------------------
 */

typedef int (*CVSpgmrPrecSetupFnB)(realtype t, N_Vector y,
                                   N_Vector yB, N_Vector fyB,
                                   booleantype jokB,
				   booleantype *jcurPtrB, realtype gammaB,
                                   void *P_dataB,
                                   N_Vector tmp1B, N_Vector tmp2B,
                                   N_Vector tmp3B);

/*
 * -----------------------------------------------------------------
 * Type : CVSpbcgPrecSolveFnB
 * -----------------------------------------------------------------
 * A preconditioner solve function psolveB for the backward
 * integration must have the prototype given below.
 * -----------------------------------------------------------------
 */

typedef int (*CVSpbcgPrecSolveFnB)(realtype t, N_Vector y,
                                   N_Vector yB, N_Vector fyB,
                                   N_Vector rB, N_Vector zB,
                                   realtype gammaB, realtype deltaB,
                                   int lrB, void *P_dataB, N_Vector tmpB);

/*
 * -----------------------------------------------------------------
 * Type : CVSpgmrPrecSolveFnB
 * -----------------------------------------------------------------
 * A preconditioner solve function psolveB for the backward
 * integration must have the prototype given below.
 * -----------------------------------------------------------------
 */

typedef int (*CVSpgmrPrecSolveFnB)(realtype t, N_Vector y,
                                   N_Vector yB, N_Vector fyB,
                                   N_Vector rB, N_Vector zB,
                                   realtype gammaB, realtype deltaB,
                                   int lrB, void *P_dataB, N_Vector tmpB);

/*
 * -----------------------------------------------------------------
 * Type : CVSpbcgJacTimesVecFnB
 * -----------------------------------------------------------------
 * A Jacobian times vector function jtimesB for the backward
 * integration must have the prototype given below.
 * -----------------------------------------------------------------
 */

typedef int (*CVSpbcgJacTimesVecFnB)(N_Vector vB, N_Vector JvB, realtype t,
                                     N_Vector y, N_Vector yB, N_Vector fyB,
                                     void *jac_dataB, N_Vector tmpB);

/*
 * -----------------------------------------------------------------
 * Type : CVSpgmrJacTimesVecFnB
 * -----------------------------------------------------------------
 * A Jacobian times vector function jtimesB for the backward
 * integration must have the prototype given below.
 * -----------------------------------------------------------------
 */

typedef int (*CVSpgmrJacTimesVecFnB)(N_Vector vB, N_Vector JvB, realtype t,
                                     N_Vector y, N_Vector yB, N_Vector fyB,
                                     void *jac_dataB, N_Vector tmpB);

/*
 * -----------------------------------------------------------------
 * Type : CVLocalFnB and CVCommFnB
 * -----------------------------------------------------------------
 * Local approximation function and inter-process communication
 * function for the BBD preconditioner on the backward phase.
 * -----------------------------------------------------------------
 */

typedef void (*CVLocalFnB)(long int NlocalB, realtype t,
                           N_Vector y, N_Vector yB, N_Vector gB,
                           void *f_dataB);

typedef void (*CVCommFnB)(long int NlocalB, realtype t,
                          N_Vector y, N_Vector yB,
                          void *f_dataB);

/*
 * -----------------------------------------------------------------
 * Function : CVadjMalloc
 * -----------------------------------------------------------------
 * CVadjMalloc space for the global CVODEA memory structure.
 * -----------------------------------------------------------------
 */

void *CVadjMalloc(void *cvode_mem, long int steps);

/*
 * -----------------------------------------------------------------
 * Function : CVodeF
 * -----------------------------------------------------------------
 * CVodeF integrates towards tout and returns solution into yout.
 * In the same time, it stores check point data every 'steps'.
 *
 * CVodeF can be called repeatedly by the user.
 *
 * ncheckPtr points to the number of check points stored so far.
 *
 * Return values:
 *    SUCCESS
 *    CVADJ_MEM_FAIL
 *    any CVode return value
 * -----------------------------------------------------------------
 */

int CVodeF(void *cvadj_mem, realtype tout, N_Vector yout,
           realtype *tret, int itask, int *ncheckPtr);

/*
 * -----------------------------------------------------------------
 * Function : CVodeCreateB, CVodeMallocB, CVodeSet*B
 * -----------------------------------------------------------------
 * These functions are just wrappers around the corresponding
 * functions in cvodes.h, with some particularizations for the
 * backward integration.
 * -----------------------------------------------------------------
 */

int CVodeCreateB(void *cvadj_mem, int lmmB, int iterB);

int CVodeSetIterTypeB(void *cvadj_mem, int iterB);

int CVodeSetFdataB(void *cvadj_mem, void *f_dataB);
int CVodeSetErrFileB(void *cvadj_mem, FILE *errfpB);
int CVodeSetMaxOrdB(void *cvadj_mem, int maxordB);
int CVodeSetMaxNumStepsB(void *cvadj_mem, long int mxstepsB);
int CVodeSetStabLimDetB(void *cvadj_mem, booleantype stldetB);
int CVodeSetInitStepB(void *cvadj_mem, realtype hinB);
int CVodeSetMinStepB(void *cvadj_mem, realtype hminB);
int CVodeSetMaxStepB(void *cvadj_mem, realtype hmaxB);

int CVodeMallocB(void *cvadj_mem, CVRhsFnB fB,
                 realtype tB0, N_Vector yB0,
                 int itolB, realtype reltolB, void *abstolB);

int CVodeReInitB(void *cvadj_mem, CVRhsFnB fB,
                 realtype tB0, N_Vector yB0,
                 int itolB, realtype reltolB, void *abstolB);

/*
 * -----------------------------------------------------------------
 * Function : CVodeSetQuad*B, CVodeQuadMallocB, CVodeQuadReInitB
 * -----------------------------------------------------------------
 */

int CVodeSetQuadFdataB(void *cvadj_mem, void *fQ_dataB);
int CVodeSetQuadErrConB(void *cvadj_mem, booleantype errconQB,
                        int itolQB, realtype reltolQB, void *abstolQB);
int CVodeQuadMallocB(void *cvadj_mem, CVQuadRhsFnB fQB, N_Vector yQB0);
int CVodeQuadReInitB(void *cvadj_mem, CVQuadRhsFnB fQB, N_Vector yQB0);

/*
 * -----------------------------------------------------------------
 * Function : CVDenseB, CVDenseSet*B
 * -----------------------------------------------------------------
 * CVDenseB links the main CVODE integrator with the CVDENSE
 * linear solver for the backward integration.
 * -----------------------------------------------------------------
 */

int CVDenseB(void *cvadj_mem, long int nB);

int CVDenseSetJacFnB(void *cvadj_mem, CVDenseJacFnB djacB, void *jac_dataB);

/*
 * -----------------------------------------------------------------
 * Function : CVDiagB
 * -----------------------------------------------------------------
 * CVDiagB links the main CVODE integrator with the CVDIAG
 * linear solver for the backward integration.
 * -----------------------------------------------------------------
 */

int CVDiagB(void *cvadj_mem);

/*
 * -----------------------------------------------------------------
 * Function : CVBandB, CVBandSet*B
 * -----------------------------------------------------------------
 * CVBandB links the main CVODE integrator with the CVBAND
 * linear solver for the backward integration.
 * -----------------------------------------------------------------
 */

int CVBandB(void *cvadj_mem, long int nB,
            long int mupperB, long int mlowerB);

int CVBandSetJacFnB(void *cvadj_mem, CVBandJacFnB bjacB, void *jac_dataB);

/*
 * -----------------------------------------------------------------
 * Function : CVSpbcgB, CVSpbcgSet*B
 * -----------------------------------------------------------------
 * CVSpbcgB links the main CVODE integrator with the CVSPBCG
 * linear solver for the backward integration.
 * -----------------------------------------------------------------
 */

int CVSpbcgB(void *cvadj_mem, int pretypeB, int maxlB);

int CVSpbcgSetPrecTypeB(void *cvadj_mem, int pretypeB);
int CVSpbcgSetDeltB(void *cvadj_mem, realtype deltB);
int CVSpbcgSetPreconditionerB(void *cvadj_mem, CVSpbcgPrecSetupFnB psetB,
                              CVSpbcgPrecSolveFnB psolveB, void *P_dataB);
int CVSpbcgSetJacTimesVecFnB(void *cvadj_mem, 
                             CVSpbcgJacTimesVecFnB jtimesB, void *jac_dataB);

/*
 * -----------------------------------------------------------------
 * Function : CVSpgmrB, CVSpgmrSet*B
 * -----------------------------------------------------------------
 * CVSpgmrB links the main CVODE integrator with the CVSPGMR
 * linear solver for the backward integration.
 * -----------------------------------------------------------------
 */

int CVSpgmrB(void *cvadj_mem, int pretypeB, int maxlB);

int CVSpgmrSetPrecTypeB(void *cvadj_mem, int pretypeB);
int CVSpgmrSetGSTypeB(void *cvadj_mem, int gstypeB);
int CVSpgmrSetDeltB(void *cvadj_mem, realtype deltB);
int CVSpgmrSetPreconditionerB(void *cvadj_mem, CVSpgmrPrecSetupFnB psetB,
                              CVSpgmrPrecSolveFnB psolveB, void *P_dataB);
int CVSpgmrSetJacTimesVecFnB(void *cvadj_mem, CVSpgmrJacTimesVecFnB jtimesB,
                             void *jac_dataB);
/*
 * -----------------------------------------------------------------
 * Function: CVBandPrecAllocB, CVBPSpgmrB, CVBPSpbcgB
 * -----------------------------------------------------------------
 * CVBandPrecAllocB interfaces to the CVBANDPRE preconditioner for
 * the backward integration. The pointer to the structure
 * returned by this routine should then be used in the call to
 * CVBPSpgmrB/CVBPSpbcgB which interfaces to CVBPSpgmr/CVBPSpbcg.
 * -----------------------------------------------------------------
 */

int CVBandPrecAllocB(void *cvadj_mem, long int nB,
                     long int muB, long int mlB);

int CVBPSpbcgB(void *cvadj_mem, int pretypeB, int maxlB);
int CVBPSpgmrB(void *cvadj_mem, int pretypeB, int maxlB);

/*
 * -----------------------------------------------------------------
 * Functions: CVBBDPrecAllocB, CVBBDSpgmrB, CVBBDSpbcgB,
 *            CVBBDPrecReInit
 * -----------------------------------------------------------------
 * Interface functions for the BBD preconditioner to be used on
 * the backward phase.
 * -----------------------------------------------------------------
 */

int CVBBDPrecAllocB(void *cvadj_mem, long int NlocalB,
                    long int mudqB, long int mldqB,
                    long int mukeepB, long int mlkeepB,
                    realtype dqrelyB,
                    CVLocalFnB glocB, CVCommFnB cfnB);

int CVBBDSpbcgB(void *cvadj_mem, int pretypeB, int maxlB);
int CVBBDSpgmrB(void *cvadj_mem, int pretypeB, int maxlB);

int CVBBDPrecReInitB(void *cvadj_mem, long int mudqB, long int mldqB,
                     realtype dqrelyB, CVLocalFnB glocB, CVCommFnB cfnB);

/*
 * -----------------------------------------------------------------
 * Function : CVodeB
 * -----------------------------------------------------------------
 * CVodeB performs the backward integration from tfinal to
 * tinitial through a sequence of forward-backward runs in
 * between consecutive check points. It returns the values of
 * the adjoint variables and any existing quadrature variables
 * at tinitial.
 * -----------------------------------------------------------------
 */

int CVodeB(void *cvadj_mem, realtype tBout, N_Vector yBout,
           realtype *tBret, int itaskB);

/*
 * -----------------------------------------------------------------
 * Function : CVodeGetQuadB
 * -----------------------------------------------------------------
 * CVodeGetQuadB extracts values for quadrature variables in
 * the N_Vector qB.
 * -----------------------------------------------------------------
 */

int CVodeGetQuadB(void *cvadj_mem, N_Vector qB);

/*
 * -----------------------------------------------------------------
 * Function : CVadjFree
 * -----------------------------------------------------------------
 * CVadjFree frees the memory allocated by CVadjMalloc.
 * -----------------------------------------------------------------
 */

void CVadjFree(void *cvadj_mem);

/*
 * -----------------------------------------------------------------
 * Function : CVadjGetCVodeBmem
 * -----------------------------------------------------------------
 * CVadjGetCVodeBmem returns a (void *) pointer to the CVODES
 * memory allocated for the backward problem. This pointer can
 * then be used to call any of the CVodeGet* CVODES routines to
 * extract optional output for the backward integration phase.
 * -----------------------------------------------------------------
 */

void *CVadjGetCVodeBmem(void *cvadj_mem);

/*
 * -----------------------------------------------------------------
 * Function : CVadjGetY
 * -----------------------------------------------------------------
 * This routine uses cubic piece-wise Hermite interpolation for
 * the forward solution vector.
 * -----------------------------------------------------------------
 */

int CVadjGetY(void *cvadj_mem, realtype t, N_Vector y);

/*
 * -----------------------------------------------------------------
 * Function : CVadjGetCheckPointsList
 * -----------------------------------------------------------------
 */

void  CVadjGetCheckPointsList(void *cvadj_mem);

/*
 * -----------------------------------------------------------------
 * Function : CVadjGetStoredData
 * -----------------------------------------------------------------
 */

void CVadjGetStoredData(void *cvadj_mem, long int which,
                        realtype *t, N_Vector yout, N_Vector ydout);

/*
 * -----------------------------------------------------------------
 * CVODEA return values
 * -----------------------------------------------------------------
 */

#define CV_ADJMEM_NULL -101
#define CV_BAD_TB0     -103
#define CV_BCKMEM_NULL -104
#define CV_REIFWD_FAIL -105
#define CV_FWD_FAIL    -106
#define CV_BAD_ITASK   -107
#define CV_BAD_TBOUT   -108
#define CV_GETY_BADT   -109

#ifdef __cplusplus
}
#endif

#endif
