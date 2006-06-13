/*
 * -----------------------------------------------------------------
 * $Revision: 1.49 $
 * $Date: 2006-06-13 01:22:12 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * This is the interface file for the CVODEA adjoint integrator.
 *
 * Function types:
 *    CVRhsFnB
 *    CVQuadRhsFnB
 * Exported functions prototypes:
 *    CVadjMalloc
 *    CVadjSetInterpType
 *    CVodeF
 *    CVodeCreateB
 *    CVodeMallocB
 *    CVodeB
 *    CVadjFree
 * -----------------------------------------------------------------
 */

#ifndef _CVODEA_H
#define _CVODEA_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  /* 
   * ===============================================================
   * INCLUDED HEADER FILES
   * ===============================================================
   */

#include <stdio.h>

#include "cvodes.h"
#include "sundials_nvector.h"

  /* 
   * ===============================================================
   * DEFINITIONS OF CVODEA INPUTS
   * ===============================================================
   */

  /*
   * -----------------------------------------------------------------
   * interp: Specifies the interpolation type used to evaluate the
   *         forward solution during the backward integration phase.
   *         CV_HERMITE specifies cubic Hermite interpolation.
   *         CV_POYNOMIAL specifies the polynomial interpolation
   * -----------------------------------------------------------------
   */
  
#define CV_HERMITE    1
#define CV_POLYNOMIAL 2

  /*
   * ===============================================================
   * CVODEA RETURN VALUES
   * ===============================================================
   */

#define CV_ADJMEM_NULL -101
#define CV_BAD_TB0     -103
#define CV_BCKMEM_NULL -104
#define CV_REIFWD_FAIL -105
#define CV_FWD_FAIL    -106
#define CV_BAD_ITASK   -107
#define CV_BAD_TBOUT   -108
#define CV_GETY_BADT   -109

  /* 
   * ===============================================================
   * FUNCTION TYPES
   * ===============================================================
   */

  /*
   * -----------------------------------------------------------------
   * CVRhsFnB
   *    The fB function which defines the right hand side of the
   *    ODE systems to be integrated backwards must have type CVRhsFnB.
   * -----------------------------------------------------------------
   * CVQuadRhsFnB
   *    The fQB function which defines the quadratures to be integrated
   *    backwards must have type CVQuadRhsFnB.
   * -----------------------------------------------------------------
   */
  
  typedef int (*CVRhsFnB)(realtype t, N_Vector y,
                          N_Vector yB, N_Vector yBdot,
                          void *f_dataB);
  
  typedef int (*CVQuadRhsFnB)(realtype t, N_Vector y,
                              N_Vector yB, N_Vector qBdot,
                              void *fQ_dataB);
     
  /* 
   * ===============================================================
   * EXPORTED FUNCTIONS
   * ===============================================================
   */

  /*
   * -----------------------------------------------------------------
   * CVadjMalloc
   * -----------------------------------------------------------------
   * CVadjMalloc specifies some parameters for the adjoint problem and
   * allocates space for the global CVODEA memory structure.
   * -----------------------------------------------------------------
   */
  
  void *CVadjMalloc(void *cvode_mem, long int steps, int interp);

  /*
   * -----------------------------------------------------------------
   * CVadjSetInterpType
   * -----------------------------------------------------------------
   * Changes the interpolation type. 
   * Must be called only after CVadjMalloc
   * -----------------------------------------------------------------
   */
  
  int CVadjSetInterpType(void *cvadj_mem, int interp);

  /*
   * -----------------------------------------------------------------
   * CVodeF
   * -----------------------------------------------------------------
   * CVodeF integrates towards tout and returns solution into yout.
   * In the same time, it stores check point data every 'steps'.
   *
   * CVodeF can be called repeatedly by the user.
   *
   * ncheckPtr points to the number of check points stored so far.
   *
   * Return values:
   *    CV_SUCCESS
   *    CVADJ_MEM_FAIL
   *    any CVode return value
   * -----------------------------------------------------------------
   */

  int CVodeF(void *cvadj_mem, realtype tout, N_Vector yout,
             realtype *tret, int itask, int *ncheckPtr);

  /*
   * -----------------------------------------------------------------
   * Interfaces to CVODES functions for setting-up the
   *  backward integration
   * -----------------------------------------------------------------
   * CVodeCreateB, CVodeMallocB, CVodeSet*B
   *    These functions are just wrappers around the corresponding
   *    functions in cvodes.h, with some particularizations for the
   *    backward integration.
   * -----------------------------------------------------------------
   * CVodeSetQuad*B, CVodeQuadMallocB, CVodeQuadReInitB
   * -----------------------------------------------------------------
   */

  int CVodeCreateB(void *cvadj_mem, int lmmB, int iterB);

  int CVodeMallocB(void *cvadj_mem, CVRhsFnB fB,
                   realtype tB0, N_Vector yB0,
                   int itolB, realtype reltolB, void *abstolB);
  
  int CVodeSetErrHandlerFnB(void *cvadj_mem, CVErrHandlerFn ehfunB, void *eh_dataB);
  int CVodeSetErrFileB(void *cvadj_mem, FILE *errfpB);
  int CVodeSetIterTypeB(void *cvadj_mem, int iterB);
  int CVodeSetFdataB(void *cvadj_mem, void *f_dataB);
  int CVodeSetMaxOrdB(void *cvadj_mem, int maxordB);
  int CVodeSetMaxNumStepsB(void *cvadj_mem, long int mxstepsB);
  int CVodeSetStabLimDetB(void *cvadj_mem, booleantype stldetB);
  int CVodeSetInitStepB(void *cvadj_mem, realtype hinB);
  int CVodeSetMinStepB(void *cvadj_mem, realtype hminB);
  int CVodeSetMaxStepB(void *cvadj_mem, realtype hmaxB);
  
  int CVodeReInitB(void *cvadj_mem, CVRhsFnB fB,
                   realtype tB0, N_Vector yB0,
                   int itolB, realtype reltolB, void *abstolB);
    
  int CVodeSetQuadFdataB(void *cvadj_mem, void *fQ_dataB);
  int CVodeSetQuadErrConB(void *cvadj_mem, booleantype errconQB,
                          int itolQB, realtype reltolQB, void *abstolQB);
  int CVodeQuadMallocB(void *cvadj_mem, CVQuadRhsFnB fQB, N_Vector yQB0);
  int CVodeQuadReInitB(void *cvadj_mem, CVQuadRhsFnB fQB, N_Vector yQB0);
    
  /*
   * -----------------------------------------------------------------
   * CVodeB
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
   * CVodeGetQuadB
   * -----------------------------------------------------------------
   * CVodeGetQuadB extracts values for quadrature variables in
   * the N_Vector qB.
   * -----------------------------------------------------------------
   */
  
  int CVodeGetQuadB(void *cvadj_mem, N_Vector qB);
  
  /*
   * -----------------------------------------------------------------
   * CVadjFree
   * -----------------------------------------------------------------
   * CVadjFree frees the memory allocated by CVadjMalloc.
   * -----------------------------------------------------------------
   */
  
  void CVadjFree(void **cvadj_mem);
  
  /*
   * -----------------------------------------------------------------
   * CVadjGetCVodeBmem
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
   * The following function returns the name of the constant 
   * associated with a CVODEA-specific return flag
   * -----------------------------------------------------------------
   */
  
  char *CVadjGetReturnFlagName(int flag);


  /*
   * -----------------------------------------------------------------
   * CVadjGetY
   *    Returns the interpolated forward solution at time t. This
   *    function is a wrapper around the interpType-dependent internal
   *    function.
   *    The calling function must allocate space for y.
   * -----------------------------------------------------------------
   */

  int CVadjGetY(void *cvadj_mem, realtype t, N_Vector y);

  /*
   * -----------------------------------------------------------------
   * CVadjGetCheckPointsInfo
   *    Loads an array of nckpnts structures of type CVadjCheckPointRec.
   *    The user must allocate space for ckpnt (ncheck+1).
   * -----------------------------------------------------------------
   */

  typedef struct {
    void *my_addr;
    void *next_addr;
    realtype t0;
    realtype t1;
    long int nstep;
    int order;
    realtype step;
  } CVadjCheckPointRec;

  int CVadjGetCheckPointsInfo(void *cvadj_mem, CVadjCheckPointRec *ckpnt);

  /*
   * -----------------------------------------------------------------
   * CVadjGetDataPointHermite
   *    Returns the 2 vectors stored for cubic Hermite interpolation 
   *    at the data point 'which'. The user must allocate space for
   *    y and yd. Returns CVADJ_MEM_NULL if cvadj_mem is NULL.
   *    Returns CV_ILL_INPUT if interpType != CV_HERMITE.
   * CVadjGetDataPointPolynomial
   *    Returns the vector stored for polynomial interpolation 
   *    at the data point 'which'. The user must allocate space for
   *    y. Returns CVADJ_MEM_NULL if cvadj_mem is NULL.
   *    Returns CV_ILL_INPUT if interpType != CV_POLYNOMIAL.
   * -----------------------------------------------------------------
   */

  int CVadjGetDataPointHermite(void *cvadj_mem, long int which,
                               realtype *t, N_Vector y, N_Vector yd);
  
  int CVadjGetDataPointPolynomial(void *cvadj_mem, long int which,
                                  realtype *t, int *order, N_Vector y);

  /* 
   * ===============================================================
   * DEVELOPMENT USER-CALLABLE FUNCTIONS
   * ===============================================================
   */

  /*
   * -----------------------------------------------------------------
   * CVadjGetCurrentCheckPoint
   *    Returns the address of the 'active' check point.
   * -----------------------------------------------------------------
   */

  int CVadjGetCurrentCheckPoint(void *cvadj_mem, void **addr);

#ifdef __cplusplus
}
#endif

#endif
