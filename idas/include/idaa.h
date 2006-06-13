/*
 * -----------------------------------------------------------------
 * $Revision: 1.19 $
 * $Date: 2006-06-13 01:33:59 $
 * ----------------------------------------------------------------- 
 * Programmers: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/idas/LICENSE
 * -----------------------------------------------------------------
 * This is the interface file for the IDAA adjoint integrator.
 * -----------------------------------------------------------------
 */

#ifndef _IDAA_H
#define _IDAA_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  /* 
   * ===============================================================
   * INCLUDED HEADER FILES
   * ===============================================================
   */

#include <stdio.h>

#include "idas.h"
#include "sundials_nvector.h"

  /* 
   * ===============================================================
   * DEFINITIONS OF IDAA INPUTS
   * ===============================================================
   */

  /*
   * -----------------------------------------------------------------
   * interp: Specifies the interpolation type used to evaluate the
   *         forward solution during the backward integration phase.
   *         IDA_HERMITE specifies cubic Hermite interpolation.
   *         IDA_POYNOMIAL specifies the polynomial interpolation
   * -----------------------------------------------------------------
   */
  
#define IDA_HERMITE    1
#define IDA_POLYNOMIAL 2


  /*
   * ===============================================================
   * IDAA RETURN VALUES
   * ===============================================================
   */

#define IDA_ADJMEM_NULL -101
#define IDA_BAD_TB0     -103
#define IDA_BCKMEM_NULL -104
#define IDA_REIFWD_FAIL -105
#define IDA_FWD_FAIL    -106
#define IDA_BAD_ITASK   -107
#define IDA_BAD_TBOUT   -108
#define IDA_GETY_BADT   -109

  /* 
   * ===============================================================
   * FUNCTION TYPES
   * ===============================================================
   */

  typedef int (*IDAResFnB)(realtype tt, 
                           N_Vector yy, N_Vector yp,
                           N_Vector yyB, N_Vector ypB, N_Vector rrB,
                           void *rdataB);

  typedef void (*IDAQuadRhsFnB)(realtype tt, 
                                N_Vector yy, N_Vector yp, 
                                N_Vector yyB, N_Vector ypB,
                                N_Vector ypQB, void *rdataQB);

  /* 
   * ===============================================================
   * EXPORTED FUNCTIONS
   * ===============================================================
   */


  /* Initialization and optional input for ADJOINT module */

  void *IDAadjMalloc(void *ida_mem, long int steps, int interp);

  int IDAadjSetInterpType(void *idaadj_mem, int interp);

  /* Forward solution function */

  int IDASolveF(void *idaadj_mem, realtype tout, realtype *tret,
                N_Vector yret, N_Vector ypret, int itask, int *ncheckPtr);

  /* Initialization and optional input for backward integration */

  int IDACreateB(void *ida_mem);
  int IDAMallocB(void *idaadj_mem, IDAResFnB resB,
                 realtype tB0, N_Vector yyB0, N_Vector ypB0, 
                 int itolB, realtype *reltolB, void *abstolB);
  int IDAReInitB(void *idaadj_mem, IDAResFnB resB,
                 realtype tB0, N_Vector yyB0, N_Vector ypB0,
                 int itolB, realtype *reltolB, void *abstolB);

  int IDASetErrHandlerFnB(void *idaadj_mem, IDAErrHandlerFn ehfunB, void *eh_dataB);
  int IDASetErrFileB(void *idaadj_mem, FILE *errfpB);
  int IDASetRdataB(void *idaadj_mem, void *res_dataB);
  int IDASetMaxOrdB(void *idaadj_mem, int maxordB);
  int IDASetMaxNumStepsB(void *idaadj_mem, long int mxstepsB);
  int IDASetInitStepB(void *idaadj_mem, realtype hinB);
  int IDASetMaxStepB(void *idaadj_mem, realtype hmaxB);
  int IDASetSuppressAlgB(void *idaadj_mem, booleantype suppressalgB);
  int IDASetIdB(void *idaadj_mem, N_Vector idB);
  int IDASetConstraintsB(void *idaadj_mem, N_Vector constraintsB);

  int IDASetQuadFdataB(void *idaadj_mem, void *rhsQ_dataB);
  int IDASetQuadErrConB(void *idaadj_mem, booleantype errconQB, 
                        int itolQB, realtype reltolQB, void *abstolQB);
  int IDAQuadMallocB(void *idaadj_mem, IDAQuadRhsFnB rhsQB, N_Vector yQB0);
  int IDAQuadReInitB(void *idaadj_mem, IDAQuadRhsFnB rhsQB, N_Vector yQB0);

  /* Backward solution function */

  int IDASolveB(void *idaadj_mem, realtype tBout, realtype *tBret,
                N_Vector yBret, N_Vector ypBret, int itaskB);

  /* Optional output from backward integration */

  int IDAGetQuadB(void *idaadj_mem, N_Vector qB);

  /* Deallocation of ADJOINT module */

  void IDAadjFree(void **idaadj_mem);

  /* Optional output from ADJOINT module */

  void *IDAadjGetIDABmem(void *idaadj_mem);

  char *IDAadjGetReturnFlagName(int flag);

  /*
   * IDAadjGetY
   *    Returns the interpolated forward solution at time t. This
   *    function is a wrapper around the interpType-dependent internal
   *    function.
   *    The calling function must allocate space for y.
   */

  int IDAadjGetY(void *idaadj_mem, realtype t, N_Vector y, N_Vector yp);

  /* 
   * ===============================================================
   * DEVELOPMENT USER-CALLABLE FUNCTIONS
   * ===============================================================
   */

  /*
  typedef struct {
    void *my_addr;
    void *next_addr;
    realtype t0;
    realtype t1;
    long int nstep;
    int order;
    realtype step;
  } IDAadjCheckPointRec;

  int IDAadjGetCheckPointsInfo(void *idaadj_mem, IDAadjCheckPointRec *ckpnt);
  int IDAadjGetCurrentCheckPoint(void *idaadj_mem, void **addr);

  int IDAadjGetDataPointHermite(void *idaadj_mem, long int which,
                                realtype *t, N_Vector y, N_Vector yd);
  
  int CVadjGetDataPointPolynomial(void *idaadj_mem, long int which,
                                  realtype *t, int *order, N_Vector y);
  
  */

#ifdef __cplusplus
}
#endif

#endif
