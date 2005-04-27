/*
 * -----------------------------------------------------------------
 * $Revision: 1.14 $
 * $Date: 2005-04-27 22:51:47 $
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

#include <stdio.h>
#include "dense.h"
#include "band.h"
#include "spgmr.h"
#include "sundialstypes.h"
#include "nvector.h"
  
  /******************************************************************
   *                                                                *
   * Type : IDAResFnB                                                  *
   *----------------------------------------------------------------*
   *                                                                *
   ******************************************************************/
  
  typedef int (*IDAResFnB)(realtype tt, 
                           N_Vector yy, N_Vector yp,
                           N_Vector yyB, N_Vector ypB, N_Vector rrB,
                           void *rdataB);

  /******************************************************************
   *                                                                *
   * Type : IDAQuadRhsFnB                                              *
   *----------------------------------------------------------------*
   * The fQB function which defines the quadratures to be integrated*
   * backwards must have type IDAQuadRhsFnB.                           *
   *                                                                *
   ******************************************************************/
  
  typedef void (*IDAQuadRhsFnB)(realtype tt, 
                                N_Vector yy, N_Vector yp, 
                                N_Vector yyB, N_Vector ypB,
                                N_Vector ypQB, void *rdataQB);
  
  /******************************************************************
   *                                                                *
   * Type : IDADenseJacFnB                                           *
   *----------------------------------------------------------------*
   * A dense Jacobian approximation function djacB for the backward *
   * integration must have the prototype given below.               *
   *                                                                *
   ******************************************************************/
  
  typedef int (*IDADenseJacFnB)(long int NeqB, realtype tt, 
                                N_Vector yy, N_Vector yp,
                                N_Vector yyB, N_Vector ypB, N_Vector rrB,
                                realtype c_jB, void *jac_dataB, 
                                DenseMat JacB, 
                                N_Vector tmp1B, N_Vector tmp2B, 
                                N_Vector tmp3B);

  /******************************************************************
   *                                                                *
   * Type : IDABandJacFnB                                            *
   *----------------------------------------------------------------*
   * A band Jacobian approximation function bjacB for the backward  *
   * integration must have the prototype given below.               *
   *                                                                *
   ******************************************************************/

  typedef int (*IDABandJacFnB)(long int NeqB, 
                               long int mupperB, long int mlowerB, 
                               realtype tt, 
                               N_Vector yy, N_Vector yp,
                               N_Vector yyB, N_Vector ypB, N_Vector rrB,
                               realtype c_jB, void *jac_dataB,
                               BandMat JacB, 
                               N_Vector tmp1B, N_Vector tmp2B, 
                               N_Vector tmp3B);

  /******************************************************************
   *                                                                *
   * Type : IDASpgmrPrecSetupFnB                                     *
   *----------------------------------------------------------------*
   * A preconditioner setup function precondB for the backward      *
   * integration must have the prototype given below.               *
   *                                                                *
   ******************************************************************/

  typedef int (*IDASpgmrPrecSetupFnB)(realtype tt, 
                                      N_Vector yy, N_Vector yp,
                                      N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                                      realtype c_jB, void *prec_dataB,
                                      N_Vector tmp1B, N_Vector tmp2B, 
                                      N_Vector tmp3B);

  /******************************************************************
   *                                                                *
   * Type : IDASpgmrPrecSolveFnB                                     *
   *----------------------------------------------------------------*
   * A preconditioner solve function psolveB for the backward       *
   * integration must have the prototype given below.               *
   *                                                                *
   ******************************************************************/

  typedef int (*IDASpgmrPrecSolveFnB)(realtype tt, 
                                      N_Vector yy, N_Vector yp,
                                      N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                                      N_Vector rvecB, N_Vector zvecB,
                                      realtype c_jB, realtype deltaB,
                                      void *prec_dataB, N_Vector tmpB);

  /******************************************************************
   *                                                                *
   * Type : IDASpgmrJacTimesVecFnB                                   *
   *----------------------------------------------------------------*
   * A Jacobian times vector function jtimesB for the backward      *
   * integration must have the prototype given below.               *
   *                                                                *
   ******************************************************************/

  typedef int (*IDASpgmrJacTimesVecFnB)(realtype t,
                                        N_Vector yy, N_Vector yp,
                                        N_Vector yyB, N_Vector ypB, N_Vector rrB,
                                        N_Vector vB, N_Vector JvB, 
                                        realtype c_jB, void *jac_dataB, 
                                        N_Vector tmp1B, N_Vector tmp2B);

  /******************************************************************
   *                                                                *
   * Function : IDAAdjMalloc                                         *
   *----------------------------------------------------------------*
   * IDAAdjMalloc space for the global IDAA memory structure.      *
   *                                                                *
   ******************************************************************/

  void *IDAAdjMalloc(void *ida_mem, long int steps);

  /******************************************************************
   *                                                                *
   * Function : IDAAdjFree                                           *
   *----------------------------------------------------------------*
   * IDAAdjFree frees the memory allocated by IDAAdjMalloc.           *
   *                                                                *
   ******************************************************************/

  void IDAAdjFree(void *idaadj_mem);

  /******************************************************************
   *                                                                *
   * Function : IDAAdjGetIDABmem                                   *
   *----------------------------------------------------------------*
   * IDAAdjGetIDABmem returns a (void *) pointer to the IDAS     *
   * memory allocated for the backward problem. This pointer can    *
   * then be used to call any of the IDAGet* IDAS routines to   *
   * extract optional output for the backward integration phase.    *
   *                                                                *
   ******************************************************************/

  void *IDAAdjGetIDABmem(void *idaadj_mem);

  /******************************************************************
   *                                                                *
   * Function : IDAAdjCheckPointsList                               *
   *----------------------------------------------------------------*
   *                                                                *
   ******************************************************************/

  void  IDAAdjCheckPointsList(void *idaadj_mem);

  /******************************************************************
   *                                                                *
   * Function : IDASolveF                                              *
   *----------------------------------------------------------------*
   * IDASolveF integrates towards tout and returns solution into yout. *
   * In the same time, it stores check point data every 'steps'.    *
   *                                                                *
   * IDASolveF can be called repeatedly by the user. The last tout     *
   * will be used as the starting time for the backward integration.*
   *                                                                *
   * ncheckPtr points to the number of check points stored so far.  *
   *                                                                *
   * Return:
   * IDA_MEM_FAIL or any IDASolve return value
   ******************************************************************/

  int IDASolveF(void *idaadj_mem, realtype tout, realtype *tret,
                N_Vector yret, N_Vector ypret, int itask, int *ncheckPtr);


  /******************************************************************
   *                                                                *
   * Function : IDACreateB, IDAMallocB, IDASet*B              *
   *----------------------------------------------------------------*
   * These functions are just wrappers around the corresponding     *
   * functions in idas.h, with some particularizations for the    *
   * backward integration.                                          *
   *                                                                *
   ******************************************************************/

  int IDACreateB(void *ida_mem);
  
  int IDASetRdataB(void *idaadj_mem, void *res_dataB);
  int IDASetErrFileB(void *idaadj_mem, FILE *errfpB);
  int IDASetMaxOrdB(void *idaadj_mem, int maxordB);
  int IDASetMaxNumStepsB(void *idaadj_mem, long int mxstepsB);
  int IDASetInitStepB(void *idaadj_mem, realtype hinB);
  int IDASetMaxStepB(void *idaadj_mem, realtype hmaxB);
  int IDASetSuppressAlgB(void *idaadj_mem, booleantype suppressalgB);
  int IDASetIdB(void *idaadj_mem, N_Vector idB);
  int IDASetConstraintsB(void *idaadj_mem, N_Vector constraintsB);

  int IDAMallocB(void *idaadj_mem, IDAResFnB resB,
                 realtype tB0, N_Vector yyB0, N_Vector ypB0, 
                 int itolB, realtype *reltolB, void *abstolB);
  
  int IDAReInitB(void *idaadj_mem, IDAResFnB resB,
                 realtype tB0, N_Vector yyB0, N_Vector ypB0,
                 int itolB, realtype *reltolB, void *abstolB);

  /******************************************************************
   *                                                                *
   * Function : IDAGetMemB                                    *
   *----------------------------------------------------------------*
   * IDAGetMemB returns a (void *) pointer to the IDAS     *
   * memory allocated for the backward problem. This pointer can    *
   * then be used to call any of the IDAGet* IDAS routines to   *
   * extract optional output for the backward integration phase.    *
   *                                                                *
   ******************************************************************/

  void *IDAGetMemB(void *idaadj_mem);


  /******************************************************************
   *                                                                *
   * Function : IDASetQuad*B, IDAQuadMallocB                    *
   *----------------------------------------------------------------*
   *                                                                *
   ******************************************************************/

  int IDASetQuadErrConB(void *idaadj_mem, booleantype errconQB);
  int IDASetQuadRdataB(void *idaadj_mem, void *rhs_dataQB);
  int IDASetQuadTolerancesB(void *idaadj_mem, int itolQB, 
                            realtype *reltolQB, void *abstolQB);

  int IDAQuadMallocB(void *idaadj_mem, IDAQuadRhsFnB rhsQB, N_Vector yQB0);

  /******************************************************************
   *                                                                *
   * Function : IDAQuadReInitB                                    *
   *----------------------------------------------------------------*
   * IDAQuadReInitB re-initializaes memory for quadrature         *
   * integration during the backward phase                          *
   *                                                                *
   ******************************************************************/

  int IDAQuadReInitB(void *idaadj_mem, IDAQuadRhsFnB rhsQB, N_Vector yQB0);

  /******************************************************************
   *                                                                *
   * Function : IDADenseB, IDADenseSet*B                              *
   *----------------------------------------------------------------*
   * IDADenseB links the main IDAS integrator with the IDADENSE      *
   * linear solver for the backward integration.                    *
   *                                                                *
   ******************************************************************/
  
  int IDADenseB(void *idaadj_mem, long int NeqB);
  
  int IDADenseSetJacFnB(void *idaadj_mem, IDADenseJacFnB djacB, void *jdataB);

  /******************************************************************
   *                                                                *
   * Function : IDABandB, IDABandSet*B                                *
   *----------------------------------------------------------------*
   * IDABandB links the main IDAS integrator with the IDABAND        *
   * linear solver for the backward integration.                    *
   *                                                                *
   ******************************************************************/

  int IDABandB(void *idaadj_mem, long int NeqB,
               long int mupperB, long int mlowerB);

  int IDABandSetJacFnB(void *idaadj_mem, IDABandJacFnB bjacB, void *jdataB);

  /******************************************************************
   *                                                                *
   * Function : IDASpgmrB, IDASpgmrSet*B                              *
   *----------------------------------------------------------------*
   * IDASpgmrB links the main IDAS integrator with the IDASPGMR      *
   * linear solver for the backward integration.                    *
   *                                                                *
   ******************************************************************/

  int IDASpgmrB(void *idaadj_mem, int maxlB);
  
  int IDASpgmrSetGSTypeB(void *idaadj_mem, int gstypeB);
  int IDASpgmrSetMaxRestartsB(void *idaadj_mem, int maxrsB);
  int IDASpgmrSetEpsLinB(void *idaadj_mem, realtype eplifacB);
  int IDASpgmrSetIncrementFactorB(void *idaadj_mem, realtype dqincfacB);
  int IDASpgmrSetPreconditionerB(void *idaadj_mem, IDASpgmrPrecSetupFnB psetB,
				 IDASpgmrPrecSolveFnB psolveB, void *pdataB);
  int IDASpgmrSetJacTimesVecFnB(void *idaadj_mem, IDASpgmrJacTimesVecFnB jtimesB, void *jdataB);

  /******************************************************************
   *                                                                *
   * Function : IDASolveB                                              *
   *----------------------------------------------------------------*
   * IDASolveB performs the backward integration from tfinal to        *
   * tinitial through a sequence of forward-backward runs in        *
   * between consecutive check points. It returns the values of     *
   * the adjoint variables and any existing quadrature variables    *
   * at tinitial.                                                   *
   *                                                                *
   ******************************************************************/

  int IDASolveB(void *idaadj_mem, N_Vector yyB, N_Vector ypB);

  /******************************************************************
   *                                                                *
   * Function : IDAGetQuadB                                       *
   *----------------------------------------------------------------*
   * IDAGetQuadB extracts values for quadrature variables in      *
   * the N_Vector qB.                                               *
   *                                                                *
   ******************************************************************/

  int IDAGetQuadB(void *idaadj_mem, N_Vector qB);

  /******************************************************************
   * Debugging routines....                                         *
   *----------------------------------------------------------------*
   *                                                                *
   ******************************************************************/
  
  int IDAAloadData(void *idaadj_mem, int which_ckpnt, long int *points);
  void IDAAgetData(void *idaadj_mem, long int which_pnt, 
                   realtype *t, N_Vector yout, N_Vector ydout);

#ifdef __cplusplus
}
#endif

#endif
