/*******************************************************************
 * File          : idaa.h                                          *
 * Programmers   : Radu Serban @ LLNL                              *
 * Version of    : 07 February 2004                                *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/cvodes/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the interface file for the IDAA adjoint integrator.     *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif
  
#ifndef _idaa_h
#define _idaa_h
  
#include <stdio.h>
#include "nvector.h"
#include "idas.h"
#include "idadense.h"
#include "idaband.h"
#include "idaspgmr.h"
  
  /******************************************************************
   *                                                                *
   * Type : ResFnB                                                  *
   *----------------------------------------------------------------*
   *                                                                *
   ******************************************************************/
  
  typedef int (*ResFnB)(realtype t, N_Vector yy, N_Vector yp,
                        N_Vector yyB, N_Vector ypB, N_Vector resvalB,
                        void *rdataB);

  /******************************************************************
   *                                                                *
   * Type : QuadRhsFnB                                              *
   *----------------------------------------------------------------*
   * The fQB function which defines the quadratures to be integrated*
   * backwards must have type QuadRhsFnB.                           *
   *                                                                *
   ******************************************************************/
  
  typedef void (*QuadRhsFnB)(realtype t, N_Vector yy, N_Vector yp, 
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
  
  typedef int (*IDADenseJacFnB)(long int NeqB, realtype t, 
                                N_Vector yy, N_Vector yp,
                                N_Vector yyB, N_Vector ypB, realtype cjB, 
                                void *jdataB, 
                                N_Vector resvecB, DenseMat JacB, 
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
                               realtype t, 
                               N_Vector yy, N_Vector yp,
                               N_Vector yyB, N_Vector ypB, realtype cjB, 
                               void *jdataB,
                               N_Vector resvecB, BandMat JacB, 
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

  typedef int (*IDASpgmrPrecSetupFnB)(realtype t, 
                                      N_Vector yy, N_Vector yp,
                                      N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                                      realtype cjB, void *pdataB,
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

  typedef int (*IDASpgmrPrecSolveFnB)(realtype t, 
                                      N_Vector yy, N_Vector yp,
                                      N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                                      N_Vector rvecB, N_Vector zvecB,
                                      realtype cjB, realtype deltaB,
                                      void *pdataB, N_Vector tmpB);

  /******************************************************************
   *                                                                *
   * Type : IDASpgmrJacTimesVecFnB                                   *
   *----------------------------------------------------------------*
   * A Jacobian times vector function jtimesB for the backward      *
   * integration must have the prototype given below.               *
   *                                                                *
   ******************************************************************/

  typedef int (*IDASpgmrJacTimesVecFnB)(N_Vector vB, N_Vector JvB, realtype t,
                                        N_Vector yy, N_Vector yp,
                                        N_Vector yyB, N_Vector ypB, N_Vector rrB,
                                        realtype cjB, void *jdataB, 
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
   ******************************************************************/

  int IDASolveF(void *idaadj_mem, realtype tout, realtype *tret,
                N_Vector yret, N_Vector ypret, int itask, int *ncheckPtr);

  /* IDASolveF return values */
#define IDASOLVEF_MEM_FAIL -20
  /* or any IDASolve return value */ 

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
  
  int IDASetRdataB(void *idaadj_mem, void *rdataB);
  int IDASetErrFileB(void *idaadj_mem, FILE *errfpB);
  int IDASetMaxOrdB(void *idaadj_mem, int maxordB);
  int IDASetMaxNumStepsB(void *idaadj_mem, long int mxstepsB);
  int IDASetInitStepB(void *idaadj_mem, realtype hinB);
  int IDASetMaxStepB(void *idaadj_mem, realtype hmaxB);
  int IDASetSuppressAlgB(void *idaadj_mem, booleantype suppressalgB);
  int IDASetIdB(void *idaadj_mem, N_Vector idB);
  int IDASetConstraintsB(void *idaadj_mem, N_Vector constraintsB);

  int IDAMallocB(void *idaadj_mem, ResFnB resB,
                 realtype tB0, N_Vector yyB0, N_Vector ypB0, 
                 int itolB, realtype *reltolB, void *abstolB,
                 NV_Spec nvspecB);
  
  int IDAReInitB(void *idaadj_mem, ResFnB resB,
                 realtype tB0, N_Vector yyB0, N_Vector ypB0,
                 int itolB, realtype *reltolB, void *abstolB);

  /* IDAMallocB return values */
  /* SUCCESS=0 */
#define IDABM_NO_MEM   -101
#define IDABM_MEM_FAIL -102
#define IDABM_BAD_TB0  -103

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
  int IDASetQuadRdataB(void *idaadj_mem, void *rdataQB);
  int IDASetQuadTolerancesB(void *idaadj_mem, int itolQB, 
                            realtype *reltolQB, void *abstolQB);

  int IDAQuadMallocB(void *idaadj_mem, QuadRhsFnB rhsQB, NV_Spec nvspecQB);

  /* IDAQuadMallocB return values */
  /* SUCCESS=0                    */
  /* IDABM_NO_MEM                 */
#define IDABM_ILL_INPUT   -104

  /******************************************************************
   *                                                                *
   * Function : IDAQuadReInitB                                    *
   *----------------------------------------------------------------*
   * IDAQuadReInitB re-initializaes memory for quadrature         *
   * integration during the backward phase                          *
   *                                                                *
   ******************************************************************/

  int IDAQuadReInitB(void *idaadj_mem, QuadRhsFnB rhsQB);

  /******************************************************************
   *                                                                *
   * Function : IDADenseB, IDADenseSet*B                              *
   *----------------------------------------------------------------*
   * IDADenseB links the main IDAS integrator with the IDADENSE      *
   * linear solver for the backward integration.                    *
   *                                                                *
   ******************************************************************/
  
  int IDADenseB(void *idaadj_mem, long int NeqB);
  
  int IDADenseSetJacFnB(void *idaadj_mem, IDADenseJacFnB djacB);
  int IDADenseSetJacDataB(void *idaadj_mem, void *jdataB);

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

  int IDABandSetJacFnB(void *idaadj_mem, IDABandJacFnB bjacB);
  int IDABandSetJacDataB(void *idaadj_mem, void *jdataB);

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
  int IDASpgmrSetPrecSetupFnB(void *idaadj_mem, IDASpgmrPrecSetupFnB psetB);
  int IDASpgmrSetPrecSolveFnB(void *idaadj_mem, IDASpgmrPrecSolveFnB psolveB);
  int IDASpgmrSetJacTimesVecFnB(void *idaadj_mem, IDASpgmrJacTimesVecFnB jtimesB);
  int IDASpgmrSetPrecDataB(void *idaadj_mem, void *pdataB);
  int IDASpgmrSetJacDataB(void *idaadj_mem, void *jdataB);

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
    N_Vector ck_phi[MXORDP1];

    /* Nordsieck History Array for quadratures */
    N_Vector ck_phiQ[MXORDP1];

    /* Do we need to carry quadratures? */
    booleantype ck_quad;

    realtype ck_psi[MXORDP1];
    realtype ck_alpha[MXORDP1];
    realtype ck_beta[MXORDP1];
    realtype ck_sigma[MXORDP1];
    realtype ck_gamma[MXORDP1];
    
    /* Step data */
    long int     ck_nst;
    long int     ck_ns;
    int          ck_kk;
    int          ck_kused;
    int          ck_knew;
    int          ck_phase;
    
    realtype     ck_hh;
    realtype     ck_hused;
    realtype     ck_rr;
    realtype     ck_cj;
    realtype     ck_cjlast;
    realtype     ck_cjold;
    realtype     ck_cjratio;
    
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
   * Types : struct IDAadjMemRec, IDAadjMem                         *
   *----------------------------------------------------------------*
   * The type IDAadjMem is type pointer to struct IDAadjMemRec.     *
   * This structure contins fields to store all information         *
   * necessary for adjoint sensitivity analysis.                    *
   *                                                                *
   ******************************************************************/

  typedef struct IDAadjMemRec {
    
    /* IDAS memory for forward runs */
    struct IDAMemRec *IDA_mem;
    
    /* IDAS memory for backward run */
    struct IDAMemRec *IDAB_mem;
    
    /* Storage for check point information */
    struct CkpntMemRec *ck_mem;
    
    /* Storage for data from forward runs */
    struct DtpntMemRec **dt_mem;
    
    /* Residual function (resB) for backward run */
    ResFnB ia_resB;
    
    /* Right hand side quadrature function (fQB) for backward run */
    QuadRhsFnB ia_rhsQB;
    
    /* Dense Jacobian function (djacB) for backward run */
    IDADenseJacFnB ia_djacB;
    
    /* Banded Jacobian function (bjacB) for backward run */
    IDABandJacFnB ia_bjacB;
    
    /* Preconditioner routines (precondB and psolveB) for backward run */
    IDASpgmrPrecSetupFnB ia_psetB;
    IDASpgmrPrecSolveFnB ia_psolveB;
    
    /* Jac times vec routine (jtimesB) for backward run */
    IDASpgmrJacTimesVecFnB ia_jtimesB;
    
    /* User rdataB */
    void *ia_rdataB;
    
    /* User rdataQB */
    void *ia_rdataQB;
    
    /* User jdataB */
    void *ia_jdataB;
    
    /* User pdataB */
    void *ia_pdataB;
    
    /* Unit roundoff */
    realtype ia_uround;
    
    /* Integration interval */
    realtype ia_tinitial, ia_tfinal;
    
    /* Time at which to extract quadratures */
    realtype ia_t_for_quad;
    
    /* Number of check points */
    int ia_nckpnts;
    
    /* Number of steps between 2 check points */
    long int ia_nsteps;
    
    /* Flag to indicate that data in dt_mem is new */
    booleantype ia_newData;
    
    /* address of the check point structure for which data is available */
    struct CkpntMemRec *ia_ckpntData;
    
    /* Actual number of data points saved in current dt_mem */
    /* Commonly, np = nsteps+1                              */
    long int ia_np;
    
    /* Temporary space used by the Hermite interpolation */
    realtype ia_dt;
    N_Vector ia_Y0, ia_Y1;
    N_Vector ia_ytmp, ia_yptmp;
    
  } *IDAadjMem;
  
#endif
  
#ifdef __cplusplus
}
#endif
