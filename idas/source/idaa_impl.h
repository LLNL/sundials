/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2004-11-05 23:55:11 $
 * ----------------------------------------------------------------- 
 * Programmers: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/idas/LICENSE
 * -----------------------------------------------------------------
 * Implementation header file for the IDAA adjoint integrator.
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif
  
#ifndef _idaa_impl_h
#define _idaa_impl_h
  
#include <stdio.h>

#include "idaa.h"

#include "idas_impl.h"
#include "idadense_impl.h"
#include "idaband_impl.h"
#include "idaspgmr_impl.h"
#include "idabbdpre_impl.h"

#include "nvector.h"
#include "sundialstypes.h"

  /*
   *                                                                
   * Types : struct CkpntMemRec, CkpntMem                           
   *----------------------------------------------------------------
   * The type CkpntMem is type pointer to struct CkpntMemRec.       
   * This structure contains fields to store all information at a   
   * check point that is needed to 'hot' start idas.  
   *                                                                
   */

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

  /*
   *                                                                
   * Types : struct DtpntMemRec, DtpntMem                           
   *----------------------------------------------------------------
   * The type DtpntMem is type pointer to struct DtpntMemRec.       
   * This structure contains fields to store all information at a   
   * data point that is needed to interpolate solution of forward   
   * simulations.                                                   
   *                                                                
   */

  typedef struct DtpntMemRec {
    
    /* time */
    realtype t;
    
    /* solution */
    N_Vector y;
    
    /* solution derivative */
    N_Vector yd;
    
  } *DtpntMem;
  
  /*
   *                                                                
   * Types : struct IDAadjMemRec, IDAadjMem                           
   *----------------------------------------------------------------
   * The type IDAadjMem is type pointer to struct IDAadjMemRec.       
   * This structure contins fields to store all information         
   * necessary for adjoint sensitivity analysis.                    
   *                                                                
   */

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
    IDAResFnB ia_resB;
    
    /* Right hand side quadrature function (fQB) for backward run */
    IDAQuadRhsFnB ia_rhsQB;
    
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

/*
 *----------------------------------------------------------------
 * IDAA Error Messages
 *----------------------------------------------------------------
 */

#define MSG_IDAAM_NO_MEM     "IDAAdjMalloc-- ida_mem = NULL illegal.\n\n"
#define MSG_IDAAM_BAD_STEPS  "IDAAdjMalloc-- steps non-positive illegal.\n\n"
#define MSG_IDAAM_MEM_FAIL   "IDAAdjMalloc-- a memory request failed.\n\n"

#define MSG_IDASOLVEF_MEM_FAIL "IDASolveF-- a memory request failed.\n\n"

#define MSG_IDABM_NO_MEM     "IDAMallocB/IDAReInitB-- idaadj_mem = NULL illegal.\n\n"
#define MSG_IDABM_BAD_TB0    "IDAMallocB/IDAReInitB-- tB0 out of range.\n\n"
#define MSG_IDABM_MEM_FAIL   "IDAMallocB/IDAReInitB-- a memory request failed.\n\n"

#define MSG_IDABQM_NO_MEM    "IDAQuadMallocB-- idaadj_mem = NULL illegal.\n\n"

#define MSG_IDASOLVEB_FWD    "IDASolveB-- an error occured during the forward phase.\n\n"
#define MSG_IDASOLVEB_BCK    "IDASolveB-- an error occured during the backward phase.\n\n"

  
#ifdef __cplusplus
}
#endif
