/*
 * -----------------------------------------------------------------
 * $Revision: 1.9 $
 * $Date: 2005-04-26 18:32:41 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Implementation header file for the CVODEA adjoint integrator.
 * -----------------------------------------------------------------
 */

#ifndef _CVODEA_IMPL_H
#define _CVODEA_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <stdio.h>

#include "cvodea.h"
#include "cvodes_impl.h"
#include "cvdense_impl.h"
#include "cvband_impl.h"
#include "cvspbcg_impl.h"
#include "cvspgmr_impl.h"
#include "cvbandpre_impl.h"
#include "cvbbdpre_impl.h"

#include "nvector.h"
#include "sundialstypes.h"

  /*
   * -----------------------------------------------------------------
   * Forward references for pointers to various structures 
   * -----------------------------------------------------------------
   */

  typedef struct CVadjMemRec *CVadjMem;
  typedef struct DtpntMemRec *DtpntMem;

  /*
   * -----------------------------------------------------------------
   * Types for functions provided by an interpolation module
   * -----------------------------------------------------------------
   * CVAGetYFn:   Function type for a function that returns the 
   *              interpolated forward solution.
   * CVAStorePnt: Function type for a function that stores a new
   *              point in the structure d
   * -----------------------------------------------------------------
   */

  typedef int (*CVAGetYFn)(CVadjMem ca_mem, realtype t, N_Vector y);
  typedef int (*CVAStorePntFn)(CVodeMem cv_mem, DtpntMem d);

  /*
   * -----------------------------------------------------------------
   * Types : struct CkpntMemRec, CkpntMem
   * -----------------------------------------------------------------
   * The type CkpntMem is type pointer to struct CkpntMemRec.
   * This structure contains fields to store all information at a
   * check point that is needed to 'hot' start cvodes.
   * -----------------------------------------------------------------
   */

  typedef struct CkpntMemRec {

    /* Integration limits */
    realtype ck_t0;
    realtype ck_t1;
    
    /* Nordsieck History Array */
    N_Vector ck_zn[L_MAX];
    
    /* Nordsieck History Array for quadratures */
    N_Vector ck_znQ[L_MAX];
    
    /* Do we need to carry quadratures? */
    booleantype ck_quadr;
    
    /* Was ck_zn[qmax] allocated?
       ck_zqm = 0    - no
       ck_zqm = qmax - yes      */
    int ck_zqm;
    
    /* Step data */
    long int ck_nst;
    realtype ck_tretlast;
    int      ck_q;
    int      ck_qprime;
    int      ck_qwait;
    int      ck_L;
    realtype ck_gammap;
    realtype ck_h;
    realtype ck_hprime;
    realtype ck_hscale;
    realtype ck_eta;
    realtype ck_etamax;
    realtype ck_tau[L_MAX+1];
    realtype ck_tq[NUM_TESTS+1];
    realtype ck_l[L_MAX];
    
    /* Saved values */
    realtype ck_saved_tq5;
    
    /* Pointer to next structure in list */
    struct CkpntMemRec *ck_next;
    
  } *CkpntMem;
  
  /*
   * -----------------------------------------------------------------
   * Type : struct DtpntMemRec
   * -----------------------------------------------------------------
   * This structure contains fields to store all information at a
   * data point that is needed to interpolate solution of forward
   * simulations. Its content field is interpType-dependent.
   * -----------------------------------------------------------------
   */
  
  struct DtpntMemRec {
    realtype t;    /* time */
    void *content; /* interpType-dependent content */
  };

  /* Data for cubic Hermite interpolation */
  typedef struct HermiteDataMemRec {
    N_Vector y;
    N_Vector yd;
  } *HermiteDataMem;

  /* Data for interpolation using the BDF interpolant */
  typedef struct BDFinterpDataMemRec {
  } *BDFinterpDataMem;

  /*
   * -----------------------------------------------------------------
   * Type : struct CVadjMemRec
   * -----------------------------------------------------------------
   * The type CVadjMem is type pointer to struct CVadjMemRec.
   * This structure contins fields to store all information
   * necessary for adjoint sensitivity analysis.
   * -----------------------------------------------------------------
   */
  
  struct CVadjMemRec {
    
    /* CVODE memory for forward runs */
    struct CVodeMemRec *cv_mem;
    
    /* CVODE memory for backward run */
    struct CVodeMemRec *cvb_mem;
    
    /* Storage for check point information */
    struct CkpntMemRec *ck_mem;

    /* Interpolation type */
    int ca_interpType;

    /* Storage for data from forward runs */
    struct DtpntMemRec **dt_mem;
    
    /* Functions set by the interpolation module */
    CVAStorePntFn ca_storePnt; /* store a new interpolation point */
    CVAGetYFn     ca_getY;     /* interpolate forward solution    */
    
    /* Right hand side function (fB) for backward run */
    CVRhsFnB ca_fB;
    
    /* Right hand side quadrature function (fQB) for backward run */
    CVQuadRhsFnB ca_fQB;
    
    /* Dense Jacobian function (djacB) for backward run */
    CVDenseJacFnB ca_djacB;
    
    /* Banded Jacobian function (bjacB) for backward run */
    CVBandJacFnB ca_bjacB;
    
    /* Jac times vec routine (jtimesB) for backward run */
    CVSpilsJacTimesVecFnB ca_jtimesB;
    
    /* Preconditioner routines (precondB and psolveB) for backward run */
    CVSpilsPrecSetupFnB ca_psetB;
    CVSpilsPrecSolveFnB ca_psolveB;
    
    /* BBD user functions (glocB and cfnB) for backward run */
    CVLocalFnB ca_glocB;
    CVCommFnB  ca_cfnB;
    
    /* User f_dataB */
    void *ca_f_dataB;
    
    /* User fQ_dataB */
    void *ca_fQ_dataB;
    
    /* User jac_dataB */
    void *ca_jac_dataB;
    
    /* User P_dataB */
    void *ca_P_dataB;
    
    /* BP prec data */
    void *ca_bp_dataB;
    
    /* BBD prec data */
    void *ca_bbd_dataB;
    
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
    
  };
  
  /* Error Messages */
  
#define _CVAM_          "CVadjMalloc-- "
#define MSGAM_NO_MEM    _CVAM_ "cvode_mem = NULL illegal.\n\n"
#define MSGAM_BAD_STEPS _CVAM_ "Steps nonpositive illegal.\n\n"
#define MSGAM_MEM_FAIL  _CVAM_ "A memory request failed.\n\n"
#define MSGAM_BAD_INTERP _CVAM_ "Illegal value for interp.\n\n"

#ifdef __cplusplus
}
#endif

#endif
