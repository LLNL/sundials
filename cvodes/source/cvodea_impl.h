/*
 * -----------------------------------------------------------------
 * $Revision: 1.18 $
 * $Date: 2006-02-02 00:32:21 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
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

#include "sundials_nvector.h"
#include "sundials_types.h"

  /*
   * -----------------------------------------------------------------
   * Forward references for pointers to various structures 
   * -----------------------------------------------------------------
   */

  typedef struct CVadjMemRec *CVadjMem;
  typedef struct CkpntMemRec *CkpntMem;
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

  struct CkpntMemRec {

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
    
  };
  
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

  /* Data for polynomial interpolation */
  typedef struct PolynomialDataMemRec {
    N_Vector y;
    int order;
  } *PolynomialDataMem;

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

    /* User f_dataB */
    void *ca_f_dataB;
    
    /* User fQ_dataB */
    void *ca_fQ_dataB;
    
    /* Memory block for a linear solver's interface to CVODEA */
    void *ca_lmemB;

    /* Function to free any memory allocated by the linear solver */
    void (*ca_lfreeB)(CVadjMem ca_mem);

    /* Memory block for a preconditioner's module interface to CVODEA */ 
    void *ca_pmemB;
    
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
    
    /* Workspace used by the Hermite interpolation */
    N_Vector ca_Y0, ca_Y1;    /* pointers to zn[0] and zn[1] */

    /* Workspace for polynomial interpolation */
    N_Vector ca_Y[L_MAX];     /* pointers to zn[i] */
    realtype ca_T[L_MAX];

    /* Workspace for wrapper functions */
    N_Vector ca_ytmp;
    
  };
  
  /* Error Messages */
  
#define MSGAM_NULL_CVMEM   "cvode_mem = NULL illegal."
#define MSGAM_NULL_CAMEM   "cvadj_mem = NULL illegal."
#define MSGAM_BAD_STEPS    "Steps nonpositive illegal."
#define MSGAM_MEM_FAIL     "A memory request failed."
#define MSGAM_BAD_INTERP   "Illegal value for interp."
#define MSGAM_BAD_ITASKB   "Illegal value for itaskB. Legal values are CV_NORMAL and CV_ONE_STEP."
#define MSGAM_BAD_TB0      "The initial time tB0 is outside the interval over which the forward problem was solved."
#define MSGAM_BAD_TBOUT    "The final time tBout is outside the interval over which the forward problem was solved."
#define MSGAM_BAD_T        "Bad t for interpolation. Abort!"
#define MSGAM_WRONG_INTERP "This function cannot be called for the specified interp type."

#ifdef __cplusplus
}
#endif

#endif
