/*
 * -----------------------------------------------------------------
 * $Revision: 1.9 $
 * $Date: 2006-06-15 15:39:20 $
 * ----------------------------------------------------------------- 
 * Programmers: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see the LICENSE file
 * -----------------------------------------------------------------
 * Implementation header file for the IDAA adjoint integrator.
 * -----------------------------------------------------------------
 */

#ifndef _IDAA_IMPL_H
#define _IDAA_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <stdio.h>

#include "idaa.h"
#include "idas_impl.h"

#include "sundials_nvector.h"
#include "sundials_types.h"


  /*
   * -----------------------------------------------------------------
   * Forward references for pointers to various structures 
   * -----------------------------------------------------------------
   */

  typedef struct IDAadjMemRec *IDAadjMem;
  typedef struct CkpntMemRec *CkpntMem;
  typedef struct DtpntMemRec *DtpntMem;

  /*
   * -----------------------------------------------------------------
   * Types for functions provided by an interpolation module
   * -----------------------------------------------------------------
   * IDAAGetYFn:   Function type for a function that returns the 
   *               interpolated forward solution.
   * IDAAStorePnt: Function type for a function that stores a new
   *               point in the structure d
   * -----------------------------------------------------------------
   */

  typedef int (*IDAAGetYFn)(IDAadjMem IDAADJ_mem, realtype t, N_Vector y, N_Vector yp);
  typedef int (*IDAAStorePntFn)(IDAMem IDA_mem, DtpntMem d);

  /*
   * -----------------------------------------------------------------
   * Types : struct CkpntMemRec, CkpntMem
   * -----------------------------------------------------------------
   * The type CkpntMem is type pointer to struct CkpntMemRec.
   * This structure contains fields to store all information at a
   * check point that is needed to 'hot' start IDAS.
   * -----------------------------------------------------------------
   */

  struct CkpntMemRec {

    /* Integration limits */
    realtype ck_t0;
    realtype ck_t1;
    
    /* Nordsieck History Array */
    N_Vector ck_phi[MXORDP1];

    /* Nordsieck History Array for quadratures */
    N_Vector ck_phiQ[MXORDP1];
    
    /* Do we need to carry quadratures? */
    booleantype ck_quadr;

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
   * Type : struct IDAadjMemRec
   * -----------------------------------------------------------------
   * The type IDAadjMem is type pointer to struct IDAadjMemRec.
   * This structure contins fields to store all information
   * necessary for adjoint sensitivity analysis.
   * -----------------------------------------------------------------
   */

  struct IDAadjMemRec {
    
    /* IDAS memory for forward runs */
    struct IDAMemRec *IDA_mem;
    
    /* IDAS memory for backward run */
    struct IDAMemRec *IDAB_mem;
    
    /* Storage for check point information */
    struct CkpntMemRec *ck_mem;

    /* Interpolation type */
    int ca_interpType;

    /* Storage for data from forward runs */
    struct DtpntMemRec **dt_mem;
    
    /* Functions set by the interpolation module */
    IDAAStorePntFn ia_storePnt; /* store a new interpolation point */
    IDAAGetYFn     ia_getY;     /* interpolate forward solution    */


    /* Residual function (resB) for backward run */
    IDAResFnB ia_resB;
    
    /* Right hand side quadrature function (fQB) for backward run */
    IDAQuadRhsFnB ia_rhsQB;

    /* User rdataB */
    void *ia_rdataB;
    
    /* User rdataQB */
    void *ia_rdataQB;

    /* Memory block for a linear solver's interface to IDAA */
    void *ia_lmemB;

    /* Function to free any memory allocated by the linear solver */
    void (*ia_lfreeB)(IDAadjMem IDAADJ_mem);

    /* Memory block for a preconditioner's module interface to IDAA */ 
    void *ia_pmemB;
    
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
    
    /* Workspace used by the Hermite interpolation */
    N_Vector ia_Y0, ia_Y1;       /* pointers to phi[0] and phi[1] */

    /* Workspace for polynomial interpolation */
    N_Vector ia_Y[MXORDP1];      /* pointers to phi[i] */
    realtype ia_T[MXORDP1];

    /* Workspace for wrapper functions */
    N_Vector ia_ytmp, ia_yptmp;
    
  };


/*
 *----------------------------------------------------------------
 * IDAA Error Messages
 *----------------------------------------------------------------
 */
#define MSGAM_NULL_IDAMEM  "ida_mem = NULL illegal."
#define MSGAM_NULL_CAMEM   "idaadj_mem = NULL illegal."
#define MSGAM_BAD_STEPS    "Steps nonpositive illegal."
#define MSGAM_MEM_FAIL     "A memory request failed."
#define MSGAM_BAD_INTERP   "Illegal value for interp."
#define MSGAM_BAD_ITASKB   "Illegal value for itaskB. Legal values are IDA_NORMAL and IDA_ONE_STEP."
#define MSGAM_BAD_TB0      "The initial time tB0 is outside the interval over which the forward problem was solved."
#define MSGAM_BAD_TBOUT    "The final time tBout is outside the interval over which the forward problem was solved."
#define MSGAM_BAD_T        "Bad t for interpolation."
#define MSGAM_WRONG_INTERP "This function cannot be called for the specified interp type."

#ifdef __cplusplus
}
#endif

#endif
