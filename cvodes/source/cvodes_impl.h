/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2004-10-26 20:13:21 $
 * ----------------------------------------------------------------- 
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, Radu Serban
 *              and Dan Shumaker @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvodes/LICENSE
 * -----------------------------------------------------------------
 * Implementation header file for the main CVODES integrator.
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvodes_impl_h
#define _cvodes_impl_h

#include <stdio.h>

#include "cvodes.h"

#include "nvector.h"
#include "sundialstypes.h"

/* 
 * =================================================================
 *   M A I N    I N T E G R A T O R    M E M O R Y    B L O C K
 * =================================================================
 */

/* Basic CVODES constants */

#define ADAMS_Q_MAX 12      /* max value of q for lmm == ADAMS    */
#define BDF_Q_MAX    5      /* max value of q for lmm == BDF      */
#define Q_MAX  ADAMS_Q_MAX  /* max value of q for either lmm      */
#define L_MAX  (Q_MAX+1)    /* max value of L for either lmm      */
#define NUM_TESTS    5      /* number of error test quantities    */

/*
 * ifS:   Type of the function returning the sensitivity right    
 *        hand side. ifS can be either CV_ALLSENS if the function    
 *        (of type CVSensRhsFn) returns right hand sides for all    
 *        sensitivity systems at once, or CV_ONESENS if the function 
 *        (of type SensRhs1Fn) returns the right hand side of one 
 *        sensitivity system at a time.                           
 *                                                                
 */
#define CV_ONESENS 1
#define CV_ALLSENS 2

/*
 * -----------------------------------------------------------------
 * Types: struct CVodeMemRec, CVodeMem
 * -----------------------------------------------------------------
 * The type CVodeMem is type pointer to struct CVodeMemRec.
 * This structure contains fields to keep track of problem state.
 * -----------------------------------------------------------------
 */
  
typedef struct CVodeMemRec {
    
  realtype cv_uround;      /* machine unit roundoff                        */   

  /*---------------------------- 
    Problem Specification Data 
  ----------------------------*/

  CVRhsFn cv_f;            /* y' = f(t,y(t))                               */
  void *cv_f_data;         /* user pointer passed to f                     */
  int cv_lmm;              /* lmm = ADAMS or BDF                           */
  int cv_iter;             /* iter = FUNCTIONAL or NEWTON                  */
  int cv_itol;             /* itol = SS or SV                              */
  realtype *cv_reltol;     /* ptr to relative tolerance                    */
  void *cv_abstol;         /* ptr to absolute tolerance                    */

  /*--------------------------
    Quadrature Related Data 
  --------------------------*/

  booleantype cv_quad;     /* TRUE if integrating quadratures              */
  CVQuadRhsFn cv_fQ;
  int cv_itolQ;
  realtype *cv_reltolQ;    /* ptr to relative tolerance for quad           */
  void *cv_abstolQ;        /* ptr to absolute tolerance for quad           */
  booleantype cv_errconQ;
  void *cv_fQ_data;        /* user pointer passed to fQ                    */

  /*--------------------------
    Sensitivity Related Data 
  --------------------------*/

  booleantype cv_sensi;    /* TRUE if computing sensitivities              */
  int cv_Ns;               /* Number of sensitivities                      */
  CVSensRhsFn cv_fS;       /* fS = (df/dy)*yS + (df/dp)                    */
  CVSensRhs1Fn cv_fS1;     /* fS1 = (df/dy)*yS_i + (df/dp)                 */
  booleantype cv_fSDQ;
  int cv_ifS;              /* ifS = ALLSENS or ONESENS                     */
  int cv_ism;              /* ism = SIMULTANEOUS or STAGGERED              */
  realtype *cv_p;          /* parameters in f(t,y,p)                       */
  realtype *cv_pbar;       /* scale factors for parameters                 */
  int *cv_plist;           /* list of sensitivities                        */
  realtype cv_rhomax;      /* cut-off value for centered/forward finite
                              differences                                  */
  booleantype cv_errconS;  /* TRUE if sensitivities are in err. control    */
  void *cv_fS_data;        /* user pointer passed to fS                    */

  int cv_itolS;
  realtype *cv_reltolS;          /* ptr to relative tolerance for sensi    */
  void *cv_abstolS;              /* ptr to absolute tolerance for sensi    */
  booleantype cv_testSensTol;    /* flag to indicate if sensi. tolerances
                                    must be checked now                    */
  booleantype cv_setSensTol;     /* flag to indicate if sensi. tolerances 
                                    must be set now                        */
  booleantype cv_atolSallocated; /* TRUE if CVODES has allocated space for
                                    sensitivity absolute tolerances        */

  /*-------------------------
    Nordsieck History Array 
  -------------------------*/

  N_Vector cv_zn[L_MAX];   /* Nordsieck array, of size N x (q+1).
                              zn[j] is a vector of length N (j=0,...,q)
                              zn[j] = [1/factorial(j)] * h^j * 
                              (jth derivative of the interpolating 
                              polynomial                                   */

  /*---------------------
    Vectors of length N 
  ---------------------*/

  N_Vector cv_ewt;         /* error weight vector                          */
  N_Vector cv_y;           /* y is used as temporary storage by the solver.
                              The memory is provided by the user to CVode 
                              where the vector is named yout.              */
  N_Vector cv_acor;        /* In the context of the solution of the
                              nonlinear equation, acor = y_n(m) - y_n(0).
                              On return, this vector is scaled to give
                              the estimated local error in y.              */
  N_Vector cv_tempv;       /* temporary storage vector                     */
  N_Vector cv_ftemp;       /* temporary storage vector                     */

  /*-----------------------------
    Quadrature Related Vectors 
  -----------------------------*/

  N_Vector cv_znQ[L_MAX];  /* Nordsieck arrays for sensitivities           */
  N_Vector cv_ewtQ;        /* error weight vector for quadratures          */
  N_Vector cv_yQ;          /* Unlike y, yQ is not allocated by the user    */
  N_Vector cv_acorQ;       /* acorQ = yQ_n(m) - yQ_n(0)                    */
  N_Vector cv_tempvQ;      /* temporary storage vector (~ tempv)           */

  /*-----------------------------
    Sensitivity Related Vectors 
  -----------------------------*/

  N_Vector *cv_znS[L_MAX]; /* Nordsieck arrays for sensitivities           */
  N_Vector *cv_ewtS;       /* error weight vectors for sensitivities       */
  N_Vector *cv_yS;         /* yS=yS0 (allocated by the user)               */
  N_Vector *cv_acorS;      /* acorS = yS_n(m) - yS_n(0)                    */
  N_Vector *cv_tempvS;     /* temporary storage vector (~ tempv)           */
  N_Vector *cv_ftempS;     /* temporary storage vector (~ ftemp)           */

  /*-------------------------------------------------
    Does CVodeSensMalloc allocate additional space?
  -------------------------------------------------*/  

  booleantype cv_stgr1alloc;     /* Are ncfS1, ncfnS1, and nniS1 allocated 
                                    by CVODES?                             */

  /*-----------------
    Tstop information
  -------------------*/
  booleantype cv_tstopset;
  realtype cv_tstop;

  /*-----------
    Step Data 
  -----------*/

  int cv_q;                    /* current order                            */
  int cv_qprime;               /* order to be used on the next step        */ 
                               /* = q-1, q, or q+1                         */
  int cv_next_q;               /* order to be used on the next step        */
  int cv_qwait;                /* number of internal steps to wait before  */
                               /* considering a change in q                */
  int cv_L;                    /* L = q + 1                                */

  realtype cv_hin;
  realtype cv_h;               /* current step size                        */
  realtype cv_hprime;          /* step size to be used on the next step    */ 
  realtype cv_next_h;          /* step size to be used on the next step    */ 
  realtype cv_eta;             /* eta = hprime / h                         */
  realtype cv_hscale;          /* value of h used in zn                    */
  realtype cv_tn;              /* current internal value of t              */

  realtype cv_tau[L_MAX+1];    /* array of previous q+1 successful step
                                  sizes indexed from 1 to q+1              */
  realtype cv_tq[NUM_TESTS+1]; /* array of test quantities indexed from
                                  1 to NUM_TESTS(=5)                       */
  realtype cv_l[L_MAX];        /* coefficients of l(x) (degree q poly)     */

  realtype cv_rl1;             /* 1 / l[1]                                 */
  realtype cv_gamma;           /* gamma = h * rl1                          */
  realtype cv_gammap;          /* gamma at the last setup call             */
  realtype cv_gamrat;          /* gamma / gammap                           */

  realtype cv_crate;           /* est. corrector conv. rate in Nls         */
  realtype cv_crateS;          /* est. corrector conv. rate in NlsStgr     */
  realtype cv_acnrm;           /* | acor |                                 */
  realtype cv_acnrmS;          /* | acorS |                                */
  realtype cv_acnrmQ;          /* | acorQ |                                */
  realtype cv_nlscoef;         /* coeficient in nonlinear convergence test */
  int  cv_mnewt;               /* Newton iteration counter                 */
  int  *cv_ncfS1;              /* Array of Ns local counters for conv.  
                                  failures (used in CVStep for STAGGERED1) */

  /*--------
    Limits 
  --------*/

  int cv_qmax;             /* q <= qmax                                    */
  long int cv_mxstep;      /* maximum number of internal steps for one 
                              user call                                    */
  int cv_maxcor;           /* maximum number of corrector iterations for 
                              the solution of the nonlinear equation       */
  int cv_maxcorS;
  int cv_mxhnil;           /* maximum number of warning messages issued to 
                              the user that t + h == t for the next 
                              internal step                                */
  int cv_maxnef;           /* maximum number of error test failures        */
  int cv_maxncf;           /* maximum number of nonlinear conv. failures   */

  realtype cv_hmin;        /* |h| >= hmin                                  */
  realtype cv_hmax_inv;    /* |h| <= 1/hmax_inv                            */
  realtype cv_etamax;      /* eta <= etamax                                */

  /*----------
    Counters 
  ----------*/

  long int cv_nst;         /* number of internal steps taken               */
  long int cv_nfe;         /* number of f calls                            */
  long int cv_nfSe;        /* number of fS calls                           */
  long int cv_nfQe;        /* number of fQ calls                           */
  long int cv_nfeS;        /* number of f calls from sensi DQ              */

  long int cv_ncfn;        /* number of corrector convergence failures     */
  long int cv_ncfnS;       /* number of total sensi. corr. conv. failures  */
  long int *cv_ncfnS1;     /* number of sensi. corrector conv. failures    */

  long int cv_nni;         /* number of nonlinear iterations performed     */
  long int cv_nniS;        /* number of total sensi. nonlinear iterations  */
  long int *cv_nniS1;      /* number of sensi. nonlinear iterations        */

  long int cv_netf;        /* number of error test failures                */
  long int cv_netfS;       /* number of sensi. error test failures         */
  long int cv_netfQ;       /* number of quadr. error test failures         */

  long int cv_nsetups;     /* number of setup calls                        */
  long int cv_nsetupsS;    /* number of setup calls due to sensitivities   */

  int cv_nhnil;            /* number of messages issued to the user that
                                   t + h == t for the next iternal step    */

  /*------------------------------- 
    Space requirements for CVODES 
  -------------------------------*/

  long int cv_lrw1;        /* no. of realtype words in 1 N_Vector y        */ 
  long int cv_liw1;        /* no. of integer words in 1 N_Vector y         */ 
  long int cv_lrw1Q;       /* no. of realtype words in 1 N_Vector yQ       */ 
  long int cv_liw1Q;       /* no. of integer words in 1 N_Vector yQ        */ 
  long int cv_lrw;         /* no. of realtype words in CVODES work vectors */
  long int cv_liw;         /* no. of integer words in CVODES work vectors  */

  /*------------------
    Step size ratios
  ------------------*/

  realtype cv_etaqm1;      /* ratio of new to old h for order q-1          */
  realtype cv_etaq;        /* ratio of new to old h for order q            */
  realtype cv_etaqp1;      /* ratio of new to old h for order q+1          */

  /*--------------------
    Linear Solver Data 
  --------------------*/

  /* Linear Solver functions to be called */

  int (*cv_linit)(struct CVodeMemRec *cv_mem);

  int (*cv_lsetup)(struct CVodeMemRec *cv_mem, int convfail, 
                   N_Vector ypred, N_Vector fpred, booleantype *jcurPtr, 
                   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3); 

  int (*cv_lsolve)(struct CVodeMemRec *cv_mem, N_Vector b, N_Vector weight,
                   N_Vector ycur, N_Vector fcur);

  void (*cv_lfree)(struct CVodeMemRec *cv_mem);

  /* Linear Solver specific memory */

  void *cv_lmem;           

  /* Flag to request a call to the setup routine */

  booleantype cv_forceSetup;

  /*--------------
    Saved Values
  --------------*/

  int cv_qu;               /* last successful q value used                 */
  long int cv_nstlp;       /* step number of last setup call               */
  realtype cv_h0u;         /* actual initial stepsize                      */
  realtype cv_hu;          /* last successful h value used                 */
  realtype cv_saved_tq5;   /* saved value of tq[5]                         */
  booleantype cv_jcur;     /* Is the Jacobian info used by
                              linear solver current?                       */
  realtype cv_tolsf;       /* tolerance scale factor                       */
  booleantype cv_setupNonNull;  /* Does setup do something?                */

  /*-------------------------------------------------------------------
    Flags turned ON by CVodeMalloc, CVodeSensMalloc, and CVodeQuadMalloc 
    and read by CVodeReInit, CVodeSensReInit, and CVodeQuadReInit
    ------------------------------------------------------------------*/

  booleantype cv_MallocDone;
  booleantype cv_sensMallocDone;
  booleantype cv_quadMallocDone;

  /*------------
    Error File 
  ------------*/

  FILE *cv_errfp;             /* CVODE error messages are sent to errfp    */

  /*-------------------------
    Stability Limit Detection
  ---------------------------*/

  booleantype cv_sldeton;     /* Is Stability Limit Detection on?          */
  realtype cv_ssdat[6][4];    /* scaled data array for STALD               */
  int cv_nscon;               /* counter for STALD method                  */
  long int cv_nor;            /* counter for number of order reductions    */

  /*----------------
    Rootfinding Data
  ------------------*/

  CVRootFn cv_gfun;     /* Function g for roots sought                     */
  int cv_nrtfn;         /* number of components of g                       */
  void *cv_g_data;      /* pointer to user data for g                      */
  int *cv_iroots;       /* int array for root information                  */
  realtype cv_tlo;      /* nearest endpoint of interval in root search     */
  realtype cv_thi;      /* farthest endpoint of interval in root search    */
  realtype cv_troot;    /* approximate root location                       */
  realtype *cv_glo;     /* saved array of g values at t = tlo              */
  realtype *cv_ghi;     /* saved array of g values at t = thi              */
  realtype *cv_groot;   /* array of g values at t = troot                  */
  realtype cv_tretlast; /* last value of t returned                        */
  realtype cv_toutc;    /* copy of tout (if NORMAL mode)                   */
  realtype cv_ttol;     /* tolerance on root location troot                */
  int cv_taskc;         /* copy of parameter task                          */
  int cv_irfnd;         /* flag showing whether last step had a root       */
  int cv_nge;           /* counter for g evaluations                       */

  /*-------------------------
    Complex step memory block
    -------------------------*/

  void *cv_csmem;

} *CVodeMem;


#endif

#ifdef __cplusplus
}
#endif
