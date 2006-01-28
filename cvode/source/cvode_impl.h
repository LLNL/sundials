/*
 * -----------------------------------------------------------------
 * $Revision: 1.17 $
 * $Date: 2006-01-28 00:47:27 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban
 *                and Dan Shumaker @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * Implementation header file for the main CVODE integrator.
 * -----------------------------------------------------------------
 */

#ifndef _CVODE_IMPL_H
#define _CVODE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <stdarg.h>

#include "cvode.h"

/*
 * =================================================================
 *   M A I N    I N T E G R A T O R    M E M O R Y    B L O C K
 * =================================================================
 */

/* Basic CVODE constants */

#define ADAMS_Q_MAX 12     /* max value of q for lmm == ADAMS     */
#define BDF_Q_MAX    5     /* max value of q for lmm == BDF       */
#define Q_MAX  ADAMS_Q_MAX /* max value of q for either lmm       */
#define L_MAX  (Q_MAX+1)   /* max value of L for either lmm       */
#define NUM_TESTS    5     /* number of error test quantities     */

#define HMIN_DEFAULT     RCONST(0.0)    /* hmin default value     */
#define HMAX_INV_DEFAULT RCONST(0.0)    /* hmax_inv default value */
#define MXHNIL_DEFAULT   10             /* mxhnil default value   */
#define MXSTEP_DEFAULT   500            /* mxstep default value   */

/*
 * -----------------------------------------------------------------
 * Types : struct CVodeMemRec, CVodeMem
 * -----------------------------------------------------------------
 * The type CVodeMem is type pointer to struct CVodeMemRec.
 * This structure contains fields to keep track of problem state.
 * -----------------------------------------------------------------
 */

typedef struct CVodeMemRec {

  realtype cv_uround;    /* machine unit roundoff */

  /*-------------------------- 
    Problem Specification Data 
    --------------------------*/

  CVRhsFn cv_f;        /* y' = f(t,y(t))                    */
  void *cv_f_data;     /* user pointer passed to f          */
  int cv_lmm;          /* lmm = CV_ADAMS or CV_BDF          */
  int cv_iter;         /* iter = CV_FUNCTIONAL or CV_NEWTON */
  int cv_itol;         /* itol = CV_SS or CV_SV             */

  realtype cv_reltol;  /* relative tolerance                */
  realtype cv_Sabstol; /* scalar absolute tolerance         */
  N_Vector cv_Vabstol; /* vector absolute tolerance         */
  CVEwtFn cv_efun;     /* function to set ewt               */
  void *cv_e_data;     /* user pointer passed to efun       */

  /*-----------------------
    Nordsieck History Array 
    -----------------------*/

  N_Vector cv_zn[L_MAX];  /* Nordsieck array, of size N x (q+1).         */
                          /* zn[j] is a vector of length N (j=0,...,q)   */
                          /* zn[j] = [1/factorial(j)] * h^j * (jth       */ 
                          /* derivative of the interpolating polynomial  */

  /*--------------------------
    other vectors of length N 
    -------------------------*/

  N_Vector cv_ewt;     /* error weight vector                          */
  N_Vector cv_y;       /* y is used as temporary storage by the solver */
                       /* The memory is provided by the user to CVode  */
                       /* where the vector is named yout.              */
  N_Vector cv_acor;    /* In the context of the solution of the        */
                       /* nonlinear equation, acor = y_n(m) - y_n(0).  */
                       /* On return, this vector is scaled to give     */
                       /* the estimated local error in y.              */
  N_Vector cv_tempv;   /* temporary storage vector                     */
  N_Vector cv_ftemp;   /* temporary storage vector                     */

  /*-----------------
    Tstop information
    -----------------*/
  booleantype cv_tstopset;
  booleantype cv_istop;
  realtype cv_tstop;

  /*---------
    Step Data 
    ---------*/  

  int cv_q;         /* current order                           */
  int cv_qprime;    /* order to be used on the next step       */ 
                    /* = q-1, q, or q+1                        */
  int cv_next_q;    /* order to be used on the next step       */
  int cv_qwait;     /* number of internal steps to wait before */
                    /* considering a change in q               */
  int cv_L;         /* L = q + 1                               */

  realtype cv_hin;    /* initial step size                      */
  realtype cv_h;      /* current step size                      */
  realtype cv_hprime; /* step size to be used on the next step  */ 
  realtype cv_next_h; /* step size to be used on the next step  */ 
  realtype cv_eta;    /* eta = hprime / h                       */
  realtype cv_hscale; /* value of h used in zn                  */
  realtype cv_tn;     /* current internal value of t            */
  realtype cv_tretlast; /* value of tret last returned by CVode */

  realtype cv_tau[L_MAX+1];    /* array of previous q+1 successful step     */
                               /* sizes indexed from 1 to q+1               */
  realtype cv_tq[NUM_TESTS+1]; /* array of test quantities indexed from     */
                               /* 1 to NUM_TESTS(=5)                        */
  realtype cv_l[L_MAX];        /* coefficients of l(x) (degree q poly)      */

  realtype cv_rl1;     /* the scalar 1/l[1]            */
  realtype cv_gamma;   /* gamma = h * rl1              */
  realtype cv_gammap;  /* gamma at the last setup call */
  realtype cv_gamrat;  /* gamma / gammap               */

  realtype cv_crate;   /* estimated corrector convergence rate     */
  realtype cv_acnrm;   /* | acor | wrms                            */
  realtype cv_nlscoef; /* coeficient in nonlinear convergence test */
  int  cv_mnewt;       /* Newton iteration counter                 */

  /*------
    Limits 
    ------*/

  int cv_qmax;        /* q <= qmax                                          */
  long int cv_mxstep; /* maximum number of internal steps for one user call */
  int cv_maxcor;      /* maximum number of corrector iterations for the     */
                      /* solution of the nonlinear equation                 */
  int cv_mxhnil;      /* maximum number of warning messages issued to the   */
                      /* user that t + h == t for the next internal step    */
  int cv_maxnef;      /* maximum number of error test failures              */
  int cv_maxncf;      /* maximum number of nonlinear convergence failures   */

  realtype cv_hmin;     /* |h| >= hmin       */
  realtype cv_hmax_inv; /* |h| <= 1/hmax_inv */
  realtype cv_etamax;   /* eta <= etamax     */

  /*--------
    Counters 
    --------*/

  long int cv_nst;              /* number of internal steps taken             */
  long int cv_nfe;              /* number of f calls                          */
  long int cv_ncfn;             /* number of corrector convergence failures   */
  long int cv_netf;             /* number of error test failures              */
  long int cv_nni;              /* number of Newton iterations performed      */
  long int cv_nsetups;          /* number of setup calls                      */
  int cv_nhnil;                 /* number of messages issued to the user that */
                                /* t + h == t for the next iternal step       */

  realtype cv_etaqm1;      /* ratio of new to old h for order q-1        */
  realtype cv_etaq;        /* ratio of new to old h for order q          */
  realtype cv_etaqp1;      /* ratio of new to old h for order q+1        */

  /*----------------------------
    Space requirements for CVODE 
    ----------------------------*/

  long int cv_lrw1;        /* no. of realtype words in 1 N_Vector         */ 
  long int cv_liw1;        /* no. of integer words in 1 N_Vector          */ 
  long int cv_lrw;         /* no. of realtype words in CVODE work vectors */
  long int cv_liw;         /* no. of integer words in CVODE work vectors  */

  /*------------------
    Linear Solver Data 
    ------------------*/

  /* Linear Solver functions to be called */

  int (*cv_linit)(struct CVodeMemRec *cv_mem);

  int (*cv_lsetup)(struct CVodeMemRec *cv_mem, int convfail, N_Vector ypred,
                   N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                   N_Vector vtemp2, N_Vector vtemp3); 

  int (*cv_lsolve)(struct CVodeMemRec *cv_mem, N_Vector b, N_Vector weight,
                   N_Vector ycur, N_Vector fcur);

  void (*cv_lfree)(struct CVodeMemRec *cv_mem);

  /* Linear Solver specific memory */

  void *cv_lmem;           

  /*------------
    Saved Values
    ------------*/

  int cv_qu;             /* last successful q value used                        */
  long int cv_nstlp;     /* step number of last setup call                      */
  realtype cv_h0u;       /* actual initial stepsize                             */
  realtype cv_hu;        /* last successful h value used                        */
  realtype cv_saved_tq5; /* saved value of tq[5]                                */
  booleantype cv_jcur;   /* Is the Jacobian info used by linear solver current? */
  realtype cv_tolsf;     /* tolerance scale factor                              */
  int cv_qmax_alloc;     /* value of qmax used when allocating memory           */
  int cv_indx_acor;      /* index of the zn vector in which acor is saved       */
  booleantype cv_setupNonNull; /* Does setup do something?                      */

  booleantype cv_VabstolMallocDone;
  booleantype cv_MallocDone;  

  /*-------------------------------------------
    Error handler function and error ouput file 
    -------------------------------------------*/

  CVErrHandlerFn cv_ehfun;    /* Error messages are handled by ehfun     */
  void *cv_eh_data;           /* user pointer passed to ehfun            */
  FILE *cv_errfp;             /* CVODE error messages are sent to errfp  */

  /*-------------------------
    Stability Limit Detection
    -------------------------*/

  booleantype cv_sldeton;     /* Is Stability Limit Detection on?         */
  realtype cv_ssdat[6][4];    /* scaled data array for STALD              */
  int cv_nscon;               /* counter for STALD method                 */
  long int cv_nor;            /* counter for number of order reductions   */

  /*----------------
    Rootfinding Data
    ----------------*/

  CVRootFn cv_gfun;     /* Function g for roots sought                     */
  int cv_nrtfn;         /* number of components of g                       */
  void *cv_g_data;      /* pointer to user data for g                      */
  int *cv_iroots;       /* int array for root information                  */
  realtype cv_tlo;      /* nearest endpoint of interval in root search     */
  realtype cv_thi;      /* farthest endpoint of interval in root search    */
  realtype cv_trout;    /* t value returned by rootfinding routine         */
  realtype *cv_glo;     /* saved array of g values at t = tlo              */
  realtype *cv_ghi;     /* saved array of g values at t = thi              */
  realtype *cv_grout;   /* array of g values at t = trout                  */
  realtype cv_toutc;    /* copy of tout (if NORMAL mode)                   */
  realtype cv_ttol;     /* tolerance on root location                      */
  int cv_taskc;         /* copy of parameter task                          */
  int cv_irfnd;         /* flag showing whether last step had a root       */
  long int cv_nge;      /* counter for g evaluations                       */

} *CVodeMem;

/*
 * =================================================================
 *   C V O D E    I N T E R N A L   F U N C T I O N S
 * =================================================================
 */


/* Prototype of internal ewtSet function */

int CVEwtSet(N_Vector ycur, N_Vector weight, void *e_data);

/* High level error handler */

void CVProcessError(CVodeMem cv_mem, 
                    int error_code, const char *module, const char *fname, 
                    const char *msgfmt, ...);

/* Prototype of internal errHandler function */

void CVErrHandler(int error_code, const char *module, const char *function, 
                  char *msg, void *eh_data);

/*
 * =================================================================
 *   C V O D E    E R R O R    M E S S A G E S
 * =================================================================
 */

#if defined(SUNDIALS_EXTENDED_PRECISION)

#define MSG_TIME      "t = %Lg"
#define MSG_TIME_H    "t = %Lg and h = %Lg"
#define MSG_TIME_INT  "t = %Lg is not between tcur - hu = %Lg and tcur = %Lg.\n\n"
#define MSG_TIME_TOUT "tout = %Lg"

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define MSG_TIME      "t = %lg"
#define MSG_TIME_H    "t = %lg and h = %lg"
#define MSG_TIME_INT  "t = %lg is not between tcur - hu = %lg and tcur = %lg.\n\n"
#define MSG_TIME_TOUT "tout = %lg"

#else

#define MSG_TIME      "t = %g"
#define MSG_TIME_H    "t = %g and h = %g"
#define MSG_TIME_INT  "t = %g is not between tcur - hu = %g and tcur = %g.\n\n"
#define MSG_TIME_TOUT "tout = %g"

#endif

/* Initialization and I/O error messages */

#define MSGCV_NO_MEM "cvode_mem = NULL illegal."
#define MSGCV_CVMEM_FAIL "Allocation of cvode_mem failed."
#define MSGCV_MEM_FAIL "A memory request failed."
#define MSGCV_BAD_LMM  "Illegal value for lmm. The legal values are CV_ADAMS and CV_BDF."
#define MSGCV_BAD_ITER  "Illegal value for iter. The legal values are CV_FUNCTIONAL and CV_NEWTON."
#define MSGCV_BAD_ITOL "Illegal value for itol. The legal values are CV_SS, CV_SV, and CV_WF."
#define MSGCV_NO_MALLOC "Attempt to call before CVodeMalloc."
#define MSGCV_NEG_MAXORD "maxord <= 0 illegal."
#define MSGCV_BAD_MAXORD  "Illegal attempt to increase maximum method order."
#define MSGCV_NEG_MXSTEPS "mxsteps < 0 illegal."
#define MSGCV_SET_SLDET  "Attempt to use stability limit detection with the CV_ADAMS method illegal."
#define MSGCV_NEG_HMIN "hmin < 0 illegal."
#define MSGCV_NEG_HMAX "hmax < 0 illegal."
#define MSGCV_BAD_HMIN_HMAX "Inconsistent step size limits: hmin > hmax."
#define MSGCV_BAD_RELTOL "reltol < 0 illegal."
#define MSGCV_BAD_ABSTOL "abstol has negative component(s) (illegal)."
#define MSGCV_NULL_ABSTOL "abstol = NULL illegal."
#define MSGCV_NULL_Y0 "y0 = NULL illegal."
#define MSGCV_NULL_F "f = NULL illegal."
#define MSGCV_NULL_G "g = NULL illegal."
#define MSGCV_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGCV_BAD_K "Illegal value for k."
#define MSGCV_NULL_DKY "dky = NULL illegal."
#define MSGCV_BAD_T "Illegal value for t." MSG_TIME_INT

/* CVode Error Messages */

#define MSGCV_LSOLVE_NULL "The linear solver's solve routine is NULL."
#define MSGCV_YOUT_NULL "yout = NULL illegal."
#define MSGCV_TRET_NULL "tret = NULL illegal."
#define MSGCV_BAD_EWT "Initial ewt has component(s) equal to zero (illegal)."
#define MSGCV_EWT_NOW_BAD "At " MSG_TIME ", a component of ewt has become <= 0."
#define MSGCV_BAD_ITASK "Illegal value for itask."
#define MSGCV_BAD_H0 "h0 and tout - t0 inconsistent."
#define MSGCV_BAD_INIT_ROOT "Root found at and very near initial t."
#define MSGCV_BAD_TOUT "Trouble interpolating at " MSG_TIME_TOUT ". tout too far back in direction of integration"
#define MSGCV_NO_EFUN "itol = CV_WF but no EwtSet function was provided."
#define MSGCV_NO_TSTOP "itask = CV_NORMAL_TSTOP or itask = CV_ONE_STEP_TSTOP but tstop was not set."
#define MSGCV_EWT_FAIL "The user-provide EwtSet function failed."
#define MSGCV_EWT_NOW_FAIL "At " MSG_TIME ", the user-provide EwtSet function failed."
#define MSGCV_LINIT_FAIL "The linear solver's init routine failed."
#define MSGCV_HNIL_DONE "The above warning has been issued mxhnil times and will not be issued again for this problem."
#define MSGCV_TOO_CLOSE "tout too close to t0 to start integration."
#define MSGCV_MAX_STEPS "At " MSG_TIME ", mxstep steps taken before reaching tout."
#define MSGCV_TOO_MUCH_ACC "At " MSG_TIME ", too much accuracy requested."
#define MSGCV_HNIL "Internal " MSG_TIME_H " are such that t + h = t on the next step. The solver will continue anyway."
#define MSGCV_ERR_FAILS "At " MSG_TIME_H ", the error test failed repeatedly or with |h| = hmin."
#define MSGCV_CONV_FAILS "At " MSG_TIME_H ", the corrector convergence test failed repeatedly or with |h| = hmin."
#define MSGCV_SETUP_FAILED "At " MSG_TIME ", the setup routine failed in an unrecoverable manner."
#define MSGCV_SOLVE_FAILED "At " MSG_TIME ", the solve routine failed in an unrecoverable manner."
#define MSGCV_CLOSE_ROOTS "Root found at and very near " MSG_TIME "."
#define MSGCV_BAD_TSTOP "tstop is behind current " MSG_TIME "in the direction of integration."


#ifdef __cplusplus
}
#endif

#endif
