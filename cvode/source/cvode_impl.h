/*
 * -----------------------------------------------------------------
 * $Revision: 1.9 $
 * $Date: 2005-04-05 01:59:46 $
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

#include <stdio.h>

#include "cvode.h"
#include "nvector.h"
#include "sundialstypes.h"

/* Prototype of internal ewtSet function */

int CVEwtSet(N_Vector ycur, N_Vector weight, void *e_data);


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

  realtype cv_hin;
  realtype cv_h;      /* current step size                     */
  realtype cv_hprime; /* step size to be used on the next step */ 
  realtype cv_next_h; /* step size to be used on the next step */ 
  realtype cv_eta;    /* eta = hprime / h                      */
  realtype cv_hscale; /* value of h used in zn                 */
  realtype cv_tn;     /* current internal value of t           */

  realtype cv_tau[L_MAX+1];    /* array of previous q+1 successful step     */
                               /* sizes indexed from 1 to q+1               */
  realtype cv_tq[NUM_TESTS+1]; /* array of test quantities indexed from     */
                               /* 1 to NUM_TESTS(=5)                        */
  realtype cv_l[L_MAX];        /* coefficients of l(x) (degree q poly)      */

  realtype cv_rl1;     /* 1 / l[1]                     */
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

  int cv_qu;             /* last successful q value used   */
  long int cv_nstlp;          /* step number of last setup call */
  realtype cv_h0u;       /* actual initial stepsize        */
  realtype cv_hu;        /* last successful h value used   */
  realtype cv_saved_tq5; /* saved value of tq[5]           */
  booleantype cv_jcur;   /* Is the Jacobian info used by   */
                         /* linear solver current?         */
  realtype cv_tolsf;     /* tolerance scale factor         */
  booleantype cv_setupNonNull; /* Does setup do something? */

  booleantype cv_VabstolMallocDone;
  booleantype cv_MallocDone;  

  /*----------
    Error File 
    ----------*/

  FILE *cv_errfp;       /* CVODE error messages are sent to errfp */

  /*-------------------------
    Stability Limit Detection
    -------------------------*/

  booleantype cv_sldeton;     /* Is Stability Limit Detection on?          */
  realtype cv_ssdat[6][4];    /* scaled data array for STALD               */
  int cv_nscon;               /* counter for STALD method                  */
  long int cv_nor;            /* counter for number of order reductions    */

  /*----------------
    Rootfinding Data
    ----------------*/

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

} *CVodeMem;

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

/* CvodeCreate Error Messages */

#define _CVC_ "CVodeCreate-- "

#define MSGCV_BAD_LMM1 _CVC_ "Illegal value for lmm.\n"
#define MSGCV_BAD_LMM2 "The legal values are CV_ADAMS and CV_BDF.\n\n"
#define MSGCV_BAD_LMM  MSGCV_BAD_LMM1 MSGCV_BAD_LMM2

#define MSGCV_BAD_ITER1 _CVC_ "Illegal value for iter.\n"
#define MSGCV_BAD_ITER2 "The legal values are CV_FUNCTIONAL "
#define MSGCV_BAD_ITER3 "and CV_NEWTON.\n\n"
#define MSGCV_BAD_ITER  MSGCV_BAD_ITER1 MSGCV_BAD_ITER2 MSGCV_BAD_ITER3

#define MSGCV_CVMEM_FAIL _CVC_ "Allocation of cv_mem failed.\n\n"

/* CVodeSet* Error Messages */

#define MSGCV_SET_NO_MEM "cvode_mem = NULL in a CVodeSet routine illegal.\n\n"

#define MSGCV_SET_BAD_ITER1 "CVodeSetIterType-- Illegal value for iter.\n"
#define MSGCV_SET_BAD_ITER2 "The legal values are CV_FUNCTIONAL "
#define MSGCV_SET_BAD_ITER3 "and CV_NEWTON.\n\n"
#define MSGCV_SET_BAD_ITER   MSGCV_SET_BAD_ITER1 MSGCV_SET_BAD_ITER2 MSGCV_SET_BAD_ITER3

#define MSGCV_SET_NEG_MAXORD "CVodeSetMaxOrd-- maxord <= 0 illegal.\n\n"

#define MSGCV_SET_BAD_MAXORD1 "CVodeSetMaxOrd-- Illegal attempt to increase "
#define MSGCV_SET_BAD_MAXORD2 "maximum method order.\n\n"
#define MSGCV_SET_BAD_MAXORD  MSGCV_SET_BAD_MAXORD1 MSGCV_SET_BAD_MAXORD2 

#define MSGCV_SET_NEG_MXSTEPS "CVodeSetMaxNumSteps-- mxsteps < 0 illegal.\n\n"

#define MSGCV_SET_SLDET1 "CVodeSetStabLimDet-- Attempt to use stability "
#define MSGCV_SET_SLDET2 "limit detection with the CV_ADAMS method illegal.\n\n"
#define MSGCV_SET_SLDET  MSGCV_SET_SLDET1 MSGCV_SET_SLDET2

#define MSGCV_SET_NEG_HMIN "CVodeSetMinStep-- hmin < 0 illegal.\n\n"

#define MSGCV_SET_NEG_HMAX "CVodeSetMaxStep-- hmax < 0 illegal.\n\n"

#define MSGCV_SET_BAD_HMM1      "CVodeSetMinStep/CVodeSetMaxStep-- Inconsistent \n"
#define MSGCV_SET_BAD_HMM2      "step size limits: hmin > hmax.\n\n"
#define MSGCV_SET_BAD_HMIN_HMAX MSGCV_SET_BAD_HMM1 MSGCV_SET_BAD_HMM2

#define _CVSET_TOL_ "CVodeSetTolerances-- "

#define MSGCV_SET_NO_MALLOC _CVSET_TOL_ "Attempt to call before CVodeMalloc.\n\n"

#define MSGCV_SET_BAD_ITOL1 _CVSET_TOL_ "Illegal value for itol.\n"
#define MSGCV_SET_BAD_ITOL2 "The legal values are CV_SS, CV_SV, CV_WF.\n\n"
#define MSGCV_SET_BAD_ITOL  MSGCV_SET_BAD_ITOL1 MSGCV_SET_BAD_ITOL2

#define MSGCV_SET_BAD_RELTOL _CVSET_TOL_ "reltol < 0 illegal.\n\n"

#define MSGCV_SET_ABSTOL_NULL _CVSET_TOL_ "abstol = NULL illegal.\n\n"

#define MSGCV_SET_BAD_ABSTOL _CVSET_TOL_ "abstol has negative component(s) (illegal).\n\n"

#define MSGCV_SET_NO_EFUN1 "CVodeSetEdata-- Attempt to set e_data before specifying efun\n"
#define MSGCV_SET_NO_EFUN2 "(through CVodeMalloc or CVodeSetTolerances) illegal.\n\n"
#define MSGCV_SET_NO_EFUN MSGCV_SET_NO_EFUN1 MSGCV_SET_NO_EFUN2

/* CVodeMalloc/CVodeReInit Error Messages */

#define _CVM_ "CVodeMalloc/CVodeReInit-- "

#define MSGCV_CVM_NO_MEM _CVM_ "cvode_mem = NULL illegal.\n\n"

#define MSGCV_Y0_NULL _CVM_ "y0 = NULL illegal.\n\n"

#define MSGCV_BAD_ITOL1 _CVM_ "Illegal value for itol.\n"
#define MSGCV_BAD_ITOL2 "The legal values are CV_SS, CV_SV, and CV_WF.\n\n"
#define MSGCV_BAD_ITOL  MSGCV_BAD_ITOL1 MSGCV_BAD_ITOL2

#define MSGCV_F_NULL _CVM_ "f = NULL illegal.\n\n"

#define MSGCV_BAD_RELTOL _CVM_ "reltol < 0 illegal.\n\n"

#define MSGCV_ABSTOL_NULL _CVM_ "abstol = NULL illegal.\n\n"

#define MSGCV_BAD_ABSTOL _CVM_ "abstol has negative component(s) (illegal).\n\n"

#define MSGCV_BAD_NVECTOR _CVM_ "A required vector operation is not implemented.\n\n"

#define MSGCV_MEM_FAIL _CVM_ "A memory request failed.\n\n"

#define MSGCV_BAD_EWT _CVM_ "Initial ewt has component(s) equal to zero (illegal).\n\n"

#define MSGCV_CVREI_NO_MALLOC "CVodeReInit-- Attempt to call before CVodeMalloc.\n\n"

/* CVodeRootInit Error Messages */

#define _CVRT_ "CVodeRootInit-- "

#define MSGCV_ROOT_NO_MEM _CVRT_ "cvode_mem = NULL illegal.\n\n"

#define MSGCV_ROOT_MEM_FAIL _CVRT_ "A memory request failed.\n\n"

#define MSGCV_ROOT_FUNC_NULL _CVRT_ "g = NULL illegal.\n\n"

/* CVode Error Messages */

#define _CVODE_ "CVode-- "

#define _NO_MEM_ "cvode_mem = NULL illegal.\n\n"

#define MSGCV_CVODE_NO_MEM _CVODE_ _NO_MEM_

#define MSGCV_CVODE_NO_MALLOC _CVODE_ "CVodeMalloc has not been called yet.\n\n"
 
#define MSGCV_LSOLVE_NULL _CVODE_ "The linear solver's solve routine is NULL.\n\n"

#define MSGCV_LINIT_FAIL _CVODE_ "The linear solver's init routine failed.\n\n"

#define MSGCV_YOUT_NULL _CVODE_ "yout = NULL illegal.\n\n"

#define MSGCV_TRET_NULL _CVODE_ "tret = NULL illegal.\n\n"

#define MSGCV_BAD_ITASK _CVODE_ "Illegal value for itask.\n"

#define MSGCV_NO_TSTOP1 _CVODE_ "itask = CV_NORMAL_TSTOP or itask = CV_ONE_STEP_TSTOP "
#define MSGCV_NO_TSTOP2 _CVODE_ "but tstop was not set.\n\n"
#define MSGCV_NO_TSTOP  MSGCV_NO_TSTOP1 MSGCV_NO_TSTOP2

#define MSGCV_BAD_H0 _CVODE_ "h0 and tout - t0 inconsistent.\n\n"

#define MSGCV_HNIL_DONE_1 _CVODE_ "The above warning has been issued mxhnil times "
#define MSGCV_HNIL_DONE_2 "and will not be\nissued again for this problem.\n\n"
#define MSGCV_HNIL_DONE   MSGCV_HNIL_DONE_1 MSGCV_HNIL_DONE_2

#define MSGCV_TOO_CLOSE_1 _CVODE_ "tout too close to t0 to start"
#define MSGCV_TOO_CLOSE_2 " integration.\n\n"
#define MSGCV_TOO_CLOSE   MSGCV_TOO_CLOSE_1 MSGCV_TOO_CLOSE_2

#define MSGCV_BAD_INIT_ROOT _CVODE_ "Root found at and very near initial t.\n\n"

#define MSGCV_BAD_TOUT_1 _CVODE_ "Trouble interpolating at " MSG_TIME_TOUT ".\n"
#define MSGCV_BAD_TOUT_2 "tout too far back in direction of integration.\n\n"
#define MSGCV_BAD_TOUT   MSGCV_BAD_TOUT_1 MSGCV_BAD_TOUT_2

#define MSGCV_MAX_STEPS_1 _CVODE_ "At " MSG_TIME ", mxstep steps taken "
#define MSGCV_MAX_STEPS_2 "before reaching tout.\n\n"
#define MSGCV_MAX_STEPS   MSGCV_MAX_STEPS_1 MSGCV_MAX_STEPS_2

#define MSGCV_EWT_NOW_BAD_1 _CVODE_ "At " MSG_TIME ", a component of ewt has become <= 0.\n\n"
#define MSGCV_EWT_NOW_BAD   MSGCV_EWT_NOW_BAD_1

#define MSGCV_TOO_MUCH_ACC _CVODE_ "At " MSG_TIME ", too much accuracy requested.\n\n"

#define MSGCV_HNIL_1 _CVODE_ "Warning: Internal " MSG_TIME_H "\n"
#define MSGCV_HNIL_2 "are such that t + h = t on the next step.\n"
#define MSGCV_HNIL_3 "The solver will continue anyway.\n\n"
#define MSGCV_HNIL   MSGCV_HNIL_1 MSGCV_HNIL_2 MSGCV_HNIL_3

#define MSGCV_ERR_FAILS_1 _CVODE_ "At " MSG_TIME_H ", the error test\n"
#define MSGCV_ERR_FAILS_2 "failed repeatedly or with |h| = hmin.\n\n"
#define MSGCV_ERR_FAILS   MSGCV_ERR_FAILS_1 MSGCV_ERR_FAILS_2

#define MSGCV_CONV_FAILS_1 _CVODE_ "At " MSG_TIME_H ", the corrector\n"
#define MSGCV_CONV_FAILS_2 "convergence failed repeatedly or "
#define MSGCV_CONV_FAILS_3 "with |h| = hmin.\n\n"
#define MSGCV_CONV_FAILS   MSGCV_CONV_FAILS_1 MSGCV_CONV_FAILS_2 MSGCV_CONV_FAILS_3

#define MSGCV_SETUP_FAILED_1 _CVODE_ "At " MSG_TIME ", the setup routine failed in an "
#define MSGCV_SETUP_FAILED_2 "unrecoverable manner.\n\n"
#define MSGCV_SETUP_FAILED   MSGCV_SETUP_FAILED_1 MSGCV_SETUP_FAILED_2

#define MSGCV_SOLVE_FAILED_1 _CVODE_ "At " MSG_TIME ", the solve routine failed in an "
#define MSGCV_SOLVE_FAILED_2 "unrecoverable manner.\n\n"
#define MSGCV_SOLVE_FAILED   MSGCV_SOLVE_FAILED_1 MSGCV_SOLVE_FAILED_2

#define MSGCV_CLOSE_ROOTS _CVODE_ "Root found at and very near " MSG_TIME ".\n\n"

#define MSGCV_BAD_TSTOP_1 _CVODE_ "tstop is behind current " MSG_TIME
#define MSGCV_BAD_TSTOP_2 "\nin the direction of integration.\n\n"
#define MSGCV_BAD_TSTOP   MSGCV_BAD_TSTOP_1 MSGCV_BAD_TSTOP_2

/* CVodeGetDky Error Messages */

#define _DKY_ "CVodeGetDky-- "

#define MSGCV_DKY_NO_MEM _DKY_ _NO_MEM_

#define MSGCV_BAD_K _DKY_ "Illegal value for k.\n\n"

#define MSGCV_BAD_DKY _DKY_ "dky = NULL illegal.\n\n"

#define MSGCV_BAD_T1 _DKY_ "Illegal value for t.\n"
#define MSGCV_BAD_T2 MSG_TIME_INT
#define MSGCV_BAD_T  MSGCV_BAD_T1 MSGCV_BAD_T2

/* CVodeGet* Error Messages */

#define MSGCV_GET_NO_MEM "cvode_mem = NULL in a CVodeGet routine illegal.\n\n"

#define MSGCV_GET_NO_SLDET1 "CVodeGetNumStabLimOrderReds-- Illegal attempt "
#define MSGCV_GET_NO_SLDET2 "to call without enabling SLDET.\n\n"
#define MSGCV_GET_NO_SLDET  MSGCV_GET_NO_SLDET1 MSGCV_GET_NO_SLDET2

#ifdef __cplusplus
}
#endif

#endif
