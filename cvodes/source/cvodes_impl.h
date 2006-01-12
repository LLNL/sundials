/*
 * -----------------------------------------------------------------
 * $Revision: 1.23 $
 * $Date: 2006-01-12 20:24:07 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Implementation header file for the main CVODES integrator.
 * -----------------------------------------------------------------
 */

#ifndef _CVODES_IMPL_H
#define _CVODES_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvodes.h"

  /* Prototype of internal ewtSet function */

  int CVEwtSet(N_Vector ycur, N_Vector weight, void *e_data);


  /* Prototypes for internal sensitivity rhs DQ functions */

  void CVSensRhsDQ(int Ns, realtype t, 
                   N_Vector y, N_Vector ydot, 
                   N_Vector *yS, N_Vector *ySdot, 
                   void *fS_data,  
                   N_Vector tempv, N_Vector ftemp);

  void CVSensRhs1DQ(int Ns, realtype t, 
                    N_Vector y, N_Vector ydot, 
                    int is, N_Vector yS, N_Vector ySdot, 
                    void *fS_data,
                    N_Vector tempv, N_Vector ftemp);

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

#define HMIN_DEFAULT     RCONST(0.0)    /* hmin default value     */
#define HMAX_INV_DEFAULT RCONST(0.0)    /* hmax_inv default value */
#define MXHNIL_DEFAULT   10             /* mxhnil default value   */
#define MXSTEP_DEFAULT   500            /* mxstep default value   */

  /*                                                                
   * ifS is the type of the function returning the sensitivity
   * right-hand side. ifS can be either CV_ALLSENS if the function    
   * (of type CVSensRhsFn) returns right hand sides for all    
   * sensitivity systems at once, or CV_ONESENS if the function 
   * (of type SensRhs1Fn) returns the right hand side of one 
   *  sensitivity system at a time.                           
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

    /*-------------------------- 
      Problem Specification Data 
      --------------------------*/

    CVRhsFn cv_f;            /* y' = f(t,y(t))                               */
    void *cv_f_data;         /* user pointer passed to f                     */
    int cv_lmm;              /* lmm = ADAMS or BDF                           */
    int cv_iter;             /* iter = FUNCTIONAL or NEWTON                  */
    int cv_itol;             /* itol = SS or SV                              */

    realtype cv_reltol;      /* relative tolerance                           */
    realtype cv_Sabstol;     /* scalar absolute tolerance                    */
    N_Vector cv_Vabstol;     /* vector absolute tolerance                    */
    void *cv_e_data;         /* user pointer passed to efun                  */

    /*-----------------------
      Quadrature Related Data 
      -----------------------*/

    booleantype cv_quadr;    /* TRUE if integrating quadratures              */

    CVQuadRhsFn cv_fQ;
    void *cv_fQ_data;        /* user pointer passed to fQ                    */
    int cv_itolQ;
    booleantype cv_errconQ;

    realtype cv_reltolQ;     /* relative tolerance for quadratures           */
    realtype cv_SabstolQ;    /* scalar absolute tolerance for quadratures    */
    N_Vector cv_VabstolQ;    /* vector absolute tolerance for quadratures    */
    CVEwtFn cv_efun;         /* function to set ewt                          */

    /*------------------------
      Sensitivity Related Data 
      ------------------------*/

    booleantype cv_sensi;    /* TRUE if computing sensitivities              */

    int cv_Ns;               /* Number of sensitivities                      */

    int cv_ism;              /* ism = SIMULTANEOUS or STAGGERED              */

    CVSensRhsFn cv_fS;       /* fS = (df/dy)*yS + (df/dp)                    */
    CVSensRhs1Fn cv_fS1;     /* fS1 = (df/dy)*yS_i + (df/dp)                 */
    void *cv_fS_data;        /* user pointer passed to fS                    */
    booleantype cv_fSDQ;
    int cv_ifS;              /* ifS = ALLSENS or ONESENS                     */

    realtype *cv_p;          /* parameters in f(t,y,p)                       */
    realtype *cv_pbar;       /* scale factors for parameters                 */
    int *cv_plist;           /* list of sensitivities                        */
    realtype cv_rhomax;      /* cut-off value for centered/forward finite
                                differences                                  */

    booleantype cv_errconS;  /* TRUE if sensitivities are in err. control    */

    int cv_itolS;
    realtype cv_reltolS;     /* relative tolerance for sensitivities         */
    realtype *cv_SabstolS;   /* scalar absolute tolerances for sensi.        */
    N_Vector *cv_VabstolS;   /* vector absolute tolerances for sensi.        */

    /*-----------------------
      Nordsieck History Array 
      -----------------------*/

    N_Vector cv_zn[L_MAX];   /* Nordsieck array, of size N x (q+1).
                                zn[j] is a vector of length N (j=0,...,q)
                                zn[j] = [1/factorial(j)] * h^j * 
                                (jth derivative of the interpolating 
                                polynomial                                   */

    /*-------------------
      Vectors of length N 
      -------------------*/

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

    /*--------------------------
      Quadrature Related Vectors 
      --------------------------*/

    N_Vector cv_znQ[L_MAX];  /* Nordsieck arrays for sensitivities           */
    N_Vector cv_ewtQ;        /* error weight vector for quadratures          */
    N_Vector cv_yQ;          /* Unlike y, yQ is not allocated by the user    */
    N_Vector cv_acorQ;       /* acorQ = yQ_n(m) - yQ_n(0)                    */
    N_Vector cv_tempvQ;      /* temporary storage vector (~ tempv)           */

    /*---------------------------
      Sensitivity Related Vectors 
      ---------------------------*/

    N_Vector *cv_znS[L_MAX]; /* Nordsieck arrays for sensitivities           */
    N_Vector *cv_ewtS;       /* error weight vectors for sensitivities       */
    N_Vector *cv_yS;         /* yS=yS0 (allocated by the user)               */
    N_Vector *cv_acorS;      /* acorS = yS_n(m) - yS_n(0)                    */
    N_Vector *cv_tempvS;     /* temporary storage vector (~ tempv)           */
    N_Vector *cv_ftempS;     /* temporary storage vector (~ ftemp)           */

    /*-----------------------------------------------
      Does CVodeSensMalloc allocate additional space?
      -----------------------------------------------*/  

    booleantype cv_stgr1alloc;     /* Are ncfS1, ncfnS1, and nniS1 allocated 
                                      by CVODES?                             */

    /*-----------------
      Tstop information
      -----------------*/
    booleantype cv_istop;
    booleantype cv_tstopset;
    realtype cv_tstop;

    /*---------
      Step Data 
      ---------*/

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
    realtype cv_tretlast;        /* last value of t returned                 */

    realtype cv_tau[L_MAX+1];    /* array of previous q+1 successful step
                                    sizes indexed from 1 to q+1              */
    realtype cv_tq[NUM_TESTS+1]; /* array of test quantities indexed from
                                    1 to NUM_TESTS(=5)                       */
    realtype cv_l[L_MAX];        /* coefficients of l(x) (degree q poly)     */

    realtype cv_rl1;             /* the scalar 1/l[1]                        */
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

    /*------
      Limits 
      ------*/

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

    /*-----------------------------
      Space requirements for CVODES 
      -----------------------------*/

    long int cv_lrw1;        /* no. of realtype words in 1 N_Vector y        */ 
    long int cv_liw1;        /* no. of integer words in 1 N_Vector y         */ 
    long int cv_lrw1Q;       /* no. of realtype words in 1 N_Vector yQ       */ 
    long int cv_liw1Q;       /* no. of integer words in 1 N_Vector yQ        */ 
    long int cv_lrw;         /* no. of realtype words in CVODES work vectors */
    long int cv_liw;         /* no. of integer words in CVODES work vectors  */

    /*----------------
      Step size ratios
      ----------------*/

    realtype cv_etaqm1;      /* ratio of new to old h for order q-1          */
    realtype cv_etaq;        /* ratio of new to old h for order q            */
    realtype cv_etaqp1;      /* ratio of new to old h for order q+1          */

    /*------------------
      Linear Solver Data 
      ------------------*/

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
    int cv_qmax_allocQ;    /* value of qmax used when allocating quad. memory     */
    int cv_qmax_allocS;    /* value of qmax used when allocating sensi. memory    */
    int cv_indx_acor;      /* index of the zn vector in which acor is saved       */
    booleantype cv_setupNonNull; /* Does setup do something?                      */

    /*--------------------------------------------------------------------
      Flags turned ON by CVodeMalloc, CVodeSensMalloc, and CVodeQuadMalloc 
      and read by CVodeReInit, CVodeSensReInit, and CVodeQuadReInit
      --------------------------------------------------------------------*/

    booleantype cv_VabstolMallocDone;
    booleantype cv_MallocDone;

    booleantype cv_VabstolQMallocDone;
    booleantype cv_quadMallocDone;

    booleantype cv_VabstolSMallocDone;
    booleantype cv_SabstolSMallocDone;
    booleantype cv_sensMallocDone;

    /*----------
      Error File 
      ----------*/

    FILE *cv_errfp;             /* CVODE error messages are sent to errfp    */

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
    realtype cv_trout;    /* t value returned by rootfinding routine         */
    realtype *cv_glo;     /* saved array of g values at t = tlo              */
    realtype *cv_ghi;     /* saved array of g values at t = thi              */
    realtype *cv_grout;   /* array of g values at t = trout                  */
    realtype cv_toutc;    /* copy of tout (if NORMAL mode)                   */
    realtype cv_ttol;     /* tolerance on root location trout                */
    int cv_taskc;         /* copy of parameter task                          */
    int cv_irfnd;         /* flag showing whether last step had a root       */
    long int cv_nge;      /* counter for g evaluations                       */

    /*-------------------------
      Complex step memory block
      -------------------------*/

    void *cv_csmem;

  } *CVodeMem;

  /* 
   * =================================================================
   *   C V O D E S    E R R O R    M E S S A G E S
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

#define MSGCVS_BAD_LMM1 _CVC_ "Illegal value for lmm.\n"
#define MSGCVS_BAD_LMM2 "The legal values are CV_ADAMS and CV_BDF.\n\n"
#define MSGCVS_BAD_LMM  MSGCVS_BAD_LMM1 MSGCVS_BAD_LMM2

#define MSGCVS_BAD_ITER1 _CVC_ "Illegal value for iter.\n"
#define MSGCVS_BAD_ITER2 "The legal values are CV_FUNCTIONAL "
#define MSGCVS_BAD_ITER3 "and CV_NEWTON.\n\n"
#define MSGCVS_BAD_ITER  MSGCVS_BAD_ITER1 MSGCVS_BAD_ITER2 MSGCVS_BAD_ITER3

#define MSGCVS_CVMEM_FAIL _CVC_ "Allocation of cv_mem failed.\n\n"

  /* CVodeSet* Error Messages */

#define MSGCVS_SET_NO_MEM "cvode_mem = NULL in a CVodeSet routine illegal.\n\n"

#define MSGCVS_SET_BAD_ITER1 "CVodeSetIterType-- Illegal value for iter.\n"
#define MSGCVS_SET_BAD_ITER2 "The legal values are CV_FUNCTIONAL "
#define MSGCVS_SET_BAD_ITER3 "and CV_NEWTON.\n\n"
#define MSGCVS_SET_BAD_ITER  MSGCVS_SET_BAD_ITER1 MSGCVS_SET_BAD_ITER2 MSGCVS_SET_BAD_ITER3

#define MSGCVS_SET_NEG_MAXORD "CVodeSetMaxOrd-- maxord <= 0 illegal.\n\n"

#define MSGCVS_SET_BAD_MAXORD1 "CVodeSetMaxOrd-- Illegal attempt to increase "
#define MSGCVS_SET_BAD_MAXORD2 "maximum method order.\n\n"
#define MSGCVS_SET_BAD_MAXORD  MSGCVS_SET_BAD_MAXORD1 MSGCVS_SET_BAD_MAXORD2 

#define MSGCVS_SET_NEG_MXSTEPS "CVodeSetMaxNumSteps-- mxsteps < 0 illegal.\n\n"

#define MSGCVS_SET_SLDET1 "CVodeSetStabLimDet-- Attempt to use stability "
#define MSGCVS_SET_SLDET2 "limit detection with the CV_ADAMS method illegal.\n\n"
#define MSGCVS_SET_SLDET  MSGCVS_SET_SLDET1 MSGCVS_SET_SLDET2

#define MSGCVS_SET_NEG_HMIN "CVodeSetMinStep-- hmin < 0 illegal.\n\n"

#define MSGCVS_SET_NEG_HMAX "CVodeSetMaxStep-- hmax < 0 illegal.\n\n"

#define MSGCVS_SET_BAD_HMM1       "CVodeSetMinStep/CVodeSetMaxStep-- Inconsistent\n"
#define MSGCVS_SET_BAD_HMM2       "step size limits: hmin > hmax.\n\n"
#define MSGCVS_SET_BAD_HMIN_HMAX  MSGCVS_SET_BAD_HMM1 MSGCVS_SET_BAD_HMM2

#define _CVSET_TOL_ "CVodeSetTolerances-- "

#define MSGCVS_SET_NO_MALLOC _CVSET_TOL_ "Attempt to call before CVodeMalloc.\n\n"

#define MSGCVS_SET_BAD_ITOL1 _CVSET_TOL_ "Illegal value for itol.\n"
#define MSGCVS_SET_BAD_ITOL2 "The legal values are CV_SS and CV_SV.\n\n"
#define MSGCVS_SET_BAD_ITOL  MSGCVS_SET_BAD_ITOL1 MSGCVS_SET_BAD_ITOL2

#define MSGCVS_SET_BAD_RELTOL _CVSET_TOL_ "reltol < 0 illegal.\n\n"

#define MSGCVS_SET_ABSTOL_NULL _CVSET_TOL_ "abstol = NULL illegal.\n\n"

#define MSGCVS_SET_BAD_ABSTOL _CVSET_TOL_ "abstol has negative component(s) (illegal).\n\n"

  /* CVodeMalloc/CVodeReInit Error Messages */

#define _CVM_ "CVodeMalloc/CVodeReInit-- "

#define MSGCVS_CVM_NO_MEM _CVM_ "cvode_mem = NULL illegal.\n\n"

#define MSGCVS_Y0_NULL _CVM_ "y0 = NULL illegal.\n\n"

#define MSGCVS_BAD_ITOL1 _CVM_ "Illegal value for itol.\n"
#define MSGCVS_BAD_ITOL2 "The legal values are CV_SS, CV_SV, and CV_WF.\n\n"
#define MSGCVS_BAD_ITOL  MSGCVS_BAD_ITOL1 MSGCVS_BAD_ITOL2

#define MSGCVS_F_NULL _CVM_ "f = NULL illegal.\n\n"

#define MSGCVS_BAD_RELTOL _CVM_ "reltol < 0 illegal.\n\n"

#define MSGCVS_ABSTOL_NULL _CVM_ "abstol = NULL illegal.\n\n"

#define MSGCVS_BAD_ABSTOL _CVM_ "abstol has negative component(s) (illegal).\n\n"

#define MSGCVS_BAD_NVECTOR _CVM_ "A required vector operation is not implemented.\n\n"

#define MSGCVS_MEM_FAIL _CVM_ "A memory request failed.\n\n"

#define MSGCVS_CVREI_NO_MALLOC "CVodeReInit-- Attempt to call before CVodeMalloc.\n\n"

  /* CVodeRootInit Error Messages */

#define _CVRT_ "CVodeRootInit-- "

#define MSGCVS_ROOT_NO_MEM _CVRT_ "cvode_mem = NULL illegal.\n\n"

#define MSGCVS_ROOT_MEM_FAIL _CVRT_ "A memory request failed.\n\n"

#define MSGCVS_ROOT_FUNC_NULL _CVRT_ "g = NULL illegal.\n\n"

  /* CVodeQuadMalloc/CVodeQuadReInit Error Messages */

#define _CVSET_QTOL_ "CVodeSetQuadErrCon-- "

#define MSGCVS_SET_NO_QUAD  _CVSET_QTOL_ "Illegal attempt to call before calling CVodeQuadMalloc.\n\n"

#define MSGCVS_BAD_ITOLQ1 "Illegal value for itolQ.\nThe legal values are CV_SS and CV_SV.\n\n"
#define MSGCVS_BAD_ITOLQ  _CVSET_QTOL_ MSGCVS_BAD_ITOLQ1

#define MSGCVS_ABSTOLQ_NULL _CVSET_QTOL_ "abstolQ = NULL illegal.\n\n"

#define MSGCVS_BAD_RELTOLQ _CVSET_QTOL_ "reltolQ < 0 illegal.\n\n"

#define MSGCVS_BAD_ABSTOLQ _CVSET_QTOL_ "abstolQ has negative component(s) (illegal).\n\n"  

#define MSGCVS_QCVM_NO_MEM "CVodeQuadMalloc/CVodeQuadReInit-- cvode_mem = NULL illegal.\n\n"

#define MSGCVS_QCVM_MEM_FAIL "CVodeQuadMalloc/CVodeQuadReInit-- A memory request failed.\n\n"

#define MSGCVS_QREI_QUAD1   "CVodeQuadReInit-- Illegal attempt to call before "
#define MSGCVS_QREI_QUAD2   "calling CVodeQuadMalloc.\n\n"
#define MSGCVS_QREI_NO_QUAD MSGCVS_QREI_QUAD1 MSGCVS_QREI_QUAD2

  /* CVodeSetSens* /CVodeSensMalloc/CVodeSensReInit Error Messages */

#define _CVSET_STOL_ "CVodeSetSensTolerances-- "

#define MSGCVS_SET_SENSI_1  "Illegal attempt to call before calling CVodeSensMalloc.\n\n"
#define MSGCVS_SET_NO_SENSI _CVSET_STOL_ MSGCVS_SET_SENSI_1

#define MSGCVS_BAD_ITOLS1 _CVSET_STOL_ "Illegal value for itolS.\n"
#define MSGCVS_BAD_ITOLS2 "The legal values are CV_SS, CV_SV, and CV_EE.\n\n"
#define MSGCVS_BAD_ITOLS  MSGCVS_BAD_ITOLS1 MSGCVS_BAD_ITOLS2

#define MSGCVS_ABSTOLS_NULL _CVSET_STOL_ "abstolS = NULL illegal.\n\n"

#define MSGCVS_BAD_RELTOLS _CVSET_STOL_ "reltolS < 0 illegal.\n\n"

#define MSGCVS_BAD_ABSTOLS _CVSET_STOL_ "abstolS has negative component(s) (illegal).\n\n"  

#define MSGCVS_BAD_PBAR "CVodeSetSensParams-- pbar has zero component(s) (illegal).\n\n"
#define MSGCVS_BAD_PLIST "CVodeSetSensParams-- plist has negative component(s) (illegal).\n\n"

#define _SCVM_ "CVodeSensMalloc/CVodeSensReInit-- "

#define MSGCVS_SCVM_NO_MEM _SCVM_ "cvode_mem = NULL illegal.\n\n"

#define MSGCVS_SCVM_MEM_FAIL _SCVM_ "A memory request failed.\n\n"

#define MSGCVS_BAD_NS _SCVM_ "NS <= 0 illegal.\n\n"


#define MSGCVS_YS0_NULL _SCVM_ "yS0 = NULL illegal.\n\n"

#define MSGCVS_BAD_ISM1 _SCVM_ "Illegal value for ism.\n"
#define MSGCVS_BAD_ISM2 "The legal values are: "
#define MSGCVS_BAD_ISM3 "CV_SIMULTANEOUS, CV_STAGGERED and CV_STAGGERED1.\n\n"
#define MSGCVS_BAD_ISM  MSGCVS_BAD_ISM1 MSGCVS_BAD_ISM2 MSGCVS_BAD_ISM3

#define MSGCVS_SREI_SENSI1 "CVodeSensReInit-- Illegal attempt to call before "
#define MSGCVS_SREI_SENSI2 "calling CVodeSensMalloc.\n\n"
#define MSGCVS_SREI_NO_SENSI MSGCVS_SREI_SENSI1 MSGCVS_SREI_SENSI2

#define MSGCVS_SCVT_NO_MEM "CVodeSensToggle-- cvode_mem = NULL illegal.\n\n"

#define MSGCVS_SCVT_SENSI1 "CVodeSensToggle-- Illegal attempt to call before "
#define MSGCVS_SCVT_SENSI2 "calling CVodeSensMalloc.\n\n"
#define MSGCVS_SCVT_NO_SENSI MSGCVS_SCVT_SENSI1 MSGCVS_SCVT_SENSI2

  /* CVode Error Messages */

#define _CVODE_ "CVode-- "
#define _CVIS_  "Initial Setup: "

#define MSGCVS_CVODE_NO_MEM _CVODE_ "cvode_mem = NULL illegal.\n\n"

#define MSGCVS_CVODE_NO_MALLOC _CVODE_ "CVodeMalloc has not been called yet.\n\n"

#define MSGCVS_NO_EFUN _CVODE_ _CVIS_ "itol = CV_WF but no EwtSet function was provided.\n\n"

#define MSGCVS_FAIL_EWT _CVODE_ _CVIS_ "The user-provide EwtSet function failed.\n\n"

#define MSGCVS_BAD_EWT _CVODE_ _CVIS_ "Initial ewt has component(s) equal to zero (illegal).\n\n"

#define MSGCVS_BAD_EWTQ _CVODE_ _CVIS_ "Initial ewtQ has component(s) equal to zero (illegal).\n\n"

#define MSGCVS_BAD_ISM_IFS _CVODE_ _CVIS_ "Illegal sens. rhs for ism = CV_STAGGERED1.\n\n"

#define MSGCVS_BAD_EWTS _CVODE_ _CVIS_ "Initial ewtS has component(s) equal to zero (illegal).\n\n"

#define MSGCVS_LSOLVE_NULL _CVODE_ _CVIS_ "The linear solver's solve routine is NULL.\n\n"

#define MSGCVS_LINIT_FAIL _CVODE_ _CVIS_ "The linear solver's init routine failed.\n\n"

#define MSGCVS_YOUT_NULL _CVODE_ "yout = NULL illegal.\n\n"

#define MSGCVS_TRET_NULL _CVODE_ "tret = NULL illegal.\n\n"

#define MSGCVS_BAD_ITASK _CVODE_ "Illegal value for itask.\n"

#define MSGCVS_NO_TSTOP1 _CVODE_ "itask = CV_NORMAL_TSTOP or itask = CV_ONE_STEP_TSTOP "
#define MSGCVS_NO_TSTOP2 "but tstop was not set.\n\n"
#define MSGCVS_NO_TSTOP  MSGCVS_NO_TSTOP1 MSGCVS_NO_TSTOP2

#define MSGCVS_BAD_H0 _CVODE_ "h0 and tout - t0 are inconsistent.\n\n"

#define MSGCVS_P_NULL _CVODE_ "p = NULL when using internal DQ for sensitivity RHS illegal.\n\n"

#define MSGCVS_HNIL_DONE_1 _CVODE_ "The above warning has been issued mxhnil times "
#define MSGCVS_HNIL_DONE_2 "and will not be\nissued again for this problem.\n\n"
#define MSGCVS_HNIL_DONE   MSGCVS_HNIL_DONE_1 MSGCVS_HNIL_DONE_2

#define MSGCVS_TOO_CLOSE _CVODE_ "tout too close to t0 to start integration.\n\n"

#define MSGCVS_BAD_INIT_ROOT _CVODE_ "Root found at and very near initial t.\n\n"

#define MSGCVS_BAD_TOUT_1 _CVODE_ "Trouble interpolating at" MSG_TIME_TOUT ".\n"
#define MSGCVS_BAD_TOUT_2 "tout too far back in direction of integration.\n\n"
#define MSGCVS_BAD_TOUT   MSGCVS_BAD_TOUT_1 MSGCVS_BAD_TOUT_2

#define MSGCVS_MAX_STEPS _CVODE_ "At " MSG_TIME ", mxstep steps taken before reaching tout.\n\n"

#define MSGCVS_EWT_NOW_FAIL _CVODE_ "At " MSG_TIME ", the user-provide EwtSet function failed.\n\n"

#define MSGCVS_EWT_NOW_BAD _CVODE_ "At " MSG_TIME ", a component of ewt has become <= 0.\n\n"

#define MSGCVS_EWTS_NOW_BAD _CVODE_ "At " MSG_TIME ", a component of ewtS has become <= 0.\n\n"

#define MSGCVS_EWTQ_NOW_BAD _CVODE_ "At " MSG_TIME ", a component of ewtQ has become <= 0.\n\n"

#define MSGCVS_TOO_MUCH_ACC _CVODE_ "At " MSG_TIME ", too much accuracy requested.\n\n"

#define MSGCVS_HNIL_1 _CVODE_ "Warning: Internal " MSG_TIME_H
#define MSGCVS_HNIL_2 "\nare such that t + h = t on the next step.\n"
#define MSGCVS_HNIL_3 "The solver will continue anyway.\n\n"
#define MSGCVS_HNIL   MSGCVS_HNIL_1 MSGCVS_HNIL_2 MSGCVS_HNIL_3

#define MSGCVS_ERR_FAILS_1 _CVODE_ "At " MSG_TIME_H ", the error test\n"
#define MSGCVS_ERR_FAILS_2 "failed repeatedly or with |h| = hmin.\n\n"
#define MSGCVS_ERR_FAILS   MSGCVS_ERR_FAILS_1 MSGCVS_ERR_FAILS_2

#define MSGCVS_CONV_FAILS_1 _CVODE_ "At " MSG_TIME_H ", the corrector\n"
#define MSGCVS_CONV_FAILS_2 "convergence failed repeatedly or with |h| = hmin.\n\n"
#define MSGCVS_CONV_FAILS   MSGCVS_CONV_FAILS_1 MSGCVS_CONV_FAILS_2

#define MSGCVS_SETUP_FAILED_1 _CVODE_ "At " MSG_TIME ", the setup routine failed "
#define MSGCVS_SETUP_FAILED_2 "in an unrecoverable manner.\n\n"
#define MSGCVS_SETUP_FAILED   MSGCVS_SETUP_FAILED_1 MSGCVS_SETUP_FAILED_2

#define MSGCVS_SOLVE_FAILED_1 _CVODE_ "At " MSG_TIME ", the solve routine failed in an "
#define MSGCVS_SOLVE_FAILED_2 "unrecoverable manner.\n\n"
#define MSGCVS_SOLVE_FAILED   MSGCVS_SOLVE_FAILED_1 MSGCVS_SOLVE_FAILED_2

#define MSGCVS_BAD_TSTOP_1 _CVODE_ "tstop is behind current " MSG_TIME
#define MSGCVS_BAD_TSTOP_2 "\nin the direction of integration.\n\n"
#define MSGCVS_BAD_TSTOP   MSGCVS_BAD_TSTOP_1 MSGCVS_BAD_TSTOP_2

#define MSGCVS_CLOSE_ROOTS _CVODE_ "Root found at and very near current " MSG_TIME ".\n\n"

  /* CVodeGetDky Error Messages */

#define _DKY_ "CVodeGetDky-- "

#define MSGCVS_DKY_NO_MEM _DKY_ "cvode_mem = NULL illegal.\n\n"

#define MSGCVS_BAD_K _DKY_ "Illegal value for k.\n\n"

#define MSGCVS_BAD_DKY _DKY_ "dky = NULL illegal.\n\n"

#define MSGCVS_BAD_T1 _DKY_ "Illegal value for t.\n"
#define MSGCVS_BAD_T2 MSG_TIME_INT
#define MSGCVS_BAD_T  MSGCVS_BAD_T1 MSGCVS_BAD_T2

  /* CVodeGetSens/CVodeGetSens1/CVodeGetSensDky1/CVodeGetSensDky Error Messages */

#define _SDKY_ "CVodeGetSens/CVodeGetSens1/CVodeGetSensDky/CVodeGetSensDky1-- "

#define MSGCVS_SDKY_NO_MEM _SDKY_ "cvode_mem = NULL illegal.\n\n"

#define MSGCVS_SDKY_SENSI_1  "Illegal attempt to call before calling CVodeSensMalloc.\n\n"
#define MSGCVS_SDKY_NO_SENSI _SDKY_ MSGCVS_SDKY_SENSI_1

#define MSGCVS_SBAD_IS _SDKY_ "Illegal value for is.\n\n"

#define MSGCVS_SBAD_K _SDKY_ "Illegal value for k.\n\n"

#define MSGCVS_SBAD_T_1 _SDKY_ "Illegal value for t.\n"
#define MSGCVS_SBAD_T_2 "t not in interval tcur - hu to tcur.\n\n"
#define MSGCVS_SBAD_T   MSGCVS_SBAD_T_1 MSGCVS_SBAD_T_2

#define MSGCVS_SBAD_DKYA _SDKY_ "dkyA = NULL illegal.\n\n"
#define MSGCVS_SBAD_DKY  _SDKY_ "dky = NULL illegal.\n\n"

  /* CVodeGetQuad/CVodeGetQuadDky Error Messages */

#define _QDKY_ "CVodeGetQuad/CVodeGetQuadDky-- "

#define MSGCVS_QDKY_NO_MEM _QDKY_ "cvode_mem = NULL illegal.\n\n"

#define MSGCVS_QDKY_QUAD_1  "Illegal attempt to call before calling CVodeQuadMalloc.\n\n"
#define MSGCVS_QDKY_NO_QUAD _QDKY_ MSGCVS_QDKY_QUAD_1

#define MSGCVS_QBAD_DKY _QDKY_ "dky = NULL illegal.\n\n"

#define MSGCVS_QBAD_K _QDKY_ "Illegal value for k.\n\n"

#define MSGCVS_QBAD_T_1 _QDKY_ "Illegal value for t.\n"
#define MSGCVS_QBAD_T_2 MSG_TIME_INT
#define MSGCVS_QBAD_T   MSGCVS_QBAD_T_1 MSGCVS_QBAD_T_2

  /* CVodeGet* Error Messages */

#define MSGCVS_GET_NO_MEM    "cvode_mem = NULL in a CVodeGet routine illegal. \n\n"

#define MSGCVS_GET_NO_QUAD1  "CVodeGetQuad*-- Illegal attempt to call before "
#define MSGCVS_GET_NO_QUAD2  "calling CVodeQuadMalloc.\n\n"
#define MSGCVS_GET_NO_QUAD   MSGCVS_GET_NO_QUAD1 MSGCVS_GET_NO_QUAD2

#define MSGCVS_GET_NO_SENSI1 "CVodeGetSens*-- Illegal attempt to call before "
#define MSGCVS_GET_NO_SENSI2 "calling CVodeSensMalloc.\n\n"
#define MSGCVS_GET_NO_SENSI  MSGCVS_GET_NO_SENSI1 MSGCVS_GET_NO_SENSI2

#define MSGCVS_GET_EWT_BAD   "CVodeGetErrWeights--  ewt has component(s) equal to zero.\n\n"

#ifdef __cplusplus
}
#endif

#endif
