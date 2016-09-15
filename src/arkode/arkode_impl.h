/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * Implementation header file for the main ARKODE integrator.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_IMPL_H
#define _ARKODE_IMPL_H

#include <stdarg.h>
#include <arkode/arkode.h>
#include <arkode/arkode_spils.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
             ARKODE Private Constants                             
===============================================================*/

/* Basic ARKODE constants */
#define Q_DEFAULT        4       /* default RK order */
#define QDENSE_DEF       3       /* default dense output order */
#define MXSTEP_DEFAULT   500     /* mxstep default value */
#define MAXNEF           7       /* maxnef default value */
#define MAXNCF           10      /* maxncf default value */
#define MXHNIL           10      /* mxhnil default value */
#define MAXCOR           3       /* maxcor default value */
#define FP_ACCEL_M       3       /* fp_m default value */

/* Numeric constants */
#define ZERO   RCONST(0.0)      /* real 0.0     */
#define TINY   RCONST(1.0e-10)  /* small number */
#define TENTH  RCONST(0.1)      /* real 0.1     */
#define POINT2 RCONST(0.2)      /* real 0.2     */
#define FOURTH RCONST(0.25)     /* real 0.25    */
#define HALF   RCONST(0.5)      /* real 0.5     */
#define ONE    RCONST(1.0)      /* real 1.0     */
#define TWO    RCONST(2.0)      /* real 2.0     */
#define THREE  RCONST(3.0)      /* real 3.0     */
#define FOUR   RCONST(4.0)      /* real 4.0     */
#define FIVE   RCONST(5.0)      /* real 5.0     */
#define SIX    RCONST(6.0)      /* real 6.0     */
#define SEVEN  RCONST(7.0)      /* real 7.0     */
#define TWELVE RCONST(12.0)     /* real 12.0    */
#define HUND   RCONST(100.0)    /* real 100.0   */

/* Time step controller default values */
#define CFLFAC    RCONST(0.5)
#define SAFETY    RCONST(0.96)  /* CVODE uses 1.0  */
#define BIAS      RCONST(1.5)   /* CVODE uses 6.0  */
#define GROWTH    RCONST(20.0)  /* CVODE uses 10.0 */
#define HFIXED_LB RCONST(1.0)   /* CVODE uses 1.0  */
#define HFIXED_UB RCONST(1.5)   /* CVODE uses 1.5  */
#define AD0_K1    RCONST(0.58)  /* PID controller constants */
#define AD0_K2    RCONST(0.21)
#define AD0_K3    RCONST(0.1)
#define AD1_K1    RCONST(0.8)   /* PI controller constants */
#define AD1_K2    RCONST(0.31)
#define AD2_K1    RCONST(1.0)   /* I controller constants */
#define AD3_K1    RCONST(0.367) /* explicit Gustafsson controller */
#define AD3_K2    RCONST(0.268)
#define AD4_K1    RCONST(0.98)  /* implicit Gustafsson controller */
#define AD4_K2    RCONST(0.95)
#define AD5_K1    RCONST(0.367) /* imex Gustafsson controller */
#define AD5_K2    RCONST(0.268)
#define AD5_K3    RCONST(0.95)

/* Default solver tolerance factor */
/* #define NLSCOEF   RCONST(0.003)   /\* Hairer & Wanner constant *\/ */
/* #define NLSCOEF   RCONST(0.2)     /\* CVODE constant *\/ */
#define NLSCOEF   RCONST(0.1)

/* Control constants for tolerances */
#define ARK_SS  0
#define ARK_SV  1
#define ARK_WF  2


/*===============================================================
             ARKODE Routine-Specific Constants                   
===============================================================*/

/*---------------------------------------------------------------
 Control constants for lower-level functions used by arkStep:
-----------------------------------------------------------------
 arkHin return values:  ARK_SUCCESS, ARK_RHSFUNC_FAIL, or 
    ARK_TOO_CLOSE

 arkStep control constants:  SOLVE_SUCCESS or PREDICT_AGAIN

 arkStep return values:  ARK_SUCCESS, ARK_LSETUP_FAIL, 
    ARK_LSOLVE_FAIL, ARK_RHSFUNC_FAIL, ARK_RTFUNC_FAIL,
    ARK_CONV_FAILURE, ARK_ERR_FAILURE or ARK_FIRST_RHSFUNC_ERR

 arkNls input nflag values:  FIRST_CALL, PREV_CONV_FAIL or 
    PREV_ERR_FAIL
    
 arkNls return values:  ARK_SUCCESS, ARK_LSETUP_FAIL,
    ARK_LSOLVE_FAIL, ARK_RHSFUNC_FAIL, CONV_FAIL or
    RHSFUNC_RECVR
 
 arkNewtonIteration return values:  ARK_SUCCESS, ARK_LSOLVE_FAIL,
    ARK_RHSFUNC_FAIL, CONV_FAIL, RHSFUNC_RECVR or TRY_AGAIN
---------------------------------------------------------------*/
#define SOLVE_SUCCESS    +2
#define PREDICT_AGAIN    +3

#define CONV_FAIL        +4 
#define TRY_AGAIN        +5

#define FIRST_CALL       +6
#define PREV_CONV_FAIL   +7
#define PREV_ERR_FAIL    +8

#define RHSFUNC_RECVR    +9


/*---------------------------------------------------------------
 Return values for lower-level rootfinding functions
-----------------------------------------------------------------
 arkRootCheck1:  ARK_SUCCESS or ARK_RTFUNC_FAIL

 arkRootCheck2:  ARK_SUCCESS, ARK_RTFUNC_FAIL, CLOSERT or RTFOUND

 arkRootCheck3:  ARK_SUCCESS, ARK_RTFUNC_FAIL or RTFOUND

 arkRootfind:  ARK_SUCCESS, ARK_RTFUNC_FAIL or RTFOUND
---------------------------------------------------------------*/
#define RTFOUND          +1
#define CLOSERT          +3


/*---------------------------------------------------------------
 Algorithmic constants
-----------------------------------------------------------------
 ARKodeGetDky and arkStep:  FUZZ_FACTOR

 arkHin:  H0_LBFACTOR, H0_UBFACTOR, H0_BIAS and H0_ITERS

 arkStep:  
    ETAMX1      maximum step size change on first step
    ETAMXF      step size reduction factor on multiple error 
                test failures (multiple implies >= SMALL_NEF)
    ETAMIN      smallest allowable step size reduction factor 
                on an error test failure
    ETACF       step size reduction factor on nonlinear 
                convergence failure
    ONEPSM      safety factor for floating point comparisons
    ONEMSM      safety factor for floating point comparisons
    SMALL_NEF   if an error failure occurs and SMALL_NEF <= nef,
                then reset  eta = MIN(eta, ETAMXF)

 arkNls:
    CRDOWN      constant used in the estimation of the 
                convergence rate (crate) of the iterates for 
                the nonlinear equation
    DGMAX       if |gamma/gammap-1| > DGMAX then call lsetup
    RDIV        declare divergence if ratio del/delp > RDIV
    MSBP        max no. of steps between lsetup calls
---------------------------------------------------------------*/
#define FUZZ_FACTOR RCONST(100.0)

#define H0_LBFACTOR RCONST(100.0)
#define H0_UBFACTOR RCONST(0.1)
#define H0_BIAS     HALF
#define H0_ITERS    4

#define ETAMX1      RCONST(10000.0)     /* default */
#define ETAMXF      RCONST(0.3)         /* default */
#define ETAMIN      RCONST(0.1)         /* default */
#define ETACF       RCONST(0.25)        /* default */
#define ONEPSM      RCONST(1.000001)
#define ONEMSM      RCONST(0.999999)
#define SMALL_NEF   2                   /* default */

#define CRDOWN      RCONST(0.3)         /* default */
#define DGMAX       RCONST(0.2)         /* default */
#define RDIV        RCONST(2.3)         /* default */
#define MSBP        20                  /* default */


/*===============================================================
  MAIN INTEGRATOR MEMORY BLOCK
===============================================================*/

/*---------------------------------------------------------------
 Types : struct ARKodeMemRec, ARKodeMem
-----------------------------------------------------------------
 The type ARKodeMem is type pointer to struct ARKodeMemRec.
 This structure contains fields to keep track of problem state.
---------------------------------------------------------------*/
typedef struct ARKodeMemRec {

  realtype ark_uround;         /* machine unit roundoff */

  /*-------------------------- 
    Problem Specification Data 
    --------------------------*/
  ARKRhsFn     ark_fe;         /* y' = fe(t,y(t)) + fi(t,y(t))          */
  ARKRhsFn     ark_fi;
  void        *ark_user_data;  /* user pointer passed to fe, fi         */
  ARKExpStabFn ark_expstab;    /* time step stability function for fe   */
  void        *ark_estab_data; /* user pointer passed to expstab        */
  int          ark_itol;       /* itol = ARK_SS (scalar, default), 
                                         ARK_SV (vector),
                                         ARK_WF (user weight function)  */
  int          ark_ritol;      /* itol = ARK_SS (scalar, default), 
                                         ARK_SV (vector),
                                         ARK_WF (user weight function)  */
  realtype     ark_reltol;     /* relative tolerance                    */
  realtype     ark_Sabstol;    /* scalar absolute solution tolerance    */
  N_Vector     ark_Vabstol;    /* vector absolute solution tolerance    */
  realtype     ark_SRabstol;   /* scalar absolute residual tolerance    */
  N_Vector     ark_VRabstol;   /* vector absolute residual tolerance    */
  booleantype  ark_user_efun;  /* TRUE if user sets efun                */
  ARKEwtFn     ark_efun;       /* function to set ewt                   */
  void        *ark_e_data;     /* user pointer passed to efun           */
  booleantype  ark_user_rfun;  /* TRUE if user sets rfun                */
  ARKRwtFn     ark_rfun;       /* function to set rwt                   */
  void        *ark_r_data;     /* user pointer passed to rfun           */
  booleantype  ark_linear;     /* TRUE if fi is linear                  */
  booleantype  ark_linear_timedep;  /* TRUE if dfi/dy depends on t      */
  booleantype  ark_explicit;   /* TRUE if fi is disabled                */
  booleantype  ark_implicit;   /* TRUE if fe is disabled                */

  /*-----------------
    Stored RHS arrays
    -----------------*/
  N_Vector ark_Fe[ARK_S_MAX];  /* explicit RHS at each stage */
  N_Vector ark_Fi[ARK_S_MAX];  /* implicit RHS at each stage */

  /*--------------------------
    other vectors of length N 
    -------------------------*/
  N_Vector ark_ewt;     /* error weight vector                               */
  N_Vector ark_rwt;     /* residual weight vector                            */
  booleantype ark_rwt_is_ewt;     /* TRUE if rwt is a pointer to ewt         */
  N_Vector ark_y;       /* y is used as temporary storage by the solver
			   The memory is provided by the user to ARKode
			   where the vector is named yout.                   */
  N_Vector ark_ycur;    /* ycur always holds the solver's current version of 
                           the solution */
  N_Vector ark_sdata;   /* Storage for old stage data in computing residual. */
  N_Vector ark_tempv;   /* temporary storage vector                          */
  N_Vector ark_acor;    /* temporary storage vector; Between steps this  
			   holds the estimaged local error                   */
  N_Vector ark_ftemp;   /* temporary storage vector                          */
  N_Vector ark_fold;    /* f(t,y) at beginning of last successful step       */
  N_Vector ark_fnew;    /* f(t,y) at end of last successful step             */
  N_Vector ark_yold;    /* y at beginning of last successful step            */
  N_Vector ark_ynew;    /* y at end of last successful step                  */
  /* N_Vector ark_fa;      /\* f at h/3 through step (high order dense output)   *\/ */
  /* N_Vector ark_fb;      /\* f at h*2/3 through step (high order dense output) *\/ */

  /*-----------------
    Tstop information
    -----------------*/
  booleantype ark_tstopset;
  realtype    ark_tstop;

  /*-----------
    Method Data 
    -----------*/  
  int ark_q;                              /* method order                   */
  int ark_p;                              /* embedding order                */
  int ark_istage;                         /* current stage                  */
  int ark_stages;                         /* number of stages               */
  int ark_dense_q;                        /* dense output polynomial order  */
  realtype ark_Ae[ARK_S_MAX*ARK_S_MAX];   /* ERK Butcher table              */
  realtype ark_Ai[ARK_S_MAX*ARK_S_MAX];   /* IRK Butcher table              */
  realtype ark_ce[ARK_S_MAX];             /* ERK method canopy nodes        */
  realtype ark_ci[ARK_S_MAX];             /* IRK method canopy nodes        */
  realtype ark_be[ARK_S_MAX];             /* ERK method solution coeffs     */
  realtype ark_bi[ARK_S_MAX];             /* IRK method solution coeffs     */
  realtype ark_b2e[ARK_S_MAX];            /* ERK method embedding coeffs    */
  realtype ark_b2i[ARK_S_MAX];            /* IRK method embedding coeffs    */

  /*---------
    Step Data 
    ---------*/  
  realtype ark_hin;             /* initial step size                        */
  realtype ark_h;               /* current step size                        */
  realtype ark_hprime;          /* step size to be used on the next step    */ 
  realtype ark_next_h;          /* step size to be used on the next step    */ 
  realtype ark_eta;             /* eta = hprime / h                         */
  realtype ark_tn;              /* current internal value of t              */
  realtype ark_tretlast;        /* value of tret last returned by ARKode    */
  realtype ark_gamma;           /* gamma = h * rl1                          */
  realtype ark_gammap;          /* gamma at the last setup call             */
  realtype ark_gamrat;          /* gamma / gammap                           */
  realtype ark_crate;           /* estimated corrector convergence rate     */
  realtype ark_eRNrm;           /* estimated residual norm, used in nonlinear 
				   and linear solver convergence tests      */
  realtype ark_nlscoef;         /* coefficient in nonlin. convergence test  */
  int      ark_mnewt;           /* Newton iteration counter                 */


  /*-------------------------
    Time Step Adaptivity Data 
    -------------------------*/
  booleantype ark_fixedstep;       /* flag to disable temporal adaptivity    */
  ARKAdaptFn ark_hadapt;           /* function to set the new time step size */
  void      *ark_hadapt_data;      /* user pointer passed to hadapt          */
  realtype   ark_hadapt_ehist[3];  /* error history for time adaptivity      */
  realtype   ark_hadapt_hhist[3];  /* step history for time adaptivity       */
  int        ark_hadapt_imethod;   /* time step adaptivity method to use:
				      -1 -> User-specified function above
				       0 -> PID controller
				       1 -> PI controller
				       2 -> I controller
				       3 -> explicit Gustafsson controller
				       4 -> implicit Gustafsson controller
				       5 -> imex Gustafsson controller       */
  realtype   ark_hadapt_cfl;       /* cfl safety factor                      */
  realtype   ark_hadapt_safety;    /* accuracy safety factor on h            */
  realtype   ark_hadapt_bias;      /* accuracy safety factor on LTE          */
  realtype   ark_hadapt_growth;    /* maximum step growth safety factor      */
  realtype   ark_hadapt_lbound;    /* eta lower bound to leave h unchanged   */
  realtype   ark_hadapt_ubound;    /* eta upper bound to leave h unchanged   */
  booleantype ark_hadapt_pq;       /* choice of using p (0) vs q (1)         */
  realtype   ark_hadapt_k1;
  realtype   ark_hadapt_k2;        /* method-specific adaptivity parameters  */
  realtype   ark_hadapt_k3;

  /*------------------------------------
    Limits and various solver parameters
    ------------------------------------*/
  long int ark_mxstep;   /* max number of internal steps for one user call */
  int      ark_maxcor;   /* max number of corrector iterations for the
			    solution of the nonlinear equation             */
  int      ark_mxhnil;   /* max number of warning messages issued to the
			    user that t+h == t for the next internal step  */
  int      ark_maxnef;   /* max number of error test failures in one step  */
  int      ark_maxncf;   /* max number of nonlin. conv. fails in one step  */
  realtype ark_hmin;     /* |h| >= hmin                                    */
  realtype ark_hmax_inv; /* |h| <= 1/hmax_inv                              */
  realtype ark_etamax;   /* eta <= etamax                                  */

  realtype ark_etamx1;   /* max step size change on first step             */
  realtype ark_etamxf;   /* h reduction factor on multiple error fails     */
  int ark_small_nef;     /* bound to determine 'multiple' above            */
  realtype ark_etacf;    /* h reduction factor on nonlinear conv fail      */

  realtype ark_crdown;   /* nonlinear convergence rate estimation constant */
  realtype ark_rdiv;     /* declare divergence if ratio del/delp > RDIV    */

  realtype ark_dgmax;    /* if |gamma/gammap-1| > DGMAX then call lsetup   */
  int ark_msbp;          /* positive => max # steps between lsetup calls
			    negative => recompute at every Newton iter.    */
  int ark_predictor;     /* choice of prediction method */

  /*--------
    Counters 
    --------*/
  long int ark_nst;          /* number of internal steps taken             */
  long int ark_nst_acc;      /* number of accuracy-limited internal steps  */
  long int ark_nst_exp;      /* number of stability-limited internal steps */
  long int ark_nst_attempts; /* number of attempted steps                  */
  long int ark_nfe;          /* number of fe calls                         */
  long int ark_nfi;          /* number of fi calls                         */
  long int ark_ncfn;         /* number of corrector convergence failures   */
  long int ark_nmassfails;   /* number of mass matrix solver failures      */
  long int ark_netf;         /* number of error test failures              */
  long int ark_nni;          /* number of Newton iterations performed      */
  long int ark_nsetups;      /* number of setup calls                      */
  int      ark_nhnil;        /* number of messages issued to the user that 
			        t+h == t for the next iternal step         */

  /*-----------------
    Diagnostic output
    -----------------*/
  booleantype ark_report;   /* flag to enable/disable diagnostic output    */
  FILE       *ark_diagfp;   /* diagnostic outputs are sent to ark_diagfp   */

  /*-----------------------------
    Space requirements for ARKODE 
    -----------------------------*/
  long int ark_lrw1;        /* no. of realtype words in 1 N_Vector          */ 
  long int ark_liw1;        /* no. of integer words in 1 N_Vector           */ 
  long int ark_lrw;         /* no. of realtype words in ARKODE work vectors */
  long int ark_liw;         /* no. of integer words in ARKODE work vectors  */


  /*-----------------------
    Fixed-point Solver Data 
    -----------------------*/
  booleantype ark_use_fp;     /* flag to enable fixed-point solver vs Newton */
  long int    ark_fp_m;       /* number of vectors to use in acceleration    */
  long int    *ark_fp_imap;   /* array of length m          */
  realtype    *ark_fp_R;      /* array of length m*m        */
  realtype    *ark_fp_gamma;  /* array of length m          */
  N_Vector    *ark_fp_df;     /* vector array of length m   */
  N_Vector    *ark_fp_dg;     /* vector array of length m   */
  N_Vector    *ark_fp_q;      /* vector array of length m   */
  N_Vector     ark_fp_fval;   /* temporary N_Vectors        */
  N_Vector     ark_fp_fold;
  N_Vector     ark_fp_gold;

  /*------------------
    Linear Solver Data 
    ------------------*/
  int (*ark_linit)(struct ARKodeMemRec *ark_mem);
  int (*ark_lsetup)(struct ARKodeMemRec *ark_mem, int convfail, N_Vector ypred,
		    N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
		    N_Vector vtemp2, N_Vector vtemp3); 
  int (*ark_lsolve)(struct ARKodeMemRec *ark_mem, N_Vector b, N_Vector weight,
		    N_Vector ycur, N_Vector fcur);
  int (*ark_lfree)(struct ARKodeMemRec *ark_mem);
  void *ark_lmem;
  int ark_lsolve_type;   /* linear solver type: 0=iterative; 1=dense; 
                                                2=band; 3=sparse; 4=custom */

  /*-----------------------
    Mass Matrix Solver Data 
    -----------------------*/
  booleantype ark_mass_matrix;   /* flag denoting use of a non-identity M  */
  long int ark_mass_solves;      /* number of mass matrix solve calls      */
  long int ark_mass_mult;        /* number of mass matrix product calls    */
  ARKSpilsMassTimesVecFn ark_mtimes;   /* mass-matrix-vector product routine */
  void *ark_mtimes_data;         /* user pointer passed to mtimes          */
  int (*ark_minit)(struct ARKodeMemRec *ark_mem);
  int (*ark_msetup)(struct ARKodeMemRec *ark_mem, N_Vector vtemp1, 
		    N_Vector vtemp2, N_Vector vtemp3); 
  int (*ark_msolve)(struct ARKodeMemRec *ark_mem, N_Vector b, N_Vector weight);
  int (*ark_mfree)(struct ARKodeMemRec *ark_mem);
  void *ark_mass_mem;
  int ark_msolve_type;   /* mass matrix type: 0=iterative; 1=dense; 
			                      2=band; 3=sparse; 4=custom */

  /*------------
    Saved Values
    ------------*/
  long int    ark_nstlp;        /* step number of last setup call             */
  realtype    ark_h0u;          /* actual initial stepsize                    */
  realtype    ark_tnew;         /* time of last successful step               */
  realtype    ark_hold;         /* last successful h value used               */
  booleantype ark_jcur;         /* is Jacobian info. for lin. solver current? */
  realtype    ark_tolsf;        /* tolerance scale factor                     */
  booleantype ark_setupNonNull; /* does ark_lsetup do anything?               */
  booleantype ark_MassSetupNonNull; /* does ark_msetup do anything?           */
  booleantype ark_VabstolMallocDone;
  booleantype ark_VRabstolMallocDone;
  booleantype ark_MallocDone;  
  booleantype ark_resized;      /* denotes first step after ARKodeResize      */
  booleantype ark_firststage;   /* denotes first stage in simulation          */

  /*-------------------------------------------
    Error handler function and error ouput file 
    -------------------------------------------*/
  ARKErrHandlerFn ark_ehfun;    /* error messages are handled by ehfun        */
  void           *ark_eh_data;  /* data pointer passed to ehfun               */
  FILE           *ark_errfp;    /* ARKODE error messages are sent to errfp    */

  /*----------------
    Rootfinding Data
    ----------------*/
  ARKRootFn    ark_gfun;        /* function g for roots sought                  */
  int          ark_nrtfn;       /* number of components of g                    */
  int         *ark_iroots;      /* array for root information                   */
  int         *ark_rootdir;     /* array specifying direction of zero-crossing  */
  realtype     ark_tlo;         /* nearest endpoint of interval in root search  */
  realtype     ark_thi;         /* farthest endpoint of interval in root search */
  realtype     ark_trout;       /* t value returned by rootfinding routine      */
  realtype    *ark_glo;         /* saved array of g values at t = tlo           */
  realtype    *ark_ghi;         /* saved array of g values at t = thi           */
  realtype    *ark_grout;       /* array of g values at t = trout               */
  realtype     ark_toutc;       /* copy of tout (if NORMAL mode)                */
  realtype     ark_ttol;        /* tolerance on root location                   */
  int          ark_taskc;       /* copy of parameter itask                      */
  int          ark_irfnd;       /* flag showing whether last step had a root    */
  long int     ark_nge;         /* counter for g evaluations                    */
  booleantype *ark_gactive;     /* array with active/inactive event functions   */
  int          ark_mxgnull;     /* num. warning messages about possible g==0    */

  /*----------------------------------------------------
    User-supplied step solution post-processing function
    ----------------------------------------------------*/
  ARKPostProcessStepFn ark_ProcessStep;

} *ARKodeMem;



/*===============================================================
     I N T E R F A C E   T O    L I N E A R   S O L V E R S
===============================================================*/
  
/*---------------------------------------------------------------
 Communication between ARKODE and a ARKODE Linear Solver
-----------------------------------------------------------------
 convfail (input to ark_lsetup)

 ARK_NO_FAILURES : Either this is the first ark_setup call for 
                   this step, or the local error test failed on 
		   the previous attempt at this step (but the 
		   Newton iteration converged).

 ARK_FAIL_BAD_J  : This value is passed to ark_lsetup if

                  (a) The previous Newton corrector iteration
                      did not converge and the linear solver's
                      setup routine indicated that its Jacobian-
                      related data is not current
                                   or
                  (b) During the previous Newton corrector
                      iteration, the linear solver's solve 
		      routine failed in a recoverable manner 
		      and the linear solver's setup routine 
		      indicated that its Jacobian-related data 
		      is not current.

 ARK_FAIL_OTHER  : During the current internal step try, the
                   previous Newton iteration failed to converge
		   even though the linear solver was using 
		   current Jacobian-related data.
---------------------------------------------------------------*/
  
/* Constants for convfail (input to ark_lsetup) */
#define ARK_NO_FAILURES 0
#define ARK_FAIL_BAD_J  1
#define ARK_FAIL_OTHER  2

/*---------------------------------------------------------------
 int (*ark_linit)(ARKodeMem ark_mem);
-----------------------------------------------------------------
 The purpose of ark_linit is to complete initializations for a
 specific linear solver, such as counters and statistics.
 An LInitFn should return 0 if it has successfully initialized 
 the ARKODE linear solver and a negative value otherwise.
 If an error does occur, an appropriate message should be sent 
 to the error handler function.
---------------------------------------------------------------*/
  
/*---------------------------------------------------------------
 int (*ark_lsetup)(ARKodeMem ark_mem, int convfail, 
                   N_Vector ypred, N_Vector fpred, 
		   booleantype *jcurPtr, N_Vector vtemp1, 
		   N_Vector vtemp2, N_Vector vtemp3);
 -----------------------------------------------------------------
 The job of ark_lsetup is to prepare the linear solver for
 subsequent calls to ark_lsolve. It may recompute Jacobian-
 related data is it deems necessary. Its parameters are as
 follows:

 ark_mem - problem memory pointer of type ARKodeMem. See the
          typedef earlier in this file.

 convfail - a flag to indicate any problem that occurred during
            the solution of the nonlinear equation on the
            current time step for which the linear solver is
            being used. This flag can be used to help decide
            whether the Jacobian data kept by a ARKODE linear
            solver needs to be updated or not.
            Its possible values have been documented above.

 ypred - the predicted y vector for the current ARKODE internal
         step.

 fpred - f(tn, ypred).

 jcurPtr - a pointer to a boolean to be filled in by ark_lsetup.
           The function should set *jcurPtr=TRUE if its Jacobian
           data is current after the call and should set
           *jcurPtr=FALSE if its Jacobian data is not current.
           Note: If ark_lsetup calls for re-evaluation of
           Jacobian data (based on convfail and ARKODE state
           data), it should return *jcurPtr=TRUE always;
           otherwise an infinite loop can result.

 vtemp1 - temporary N_Vector provided for use by ark_lsetup.

 vtemp3 - temporary N_Vector provided for use by ark_lsetup.

 vtemp3 - temporary N_Vector provided for use by ark_lsetup.

 The ark_lsetup routine should return 0 if successful, a positive
 value for a recoverable error, and a negative value for an
 unrecoverable error.
---------------------------------------------------------------*/

/*---------------------------------------------------------------
 int (*ark_lsolve)(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
                   N_Vector ycur, N_Vector fcur);
-----------------------------------------------------------------
 ark_lsolve must solve the linear equation P x = b, where
 P is some approximation to (M - gamma J), M is the system mass
 matrix, J = (df/dy)(tn,ycur), and the RHS vector b is input. The 
 N-vector ycur contains the solver's current approximation to 
 y(tn) and the vector fcur contains the N_Vector f(tn,ycur). The 
 solution is to be returned in the vector b. ark_lsolve returns 
 a positive value for a recoverable error and a negative value 
 for an unrecoverable error. Success is indicated by a 0 return 
 value.
---------------------------------------------------------------*/

/*---------------------------------------------------------------
 int (*ark_lfree)(ARKodeMem ark_mem);
-----------------------------------------------------------------
 ark_lfree should free up any memory allocated by the linear
 solver. This routine is called once a problem has been
 completed and the linear solver is no longer needed.  It should 
 return 0 upon success, or a nonzero on failure.
---------------------------------------------------------------*/



/*---------------------------------------------------------------
 int (*ark_minit)(ARKodeMem ark_mem);
-----------------------------------------------------------------
 The purpose of ark_minit is to complete initializations for a
 specific mass matrix linear solver, such as counters and 
 statistics. An function of this type should return 0 if it has 
 successfully initialized the mass matrix linear solver and a 
 negative value otherwise.  If an error does occur, an 
 appropriate message should be sent to the error handler function.
---------------------------------------------------------------*/
  
/*---------------------------------------------------------------
 int (*ark_msetup)(ARKodeMem ark_mem, N_Vector vtemp1, 
                   N_Vector vtemp2, N_Vector vtemp3);
 -----------------------------------------------------------------
 The job of ark_msetup is to prepare the mass matrix solver for
 subsequent calls to ark_msolve. It may recompute mass matrix
 related data is it deems necessary. Its parameters are as
 follows:

 ark_mem - problem memory pointer of type ARKodeMem. See the
          typedef earlier in this file.

 vtemp1 - temporary N_Vector provided for use by ark_lsetup.

 vtemp3 - temporary N_Vector provided for use by ark_lsetup.

 vtemp3 - temporary N_Vector provided for use by ark_lsetup.

 The ark_msetup routine should return 0 if successful, and a 
 negative value for an unrecoverable error.
---------------------------------------------------------------*/

/*---------------------------------------------------------------
 int (*ark_msolve)(ARKodeMem ark_mem, N_Vector b, N_Vector weight);
-----------------------------------------------------------------
 ark_msolve must solve the linear equation M x = b, where
 M is the system mass matrix, and the RHS vector b is input. The
 solution is to be returned in the vector b.  The ark_msolve
 routine returns a positive value for a recoverable error and 
 a negative value for an unrecoverable error. Success is 
 indicated by a 0 return value.
---------------------------------------------------------------*/

/*---------------------------------------------------------------
 int (*ark_mfree)(ARKodeMem ark_mem);
-----------------------------------------------------------------
 ark_mfree should free up any memory allocated by the mass matrix
 solver. This routine is called once a problem has been
 completed and the solver is no longer needed.  It should return
 0 upon success, or a nonzero on failure.
---------------------------------------------------------------*/

  
/*===============================================================
   ARKODE PROTOTYPE FUNCTIONS (MAY BE REPLACED BY USER)
===============================================================*/

/* Prototype of internal ewtSet function */
int arkEwtSet(N_Vector ycur, N_Vector weight, void *data);

/* Prototype of internal rwtSet function */
int arkRwtSet(N_Vector ycur, N_Vector weight, void *data);

/* Prototype of internal errHandler function */
void arkErrHandler(int error_code, const char *module, 
		   const char *function, char *msg, void *data);

/* Prototype of internal explicit stability estimation function */
int arkExpStab(N_Vector y, realtype t, realtype *hstab, void *user_data);

/*===============================================================
   HIGH LEVEL ERROR HANDLER, USED THROUGHOUT ARKODE
===============================================================*/

void arkProcessError(ARKodeMem ark_mem, int error_code, 
		     const char *module, const char *fname, 
		     const char *msgfmt, ...);

/*===============================================================
   ARKODE ERROR MESSAGES
===============================================================*/

#if defined(SUNDIALS_EXTENDED_PRECISION)

#define MSG_TIME        "t = %Lg"
#define MSG_TIME_H      "t = %Lg and h = %Lg"
#define MSG_TIME_INT    "t = %Lg is not between tcur - hold = %Lg and tcur = %Lg."
#define MSG_TIME_TOUT   "tout = %Lg"
#define MSG_TIME_TSTOP  "tstop = %Lg"

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define MSG_TIME        "t = %lg"
#define MSG_TIME_H      "t = %lg and h = %lg"
#define MSG_TIME_INT    "t = %lg is not between tcur - hold = %lg and tcur = %lg."
#define MSG_TIME_TOUT   "tout = %lg"
#define MSG_TIME_TSTOP  "tstop = %lg"

#else

#define MSG_TIME        "t = %g"
#define MSG_TIME_H      "t = %g and h = %g"
#define MSG_TIME_INT    "t = %g is not between tcur - hold = %g and tcur = %g."
#define MSG_TIME_TOUT   "tout = %g"
#define MSG_TIME_TSTOP  "tstop = %g"

#endif

/* Initialization and I/O error messages */
#define MSGARK_NO_MEM        "arkode_mem = NULL illegal."
#define MSGARK_ARKMEM_FAIL   "Allocation of arkode_mem failed."
#define MSGARK_MEM_FAIL      "A memory request failed."
#define MSGARK_NO_MALLOC     "Attempt to call before ARKodeInit."
#define MSGARK_NEG_MAXORD    "maxord <= 0 illegal."
#define MSGARK_BAD_MAXORD    "Illegal attempt to increase maximum method order."
#define MSGARK_NEG_HMIN      "hmin < 0 illegal."
#define MSGARK_NEG_HMAX      "hmax < 0 illegal."
#define MSGARK_BAD_HMIN_HMAX "Inconsistent step size limits: hmin > hmax."
#define MSGARK_BAD_RELTOL    "reltol < 0 illegal."
#define MSGARK_BAD_ABSTOL    "abstol has negative component(s) (illegal)."
#define MSGARK_NULL_ABSTOL   "abstol = NULL illegal."
#define MSGARK_BAD_RABSTOL   "rabstol has negative component(s) (illegal)."
#define MSGARK_NULL_RABSTOL  "rabstol = NULL illegal."
#define MSGARK_NULL_Y0       "y0 = NULL illegal."
#define MSGARK_NULL_F        "Must specify at least one of fe, fi (both NULL)."
#define MSGARK_NULL_G        "g = NULL illegal."
#define MSGARK_BAD_NVECTOR   "A required vector operation is not implemented."
#define MSGARK_BAD_K         "Illegal value for k."
#define MSGARK_NULL_DKY      "dky = NULL illegal."
#define MSGARK_BAD_T         "Illegal value for t." MSG_TIME_INT
#define MSGARK_NO_ROOT       "Rootfinding was not initialized."

/* ARKode Error Messages */
#define MSGARK_LSOLVE_NULL    "The linear solver's solve routine is NULL."
#define MSGARK_YOUT_NULL      "yout = NULL illegal."
#define MSGARK_TRET_NULL      "tret = NULL illegal."
#define MSGARK_BAD_EWT        "Initial ewt has component(s) equal to zero (illegal)."
#define MSGARK_EWT_NOW_BAD    "At " MSG_TIME ", a component of ewt has become <= 0."
#define MSGARK_BAD_RWT        "Initial rwt has component(s) equal to zero (illegal)."
#define MSGARK_RWT_NOW_BAD    "At " MSG_TIME ", a component of rwt has become <= 0."
#define MSGARK_BAD_ITASK      "Illegal value for itask."
#define MSGARK_BAD_H0         "h0 and tout - t0 inconsistent."
#define MSGARK_BAD_TOUT       "Trouble interpolating at " MSG_TIME_TOUT ". tout too far back in direction of integration"
#define MSGARK_EWT_FAIL       "The user-provide EwtSet function failed."
#define MSGARK_EWT_NOW_FAIL   "At " MSG_TIME ", the user-provide EwtSet function failed."
#define MSGARK_RWT_FAIL       "The user-provide RwtSet function failed."
#define MSGARK_RWT_NOW_FAIL   "At " MSG_TIME ", the user-provide RwtSet function failed."
#define MSGARK_LINIT_FAIL     "The linear solver's init routine failed."
#define MSGARK_LFREE_FAIL     "The linear solver's free routine failed."
#define MSGARK_HNIL_DONE      "The above warning has been issued mxhnil times and will not be issued again for this problem."
#define MSGARK_TOO_CLOSE      "tout too close to t0 to start integration."
#define MSGARK_MAX_STEPS      "At " MSG_TIME ", mxstep steps taken before reaching tout."
#define MSGARK_TOO_MUCH_ACC   "At " MSG_TIME ", too much accuracy requested."
#define MSGARK_HNIL           "Internal " MSG_TIME_H " are such that t + h = t on the next step. The solver will continue anyway."
#define MSGARK_ERR_FAILS      "At " MSG_TIME_H ", the error test failed repeatedly or with |h| = hmin."
#define MSGARK_CONV_FAILS     "At " MSG_TIME_H ", the solver convergence test failed repeatedly or with |h| = hmin."
#define MSGARK_SETUP_FAILED   "At " MSG_TIME ", the setup routine failed in an unrecoverable manner."
#define MSGARK_SOLVE_FAILED   "At " MSG_TIME ", the solve routine failed in an unrecoverable manner."
#define MSGARK_RHSFUNC_FAILED "At " MSG_TIME ", the right-hand side routine failed in an unrecoverable manner."
#define MSGARK_RHSFUNC_UNREC  "At " MSG_TIME ", the right-hand side failed in a recoverable manner, but no recovery is possible."
#define MSGARK_RHSFUNC_REPTD  "At " MSG_TIME " repeated recoverable right-hand side function errors."
#define MSGARK_RHSFUNC_FIRST  "The right-hand side routine failed at the first call."
#define MSGARK_RTFUNC_FAILED  "At " MSG_TIME ", the rootfinding routine failed in an unrecoverable manner."
#define MSGARK_CLOSE_ROOTS    "Root found at and very near " MSG_TIME "."
#define MSGARK_BAD_TSTOP      "The value " MSG_TIME_TSTOP " is behind current " MSG_TIME " in the direction of integration."
#define MSGARK_INACTIVE_ROOTS "At the end of the first step, there are still some root functions identically 0. This warning will not be issued again."
#define MSGARK_MISSING_FE     "Cannot specify that method is explicit without providing a function pointer to fe(t,y)."
#define MSGARK_MISSING_FI     "Cannot specify that method is explicit without providing a function pointer to fe(t,y)."
#define MSGARK_MISSING_F      "Cannot specify that method is ImEx without providing function pointers to fi(t,y) and fe(t,y)."
#define MSGARK_RESIZE_FAIL    "Error in user-supplied resize() function."
#define MSGARK_MASSINIT_FAIL  "The mass matrix solver's init routine failed."
#define MSGARK_MASSSETUP_FAIL "The mass matrix solver's setup routine failed."
#define MSGARK_MASSSOLVE_NULL "The mass matrix solver's solve routine is NULL."
#define MSGARK_MASSSOLVE_FAIL "The mass matrix solver failed."
#define MSGARK_MASSFREE_FAIL  "The mass matrixsolver's free routine failed."

#ifdef __cplusplus
}
#endif

#endif
