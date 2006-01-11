/*
 * -----------------------------------------------------------------
 * $Revision: 1.18 $
 * $Date: 2006-01-11 21:13:54 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Allan G. Taylor, Alan C. Hindmarsh, Radu Serban,
 *                and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/ida/LICENSE.
 * -----------------------------------------------------------------
 * This is the header file (private version) for the main IDA solver.
 * -----------------------------------------------------------------
 */

#ifndef _IDA_IMPL_H
#define _IDA_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "ida.h"

/* Prototype of internal ewtSet function */

int IDAEwtSet(N_Vector ycur, N_Vector weight, void *e_data);

/* Basic IDA constants */

#define HMAX_INV_DEFAULT RCONST(0.0) /* hmax_inv default value          */
#define MAXORD_DEFAULT   5           /* maxord default value            */
#define MXORDP1          6           /* max. number of N_Vectors in phi */
#define MXSTEP_DEFAULT   500         /* mxstep default value            */

/******************************************************************
 * Types : struct IDAMemRec, IDAMem                               *
 *----------------------------------------------------------------*
 * The type IDAMem is type pointer to struct IDAMemRec. This      *
 * structure contains fields to keep track of problem state.      *
 *                                                                *
 ******************************************************************/

typedef struct IDAMemRec {

  realtype ida_uround;    /* machine unit roundoff */

  /* Problem Specification Data */

  IDAResFn       ida_res;            /* F(t,y(t),y'(t))=0; the function F  */
  void          *ida_rdata;          /* user pointer passed to res         */

  int            ida_itol;           /* itol = IDA_SS, IDA_SV or IDA_WF    */
  realtype       ida_rtol;           /* relative tolerance                 */
  realtype       ida_Satol;          /* scalar absolute tolerance          */  
  N_Vector       ida_Vatol;          /* vector absolute tolerance          */  
  IDAEwtFn       ida_efun;           /* function to set ewt                */
  void          *ida_edata;          /* user pointer passed to efun        */
  

  booleantype    ida_setupNonNull;   /* Does setup do something?           */
  booleantype    ida_constraintsSet; /* constraints vector present: 
                                        do constraints calc                */
  booleantype    ida_suppressalg;    /* true means suppress algebraic vars
                                        in local error tests               */

  /* Divided differences array and associated minor arrays */

  N_Vector ida_phi[MXORDP1];   /* phi = (maxord+1) arrays of divided differences */

  realtype ida_psi[MXORDP1];   /* differences in t (sums of recent step sizes)   */
  realtype ida_alpha[MXORDP1]; /* ratios of current stepsize to psi values       */
  realtype ida_beta[MXORDP1];  /* ratios of current to previous product of psi's */
  realtype ida_sigma[MXORDP1]; /* product successive alpha values and factorial  */
  realtype ida_gamma[MXORDP1]; /* sum of reciprocals of psi values               */

  /* N_Vectors */

  N_Vector ida_ewt;         /* error weight vector                           */
  N_Vector ida_yy;          /* work space for y vector (= user's yret)       */
  N_Vector ida_yp;          /* work space for y' vector (= user's ypret)     */
  N_Vector ida_delta;       /* residual vector                               */
  N_Vector ida_id;          /* bit vector for diff./algebraic components     */
  N_Vector ida_constraints; /* vector of inequality constraint options       */
  N_Vector ida_savres;      /* saved residual vector (= tempv1)              */
  N_Vector ida_ee;          /* accumulated corrections to y vector, but
			       set equal to estimated local errors upon
			       successful return                             */
  N_Vector ida_mm;          /* mask vector in constraints tests (= tempv2)   */
  N_Vector ida_tempv1;      /* work space vector                             */
  N_Vector ida_tempv2;      /* work space vector                             */
  N_Vector ida_ynew;        /* work vector for y in IDACalcIC (= tempv2)     */
  N_Vector ida_ypnew;       /* work vector for yp in IDACalcIC (= ee)        */
  N_Vector ida_delnew;      /* work vector for delta in IDACalcIC (= phi[2]) */
  N_Vector ida_dtemp;       /* work vector in IDACalcIC (= phi[3])           */

  /* Variables for use by IDACalcIC*/

  realtype ida_t0;          /* initial t                                     */
  N_Vector ida_yy0;         /* initial y vector (user-supplied).             */
  N_Vector ida_yp0;         /* initial y' vector (user-supplied).            */

  int ida_icopt;            /* IC calculation user option                    */
  booleantype ida_lsoff;    /* IC calculation linesearch turnoff option      */
  int ida_maxnh;            /* max. number of h tries in IC calculation      */
  int ida_maxnj;            /* max. number of J tries in IC calculation      */
  int ida_maxnit;           /* max. number of Netwon iterations in IC calc.  */
  int ida_nbacktr;          /* number of IC linesearch backtrack operations  */
  int ida_sysindex;         /* computed system index (0 or 1)                */
  realtype ida_epiccon;     /* IC nonlinear convergence test constant        */
  realtype ida_steptol;     /* minimum Newton step size in IC calculation    */
  realtype ida_tscale;      /* time scale factor = abs(tout1 - t0)           */

  /* Tstop information */

  booleantype ida_istop;
  booleantype ida_tstopset;
  realtype ida_tstop;

  /* Step Data */

  int ida_kk;        /* current BDF method order                              */
  int ida_kused;     /* method order used on last successful step             */
  int ida_knew;      /* order for next step from order decrease decision      */
  int ida_phase;     /* flag to trigger step doubling in first few steps      */
  int ida_ns;        /* counts steps at fixed stepsize and order              */

  realtype ida_hin;      /* initial step                                      */
  realtype ida_h0u;      /* actual initial stepsize                           */
  realtype ida_hh;       /* current step size h                               */
  realtype ida_hused;    /* step size used on last successful step            */
  realtype ida_rr;       /* rr = hnext / hused                                */
  realtype ida_tn;       /* current internal value of t                       */
  realtype ida_tretlast; /* value of tret previously returned by IDASolve     */
  realtype ida_cj;       /* current value of scalar (-alphas/hh) in Jacobian  */
  realtype ida_cjlast;   /* cj value saved from last successful step          */
  realtype ida_cjold;    /* cj value saved from last call to lsetup           */
  realtype ida_cjratio;  /* ratio of cj values: cj/cjold                      */
  realtype ida_ss;       /* scalar used in Newton iteration convergence test  */
  realtype ida_epsNewt;  /* test constant in Newton convergence test          */
  realtype ida_epcon;    /* coeficient of the Newton covergence test          */
  realtype ida_toldel;   /* tolerance in direct test on Newton corrections    */

 /* Limits */

  int ida_maxncf;        /* max numer of convergence failures                 */
  int ida_maxcor;        /* max number of Newton corrections                  */
  int ida_maxnef;        /* max number of error test failures                 */

  int ida_maxord;        /* max value of method order k:                      */
  int ida_maxord_alloc;  /* value of maxord used when allocating memory       */
  long int ida_mxstep;   /* max number of internal steps for one user call    */
  realtype ida_hmax_inv; /* inverse of max. step size hmax (default = 0.0)    */

  /* Counters */

  long int ida_nst;      /* number of internal steps taken                    */
  long int ida_nre;      /* number of function (res) calls                    */
  long int ida_ncfn;     /* number of corrector convergence failures          */
  long int ida_netf;     /* number of error test failures                     */
  long int ida_nni;      /* number of Newton iterations performed             */
  long int ida_nsetups;  /* number of lsetup calls                            */

  /* Space requirements for IDA */

  long int ida_lrw1;     /* no. of realtype words in 1 N_Vector               */
  long int ida_liw1;     /* no. of integer words in 1 N_Vector                */
  long int ida_lrw;      /* number of realtype words in IDA work vectors      */
  long int ida_liw;      /* no. of integer words in IDA work vectors          */

  realtype ida_tolsf;    /* tolerance scale factor (saved value)              */

  FILE *ida_errfp;       /* IDA error messages are sent to errfp              */

  /* Flags to verify correct calling sequence */

  booleantype ida_SetupDone;     /* set to FALSE by IDAMalloc and IDAReInit   */
                                 /* set to TRUE by IDACalcIC or IDASolve      */

  booleantype ida_VatolMallocDone;
  booleantype ida_constraintsMallocDone;
  booleantype ida_idMallocDone;

  booleantype ida_MallocDone;    /* set to FALSE by IDACreate                 */
                                 /* set to TRUE by IDAMAlloc                  */
                                 /* tested by IDAReInit and IDASolve          */

  /* Linear Solver Data */

  /* Linear Solver functions to be called */

  int (*ida_linit)(struct IDAMemRec *idamem);

  int (*ida_lsetup)(struct IDAMemRec *idamem, N_Vector yyp, 
                    N_Vector ypp, N_Vector resp, 
                    N_Vector tempv1, N_Vector tempv2, N_Vector tempv3); 

  int (*ida_lsolve)(struct IDAMemRec *idamem, N_Vector b, N_Vector weight,
                    N_Vector ycur, N_Vector ypcur, N_Vector rescur);

  int (*ida_lperf)(struct IDAMemRec *idamem, int perftask);

  int (*ida_lfree)(struct IDAMemRec *idamem);

  /* Linear Solver specific memory */

  void *ida_lmem;           

  /* Flag to indicate successful ida_linit call */

  booleantype ida_linitOK;

  /* Rootfinding Data */

  IDARootFn ida_gfun;    /* Function g for roots sought                     */
  int ida_nrtfn;         /* number of components of g                       */
  void *ida_g_data;      /* pointer to user data for g                      */
  int *ida_iroots;       /* int array for root information                  */
  realtype ida_tlo;      /* nearest endpoint of interval in root search     */
  realtype ida_thi;      /* farthest endpoint of interval in root search    */
  realtype ida_trout;    /* t return value from rootfinder routine          */
  realtype *ida_glo;     /* saved array of g values at t = tlo              */
  realtype *ida_ghi;     /* saved array of g values at t = thi              */
  realtype *ida_grout;   /* array of g values at t = trout                  */
  realtype ida_toutc;    /* copy of tout (if NORMAL mode)                   */
  realtype ida_ttol;     /* tolerance on root location                      */
  int ida_taskc;         /* copy of parameter task                          */
  int ida_irfnd;         /* flag showing whether last step had a root       */
  long int ida_nge;      /* counter for g evaluations                       */


} *IDAMem;


/*
 *----------------------------------------------------------------
 * IDA Error Messages
 *----------------------------------------------------------------
 */

#if defined(SUNDIALS_EXTENDED_PRECISION)

#define MSG_TIME "at t = %Lg, "
#define MSG_TIME_H "at t = %Lg and h = %Lg, "
#define MSG_TIME_INT "t is not between tcur - hu = %Lg and tcur = %Lg.\n\n"
#define MSG_TIME_TOUT "tout = %Lg"

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define MSG_TIME "at t = %lg, "
#define MSG_TIME_H "at t = %lg and h = %lg, "
#define MSG_TIME_INT "t is not between tcur - hu = %lg and tcur = %lg.\n\n"
#define MSG_TIME_TOUT "tout = %lg"

#else

#define MSG_TIME "at t = %g, "
#define MSG_TIME_H "at t = %g and h = %g, "
#define MSG_TIME_INT "t is not between tcur - hu = %g and tcur = %g.\n\n"
#define MSG_TIME_TOUT "tout = %g"

#endif

/* IDACreate error messages */

#define MSG_IDAMEM_FAIL    "IDACreate-- allocation of ida_mem failed. \n\n"

/* IDAMalloc/IDAReInit error messages */

#define _IDAM_             "IDAMalloc/IDAReInit-- "

#define MSG_IDAM_NO_MEM    _IDAM_ "ida_mem = NULL illegal.\n\n"

#define MSG_Y0_NULL        _IDAM_ "y0 = NULL illegal.\n\n"
#define MSG_YP0_NULL       _IDAM_ "yp0 = NULL illegal.\n\n"

#define MSG_BAD_ITOL       _IDAM_ "itol has an illegal value.\n"

#define MSG_RES_NULL       _IDAM_ "res = NULL illegal.\n\n"

#define MSG_BAD_RTOL       _IDAM_ "*reltol < 0 illegal.\n\n"

#define MSG_ATOL_NULL      _IDAM_ "abstol = NULL illegal.\n\n"

#define MSG_BAD_ATOL       _IDAM_ "some abstol component < 0.0 illegal.\n\n"

#define MSG_BAD_NVECTOR    _IDAM_ "a required vector operation is not implemented.\n\n"

#define MSG_MEM_FAIL       _IDAM_ "a memory request failed.\n\n"

#define MSG_REI_NO_MALLOC  "IDAReInit-- attempt to call before IDAMalloc. \n\n"

/* IDAInitialSetup error messages -- called from IDACalcIC or IDASolve */

#define _IDAIS_              "Initial setup-- "

#define MSG_MISSING_ID      _IDAIS_ "id = NULL but suppressalg option on.\n\n"

#define MSG_NO_EFUN         _IDAIS_ "itol = IDA_WF but no EwtSet function was provided.\n\n"

#define MSG_FAIL_EWT        _IDAIS_ "The user-provide EwtSet function failed.\n\n"

#define MSG_BAD_EWT         _IDAIS_ "some initial ewt component = 0.0 illegal.\n\n"

#define MSG_Y0_FAIL_CONSTR  _IDAIS_ "y0 fails to satisfy constraints.\n\n"

#define MSG_LSOLVE_NULL     _IDAIS_ "the linear solver's solve routine is NULL.\n\n"

#define MSG_LINIT_FAIL      _IDAIS_ "the linear solver's init routine failed.\n\n"

/* IDARootInit Error Messages */

#define _IDART_ "IDARootInit-- "

#define MSG_ROOT_NO_MEM _IDART_ "ida_mem = NULL illegal.\n\n"

#define MSG_ROOT_MEM_FAIL _IDART_ "A memory request failed.\n\n"

#define MSG_ROOT_FUNC_NULL _IDART_ "g = NULL illegal.\n\n"

/* IDASolve error messages */

#define _IDASLV_             "IDASolve-- "

#define MSG_IDA_NO_MEM     _IDASLV_ "ida_mem = NULL illegal.\n\n"

#define MSG_NO_MALLOC      _IDASLV_ "attempt to call before IDAMalloc. \n\n"
 
#define MSG_BAD_HINIT      _IDASLV_ "initial step is not towards tout.\n\n"

#define MSG_BAD_TOUT1      _IDASLV_ "trouble interpolating at " MSG_TIME_TOUT ".\n"
#define MSG_BAD_TOUT2      "tout too far back in direction of integration.\n\n"
#define MSG_BAD_TOUT       MSG_BAD_TOUT1 MSG_BAD_TOUT2

#define MSG_BAD_TSTOP      _IDASLV_ MSG_TIME "tstop is behind.\n\n"

#define MSG_MAX_STEPS      _IDASLV_ MSG_TIME "maximum number of steps reached.\n\n"

#define MSG_EWT_NOW_FAIL   _IDASLV_ MSG_TIME "the user-provide EwtSet function failed.\n\n"

#define MSG_EWT_NOW_BAD    _IDASLV_ MSG_TIME "some ewt component has become <= 0.0.\n\n"

#define MSG_TOO_MUCH_ACC   _IDASLV_ MSG_TIME "too much accuracy requested.\n\n"

#define MSG_ERR_FAILS1     "the error test\nfailed repeatedly or with |h| = hmin.\n\n"
#define MSG_ERR_FAILS      _IDASLV_ MSG_TIME_H MSG_ERR_FAILS1

#define MSG_CONV_FAILS1    "the corrector convergence\nfailed repeatedly.\n\n"
#define MSG_CONV_FAILS     _IDASLV_ MSG_TIME_H MSG_CONV_FAILS1

#define MSG_SETUP_FAILED1  "the linear solver setup failed unrecoverably.\n\n"
#define MSG_SETUP_FAILED   _IDASLV_ MSG_TIME MSG_SETUP_FAILED1

#define MSG_SOLVE_FAILED1  "the linear solver solve failed unrecoverably.\n\n"
#define MSG_SOLVE_FAILED   _IDASLV_ MSG_TIME MSG_SOLVE_FAILED1

#define MSG_TOO_CLOSE      _IDASLV_ "tout too close to t0 to start integration.\n\n"

#define MSG_BAD_INIT_ROOT  _IDASLV_ "Root found at and very near initial t.\n\n"

#define MSG_YRET_NULL      _IDASLV_ "yret = NULL illegal.\n\n"
#define MSG_YPRET_NULL     _IDASLV_ "ypret = NULL illegal.\n\n"
#define MSG_TRET_NULL      _IDASLV_ "tret = NULL illegal.\n\n"

#define MSG_BAD_ITASK      _IDASLV_ "itask has an illegal value.\n\n"

#define MSG_NO_TSTOP       _IDASLV_ "tstop not set for this itask. \n\n"

#define MSG_CLOSE_ROOTS    _IDASLV_ "Root found at and very near " MSG_TIME ".\n\n"

#define MSG_REP_RES_ERR1   "repeated recoverable residual errors.\n\n"
#define MSG_REP_RES_ERR    _IDASLV_ MSG_TIME MSG_REP_RES_ERR1

#define MSG_RES_NONRECOV1  "the residual function failed unrecoverably. \n\n"
#define MSG_RES_NONRECOV   _IDASLV_ MSG_TIME MSG_RES_NONRECOV1

#define MSG_FAILED_CONSTR1 "unable to satisfy inequality constraints. \n\n"
#define MSG_FAILED_CONSTR  _IDASLV_ MSG_TIME MSG_FAILED_CONSTR1

/* IDACalcIC error messages */

#define _IDAIC_            "IDACalcIC-- "

#define MSG_IC_NO_MEM      _IDAIC_ "IDA_mem = NULL illegal.\n\n"

#define MSG_IC_NO_MALLOC   _IDAIC_ "attempt to call before IDAMalloc. \n\n"
 
#define MSG_IC_Y0_NULL     _IDAIC_ "y0 = NULL illegal.\n\n"
#define MSG_IC_YP0_NULL    _IDAIC_ "yp0 = NULL illegal.\n\n"

#define MSG_IC_BAD_ICOPT   _IDAIC_ "icopt has an illegal value.\n\n"

#define MSG_IC_MISSING_ID  _IDAIC_ "id = NULL conflicts with icopt.\n\n"

#define MSG_IC_BAD_ID      _IDAIC_ "id has illegal values.\n\n"

#define MSG_IC_TOO_CLOSE1  _IDAIC_ "tout1 too close to t0 to attempt "
#define MSG_IC_TOO_CLOSE2  "initial condition calculation.\n\n"
#define MSG_IC_TOO_CLOSE   MSG_IC_TOO_CLOSE1 MSG_IC_TOO_CLOSE2

#define MSG_IC_BAD_EWT     _IDAIC_ "some ewt component = 0.0 illegal.\n\n"

#define MSG_IC_RES_NONR1   "the residual function failed unrecoverably. \n\n"
#define MSG_IC_RES_NONREC  _IDAIC_ MSG_IC_RES_NONR1

#define MSG_IC_RES_FAIL1   "the residual function failed at the first call. \n\n"
#define MSG_IC_RES_FAIL    _IDAIC_ MSG_IC_RES_FAIL1

#define MSG_IC_SETUP_FAIL1 "the linear solver setup failed unrecoverably.\n\n"
#define MSG_IC_SETUP_FAIL  _IDAIC_ MSG_IC_SETUP_FAIL1

#define MSG_IC_SOLVE_FAIL1 "the linear solver solve failed unrecoverably.\n\n"
#define MSG_IC_SOLVE_FAIL  _IDAIC_ MSG_IC_SOLVE_FAIL1

#define MSG_IC_NO_RECOV1   _IDAIC_ "The residual routine or the linear"
#define MSG_IC_NO_RECOV2   " setup or solve routine had a recoverable"
#define MSG_IC_NO_RECOV3   " error, but IDACalcIC was unable to recover.\n\n"
#define MSG_IC_NO_RECOVERY MSG_IC_NO_RECOV1 MSG_IC_NO_RECOV2 MSG_IC_NO_RECOV3

#define MSG_IC_FAIL_CON1   "Unable to satisfy the inequality constraints.\n\n"
#define MSG_IC_FAIL_CONSTR _IDAIC_ MSG_IC_FAIL_CON1

#define MSG_IC_FAILED_LS1  "the linesearch algorithm failed with too small a step.\n\n"
#define MSG_IC_FAILED_LINS _IDAIC_ MSG_IC_FAILED_LS1

#define MSG_IC_CONV_FAIL1  "Newton/Linesearch algorithm failed to converge.\n\n"
#define MSG_IC_CONV_FAILED _IDAIC_ MSG_IC_CONV_FAIL1

/* IDASet* error messages */

#define MSG_IDAS_NO_MEM      "IDASet*-- ida_mem = NULL illegal. \n\n"

#define MSG_IDAS_NEG_MAXORD  "IDASetMaxOrd-- maxord<=0 illegal. \n\n"

#define MSG_IDAS_BAD_MAXORD  "IDASetMaxOrd-- illegal to increase maximum order.\n\n"

#define MSG_IDAS_NEG_MXSTEPS "IDASetMaxNumSteps-- mxsteps < 0 illegal. \n\n"

#define MSG_IDAS_NEG_HMAX    "IDASetMaxStep-- hmax < 0 illegal. \n\n"

#define MSG_IDAS_NEG_EPCON   "IDASetNonlinConvCoef-- epcon < 0.0 illegal. \n\n"

#define MSG_IDAS_BAD_NVECTOR  "IDASetConstraints-- a required vector operation is not implemented.\n\n"

#define MSG_IDAS_BAD_CONSTRAINTS "IDASetConstraints-- illegal values in constraints vector.\n\n"

#define MSG_IDAS_NO_MALLOC   "IDASetTolerances-- Attemp to call before IDAMalloc.\n\n"

#define MSG_IDAS_BAD_ITOL    "IDASetTolerances-- itol has an illegal value.\n\n"

#define MSG_IDAS_BAD_RTOL    "IDASetTolerances-- *rtol < 0 illegal.\n\n"

#define MSG_IDAS_ATOL_NULL   "IDASetTolerances-- atol = NULL illegal.\n\n"

#define MSG_IDAS_BAD_ATOL    "IDASetTolerances-- some atol component < 0.0 illegal.\n\n"

#define MSG_IDAS_BAD_EPICCON "IDASetNonlinConvCoefIC-- epiccon < 0.0 illegal.\n\n"

#define MSG_IDAS_BAD_MAXNH   "IDASetMaxNumStepsIC-- maxnh < 0 illegal.\n\n"

#define MSG_IDAS_BAD_MAXNJ   "IDASetMaxNumJacsIC-- maxnj < 0 illegal.\n\n"

#define MSG_IDAS_BAD_MAXNIT  "IDASetMaxNumItersIC-- maxnit < 0 illegal.\n\n"

#define MSG_IDAS_BAD_STEPTOL "IDASetLineSearchOffIC-- steptol < 0.0 illegal.\n\n"

/* IDAGet* Error Messages */

#define MSG_IDAG_NO_MEM    "IDAGet*-- ida_mem = NULL illegal. \n\n"

#define MSG_IDAG_BAD_T1    "IDAGetSolution-- "
#define MSG_IDAG_BAD_T     MSG_IDAG_BAD_T1 MSG_TIME MSG_TIME_INT

#define MSG_IDAG_EWT_BAD   "IDAGetErrWeights--  ewt has component(s) equal to zero.\n\n"

#ifdef __cplusplus
}
#endif

#endif
