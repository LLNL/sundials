/*
 * -----------------------------------------------------------------
 * $Revision: 1.9 $
 * $Date: 2004-11-24 23:03:44 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/idas/LICENSE
 * -----------------------------------------------------------------
 * Implementation header file for the main IDAS integrator.
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _IDAS_IMPL_H
#define _IDAS_IMPL_H

#include <stdio.h>

#include "idas.h"

#include "nvector.h"
#include "sundialstypes.h"

/* 
 * =================================================================
 *   M A I N    I N T E G R A T O R    M E M O R Y    B L O C K
 * =================================================================
 */

/* Basic IDAS constants */

#define MXORDP1 6 /* max. number of vectors kept in the phi array */

/*
 * ifS:   Type of the function returning the sensitivity residual    
 *        ifS can be either IDA_ALLSENS if the function    
 *        (of type SensResFn) returns residuals for all    
 *        sensitivity systems at once, or IDA_ONESENS if the function 
 *        (of type SensRes1Fn) returns the residual of one 
 *        sensitivity system at a time.                           
 *                                                                
 */

#define IDA_ONESENS 1
#define IDA_ALLSENS 2

/*
 *----------------------------------------------------------------
 * Types : struct IDAMemRec, IDAMem                               
 *----------------------------------------------------------------
 * The type IDAMem is type pointer to struct IDAMemRec. This      
 * structure contains fields to keep track of problem state.      
 *                                                                
 *----------------------------------------------------------------
 */

typedef struct IDAMemRec {

  realtype ida_uround;               /* machine unit roundoff              */

  /*--------------------------
    Problem Specification Data 
    --------------------------*/

  IDAResFn       ida_res;            /* F(t,y(t),y'(t))=0; the function F  */
  void          *ida_rdata;          /* user pointer passed to res         */
  int            ida_itol;           /* itol = SS or SV                    */
  realtype      *ida_reltol;         /* ptr to relative tolerance          */
  void          *ida_abstol;         /* ptr to absolute tolerance          */  
  booleantype    ida_setupNonNull;   /* Does setup do something?           */
  booleantype    ida_constraintsSet; /* constraints vector present: 
                                        do constraints calc                */
  booleantype    ida_suppressalg;    /* true means suppress algebraic vars
                                        in local error tests               */

  /*-----------------------
    Quadrature Related Data 
    -----------------------*/

  booleantype    ida_quad;
  IDAQuadRhsFn   ida_rhsQ;
  int            ida_itolQ;
  realtype      *ida_reltolQ;
  void          *ida_abstolQ;
  booleantype    ida_errconQ;
  void          *ida_rdataQ;

  /*------------------------
    Sensitivity Related Data
    ------------------------*/

  booleantype    ida_sensi;
  int            ida_Ns;
  IDASensResFn   ida_resS;
  IDASensRes1Fn  ida_resS1;
  booleantype    ida_resSDQ;
  int            ida_iresS;
  int            ida_ism;
  realtype      *ida_p;
  realtype      *ida_pbar;
  int           *ida_plist;
  realtype       ida_rhomax;
  booleantype    ida_errconS;
  void          *ida_rdataS;

  int            ida_itolS;
  realtype      *ida_reltolS;
  void          *ida_abstolS;
  booleantype    ida_testSensTol;    /* flag to indicate if sensi. tolerances
                                        must be checked now                    */
  booleantype    ida_setSensTol;     /* flag to indicate if sensi. tolerances 
                                        must be set now                        */
  booleantype    ida_atolSallocated; /* TRUE if IDAS has allocated space for
                                        sensitivity absolute tolerances        */

  booleantype    ida_stgr1alloc;     /* ncfS1,ncfnS1,and nniS1 allocated by IDAS? */

  /*-----------------------------------------------
    Divided differences array and associated arrays
    -----------------------------------------------*/

  int ida_maxcol;              /* Actual number of phi arrays allocated          */
  N_Vector ida_phi[MXORDP1];   /* phi = (maxord+1) arrays of divided differences */

  realtype ida_psi[MXORDP1];   /* differences in t (sums of recent step sizes)   */
  realtype ida_alpha[MXORDP1]; /* ratios of current stepsize to psi values       */
  realtype ida_beta[MXORDP1];  /* ratios of current to previous product of psi's */
  realtype ida_sigma[MXORDP1]; /* product successive alpha values and factorial  */
  realtype ida_gamma[MXORDP1]; /* sum of reciprocals of psi values               */

  /*-------------------------
    N_Vectors for integration
    -------------------------*/

  N_Vector ida_ewt;         /* error weight vector                           */
  N_Vector ida_y0;          /* initial y vector (user-supplied)              */
  N_Vector ida_yp0;         /* initial y' vector (user-supplied)             */
  N_Vector ida_yy;          /* work space for y vector (= user's yret)       */
  N_Vector ida_yp;          /* work space for y' vector (= user's ypret)     */
  N_Vector ida_delta;       /* residual vector                               */
  N_Vector ida_id;          /* bit vector for diff./algebraic components     */
  N_Vector ida_constraints; /* vector of inequality constraint options       */
  N_Vector ida_savres;      /* saved residual vector (= tempv1)              */
  N_Vector ida_ee;          /* accumulated corrections to y                  */
  N_Vector ida_mm;          /* mask vector in constraints tests (= tempv2)   */
  N_Vector ida_tempv1;      /* work space vector                             */
  N_Vector ida_tempv2;      /* work space vector                             */
  N_Vector ida_ynew;        /* work vector for y in IDACalcIC (= tempv2)     */
  N_Vector ida_ypnew;       /* work vector for yp in IDACalcIC (= ee)        */
  N_Vector ida_delnew;      /* work vector for delta in IDACalcIC (= phi[2]) */
  N_Vector ida_dtemp;       /* work vector in IDACalcIC (= phi[3])           */

  /*----------------------------
    Quadrature Related N_Vectors 
    ----------------------------*/

  N_Vector ida_phiQ[MXORDP1];
  N_Vector ida_yyQ;
  N_Vector ida_ypQ;
  N_Vector ida_ewtQ;
  N_Vector ida_eeQ;

  /*---------------------------
    Sensitivity Related Vectors 
    ---------------------------*/

  N_Vector *ida_phiS[MXORDP1];
  N_Vector *ida_ewtS;

  N_Vector *ida_yS0;        /* initial yS vector (user-supplied)            */
  N_Vector *ida_ypS0;       /* initial yS' vector (user-supplied)           */

  N_Vector *ida_eeS;        /* cumulative sensitivity corrections           */

  N_Vector *ida_yyS;        /* allocated and used for:                      */
  N_Vector *ida_ypS;        /*                 ism = SIMULTANEOUS           */
  N_Vector *ida_deltaS;     /*                 ism = STAGGERED              */

  N_Vector ida_yyS1;        /* allocated and used for:                      */
  N_Vector ida_ypS1;        /*                                              */
  N_Vector ida_deltaS1;     /*                 ism = STAGGERED1             */

  N_Vector ida_tmpS1;       /* work space vectors: tmpS1 = tempv1           */
  N_Vector ida_tmpS2;       /*                     tmpS2 = tempv2           */
  N_Vector ida_tmpS3;       /*                     tmpS3 = allocated        */    
  

  /*----------------------------
    Scalars for use by IDACalcIC
    ----------------------------*/

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



  /*-----------------
    Tstop information
  -------------------*/
  booleantype ida_tstopset;
  realtype ida_tstop;

  /*---------
    Step Data
    ---------*/

  int ida_kk;        /* current BDF method order                          */
  int ida_kused;     /* method order used on last successful step         */
  int ida_knew;      /* order for next step from order decrease decision  */
  int ida_phase;     /* flag to trigger step doubling in first few steps  */
  long int ida_ns;   /* counts steps at fixed stepsize and order          */

  realtype ida_hin;      /* initial step                                      */
  realtype ida_h0u;      /* actual initial stepsize                           */
  realtype ida_hh;       /* current step size h                               */
  realtype ida_hused;    /* step size used on last successful step            */
  realtype ida_rr;       /* rr = hnext / hused                                */
  realtype ida_tn;       /* current internal value of t                       */
  realtype ida_tretp;    /* value of tret previously returned by IDASolve     */
  realtype ida_cj;       /* current value of scalar (-alphas/hh) in Jacobian  */
  realtype ida_cjlast;   /* cj value saved from last successful step          */
  realtype ida_cjold;    /* cj value saved from last call to lsetup           */
  realtype ida_cjratio;  /* ratio of cj values: cj/cjold                      */
  realtype ida_ss;       /* scalar used in Newton iteration convergence test  */
  realtype ida_epsNewt;  /* test constant in Newton convergence test          */
  realtype ida_epcon;    /* Newton convergence test constant                  */
  realtype ida_toldel;   /* tolerance in direct test on Newton corrections    */
  
  realtype ida_ssS;      /* scalar ss for sensitivity variables (STAGGERED)   */
  realtype *ida_ssS1;    /* scalars ss for sensitivity variables (STAGGERED1) */

  /*------
    Limits
    ------*/

  int ida_maxncf;        /* max numer of convergence failures                 */
  int ida_maxcor;        /* max number of Newton corrections                  */
  int ida_maxnef;        /* max number of error test failures                 */

  int ida_maxord;        /* max value of method order k:                      */
  long int ida_mxstep;   /* max number of internal steps for one user call    */
  realtype ida_hmax_inv; /* inverse of max. step size hmax (default = 0.0)    */

  int ida_maxcorS;       /* max number of Newton corrections for sensitivity
                            systems (staggered method)                        */

  /*--------
    Counters
    --------*/

  long int ida_nst;      /* number of internal steps taken                    */

  long int ida_nre;      /* number of user function calls                     */
  long int ida_nrQe;
  long int ida_nrSe;
  long int ida_nreS;

  long int ida_ncfn;     /* number of corrector convergence failures          */
  long int ida_ncfnS;
  long int *ida_ncfnS1;

  long int ida_nni;      /* number of Newton iterations performed             */
  long int ida_nniS;
  long int *ida_nniS1;

  long int ida_netf;     /* number of error test failures                     */
  long int ida_netfQ;
  long int ida_netfS;
  long int *ida_netfS1;

  long int ida_nsetups;  /* number of lsetup calls                            */
  long int ida_nsetupsS;
  
  /*---------------------------
    Space requirements for IDAS
    ---------------------------*/

  long int ida_lrw1;     /* no. of realtype words in 1 N_Vector               */
  long int ida_liw1;     /* no. of integer words in 1 N_Vector                */
  long int ida_lrw1Q;
  long int ida_liw1Q;
  long int ida_lrw;      /* number of realtype words in IDAS work vectors     */
  long int ida_liw;      /* no. of integer words in IDAS work vectors         */

  realtype ida_tolsf;    /* tolerance scale factor (saved value)              */

  /*-------------------------------------------------------
    Flags to verify correct calling sequence
    Turned ON by IDAMalloc, IDASensMalloc, and IDAQuadMalloc 
    and read by IDAMalloc, IDASensReInit, and IDAQuadReInit
    --------------------------------------------------------*/

  booleantype ida_SetupDone;     /* set to FALSE by IDAMalloc and IDAReInit */
                                 /* set to TRUE by IDACalcIC or IDASolve    */

  booleantype ida_MallocDone;    /* set to FALSE by IDACreate               */
                                 /* set to TRUE by IDAMAlloc                */
                                 /* tested by IDAReInit and IDASolve        */

  booleantype ida_sensMallocDone;

  booleantype ida_quadMallocDone;
  
  /*-------------------------------------
    IDAS error messages are sent to errfp
    -------------------------------------*/

  FILE *ida_errfp;

  /*------------------
    Linear Solver Data
    ------------------*/

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
  
  /* Flag to request a call to the setup routine */
  
  booleantype ida_forceSetup;

  /* Flag to indicate successful ida_linit call */

  booleantype ida_linitOK;

} *IDAMem;


/*
 *----------------------------------------------------------------
 * IDAS Error Messages
 *----------------------------------------------------------------
 */

#ifdef SUNDIALS_EXTENDED_PRECISION

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

#define MSG_RELTOL_NULL    _IDAM_ "reltol = NULL illegal.\n\n"

#define MSG_BAD_RTOL       _IDAM_ "*reltol < 0 illegal.\n\n"

#define MSG_ATOL_NULL      _IDAM_ "abstol = NULL illegal.\n\n"

#define MSG_BAD_ATOL       _IDAM_ "some abstol component < 0.0 illegal.\n\n"

#define MSG_BAD_NVECTOR    _IDAM_ "a required vector operation is not implemented.\n\n"

#define MSG_MEM_FAIL       _IDAM_ "a memory request failed.\n\n"

#define MSG_REI_NO_MALLOC  "IDAReInit-- attempt to call before IDAMalloc. \n\n"

/* IDAQuadMalloc/IDAQuadReInit error messages */

#define _QIDAM_            "IDAQuadMalloc/IDAQuadReInit-- "

#define MSG_QIDAM_NO_MEM   _QIDAM_ "ida_mem = NULL illegal.\n\n"

#define MSG_QIDAM_MEM_FAIL _QIDAM_ "a memory request failed.\n\n"

#define MSG_BAD_RHSQ       _QIDAM_ "rhsQ = NULL illegal.\n\n"

#define MSG_QREI_QUAD1     "IDAQuadReInit-- illegal attempt to call before "
#define MSG_QREI_QUAD2     "calling IDAQuadMalloc.\n\n"
#define MSG_QREI_NO_QUAD   MSG_QREI_QUAD1 MSG_QREI_QUAD2

/* IDASensMalloc/ IDASensReInit error messages */

#define _SIDAM_            "IDASensMalloc/IDASensReInit-- "

#define MSG_SIDAM_NO_MEM   _SIDAM_ "ida_mem = NULL illegal.\n\n"

#define MSG_SIDAM_MEM_FAIL _SIDAM_ "a memory request failed.\n\n"

#define MSG_BAD_NS         _SIDAM_ "NS <= 0 illegal.\n\n"

#define MSG_P_NULL         _SIDAM_ "p = NULL illegal.\n\n"

#define MSG_YS0_NULL       _SIDAM_ "yS0 = NULL illegal.\n\n"
#define MSG_YPS0_NULL      _SIDAM_ "ypS0 = NULL illegal.\n\n"

#define MSG_BAD_ISM        _SIDAM_ "ism has an illegal value.\n"

#define MSG_SREI_SENSI1    "IDASensReInit-- illegal attempt to call before "
#define MSG_SREI_SENSI2    "calling IDASensMalloc.\n\n"
#define MSG_SREI_NO_SENSI  MSG_SREI_SENSI1 MSG_SREI_SENSI2

/* IDAInitialSetup error messages -- called from IDACalcIC or IDASolve */

#define _IDAIS_              "Initial setup-- "

#define MSG_MISSING_ID      _IDAIS_ "id = NULL but suppressalg option on.\n\n"

#define MSG_BAD_EWT         _IDAIS_ "some initial ewt component = 0.0 illegal.\n\n"

#define MSG_BAD_CONSTRAINTS _IDAIS_ "illegal values in constraints vector.\n\n"

#define MSG_Y0_FAIL_CONSTR  _IDAIS_ "y0 fails to satisfy constraints.\n\n"

#define MSG_NO_QUADTOL      _IDAIS_ "no quad tolerances set. Illegal for errconQ=TRUE.\n\n"

#define MSG_BAD_EWTQ        _IDAIS_ "some initial ewtQ component = 0.0 illegal.\n\n"

#define MSG_BAD_ISM_IRESS1  "illegal use of ism=IDA_STAGGERED1 "
#define MSG_BAD_ISM_IRESS2  "with the provided sensitivity residual function.\n\n"
#define MSG_BAD_ISM_IRESS   _IDAIS_ MSG_BAD_ISM_IRESS1 MSG_BAD_ISM_IRESS2 

#define MSG_BAD_EWTS        _IDAIS_ "some initial ewtS component = 0.0 illegal.\n\n"

#define MSG_LSOLVE_NULL     _IDAIS_ "the linear solver's solve routine is NULL.\n\n"

#define MSG_LINIT_FAIL      _IDAIS_ "the linear solver's init routine failed.\n\n"

/* IDASolve error messages */

#define _IDASLV_             "IDASolve-- "

#define MSG_BAD_RTOLS      _IDASLV_ "*reltolS < 0.0 illegal.\n\n"

#define MSG_BAD_ATOLS      _IDASLV_ "some abstolS component < 0.0 illegal.\n\n"  

#define MSG_BAD_PBAR       _IDASLV_ "some pbar component = 0.0 illegal.\n\n"

#define MSG_ATOLS_MEM_FAIL _IDASLV_ "a memory request failed (abstolS).\n\n"

#define MSG_IDA_NO_MEM     _IDASLV_ "ida_mem = NULL illegal.\n\n"

#define MSG_NO_MALLOC      _IDASLV_ "attempt to call before IDAMalloc.\n\n"
 
#define MSG_BAD_HINIT      _IDASLV_ "initial step is not towards tout.\n\n"

#define MSG_BAD_TOUT1      _IDASLV_ "trouble interpolating at " MSG_TIME_TOUT ".\n"
#define MSG_BAD_TOUT2      "tout too far back in direction of integration.\n\n"
#define MSG_BAD_TOUT       MSG_BAD_TOUT1 MSG_BAD_TOUT2

#define MSG_BAD_TSTOP      _IDASLV_ MSG_TIME "tstop is behind.\n\n"

#define MSG_MAX_STEPS      _IDASLV_ MSG_TIME "maximum number of steps reached.\n\n"

#define MSG_EWT_NOW_BAD    _IDASLV_ MSG_TIME "some ewt component has become <= 0.0.\n\n"

#define MSG_EWTQ_NOW_BAD   _IDASLV_ MSG_TIME "some ewtQ component has become <= 0.0.\n\n"

#define MSG_EWTS_NOW_BAD   _IDASLV_ MSG_TIME "some ewtS component has become <= 0.0.\n\n"

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

#define MSG_YRET_NULL      _IDASLV_ "yret = NULL illegal.\n\n"
#define MSG_YPRET_NULL     _IDASLV_ "ypret = NULL illegal.\n\n"
#define MSG_TRET_NULL      _IDASLV_ "tret = NULL illegal.\n\n"

#define MSG_BAD_ITASK      _IDASLV_ "itask has an illegal value.\n\n"

#define MSG_NO_TSTOP       _IDASLV_ "tstop not set for this itask. \n\n"

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

#define MSG_IDAS_NEG_MXSTEPS "IDASetMaxNumSteps-- mxsteps <= 0 illegal. \n\n"

#define MSG_IDAS_NEG_HMAX    "IDASetMaxStep-- hmax <= 0 illegal. \n\n"

#define MSG_IDAS_NEG_EPCON   "IDASetNonlinConvCoef-- epcon < 0.0 illegal. \n\n"

#define MSG_IDAS_BAD_ITOL    "IDASetTolerances-- itol has an illegal value.\n\n"

#define MSG_IDAS_RTOL_NULL   "IDASetTolerances-- rtol = NULL illegal.\n\n"

#define MSG_IDAS_BAD_RTOL    "IDASetTolerances-- *rtol < 0 illegal.\n\n"

#define MSG_IDAS_ATOL_NULL   "IDASetTolerances-- atol = NULL illegal.\n\n"

#define MSG_IDAS_BAD_ATOL    "IDASetTolerances-- some atol component < 0.0 illegal.\n\n"

#define MSG_IDAS_BAD_EPICCON "IDASetNonlinConvCoefIC-- epiccon < 0.0 illegal.\n\n"

#define MSG_IDAS_BAD_MAXNH   "IDASetMaxNumStepsIC-- maxnh < 0 illegal.\n\n"

#define MSG_IDAS_BAD_MAXNJ   "IDASetMaxNumJacsIC-- maxnj < 0 illegal.\n\n"

#define MSG_IDAS_BAD_MAXNIT  "IDASetMaxNumItersIC-- maxnit < 0 illegal.\n\n"

#define MSG_IDAS_BAD_STEPTOL "IDASetLineSearchOffIC-- steptol < 0.0 illegal.\n\n"

#define MSG_IDAS_BAD_ITOLQ   "IDASetQuadTolerances-- itolQ has an illegal value.\n"

#define MSG_IDAS_RTOLQ_NULL  "IDASetQuadTolerances-- reltolQ = NULL illegal.\n\n"
 
#define MSG_IDAS_ATOLQ_NULL  "IDASetQuadTolerances-- abstolQ = NULL illegal.\n\n"

#define MSG_IDAS_BAD_RTOLQ   "IDASetQuadTolerances-- *reltolQ < 0.0 illegal.\n\n"

#define MSG_IDAS_BAD_ATOLQ   "IDASetQuadTolerances-- some abstolQ component < 0.0 illegal.\n\n"  

#define MSG_IDAS_BAD_ITOLS   "IDASetSensTolerances-- itolS has an illegal value.\n\n"

#define MSG_IDAS_RTOLS_NULL  "IDASetSensTolerances-- reltolS = NULL illegal.\n\n"
 
#define MSG_IDAS_ATOLS_NULL  "IDASetSensTolerances-- abstolS = NULL illegal.\n\n"

/* IDAGet* Error Messages */

#define MSG_IDAG_NO_MEM    "IDAGet*-- ida_mem = NULL illegal. \n\n"

#define MSG_IDAG_BAD_T1    "IDAGetSolution/IDAGetQuad/IDAGetSens-- "
#define MSG_IDAG_BAD_T     MSG_IDAG_BAD_T1 MSG_TIME MSG_TIME_INT

#define MSG_IDAG_NO_QUAD   "IDAGetQuad*-- quadrature were not computed.\n\n"

#define MSG_IDAG_NO_SENS   "IDAGetSens*-- sensitivities were not computed.\n\n"

#define MSG_IDAG_BAD_IS    "IDAGetSens1-- is has an illegal value.\n\n"

#define MSG_IDAG_NO_STGR1  "IDAGetSensStgr*-- STAGGERED1 method was not used.\n\n"


#endif

#ifdef __cplusplus
}
#endif
