/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2004-10-21 15:58:17 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Allan G. Taylor, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California  
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/ida/LICENSE
 * -----------------------------------------------------------------
 * This is the header file (private version) for the main IDA solver.
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _ida_impl_h
#define _ida_impl_h

#include <stdio.h>
#include "sundialstypes.h"
#include "nvector.h"
#include "ida.h"

/* Basic IDA constants */

#define MXORDP1   6     /* max. number of N_Vectors kept in the phi array */

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
  int            ida_itol;           /* itol = SS or SV                    */
  realtype      *ida_rtol;           /* ptr to relative tolerance          */
  void          *ida_atol;           /* ptr to absolute tolerance          */  
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

  /* Scalars for use by IDACalcIC*/

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

  booleantype ida_tstopset;
  realtype ida_tstop;

  /* Step Data */

  int ida_kk;        /* current BDF method order                              */
  int ida_kused;     /* method order used on last successful step             */
  int ida_knew;      /* order for next step from order decrease decision      */
  int ida_phase;     /* flag to trigger step doubling in first few steps      */
  long int ida_ns;   /* counts steps at fixed stepsize and order              */

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
  realtype ida_epcon;    /* coeficient of the Newton covergence test          */
  realtype ida_toldel;   /* tolerance in direct test on Newton corrections    */

 /* Limits */

  int ida_maxncf;        /* max numer of convergence failures                 */
  int ida_maxcor;        /* max number of Newton corrections                  */
  int ida_maxnef;        /* max number of error test failures                 */

  int ida_maxord;        /* max value of method order k:                      */
  int ida_mxstep;        /* max number of internal steps for one user call    */
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

} *IDAMem;

#endif

#ifdef __cplusplus
}
#endif
