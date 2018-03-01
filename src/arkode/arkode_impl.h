/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
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
#include <arkode/arkode_butcher.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include "arkode_adapt_impl.h"
#include "arkode_interp_impl.h"
#include "arkode_root_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
  ARKODE Private Constants                             
  ===============================================================*/

/* Basic ARKODE constants */
#define Q_DEFAULT        4       /* default RK order */
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
  ---------------------------------------------------------------
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
  ---------------------------------------------------------------
  arkRootCheck1:  ARK_SUCCESS or ARK_RTFUNC_FAIL

  arkRootCheck2:  ARK_SUCCESS, ARK_RTFUNC_FAIL, CLOSERT or RTFOUND

  arkRootCheck3:  ARK_SUCCESS, ARK_RTFUNC_FAIL or RTFOUND

  arkRootfind:  ARK_SUCCESS, ARK_RTFUNC_FAIL or RTFOUND
  ---------------------------------------------------------------*/
#define RTFOUND          +1
#define CLOSERT          +3


/*---------------------------------------------------------------
  Algorithmic constants
  ---------------------------------------------------------------
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
  ARKODE Interface function definitions
  ===============================================================*/

/* NOTE: documentation for the purpose of these functions is 
   located at the end of this file */
  
/* linear solver interface functions */  
typedef int (*ARKLinsolInitFn)(void* arkode_mem);
typedef int (*ARKLinsolSetupFn)(void* arkode_mem, int convfail, 
                                realtype tpred, N_Vector ypred,
                                N_Vector fpred,
                                booleantype *jcurPtr,
                                N_Vector vtemp1,
                                N_Vector vtemp2, N_Vector vtemp3); 
typedef int (*ARKLinsolSolveFn)(void* arkode_mem, N_Vector b, 
                                realtype tcur, N_Vector ycur,
                                N_Vector fcur, realtype client_tol,
                                int mnewt);
typedef int (*ARKLinsolFreeFn)(void* arkode_mem);
  
/* mass-matrix solver interface functions */  
typedef int (*ARKMassInitFn)(void *arkode_mem);
typedef int (*ARKMassSetupFn)(void *arkode_mem, N_Vector vtemp1,
                              N_Vector vtemp2, N_Vector vtemp3); 
typedef int (*ARKMassMultFn)(void *arkode_mem, N_Vector v,
                             N_Vector Mv);
typedef int (*ARKMassSolveFn)(void *arkode_mem, N_Vector b,
                              realtype client_tol);
typedef int (*ARKMassFreeFn)(void *arkode_mem);

/* time stepper interface functions */
typedef int (*ARKTimestepInitFn)(void* arkode_mem);
typedef int (*ARKTimestepAttachLinsolFn)(void* arkode_mem,
                                         ARKLinsolInitFn linit,
                                         ARKLinsolSetupFn lsetup,
                                         ARKLinsolSolveFn lsolve,
                                         ARKLinsolFreeFn lfree,
                                         int lsolve_type,
                                         void *lmem);
typedef int (*ARKTimestepAttachMasssolFn)(void* arkode_mem,
                                          ARKMassInitFn minit,
                                          ARKMassSetupFn msetup,
                                          ARKMassMultFn mmult,
                                          ARKMassSolveFn msolve,
                                          ARKMassFreeFn mfree,
                                          int msolve_type,
                                          void *mass_mem);
typedef void (*ARKTimestepDisableLSetup)(void* arkode_mem);
typedef void (*ARKTimestepDisableMSetup)(void* arkode_mem);
typedef void* (*ARKTimestepGetLinMemFn)(void* arkode_mem);
typedef void* (*ARKTimestepGetMassMemFn)(void* arkode_mem);
typedef ARKRhsFn (*ARKTimestepGetImplicitRHSFn)(void* arkode_mem);
typedef int (*ARKTimestepGetGammasFn)(void* arkode_mem,
                                      realtype *gamma,
                                      realtype *gamrat,
                                      booleantype **jcur,
                                      booleantype *dgamma_fail);
typedef int (*ARKTimestepFullRHSFn)(void* arkode_mem, realtype t, 
                                    N_Vector y, N_Vector f, int mode);
typedef int (*ARKTimestepResizeFn)(void* arkode_mem,
                                   ARKVecResizeFn resize,
                                   void *resize_data,
                                   sunindextype lrw_diff,
                                   sunindextype liw_diff,
                                   N_Vector tmpl);
typedef void (*ARKTimestepPrintFn)(void* arkode_mem, FILE* outfile);
typedef int (*ARKTimestepStepFn)(void* arkode_mem);
typedef int (*ARKTimestepFreeFn)(void* arkode_mem);

  
/*===============================================================
  ARKODE DATA STRUCTURES
  ===============================================================*/

/*---------------------------------------------------------------
  Types : struct ARKodeMassMemRec, ARKodeMassMem
  ---------------------------------------------------------------
  The type ARKodeMassMem is type pointer to struct 
  ARKodeMassMemRec.  This structure contains data pertaining to 
  the use of a non-identity mass matrix.
  ---------------------------------------------------------------*/
typedef struct ARKodeMassMemRec {

  /* mass matrix linear solver interface function pointers */
  ARKMassInitFn   minit;
  ARKMassSetupFn  msetup;
  ARKMassMultFn   mmult;
  ARKMassSolveFn  msolve;
  ARKMassFreeFn   mfree;
  void*           sol_mem;     /* mass matrix solver interface data */
  int             msolve_type; /* mass matrix interface type: 
                                  0=iterative; 1=direct; 2=custom */

} *ARKodeMassMem;


/*---------------------------------------------------------------
  Types : struct ARKodeFPMemRec, ARKodeFPMem
  ---------------------------------------------------------------
  The type ARKodeFPMem is type pointer to struct ARKodeFPMemRec.  
  This structure contains data pertaining to the use of the 
  [accelerated] fixed-point nonlinear solver.
  ---------------------------------------------------------------*/
typedef struct ARKodeFPMemRec {

  long int  m;      /* number of acceleration vectors to use */
  long int *imap;   /* array of length m                     */
  realtype *R;      /* array of length m*m                   */
  realtype *gamma;  /* array of length m                     */
  N_Vector *df;     /* vector array of length m              */
  N_Vector *dg;     /* vector array of length m              */
  N_Vector *q;      /* vector array of length m              */
  N_Vector  fval;   /* temporary N_Vectors                   */
  N_Vector  fold;
  N_Vector  gold;

} *ARKodeFPMem;


/*---------------------------------------------------------------
  Types : struct ARKodeMemRec, ARKodeMem
  ---------------------------------------------------------------
  The type ARKodeMem is type pointer to struct ARKodeMemRec.
  This structure contains fields to keep track of problem state.
  ---------------------------------------------------------------*/
typedef struct ARKodeMemRec {

  realtype uround;         /* machine unit roundoff */

  /* Problem specification data */
  void        *user_data;  /* user ptr passed to supplied functions */
  int          itol;       /* itol = ARK_SS (scalar, default), 
                                     ARK_SV (vector),
                                     ARK_WF (user weight function)  */
  int          ritol;      /* itol = ARK_SS (scalar, default), 
                                     ARK_SV (vector),
                                     ARK_WF (user weight function)  */
  realtype     reltol;     /* relative tolerance                    */
  realtype     Sabstol;    /* scalar absolute solution tolerance    */
  N_Vector     Vabstol;    /* vector absolute solution tolerance    */
  realtype     SRabstol;   /* scalar absolute residual tolerance    */
  N_Vector     VRabstol;   /* vector absolute residual tolerance    */
  booleantype  user_efun;  /* SUNTRUE if user sets efun             */
  ARKEwtFn     efun;       /* function to set ewt                   */
  void        *e_data;     /* user pointer passed to efun           */
  booleantype  user_rfun;  /* SUNTRUE if user sets rfun             */
  ARKRwtFn     rfun;       /* function to set rwt                   */
  void        *r_data;     /* user pointer passed to rfun           */

  /* Time stepper module */
  ARKTimestepAttachLinsolFn   step_attachlinsol;
  ARKTimestepAttachMasssolFn  step_attachmasssol;
  ARKTimestepDisableLSetup    step_disablelsetup;
  ARKTimestepDisableMSetup    step_disablemsetup;
  ARKTimestepGetLinMemFn      step_getlinmem;
  ARKTimestepGetMassMemFn     step_getmassmem;
  ARKTimestepGetImplicitRHSFn step_getimplicitrhs;
  ARKMassMultFn               step_mmult;
  ARKTimestepGetGammasFn      step_getgammas;
  ARKTimestepInitFn           step_init;
  ARKTimestepFullRHSFn        step_fullrhs;
  ARKTimestepStepFn           step;
  ARKTimestepResizeFn         step_resize;
  ARKTimestepPrintFn          step_print;
  ARKTimestepFreeFn           step_free;
  void                       *step_mem;

  /* N_Vector storage */
  N_Vector ewt;     /* error weight vector                               */
  N_Vector rwt;     /* residual weight vector                            */
  booleantype rwt_is_ewt;     /* SUNTRUE if rwt is a pointer to ewt      */
  N_Vector ycur;    /* pointer to user-provided solution memory; used as 
                       evolving solution by the timestepper modules      */
  N_Vector yn;      /* solution from the last successful step            */
  N_Vector tempv1;  /* temporary storage vector                          */
  N_Vector tempv2;

  /* Temporal interpolation module */
  ARKodeInterpMem interp;  
  int dense_q;      /* interpolation order (user request)      */

  /* Tstop information */
  booleantype tstopset;
  realtype    tstop;

  /* Time step data */  
  realtype hin;             /* initial step size                        */
  realtype h;               /* current step size                        */
  realtype hmin;            /* |h| >= hmin                              */
  realtype hmax_inv;        /* |h| <= 1/hmax_inv                        */
  realtype hprime;          /* step size to be used on the next step    */ 
  realtype next_h;          /* step size to be used on the next step    */ 
  realtype eta;             /* eta = hprime / h                         */
  realtype tcur;            /* current internal value of t 
                               (changes with each stage)                */
  realtype tretlast;        /* value of tret last returned by ARKode    */
  booleantype fixedstep;    /* flag to disable temporal adaptivity      */

  
  /* Limits and various solver parameters */
  long int mxstep;   /* max number of internal steps for one user call */
  int      mxhnil;   /* max number of warning messages issued to the
                        user that t+h == t for the next internal step  */

  /* Counters */
  long int nst;          /* number of internal steps taken             */
  int      nhnil;        /* number of messages issued to the user that 
                            t+h == t for the next iternal step         */

  /* Diagnostic output */
  booleantype report;   /* flag to enable/disable diagnostic output    */
  FILE       *diagfp;   /* diagnostic outputs are sent to diagfp   */

  /* Space requirements for ARKODE */
  sunindextype lrw1;        /* no. of realtype words in 1 N_Vector          */ 
  sunindextype liw1;        /* no. of integer words in 1 N_Vector           */ 
  long int lrw;             /* no. of realtype words in ARKODE work vectors */
  long int liw;             /* no. of integer words in ARKODE work vectors  */

  /* Saved Values */
  realtype    h0u;          /* actual initial stepsize                     */
  realtype    tn;           /* time of last successful step                */
  realtype    hold;         /* last successful h value used                */
  realtype    tolsf;        /* tolerance scale factor (suggestion to user) */
  booleantype VabstolMallocDone;
  booleantype VRabstolMallocDone;
  booleantype MallocDone;  
  booleantype resized;      /* denotes first step after ARKodeResize      */
  booleantype firststage;   /* denotes first stage in simulation          */

  /* Error handler function and error ouput file */
  ARKErrHandlerFn ehfun;    /* error messages are handled by ehfun        */
  void           *eh_data;  /* data pointer passed to ehfun               */
  FILE           *errfp;    /* ARKODE error messages are sent to errfp    */

  /* Rootfinding Data */
  ARKodeRootMem root_mem;          /* root-finding structure */

  /* User-supplied step solution post-processing function */
  ARKPostProcessStepFn ProcessStep;

} *ARKodeMem;



/*===============================================================
  Interface To Linear Solvers
  ===============================================================*/
  
/*---------------------------------------------------------------
  Communication between ARKODE and a ARKODE Linear Solver
  -----------------------------------------------------------------
  convfail (input to lsetup)

  ARK_NO_FAILURES : Either this is the first lsetup call for 
                    this step, or the local error test failed on 
                    the previous attempt at this step (but the 
                    Newton iteration converged).

  ARK_FAIL_BAD_J  : This value is passed to lsetup if
 
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
  --------------------------------------------------------------*/
  
/* Constants for convfail (input to lsetup) */
#define ARK_NO_FAILURES 0
#define ARK_FAIL_BAD_J  1
#define ARK_FAIL_OTHER  2

/*---------------------------------------------------------------
  ARKLinsolInitFn
  ---------------------------------------------------------------
  This function should complete initializations for a specific 
  ARKode linear solver interface, such as counters and statistics.
  This should return 0 if it has successfully initialized the 
  ARKODE linear solver interface and a negative value otherwise. 
  If an error does occur, an appropriate message should be sent 
  to the error handler function.
  ---------------------------------------------------------------*/
  
/*---------------------------------------------------------------
  ARKLinsolSetupFn
  ---------------------------------------------------------------
  This function should prepare the linear solver interface for 
  subsequent calls to the ARKLinsolSolveFn routine. It may 
  recompute Jacobian-related data is it deems necessary. Its 
  parameters are as follows:

  arkode_mem - void* problem memory pointer of type ARKodeMem. See 
           the typedef earlier in this file.

  convfail - a flag to indicate any problem that occurred during
             the solution of the nonlinear equation on the
             current time step for which the linear solver is
             being used. This flag can be used to help decide
             whether the Jacobian data kept by a ARKODE linear
             solver needs to be updated or not.
             Its possible values have been documented above.

  tpred - the time for the current ARKODE internal step.

  ypred - the predicted y vector for the current ARKODE internal
          step.

  fpred - f(tpred, ypred).

  jcurPtr - a pointer to a boolean to be filled in by lsetup.
            The function should set *jcurPtr=SUNTRUE if its Jacobian
            data is current after the call and should set
            *jcurPtr=SUNFALSE if its Jacobian data is not current.
            Note: If lsetup calls for re-evaluation of
            Jacobian data (based on convfail and ARKODE state
            data), it should return *jcurPtr=SUNTRUE always;
            otherwise an infinite loop can result.

  vtemp1 - temporary N_Vector provided for use by lsetup.

  vtemp3 - temporary N_Vector provided for use by lsetup.

  vtemp3 - temporary N_Vector provided for use by lsetup.

  This routine should return 0 if successful, a positive value 
  for a recoverable error, and a negative value for an 
  unrecoverable error.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKLinsolSolveFn
  ---------------------------------------------------------------
  This routine must solve the linear equation P x = b, where
  P is some approximation to (M - gamma J), M is the system mass
  matrix, J = (df/dy)(tcur,ycur), and the RHS vector b is input. The 
  N-vector ycur contains the solver's current approximation to 
  y(tcur) and the vector fcur contains the N_Vector f(tcur,ycur).  
  The input client_tol contains the desired accuracy (in the wrms 
  norm) of the routine calling the solver; the ARKDLS solver 
  ignores this value and the ARKSPILS solver tightens it by the 
  factor eplifac.  The input mnewt is the current nonlinear 
  iteration index (ignored by ARKDLS, used by ARKSPILS).
  
  Additional vectors that are set within the ARKode memory 
  structure, and that may be of use within an iterative linear 
  solver, include:
  
  ewt - the error weight vector (scaling for solution vector)
  
  rwt - the residual weight vector (scaling for rhs vector)
  
  The solution is to be returned in the vector b.  This should
  return a positive value for a recoverable error and a 
  negative value for an unrecoverable error. Success is 
  indicated by a 0 return value.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKLinsolFreeFn
  ---------------------------------------------------------------
  This should free up any memory allocated by the linear solver 
  interface. This routine is called once a problem has been
  completed and the linear solver is no longer needed.  It should 
  return 0 upon success, or a nonzero on failure.
  ---------------------------------------------------------------*/



/*---------------------------------------------------------------
  ARKMassInitFn
  ---------------------------------------------------------------
  This function should complete initializations for a specific 
  mass matrix linear solver interface, such as counters and  
  statistics. A function of this type should return 0 if it 
  has successfully initialized the mass matrix linear solver and 
  a negative value otherwise.  If an error does occur, an 
  appropriate message should be sent to the error handler function.
  ---------------------------------------------------------------*/
  
/*---------------------------------------------------------------
  ARKMassSetupFn
  ---------------------------------------------------------------
  This should prepare the mass matrix solver interface for 
  subsequent calls to the ARKMassMultFn and ARKMassSolveFn 
  routines. It may recompute mass matrix related data is it deems
  necessary. Its parameters are as follows:
  
  arkode_mem - void* problem memory pointer of type ARKodeMem. See 
           the typedef earlier in this file.

  vtemp1, vtemp2, vtemp3 - temporary N_Vectors

  This routine should return 0 if successful, and a negative 
  value for an unrecoverable error.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKMassMultFn
  ---------------------------------------------------------------
  This must compute the matrix-vector product, z = M*v, where M is
  the system mass matrix the vector v is input, and the vector z
  is output. The mmult routine returns a positive value for a
  recoverable error and a negative value for an unrecoverable 
  error. Success is indicated by a 0 return value.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKMassSolveFn
  ---------------------------------------------------------------
  This must solve the linear equation M x = b, where M is the 
  system mass matrix, and the RHS vector b is input. The 
  realtype client_tol contains the desired accuracy (in the wrms 
  norm) of the routine calling the solver; the ARKDLS solver 
  ignore this value and the ARKSPILS solver tightens it by the 
  factor eplifac.  The solution is to be returned in the vector b.  
  
  Additional vectors that are set within the ARKode memory 
  structure, and that may be of use within an iterative linear 
  solver, include:
  
  ewt - the error weight vector (scaling for solution vector)
  
  rwt - the residual weight vector (scaling for rhs vector)
  
  This routine should return a positive value for a recoverable  
  error and a negative value for an unrecoverable error. Success 
  is indicated by a 0 return value.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKMassFreeFn
  ---------------------------------------------------------------
  This should free up any memory allocated by the mass matrix
  solver interface. This routine is called once a problem has been
  completed and the solver is no longer needed.  It should return
  0 upon success, or a nonzero on failure.
  ---------------------------------------------------------------*/


  
  
/*===============================================================
  Interface to Time Steppers
  ===============================================================*/
  
/*---------------------------------------------------------------
  ARKTimestepAttachLinsolFn
  ---------------------------------------------------------------
  This routine should attach the various set of system linear 
  solver interface routines, linear solver interface data 
  structure, and system linear solver type to the ARKode time 
  stepping module pointed to in ark_mem->step_mem.  This will 
  be called by one of the various ARKode linear solver interfaces 
  (SPILS, DLS).
  
  This routine should return 0 if it has successfully attached 
  these items and a negative value otherwise.  If an error does 
  occur, an appropriate message should be sent to the ARKode 
  error handler function.
  ---------------------------------------------------------------*/
  
/*---------------------------------------------------------------
  ARKTimestepAttachMasssolFn
  ---------------------------------------------------------------
  This routine should attach the various set of mass matrix linear 
  solver interface routines, data structure, and solver type to 
  the ARKode time stepping module pointed to in 
  ark_mem->step_mem.  This will be called by one of the 
  various ARKode linear solver interfaces (SPILS, DLS).

  This routine should return 0 if it has successfully attached 
  these items, and a negative value otherwise.  If an error does 
  occur, an appropriate message should be sent to the ARKode 
  error handler function.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepDisableLSetup
  ---------------------------------------------------------------
  This routine should NULLify any ARKLinsolSetupFn function 
  pointer stored in the ARKode time stepping module (initially set 
  in a call to ARKTimestepAttachLinsolFn).  This can be called by 
  ARKSPILS when preconditioning is disabled.  

  This routine has no return value.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepDisableMSetup
  ---------------------------------------------------------------
  This routine should NULLify any ARKMassSetupFn function pointer
  stored in the ARKode time stepping module (initially set in a 
  call to ARKTimestepAttachMasssolFn).  This can be called by 
  ARKSPILS when preconditioning is disabled.  

  This routine has no return value.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepGetLinMemFn
  ---------------------------------------------------------------
  This routine should return the linear solver memory structure
  used by the ARKode time stepping module pointed to in 
  ark_mem->step_mem.  This will be called by one of the 
  various ARKode linear solver interfaces (SPILS, DLS).
  
  This routine should return NULL if no linear solver memory 
  structure is attached.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepGetMassMemFn
  ---------------------------------------------------------------
  This routine should return the mass matrix linear solver memory 
  structure used by the ARKode time stepping module pointed to in 
  ark_mem->step_mem.  This will be called by one of the 
  various ARKode mass matrix solver interfaces (SPILS, DLS).
  
  This routine should return NULL if no mass matrix solver memory 
  structure is attached.
  ---------------------------------------------------------------*/
  
/*---------------------------------------------------------------
  ARKTimestepGetImplicitRHSFn
  ---------------------------------------------------------------
  This routine should return the implicit RHS function pointer for 
  the current nonlinear solve (if there are multiple); it is used 
  inside the linear solver interfaces for approximation of 
  Jacobian matrix elements and/or matrix-vector products.
  
  This routine should return NULL if no implicit RHS function is 
  active.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepGetGammasFn
  ---------------------------------------------------------------
  This routine should fill the current value of gamma, the ratio 
  of the current gamma value to the gamma value when the 
  Jacobian/preconditioner was last updated, a pointer to the 
  time step module internal booleantype variable indicating 
  whether the preconditioner is current, and a logic value 
  indicating whether the gamma value is sufficiently stale 
  to cause recomputation of Jacobian/preconditioner data.  Here, 
  gamma is the coefficient preceding the RHS Jacobian 
  matrix, J, in the full nonlinear system Jacobian, 
  A = M - gamma*J.  

  The time step module must contain a booleantype variable to 
  provide for the boolentype pointer (jcur).  This is only used
  by the ARKSPILS interface, so could be NULL for time step 
  modules that only work with ARKDLS.  Optionally, the value of 
  this parameter could be set to SUNFALSE prior to return from
  the ARKTimestepGetGammasFn to force recalculation of 
  preconditioner information.

  The value of the logic flag is used as follows:  if a previous 
  Newton iteration failed due to a bad Jacobian/preconditioner, 
  and this flag is SUNFALSE, this will trigger recalculation of 
  the Jacobian/preconditioner.

  This routine should return 0 if it has successfully attached 
  these items, and a negative value otherwise.  If an error does 
  occur, an appropriate message should be sent to the ARKode 
  error handler function.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepInitFn
  ---------------------------------------------------------------
  This routine should complete initializations for a specific 
  ARKode time stepping module, such as (non)linear solver data, 
  Butcher table(s), counters and statistics.  

  This routine should return 0 if it has successfully initialized 
  the ARKODE time stepper module and a negative value otherwise. 
  If an error does occur, an appropriate message should be sent 
  to the error handler function.
  ---------------------------------------------------------------*/
  
/*---------------------------------------------------------------
  ARKTimestepFullRHSFn
  ---------------------------------------------------------------
  This routine must compute the full ODE right-hand side function 
  at the inputs (t,y), and store the result in the N_Vector f. 
  Depending on the type of stepper, this may be just the single 
  ODE RHS function supplied (e.g. ERK, DIRK, IRK), or it may be 
  the sum of many ODE RHS functions (e.g. ARK, RMIS).  The 'mode' 
  flag indicates where this routine is called:
     0 -> called at the beginning of a simulation
     1 -> called at the end of a successful step
     2 -> called elsewhere (e.g. for dense output)
  It is recommended that the stepper use the mode information to
  maximize reuse between calls to this function and RHS 
  evaluations inside the stepper itself.

  This routine should return 0 if successful, and a negative value
  otherwise.  If an error does occur, an appropriate message 
  should be sent to the error handler function.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepStepFn
  ---------------------------------------------------------------
  This is the primary computational routine for any ARKode 
  time-stepping module.  It must attempt to advance the solution 
  vector one internal time step, from tn to tn+h.  The routine may 
  internally adjust the value of h if needed to complete the step 
  successfully (e.g. to meet error goals or to converge the 
  (non)linear solution).

  It is assumed that this routine uses/modifies general problem 
  data directly out of the main ARKodeMem structure, but that all
  method-specific data be stored in the step-module-specific data
  structure.  Relevant items in the ARKodeMem structure for this
  purpose include:
  - tcur -- the current "t" value
  - ycur -- the current "y" value on input; should hold the 
    time-evolved solution on output
  - h -- the suggested/maximum "h" value to use; if the step 
    eventually completes with a smaller "h" value, then that 
    should be stored here
  - tn -- "t" value at end of the last successful step
  - nst -- the counter for overall successful steps
  - nst_attempts -- the counter for overall step attempts
  - user_data -- the (void *) pointer returned to user for 
    RHS calls
  - report / diagfp -- if any diagnostic information is 
    to be saved to disk, the report flag indicates whether 
    this is enabled, and diagfp provides the file pointer 
    where this information should be written

  Possible return values (not all must be used by the stepper):
  - ARK_SUCCESS -- the step completed successfully (although 
    perhaps with a shortened h)
  - ARK_ERR_FAILURE -- the error test failed repeatedly or 
    with |h| = hmin
  - ARK_CONV_FAILURE -- the solver convergence test failed 
    repeatedly or with |h| = hmin
  - ARK_LSETUP_FAIL -- the linear solver setup routine failed 
    in an unrecoverable manner
  - ARK_LSOLVE_FAIL -- the linear solve routine failed in an 
    unrecoverable manner
  - ARK_RHSFUNC_FAIL -- the ODE right-hand side routine failed 
    in an unrecoverable manner
  - ARK_UNREC_RHSFUNC_ERR -- the right-hand side failed in a 
    recoverable manner, but no recovery is possible
  - ARK_REPTD_RHSFUNC_ERR -- repeated recoverable right-hand 
    side function errors
  - ARK_RTFUNC_FAIL -- the rootfinding routine failed in an 
    unrecoverable manner    
  - ARK_TOO_CLOSE -- tout too close to t0 to start integration
  - ARK_MASSSOLVE_FAIL -- the mass matrix solver failed

  If additional failure modes need to be added for future 
  steppers, the relevant flags should be added to 
  include/arkode/arkode.h, and the routine arkHandleFailure (in 
  src/arkode/arkode.c) must be modified to handle the new flags.
  ---------------------------------------------------------------*/
  
/*---------------------------------------------------------------
  ARKTimestepResizeFn
  ---------------------------------------------------------------
  This routine must resize any N_Vectors in the time step module
  based on a call to ARKodeResize.  

  This routine should return 0 if successful, and a negative value
  otherwise.  If an error does occur, an appropriate message 
  should be sent to the error handler function.
  ---------------------------------------------------------------*/
  
/*---------------------------------------------------------------
  ARKTimestepPrintFn
  ---------------------------------------------------------------
  If provided (non-NULL), this routine should output the time 
  step memory structure data to a specified file pointer.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepFreeFn
  ---------------------------------------------------------------
  This should free up any memory allocated by the time stepper 
  module. This routine is called once a problem has been completed 
  and the time stepper module is no longer needed.  It should 
  return 0 upon success, or a nonzero on failure.
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
  ARKODE PRIVATE FUNCTION PROTOTYPES
  ===============================================================*/
#ifdef __GNUC__
#define SUNDIALS_UNUSED __attribute__ ((unused))
#else
#define SUNDIALS_UNUSED
#endif

int arkodeInit(ARKodeMem ark_mem, realtype t0, N_Vector y0);
int arkodeReInit(ARKodeMem ark_mem, realtype t0, N_Vector y0);
booleantype arkAllocVec(ARKodeMem ark_mem, 
                        N_Vector tmpl,
                        N_Vector *v);
void arkFreeVec(ARKodeMem ark_mem, N_Vector *v);
int arkResizeVec(ARKodeMem ark_mem, 
                 ARKVecResizeFn resize,
                 void *resize_data,
                 sunindextype lrw_diff,
                 sunindextype liw_diff,
                 N_Vector tmpl,
                 N_Vector *v);
void arkPrintMem(ARKodeMem ark_mem, FILE *outfile) SUNDIALS_UNUSED;
booleantype arkCheckTimestepper(ARKodeMem ark_mem);
booleantype arkCheckNvector(N_Vector tmpl);
booleantype arkAllocVectors(ARKodeMem ark_mem, 
                            N_Vector tmpl);
/* int arkAllocRKVectors(ARKodeMem ark_mem); */
void arkFreeVectors(ARKodeMem ark_mem);

ARKodeFPMem arkAllocFPData(ARKodeMem ark_mem, long int m);
int arkResizeFPData(ARKodeMem ark_mem,
                    ARKodeFPMem fp_mem,
                    ARKVecResizeFn resize,
                    void *resize_data,
                    sunindextype lrw_diff,
                    sunindextype liw_diff,
                    N_Vector tmpl);
void arkFreeFPData(ARKodeMem ark_mem, ARKodeFPMem fp_mem);

int arkInitialSetup(ARKodeMem ark_mem, realtype tout);
int arkPostResizeSetup(ARKodeMem ark_mem);
int arkStopTests(ARKodeMem ark_mem, realtype tout, N_Vector yout, 
                 realtype *tret, int itask, int *ier);
int arkHin(ARKodeMem ark_mem, realtype tout);
realtype arkUpperBoundH0(ARKodeMem ark_mem, 
                         realtype tdist);
int arkYddNorm(ARKodeMem ark_mem, realtype hg, 
               realtype *yddnrm);

int arkCompleteStep(ARKodeMem ark_mem);
int arkHandleFailure(ARKodeMem ark_mem,int flag);

int arkEwtSetSS(ARKodeMem ark_mem, N_Vector ycur, 
                N_Vector weight);
int arkEwtSetSV(ARKodeMem ark_mem, N_Vector ycur, 
                N_Vector weight);
int arkRwtSetSS(ARKodeMem ark_mem, N_Vector My, 
                N_Vector weight);
int arkRwtSetSV(ARKodeMem ark_mem, N_Vector My, 
                N_Vector weight);


/*===============================================================
  Reusable ARKode Error Messages
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
#define MSGARK_LSOLVE_NULL    "The linear solver object is NULL."
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
#define MSGARK_MASSSOLVE_NULL "The mass matrix linear solver object is NULL."
#define MSGARK_MASSINIT_FAIL  "The mass matrix solver's init routine failed."
#define MSGARK_MASSSETUP_FAIL "The mass matrix solver's setup routine failed."
#define MSGARK_MASSSOLVE_FAIL "The mass matrix solver failed."
#define MSGARK_MASSFREE_FAIL  "The mass matrixsolver's free routine failed."

#ifdef __cplusplus
}
#endif

#endif
