/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Implementation header file for the main ARKODE integrator.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_IMPL_H
#define _ARKODE_IMPL_H

#include <stdarg.h>

#include <arkode/arkode.h>
#include <arkode/arkode_butcher.h>
#include <arkode/arkode_butcher_dirk.h>
#include <arkode/arkode_butcher_erk.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_adaptcontroller.h>

#include <sundials/priv/sundials_context_impl.h>
#include <sundials/priv/sundials_errors_impl.h>
#include "sundials_logger_impl.h"
#include "arkode_types_impl.h"
#include "arkode_adapt_impl.h"
#include "arkode_root_impl.h"
#include "arkode_relaxation_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM  ".32Lg"
#define RSYMW "41.32Lg"
#else
#define RSYM  ".16g"
#define RSYMW "23.16g"
#endif

/*===============================================================
  SHORTCUTS
  ===============================================================*/

#define ARK_PROFILER ark_mem->sunctx->profiler
#define ARK_LOGGER ark_mem->sunctx->logger

/*===============================================================
  MACROS
  ===============================================================*/

/* TODO(DJG): replace with signbit when C99+ is required */
#define DIFFERENT_SIGN(a,b) ( ( (a) < 0 && (b) > 0 ) || ( (a) > 0 && (b) < 0 ) )
#define SAME_SIGN(a,b) ( ( (a) > 0 && (b) > 0 ) || ( (a) < 0 && (b) < 0 ) )

/*===============================================================
  ARKODE Private Constants
  ===============================================================*/

/* Basic ARKODE defaults */
#define Q_DEFAULT        4      /* method order                       */
#define MXSTEP_DEFAULT   500    /* max steps between returns          */
#define MAXNEF           7      /* max number of error failures       */
#define MAXNCF           10     /* max number of convergence failures */
#define MAXCONSTRFAILS   10     /* max number of constraint failures  */
#define MXHNIL           10     /* max number of t+h==h warnings      */

/* Numeric constants */
#define ZERO   SUN_RCONST(0.0)      /* real 0.0     */
#define TINY   SUN_RCONST(1.0e-10)  /* small number */
#define TENTH  SUN_RCONST(0.1)      /* real 0.1     */
#define HALF   SUN_RCONST(0.5)      /* real 0.5     */
#define ONE    SUN_RCONST(1.0)      /* real 1.0     */
#define TWO    SUN_RCONST(2.0)      /* real 2.0     */
#define THREE  SUN_RCONST(3.0)      /* real 3.0     */
#define FOUR   SUN_RCONST(4.0)      /* real 4.0     */
#define FIVE   SUN_RCONST(5.0)      /* real 5.0     */

/* Control constants for tolerances */
#define ARK_SS  0
#define ARK_SV  1
#define ARK_WF  2


/*===============================================================
  ARKODE Routine-Specific Constants
  ===============================================================*/

/*---------------------------------------------------------------
  Initialization types
  ---------------------------------------------------------------*/
#define FIRST_INIT   0  /* first step (re-)initialization */
#define RESET_INIT   1  /* reset initialization           */
#define RESIZE_INIT  2  /* resize initialization          */

/*---------------------------------------------------------------
  Control constants for lower-level time-stepping functions
  ---------------------------------------------------------------*/
#define PREDICT_AGAIN    +3
#define CONV_FAIL        +4
#define TRY_AGAIN        +5
#define FIRST_CALL       +6
#define PREV_CONV_FAIL   +7
#define PREV_ERR_FAIL    +8
#define RHSFUNC_RECVR    +9
#define CONSTR_RECVR     +10

/*---------------------------------------------------------------
  Return values for lower-level rootfinding functions
  ---------------------------------------------------------------*/
#define RTFOUND          +1
#define CLOSERT          +3


/*---------------------------------------------------------------
  Algorithmic constants
  ---------------------------------------------------------------
  ARKodeGetDky and arkStep:  FUZZ_FACTOR

  arkHin:  H0_LBFACTOR, H0_UBFACTOR, H0_BIAS and H0_ITERS

  time comparison factors:
     ONEPSM      safety factor for floating point comparisons
     ONEMSM      safety factor for floating point comparisons
  ---------------------------------------------------------------*/
#define FUZZ_FACTOR SUN_RCONST(100.0)

#define H0_LBFACTOR SUN_RCONST(100.0)
#define H0_UBFACTOR SUN_RCONST(0.1)
#define H0_BIAS     HALF
#define H0_ITERS    4

#define ONEPSM      SUN_RCONST(1.000001)
#define ONEMSM      SUN_RCONST(0.999999)


/*===============================================================
  ARKODE Interface function definitions
  ===============================================================*/

/* NOTE: documentation for the purpose of these functions is
   located at the end of this file */

/* linear solver interface functions */
typedef int (*ARKLinsolInitFn)(void* arkode_mem);
typedef int (*ARKLinsolSetupFn)(void* arkode_mem, int convfail,
                                sunrealtype tpred, N_Vector ypred,
                                N_Vector fpred,
                                sunbooleantype *jcurPtr,
                                N_Vector vtemp1,
                                N_Vector vtemp2, N_Vector vtemp3);
typedef int (*ARKLinsolSolveFn)(void* arkode_mem, N_Vector b,
                                sunrealtype tcur, N_Vector ycur,
                                N_Vector fcur, sunrealtype client_tol,
                                int mnewt);
typedef int (*ARKLinsolFreeFn)(void* arkode_mem);

/* mass-matrix solver interface functions */
typedef int (*ARKMassInitFn)(void *arkode_mem);
typedef int (*ARKMassSetupFn)(void *arkode_mem, sunrealtype t,
                              N_Vector vtemp1, N_Vector vtemp2,
                              N_Vector vtemp3);
typedef int (*ARKMassMultFn)(void *arkode_mem, N_Vector v,
                             N_Vector Mv);
typedef int (*ARKMassSolveFn)(void *arkode_mem, N_Vector b,
                              sunrealtype client_tol);
typedef int (*ARKMassFreeFn)(void *arkode_mem);

/* time stepper interface functions */
typedef int (*ARKTimestepInitFn)(void* arkode_mem, int init_type);
typedef int (*ARKTimestepAttachLinsolFn)(void* arkode_mem,
                                         ARKLinsolInitFn linit,
                                         ARKLinsolSetupFn lsetup,
                                         ARKLinsolSolveFn lsolve,
                                         ARKLinsolFreeFn lfree,
                                         SUNLinearSolver_Type lsolve_type,
                                         void *lmem);
typedef int (*ARKTimestepAttachMasssolFn)(void* arkode_mem,
                                          ARKMassInitFn minit,
                                          ARKMassSetupFn msetup,
                                          ARKMassMultFn mmult,
                                          ARKMassSolveFn msolve,
                                          ARKMassFreeFn mfree,
                                          sunbooleantype time_dep,
                                          SUNLinearSolver_Type msolve_type,
                                          void *mass_mem);
typedef void (*ARKTimestepDisableLSetup)(void* arkode_mem);
typedef void (*ARKTimestepDisableMSetup)(void* arkode_mem);
typedef void* (*ARKTimestepGetLinMemFn)(void* arkode_mem);
typedef void* (*ARKTimestepGetMassMemFn)(void* arkode_mem);
typedef ARKRhsFn (*ARKTimestepGetImplicitRHSFn)(void* arkode_mem);
typedef int (*ARKTimestepGetGammasFn)(void* arkode_mem,
                                      sunrealtype *gamma,
                                      sunrealtype *gamrat,
                                      sunbooleantype **jcur,
                                      sunbooleantype *dgamma_fail);
typedef int (*ARKTimestepFullRHSFn)(void* arkode_mem, sunrealtype t,
                                    N_Vector y, N_Vector f, int mode);
typedef int (*ARKTimestepStepFn)(void* arkode_mem, sunrealtype *dsm,
                                 int *nflag);


/*===============================================================
  ARKODE interpolation module definition
  ===============================================================*/

/* Forward reference for pointer to ARKInterp_Ops object */
typedef struct _generic_ARKInterpOps *ARKInterpOps;

/* Forward reference for pointer to ARKInterp object */
typedef struct _generic_ARKInterp *ARKInterp;

/* Structure containing function pointers to interpolation operations  */
struct _generic_ARKInterpOps {
  int (*resize)(void* arkode_mem, ARKInterp interp,
                ARKVecResizeFn resize, void *resize_data,
                sunindextype lrw_diff, sunindextype liw_diff,
                N_Vector tmpl);
  void (*free)(void* arkode_mem, ARKInterp interp);
  void (*print)(ARKInterp interp, FILE *outfile);
  int (*setdegree)(void *arkode_mem, ARKInterp interp, int degree);
  int (*init)(void* arkode_mem, ARKInterp interp, sunrealtype tnew);
  int (*update)(void* arkode_mem, ARKInterp interp, sunrealtype tnew);
  int (*evaluate)(void* arkode_mem, ARKInterp interp,
                  sunrealtype tau, int d, int order, N_Vector yout);
};

/* An interpolation module consists of an implementation-dependent 'content'
   structure, and a pointer to a structure of implementation-dependent operations. */
struct _generic_ARKInterp {
  void *content;
  ARKInterpOps ops;
};

/* ARKInterp module functions */
int arkInterpResize(void* arkode_mem, ARKInterp interp,
                    ARKVecResizeFn resize, void *resize_data,
                    sunindextype lrw_diff, sunindextype liw_diff,
                    N_Vector tmpl);
void arkInterpFree(void* arkode_mem, ARKInterp interp);
void arkInterpPrintMem(ARKInterp interp, FILE *outfile);
int arkInterpSetDegree(void *arkode_mem, ARKInterp interp, int degree);
int arkInterpInit(void* arkode_mem, ARKInterp interp, sunrealtype tnew);
int arkInterpUpdate(void* arkode_mem, ARKInterp interp, sunrealtype tnew);
int arkInterpEvaluate(void* arkode_mem, ARKInterp interp,
                      sunrealtype tau, int d, int order, N_Vector yout);


/*===============================================================
  ARKODE data structures
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
  Types : struct ARKodeMemRec, ARKodeMem
  ---------------------------------------------------------------
  The type ARKodeMem is type pointer to struct ARKodeMemRec.
  This structure contains fields to keep track of problem state.
  ---------------------------------------------------------------*/
struct ARKodeMemRec
{
  SUNContext sunctx;

  sunrealtype uround;             /* machine unit roundoff */

  /* Problem specification data */
  void        *user_data;      /* user ptr passed to supplied functions */
  int          itol;           /* itol = ARK_SS (scalar, default),
                                         ARK_SV (vector),
                                         ARK_WF (user weight function)  */
  int          ritol;          /* itol = ARK_SS (scalar, default),
                                         ARK_SV (vector),
                                         ARK_WF (user weight function)  */
  sunrealtype     reltol;         /* relative tolerance                    */
  sunrealtype     Sabstol;        /* scalar absolute solution tolerance    */
  N_Vector     Vabstol;        /* vector absolute solution tolerance    */
  sunbooleantype  atolmin0;       /* flag indicating that min(abstol) = 0  */
  sunrealtype     SRabstol;       /* scalar absolute residual tolerance    */
  N_Vector     VRabstol;       /* vector absolute residual tolerance    */
  sunbooleantype  Ratolmin0;      /* flag indicating that min(Rabstol) = 0 */
  sunbooleantype  user_efun;      /* SUNTRUE if user sets efun             */
  ARKEwtFn     efun;           /* function to set ewt                   */
  void        *e_data;         /* user pointer passed to efun           */
  sunbooleantype  user_rfun;      /* SUNTRUE if user sets rfun             */
  ARKRwtFn     rfun;           /* function to set rwt                   */
  void        *r_data;         /* user pointer passed to rfun           */
  sunbooleantype  constraintsSet; /* check inequality constraints          */

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
  void                       *step_mem;

  /* N_Vector storage */
  N_Vector ewt;                 /* error weight vector                        */
  N_Vector rwt;                 /* residual weight vector                     */
  sunbooleantype rwt_is_ewt;       /* SUNTRUE if rwt is a pointer to ewt         */
  N_Vector ycur;                /* pointer to user-provided solution memory;
                                   used as evolving solution by the time stepper
                                   modules */
  N_Vector yn;                  /* solution from the last successful step     */
  N_Vector fn;                  /* full IVP right-hand side from last step    */
  sunbooleantype fn_is_current; /* SUNTRUE if fn has been evaluated at yn     */
  N_Vector tempv1;              /* temporary storage vectors (for local use   */
  N_Vector tempv2;              /* and by time-stepping modules)              */
  N_Vector tempv3;
  N_Vector tempv4;

  N_Vector constraints;   /* vector of inequality constraint options         */

  /* Temporal interpolation module */
  ARKInterp interp;
  int interp_type;

  /* Tstop information */
  sunbooleantype tstopset;
  sunbooleantype tstopinterp;
  sunrealtype    tstop;

  /* Time step data */
  sunrealtype hin;                /* initial step size                        */
  sunrealtype h;                  /* current step size                        */
  sunrealtype hmin;               /* |h| >= hmin                              */
  sunrealtype hmax_inv;           /* |h| <= 1/hmax_inv                        */
  sunrealtype hprime;             /* next actual step size to be used         */
  sunrealtype next_h;             /* next dynamical step size (only used in
                                  getCurrentStep); note that this could
                                  overtake tstop */
  sunrealtype eta;                /* eta = hprime / h                         */
  sunrealtype tcur;               /* current internal value of t
                                  (changes with each stage)                */
  sunrealtype tretlast;           /* value of tret last returned by ARKODE    */
  sunbooleantype fixedstep;       /* flag to disable temporal adaptivity      */
  ARKodeHAdaptMem hadapt_mem;  /* time step adaptivity structure           */


  /* Limits and various solver parameters */
  long int mxstep;         /* max number of internal steps for one user call */
  int      mxhnil;         /* max number of warning messages issued to the
                              user that t+h == t for the next internal step  */
  int      maxconstrfails; /* max number of constraint check failures        */
  int      maxnef;         /* max error test fails in one step               */
  int      maxncf;         /* max num alg. solver conv. fails in one step    */

  /* Counters */
  long int nst_attempts;  /* number of attempted steps                  */
  long int nst;           /* number of internal steps taken             */
  int      nhnil;         /* number of messages issued to the user that
                             t+h == t for the next iternal step         */
  long int ncfn;          /* num corrector convergence failures         */
  long int netf;          /* num error test failures                    */
  long int nconstrfails;  /* number of constraint failures              */

  /* Space requirements for ARKODE */
  sunindextype lrw1;        /* no. of sunrealtype words in 1 N_Vector          */
  sunindextype liw1;        /* no. of integer words in 1 N_Vector           */
  long int lrw;             /* no. of sunrealtype words in ARKODE work vectors */
  long int liw;             /* no. of integer words in ARKODE work vectors  */

  /* Saved Values */
  sunrealtype    h0u;          /* actual initial stepsize                     */
  sunrealtype    tn;           /* time of last successful step                */
  sunrealtype    terr;         /* error in tn for compensated sums            */
  sunrealtype    hold;         /* last successful h value used                */
  sunrealtype    tolsf;        /* tolerance scale factor (suggestion to user) */
  sunbooleantype VabstolMallocDone;
  sunbooleantype VRabstolMallocDone;
  sunbooleantype MallocDone;
  sunbooleantype initsetup;    /* denotes a call to InitialSetup is needed   */
  int         init_type;    /* initialization type (see constants above)  */
  sunbooleantype firststage;   /* denotes first stage in simulation          */
  sunbooleantype initialized;  /* denotes arkInitialSetup has been done      */
  sunbooleantype call_fullrhs; /* denotes the full RHS fn will be called     */

  /* Error handler function and error ouput file */
  ARKErrHandlerFn ehfun;    /* error messages are handled by ehfun        */
  void           *eh_data;  /* data pointer passed to ehfun               */
  FILE           *errfp;    /* ARKODE error messages are sent to errfp    */

  /* Rootfinding Data */
  ARKodeRootMem root_mem;          /* root-finding structure */

  /* Relaxation Data */
  sunbooleantype relax_enabled;    /* is relaxation enabled?    */
  ARKodeRelaxMem relax_mem;        /* relaxation data structure */

  /* User-supplied step solution post-processing function */
  ARKPostProcessFn ProcessStep;
  void*                ps_data; /* pointer to user_data */

  /* User-supplied stage solution post-processing function */
  ARKPostProcessFn ProcessStage;

  sunbooleantype use_compensated_sums;

  /* XBraid interface variables */
  sunbooleantype force_pass;  /* when true the step attempt loop will ignore the
                              return value (kflag) from arkCheckTemporalError
                              and set kflag = ARK_SUCCESS to force the step
                              attempt to always pass (if a solver failure did
                              not occur before the error test). */
  int         last_kflag;  /* last value of the return flag (kflag) from a call
                              to arkCheckTemporalError. This is only set when
                              force_pass is true and is used by the XBraid
                              interface to determine if a time step passed or
                              failed the time step error test.  */
};



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
  ARKODE linear solver interface, such as counters and statistics.
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

  Additional vectors that are set within the ARKODE memory
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
  t - the 'time' at which to setup the mass matrix
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
  sunrealtype client_tol contains the desired accuracy (in the wrms
  norm) of the routine calling the solver; the ARKDLS solver
  ignore this value and the ARKSPILS solver tightens it by the
  factor eplifac.  The solution is to be returned in the vector b.

  Additional vectors that are set within the ARKODE memory
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
  structure, and system linear solver type to the ARKODE time
  stepping module pointed to in ark_mem->step_mem.  This will
  be called by the ARKODE linear solver interface.

  This routine should return 0 if it has successfully attached
  these items and a negative value otherwise.  If an error does
  occur, an appropriate message should be sent to the ARKODE
  error handler function.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepAttachMasssolFn
  ---------------------------------------------------------------
  This routine should attach the various set of mass matrix
  linear solver interface routines, data structure, mass matrix
  type, and solver type to the ARKODE time stepping module
  pointed to in ark_mem->step_mem.  This will be called by the
  ARKODE linear solver interface.

  This routine should return 0 if it has successfully attached
  these items, and a negative value otherwise.  If an error does
  occur, an appropriate message should be sent to the ARKODE
  error handler function.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepDisableLSetup
  ---------------------------------------------------------------
  This routine should NULLify any ARKLinsolSetupFn function
  pointer stored in the ARKODE time stepping module (initially set
  in a call to ARKTimestepAttachLinsolFn).

  This routine has no return value.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepDisableMSetup
  ---------------------------------------------------------------
  This routine should NULLify any ARKMassSetupFn function pointer
  stored in the ARKODE time stepping module (initially set in a
  call to ARKTimestepAttachMasssolFn).

  This routine has no return value.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepGetLinMemFn
  ---------------------------------------------------------------
  This routine should return the linear solver memory structure
  used by the ARKODE time stepping module pointed to in
  ark_mem->step_mem.  This will be called by the ARKODE linear
  solver interface.

  This routine should return NULL if no linear solver memory
  structure is attached.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepGetMassMemFn
  ---------------------------------------------------------------
  This routine should return the mass matrix linear solver memory
  structure used by the ARKODE time stepping module pointed to in
  ark_mem->step_mem.  This will be called the ARKODE mass matrix
  solver interface.

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
  time step module internal sunbooleantype variable indicating
  whether the preconditioner is current, and a logic value
  indicating whether the gamma value is sufficiently stale
  to cause recomputation of Jacobian/preconditioner data.  Here,
  gamma is the coefficient preceding the RHS Jacobian
  matrix, J, in the full nonlinear system Jacobian,
  A = M - gamma*J.

  The time step module must contain a sunbooleantype variable to
  provide for the boolentype pointer (jcur).  This is only used
  by iterative linear solvers, so could be NULL for time step
  modules that only work with direct linear solvers.  Optionally,
  the value of this parameter could be set to SUNFALSE prior to
  return from the ARKTimestepGetGammasFn to force recalculation
  of preconditioner information.

  The value of the logic flag is used as follows:  if a previous
  Newton iteration failed due to a bad Jacobian/preconditioner,
  and this flag is SUNFALSE, this will trigger recalculation of
  the Jacobian/preconditioner.

  This routine should return 0 if it has successfully attached
  these items, and a negative value otherwise.  If an error does
  occur, an appropriate message should be sent to the ARKODE
  error handler function.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepInitFn
  ---------------------------------------------------------------
  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup. It should complete initializations for
  a specific ARKODE time stepping module, such as verifying
  compatibility of user-specified linear and nonlinear solver
  objects. The input init_type flag indicates if the call is
  for (re-)initializing, resizing, or resetting the problem.

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
  the sum of many ODE RHS functions (e.g. ARK, MRI).  The 'mode'
  indicates where this routine is called:

     ARK_FULLRHS_START -> called at the beginning of a simulation
                          i.e., at (tn, yn) = (t0, y0) or (tR, yR)

     ARK_FULLRHS_END   -> called at the end of a successful step i.e,
                          at (tcur, ycur) or the start of the subsequent
                          step i.e., at (tn, yn) = (tcur, ycur) from
                          the end of the last step

     ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  It is recommended that the stepper use the mode information to
  maximize reuse between calls to this function and RHS
  evaluations inside the stepper itself.

  This routine is only required to be supplied to ARKODE if:
  * ARKODE's initial time step selection algorithm is used,
  * the user requests temporal root-finding,
  * the Hermite interpolation module is used, or
  * the user requests the "bootstrap" implicit predictor.
  Note that any stepper can itself require that this routine
  exist for its own internal business (e.g., ERKStep).

  This routine should return 0 if successful, and a negative value
  otherwise.  If an error does occur, an appropriate message
  should be sent to the error handler function.
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  ARKTimestepStepFn
  ---------------------------------------------------------------
  This routine serves the primary purpose of any ARKODE
  time-stepping module: it performs a single time step of the
  method (with embedding, if possible).

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
  - user_data -- the (void *) pointer returned to user for
    RHS calls
  - report / diagfp -- if any diagnostic information is
    to be saved to disk, the report flag indicates whether
    this is enabled, and diagfp provides the file pointer
    where this information should be written

  The output variable dsmPtr should contain estimate of the
  weighted local error if an embedding is present; otherwise it
  should be 0.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  At the start of a new
  time step, this will initially have the value FIRST_CALL.  On
  return from this function, nflagPtr should have a value:
            0 => algebraic solve completed successfully
           >0 => solve did not converge at this step size
                 (but may with a smaller stepsize)
           <0 => solve encountered an unrecoverable failure

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/


/*===============================================================
  ARKODE PROTOTYPE FUNCTIONS (MAY BE REPLACED BY USER)
  ===============================================================*/

/* Prototype of internal rwtSet function */
int arkRwtSet(N_Vector ycur, N_Vector weight, void *data);

/* Prototype of internal errHandler function */
void arkErrHandler(int error_code, const char *module,
                   const char *function, char *msg, void *data);

/* Prototype of internal explicit stability estimation function */
int arkExpStab(N_Vector y, sunrealtype t, sunrealtype *hstab, void *user_data);

/*===============================================================
  HIGH LEVEL ERROR HANDLER, USED THROUGHOUT ARKODE
  ===============================================================*/

void arkProcessError(ARKodeMem ark_mem, int error_code, int line, const char *func,
                     const char* file, const char *msgfmt, ...);

/*===============================================================
  ARKODE PRIVATE FUNCTION PROTOTYPES
  ===============================================================*/
#ifdef __GNUC__
#define SUNDIALS_UNUSED __attribute__ ((unused))
#else
#define SUNDIALS_UNUSED
#endif

int arkInit(ARKodeMem ark_mem, sunrealtype t0, N_Vector y0, int init_type);
sunbooleantype arkAllocVec(ARKodeMem ark_mem, N_Vector tmpl, N_Vector *v);
sunbooleantype arkAllocVecArray(int count, N_Vector tmpl, N_Vector **v,
                             sunindextype lrw1, long int *lrw,
                             sunindextype liw1, long int *liw);
void arkFreeVec(ARKodeMem ark_mem, N_Vector *v);
void arkFreeVecArray(int count, N_Vector **v,
                     sunindextype lrw1, long int *lrw,
                     sunindextype liw1, long int *liw);
sunbooleantype arkResizeVec(ARKodeMem ark_mem,
                         ARKVecResizeFn resize,
                         void *resize_data,
                         sunindextype lrw_diff,
                         sunindextype liw_diff,
                         N_Vector tmpl,
                         N_Vector *v);
sunbooleantype arkResizeVecArray(ARKVecResizeFn resize, void *resize_data,
                              int count, N_Vector tmpl, N_Vector **v,
                              sunindextype lrw_diff, long int *lrw,
                              sunindextype liw_diff, long int *liw);
void arkPrintMem(ARKodeMem ark_mem, FILE *outfile);
sunbooleantype arkCheckTimestepper(ARKodeMem ark_mem);
sunbooleantype arkCheckNvector(N_Vector tmpl);
sunbooleantype arkAllocVectors(ARKodeMem ark_mem,
                            N_Vector tmpl);
sunbooleantype arkResizeVectors(ARKodeMem ark_mem,
                             ARKVecResizeFn resize,
                             void *resize_data,
                             sunindextype lrw_diff,
                             sunindextype liw_diff,
                             N_Vector tmpl);
void arkFreeVectors(ARKodeMem ark_mem);

int arkInitialSetup(ARKodeMem ark_mem, sunrealtype tout);
int arkStopTests(ARKodeMem ark_mem, sunrealtype tout, N_Vector yout,
                 sunrealtype *tret, int itask, int *ier);
int arkHin(ARKodeMem ark_mem, sunrealtype tout);
sunrealtype arkUpperBoundH0(ARKodeMem ark_mem,
                         sunrealtype tdist);
int arkYddNorm(ARKodeMem ark_mem, sunrealtype hg,
               sunrealtype *yddnrm);

int arkCompleteStep(ARKodeMem ark_mem, sunrealtype dsm);
int arkHandleFailure(ARKodeMem ark_mem,int flag);

int arkEwtSetSS(N_Vector ycur, N_Vector weight, void* arkode_mem);
int arkEwtSetSV(N_Vector ycur, N_Vector weight, void* arkode_mem);
int arkEwtSetSmallReal(N_Vector ycur, N_Vector weight, void* arkode_mem);
int arkRwtSetSS(ARKodeMem ark_mem, N_Vector My,
                N_Vector weight);
int arkRwtSetSV(ARKodeMem ark_mem, N_Vector My,
                N_Vector weight);

ARKodeMem arkCreate(SUNContext sunctx);
int arkResize(ARKodeMem ark_mem, N_Vector ynew, sunrealtype hscale,
              sunrealtype t0, ARKVecResizeFn resize, void *resize_data);
int arkSStolerances(ARKodeMem ark_mem, sunrealtype reltol, sunrealtype abstol);
int arkSVtolerances(ARKodeMem ark_mem, sunrealtype reltol, N_Vector abstol);
int arkWFtolerances(ARKodeMem ark_mem, ARKEwtFn efun);
int arkResStolerance(ARKodeMem ark_mem, sunrealtype rabstol);
int arkResVtolerance(ARKodeMem ark_mem, N_Vector rabstol);
int arkResFtolerance(ARKodeMem ark_mem, ARKRwtFn rfun);
int arkRootInit(ARKodeMem ark_mem, int nrtfn, ARKRootFn g);
int arkEvolve(ARKodeMem ark_mem, sunrealtype tout, N_Vector yout,
              sunrealtype *tret, int itask);
int arkGetDky(ARKodeMem ark_mem, sunrealtype t, int k, N_Vector dky);
void arkFree(void **arkode_mem);

int arkWriteParameters(ARKodeMem ark_mem, FILE *fp);
int arkPredict_MaximumOrder(ARKodeMem ark_mem, sunrealtype tau,
                            N_Vector yguess);
int arkPredict_VariableOrder(ARKodeMem ark_mem, sunrealtype tau,
                             N_Vector yguess);
int arkPredict_CutoffOrder(ARKodeMem ark_mem, sunrealtype tau,
                           N_Vector yguess);
int arkPredict_Bootstrap(ARKodeMem ark_mem, sunrealtype hj,
                         sunrealtype tau, int nvec, sunrealtype *cvals,
                         N_Vector *Xvecs, N_Vector yguess);
int arkCheckConvergence(ARKodeMem ark_mem, int *nflagPtr, int *ncfPtr);
int arkCheckConstraints(ARKodeMem ark_mem, int *nflag, int *constrfails);
int arkCheckTemporalError(ARKodeMem ark_mem, int *nflagPtr, int *nefPtr,
                          sunrealtype dsm);
int arkAccessHAdaptMem(void* arkode_mem, const char *fname,
                       ARKodeMem *ark_mem, ARKodeHAdaptMem *hadapt_mem);

int arkSetAdaptController(void *arkode_mem, SUNAdaptController C);
int arkSetDefaults(void *arkode_mem);
int arkSetDenseOrder(void *arkode_mem, int dord);
int arkSetInterpolantType(void *arkode_mem, int itype);
int arkSetInterpolantDegree(void *arkode_mem, int degree);
int arkSetErrHandlerFn(void *arkode_mem,
                       ARKErrHandlerFn ehfun,
                       void *eh_data);
int arkSetErrFile(void *arkode_mem, FILE *errfp);
int arkSetUserData(void *arkode_mem, void *user_data);
int arkSetMaxNumSteps(void *arkode_mem, long int mxsteps);
int arkSetMaxHnilWarns(void *arkode_mem, int mxhnil);
int arkSetInitStep(void *arkode_mem, sunrealtype hin);
int arkSetMinStep(void *arkode_mem, sunrealtype hmin);
int arkSetMaxStep(void *arkode_mem, sunrealtype hmax);
int arkSetStopTime(void *arkode_mem, sunrealtype tstop);
int arkSetInterpolateStopTime(void *arkode_mem, sunbooleantype interp);
int arkClearStopTime(void *arkode_mem);
int arkSetFixedStep(void *arkode_mem, sunrealtype hfixed);
int arkSetRootDirection(void *arkode_mem, int *rootdir);
int arkSetNoInactiveRootWarn(void *arkode_mem);
int arkSetPostprocessStepFn(void *arkode_mem,
                            ARKPostProcessFn ProcessStep);
int arkSetPostprocessStageFn(void *arkode_mem,
                             ARKPostProcessFn ProcessStage);
int arkSetConstraints(void *arkode_mem, N_Vector constraints);
int arkSetMaxNumConstrFails(void *arkode_mem, int maxfails);
int arkSetAdaptivityAdjustment(void *arkode_mem, int adjust);
int arkSetCFLFraction(void *arkode_mem, sunrealtype cfl_frac);
int arkSetSafetyFactor(void *arkode_mem, sunrealtype safety);
int arkSetErrorBias(void *arkode_mem, sunrealtype bias);
int arkSetMaxGrowth(void *arkode_mem, sunrealtype mx_growth);
int arkSetMinReduction(void *arkode_mem, sunrealtype eta_min);
int arkSetFixedStepBounds(void *arkode_mem, sunrealtype lb, sunrealtype ub);
int arkSetAdaptivityMethod(void *arkode_mem, int imethod, int idefault,
                           int pq, sunrealtype adapt_params[3]);
int arkSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun, void *h_data);
int arkSetMaxFirstGrowth(void *arkode_mem, sunrealtype etamx1);
int arkSetMaxEFailGrowth(void *arkode_mem, sunrealtype etamxf);
int arkSetSmallNumEFails(void *arkode_mem, int small_nef);
int arkSetMaxCFailGrowth(void *arkode_mem, sunrealtype etacf);
int arkSetStabilityFn(void *arkode_mem, ARKExpStabFn EStab, void *estab_data);
int arkSetMaxErrTestFails(void *arkode_mem, int maxnef);
int arkSetMaxConvFails(void *arkode_mem, int maxncf);
int arkSetUseCompensatedSums(void *arkode_mem, sunbooleantype onoff);
int arkGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw);
int arkGetNumStepAttempts(void *arkode_mem, long int *nstep_attempts);
int arkGetNumSteps(void *arkode_mem, long int *nsteps);
int arkGetActualInitStep(void *arkode_mem, sunrealtype *hinused);
int arkGetLastStep(void *arkode_mem, sunrealtype *hlast);
int arkGetCurrentStep(void *arkode_mem, sunrealtype *hcur);
int arkGetCurrentState(void *arkode_mem, N_Vector *ycur);
int arkGetCurrentTime(void *arkode_mem, sunrealtype *tcur);
int arkGetTolScaleFactor(void *arkode_mem, sunrealtype *tolsfac);
int arkGetErrWeights(void *arkode_mem, N_Vector eweight);
int arkGetResWeights(void *arkode_mem, N_Vector rweight);
int arkGetNumGEvals(void *arkode_mem, long int *ngevals);
int arkGetRootInfo(void *arkode_mem, int *rootsfound);
int arkGetNumConstrFails(void *arkode_mem, long int *nconstrfails);
int arkGetNumExpSteps(void *arkode_mem, long int *nsteps);
int arkGetNumAccSteps(void *arkode_mem, long int *nsteps);
int arkGetNumErrTestFails(void *arkode_mem, long int *netfails);
int arkGetNumStepSolveFails(void *arkode_mem, long int *nncfails);
int arkGetStepStats(void *arkode_mem, long int *nsteps,
                    sunrealtype *hinused, sunrealtype *hlast,
                    sunrealtype *hcur, sunrealtype *tcur);
int arkGetUserData(void *arkode_mem, void** user_data);
int arkPrintAllStats(void *arkode_mem, FILE *outfile,
                     SUNOutputFormat fmt);
char *arkGetReturnFlagName(long int flag);

ARKODE_DIRKTableID arkButcherTableDIRKNameToID(const char *imethod);
ARKODE_ERKTableID arkButcherTableERKNameToID(const char *emethod);

/* XBraid interface functions */
int arkSetForcePass(void *arkode_mem, sunbooleantype force_pass);
int arkGetLastKFlag(void *arkode_mem, int *last_kflag);


/*===============================================================
  Reusable ARKODE Error Messages
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
#define MSG_ARK_NO_MEM         "arkode_mem = NULL illegal."
#define MSG_ARK_ARKMEM_FAIL    "Allocation of arkode_mem failed."
#define MSG_ARK_MEM_FAIL       "A memory request failed."
#define MSG_ARK_NO_MALLOC      "Attempt to call before ARKodeInit."
#define MSG_ARK_BAD_HMIN_HMAX  "Inconsistent step size limits: hmin > hmax."
#define MSG_ARK_BAD_RELTOL     "reltol < 0 illegal."
#define MSG_ARK_BAD_ABSTOL     "abstol has negative component(s) (illegal)."
#define MSG_ARK_NULL_ABSTOL    "abstol = NULL illegal."
#define MSG_ARK_BAD_RABSTOL    "rabstol has negative component(s) (illegal)."
#define MSG_ARK_NULL_RABSTOL   "rabstol = NULL illegal."
#define MSG_ARK_NULL_Y0        "y0 = NULL illegal."
#define MSG_ARK_Y0_FAIL_CONSTR "y0 fails to satisfy constraints."
#define MSG_ARK_NULL_F         "Must specify at least one of fe, fi (both NULL)."
#define MSG_ARK_NULL_G         "g = NULL illegal."
#define MSG_ARK_BAD_NVECTOR    "A required vector operation is not implemented."
#define MSG_ARK_BAD_CONSTR     "Illegal values in constraints vector."
#define MSG_ARK_NULL_DKY       "dky = NULL illegal."
#define MSG_ARK_BAD_T          "Illegal value for t. " MSG_TIME_INT
#define MSG_ARK_NO_ROOT        "Rootfinding was not initialized."

/* ARKODE Error Messages */
#define MSG_ARK_YOUT_NULL      "yout = NULL illegal."
#define MSG_ARK_TRET_NULL      "tret = NULL illegal."
#define MSG_ARK_BAD_EWT        "Initial ewt has component(s) equal to zero (illegal)."
#define MSG_ARK_EWT_NOW_BAD    "At " MSG_TIME ", a component of ewt has become <= 0."
#define MSG_ARK_BAD_RWT        "Initial rwt has component(s) equal to zero (illegal)."
#define MSG_ARK_RWT_NOW_BAD    "At " MSG_TIME ", a component of rwt has become <= 0."
#define MSG_ARK_BAD_ITASK      "Illegal value for itask."
#define MSG_ARK_BAD_H0         "h0 and tout - t0 inconsistent."
#define MSG_ARK_BAD_TOUT       "Trouble interpolating at " MSG_TIME_TOUT ". tout too far back in direction of integration"
#define MSG_ARK_EWT_FAIL       "The user-provide EwtSet function failed."
#define MSG_ARK_EWT_NOW_FAIL   "At " MSG_TIME ", the user-provide EwtSet function failed."
#define MSG_ARK_RWT_FAIL       "The user-provide RwtSet function failed."
#define MSG_ARK_RWT_NOW_FAIL   "At " MSG_TIME ", the user-provide RwtSet function failed."
#define MSG_ARK_LINIT_FAIL     "The linear solver's init routine failed."
#define MSG_ARK_HNIL_DONE      "The above warning has been issued mxhnil times and will not be issued again for this problem."
#define MSG_ARK_TOO_CLOSE      "tout too close to t0 to start integration."
#define MSG_ARK_MAX_STEPS      "At " MSG_TIME ", mxstep steps taken before reaching tout."
#define MSG_ARK_TOO_MUCH_ACC   "At " MSG_TIME ", too much accuracy requested."
#define MSG_ARK_HNIL           "Internal " MSG_TIME_H " are such that t + h = t on the next step. The solver will continue anyway."
#define MSG_ARK_ERR_FAILS      "At " MSG_TIME_H ", the error test failed repeatedly or with |h| = hmin."
#define MSG_ARK_CONV_FAILS     "At " MSG_TIME_H ", the solver convergence test failed repeatedly or with |h| = hmin."
#define MSG_ARK_SETUP_FAILED   "At " MSG_TIME ", the setup routine failed in an unrecoverable manner."
#define MSG_ARK_SOLVE_FAILED   "At " MSG_TIME ", the solve routine failed in an unrecoverable manner."
#define MSG_ARK_FAILED_CONSTR  "At " MSG_TIME ", unable to satisfy inequality constraints."
#define MSG_ARK_RHSFUNC_FAILED "At " MSG_TIME ", the right-hand side routine failed in an unrecoverable manner."
#define MSG_ARK_RHSFUNC_UNREC  "At " MSG_TIME ", the right-hand side failed in a recoverable manner, but no recovery is possible."
#define MSG_ARK_RHSFUNC_REPTD  "At " MSG_TIME " repeated recoverable right-hand side function errors."
#define MSG_ARK_RTFUNC_FAILED  "At " MSG_TIME ", the rootfinding routine failed in an unrecoverable manner."
#define MSG_ARK_CLOSE_ROOTS    "Root found at and very near " MSG_TIME "."
#define MSG_ARK_BAD_TSTOP      "The value " MSG_TIME_TSTOP " is behind current " MSG_TIME " in the direction of integration."
#define MSG_ARK_INACTIVE_ROOTS "At the end of the first step, there are still some root functions identically 0. This warning will not be issued again."
#define MSG_ARK_RESIZE_FAIL    "Error in user-supplied resize() function."
#define MSG_ARK_MASSINIT_FAIL  "The mass matrix solver's init routine failed."
#define MSG_ARK_MASSSETUP_FAIL "The mass matrix solver's setup routine failed."
#define MSG_ARK_MASSSOLVE_FAIL "The mass matrix solver failed."
#define MSG_ARK_NLS_FAIL       "At " MSG_TIME " the nonlinear solver failed in an unrecoverable manner."
#define MSG_ARK_USER_PREDICT_FAIL "At " MSG_TIME " the user-supplied predictor failed in an unrecoverable manner."
#define MSG_ARKADAPT_NO_MEM    "Adaptivity memory structure not allocated."
#define MSG_ARK_VECTOROP_ERR      "At " MSG_TIME ", a vector operation failed."
#define MSG_ARK_INNERSTEP_FAILED  "At " MSG_TIME ", the inner stepper failed in an unrecoverable manner."
#define MSG_ARK_POSTPROCESS_STEP_FAIL "At " MSG_TIME ", the step postprocessing routine failed in an unrecoverable manner."
#define MSG_ARK_POSTPROCESS_STAGE_FAIL "At " MSG_TIME ", the stage postprocessing routine failed in an unrecoverable manner."
#define MSG_ARK_NULL_SUNCTX "sunctx = NULL illegal."
#define MSG_ARK_CONTEXT_MISMATCH "Outer and inner steppers have different contexts."
#define MSG_ARK_MISSING_FULLRHS "Time-stepping module missing fullrhs routine (required by requested solver configuration)."
#define MSG_ARK_INTERPOLATION_FAIL "At " MSG_TIME ", interpolating the solution failed."

#ifdef __cplusplus
}
#endif

#endif
