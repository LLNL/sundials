/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Implementation header file for ARKode's Implicit-Explicit (IMEX) Generalized
 * Additive Runge-Kutta time stepper module.
 * ---------------------------------------------------------------------------*/

#ifndef _ARKODE_IMEXGARKSTEP_IMPL_H
#define _ARKODE_IMEXGARKSTEP_IMPL_H

#include "arkode/arkode_imexgarkstep.h"
#include "arkode_impl.h"
#include "arkode_ls_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* =============================================================================
 * IMEXGARK time step module data structure
 * ===========================================================================*/

/* -----------------------------------------------------------------------------
 * Types : struct ARKodeIMEXGARKStepMemRec, ARKodeIMEXGARKStepMem
 * -----------------------------------------------------------------------------
 * The type ARKodeARKStepMem is type pointer to struct
 * ARKodeARKStepMemRec.  This structure contains fields to
 * perform an additive Runge-Kutta time step.
 * ---------------------------------------------------------------------------*/

typedef struct ARKodeIMEXGARKStepMemRec {

  /* Problem specification: My' = fe(t,y) + fi(t,y) */
  ARKRhsFn    fe;              /* Explicit Rhs function          */
  ARKRhsFn    fi;              /* Implicit Rhs function          */
  booleantype linear;          /* SUNTRUE if fi is linear        */
  booleantype linear_timedep;  /* SUNTRUE if dfi/dy depends on t */

  /* Method storage and parameters */
  N_Vector *Fe;           /* explicit RHS at each stage */
  N_Vector *Fi;           /* implicit RHS at each stage */
  N_Vector sdata;         /* old stage data in residual */
  N_Vector zpred;         /* predicted stage solution   */
  N_Vector zcor;          /* stage correction           */
  int q;                  /* method order               */
  int p;                  /* embedding order            */
  int istage;             /* current stage              */
  int stages;             /* number of stages           */
  ARKodeButcherTable Bee; /* Explicit Butcher table                */
  ARKodeButcherTable Bei; /* Coupling table explicit with implicit */
  ARKodeButcherTable Bie; /* Coupling table implicit with explicit */
  ARKodeButcherTable Bii; /* Implicit Butcher table                */

  /* Time step adaptivity data */
  ARKodeHAdaptMem hadapt_mem;  /* time step adaptivity structure   */
  booleantype     hadapt_pq;   /* choice of using p (0) vs q (1)   */
  int             maxnef;      /* max error test fails in one step */

  /* (Non)Linear solver parameters & data */
  SUNNonlinearSolver NLS;      /* generic SUNNonlinearSolver object     */
  booleantype        ownNLS;   /* flag indicating ownership of NLS      */

  realtype gamma;    /* gamma = h * A(i,i)                       */
  realtype gammap;   /* gamma at the last setup call             */
  realtype gamrat;   /* gamma / gammap                           */
  realtype dgmax;    /* call lsetup if |gamma/gammap-1| >= dgmax */

  int      predictor;    /* implicit prediction method to use        */
  realtype crdown;       /* nonlinear conv rate estimation constant  */
  realtype rdiv;         /* nonlin divergence if del/delp > rdiv     */
  realtype crate;        /* estimated nonlin convergence rate        */
  realtype delp;         /* norm of previous nonlinear solver update */
  realtype eRNrm;        /* estimated residual norm, used in nonlin
                            and linear solver convergence tests      */
  realtype nlscoef;      /* coefficient in nonlin. convergence test  */
  int      mnewt;        /* internal Newton iteration counter        */

  int      msbp;         /* positive => max # steps between lsetup
                            negative => call at each Newton iter     */
  long int nstlp;        /* step number of last setup call           */

  int      maxcor;       /* max num iterations for solving the
                            nonlinear equation                       */
  int      maxncf;       /* max num nonlin. conv. fails in one step  */

  int      convfail;     /* NLS fail flag (for interface routines)   */
  booleantype jcur;      /* is Jacobian info for lin solver current? */

  /* Linear Solver Data */
  ARKLinsolInitFn  linit;
  ARKLinsolSetupFn lsetup;
  ARKLinsolSolveFn lsolve;
  ARKLinsolFreeFn  lfree;
  void             *lmem;
  int              lsolve_type;  /* interface type: 0=iterative; 1=direct; 2=custom */

  /* Mass matrix solver data */
  ARKMassInitFn   minit;
  ARKMassSetupFn  msetup;
  ARKMassMultFn   mmult;
  ARKMassSolveFn  msolve;
  ARKMassFreeFn   mfree;
  void*           mass_mem;
  realtype        msetuptime;  /* "t" value at last msetup call */
  int             msolve_type; /* interface type: 0=iterative; 1=direct; 2=custom */

  /* Counters */
  long int nst_attempts;  /* num attempted steps                */
  long int nfe;           /* num fe calls                       */
  long int nfi;           /* num fi calls                       */
  long int ncfn;          /* num corrector convergence failures */
  long int netf;          /* num error test failures            */
  long int nsetups;       /* num setup calls                    */

  /* Reusable arrays for fused vector operations */
  realtype *cvals;
  N_Vector *Xvecs;

} *ARKodeIMEXGARKStepMem;


/* ===============================================================
 * ARK time step module private function prototypes
 * ===============================================================*/

/* Interface routines supplied to ARKode */
int imexgarkStep_AttachLinsol(void* arkode_mem, ARKLinsolInitFn linit,
                              ARKLinsolSetupFn lsetup,
                              ARKLinsolSolveFn lsolve,
                              ARKLinsolFreeFn lfree,
                              int lsolve_type, void *lmem);
int imexgarkStep_AttachMasssol(void* arkode_mem, ARKMassInitFn minit,
                               ARKMassSetupFn msetup,
                               ARKMassMultFn mmult,
                               ARKMassSolveFn msolve,
                               ARKMassFreeFn lfree,
                               int msolve_type, void *mass_mem);
void imexgarkStep_DisableLSetup(void* arkode_mem);
void imexgarkStep_DisableMSetup(void* arkode_mem);
int imexgarkStep_Init(void* arkode_mem, int init_type);
void* imexgarkStep_GetLmem(void* arkode_mem);
void* imexgarkStep_GetMassMem(void* arkode_mem);
ARKRhsFn imexgarkStep_GetImplicitRHS(void* arkode_mem);
int imexgarkStep_GetGammas(void* arkode_mem, realtype *gamma,
                           realtype *gamrat, booleantype **jcur,
                           booleantype *dgamma_fail);
int imexgarkStep_FullRHS(void* arkode_mem, realtype t,
                         N_Vector y, N_Vector f, int mode);
int imexgarkStep_TakeStep(void* arkode_mem);

/* Internal utility routines */
int imexgarkStep_AccessStepMem(void* arkode_mem, const char *fname,
                               ARKodeMem *ark_mem, ARKodeIMEXGARKStepMem *step_mem);
booleantype imexgarkStep_CheckNVector(N_Vector tmpl);
int imexgarkStep_SetButcherTables(ARKodeMem ark_mem);
int imexgarkStep_CheckButcherTables(ARKodeMem ark_mem);
int imexgarkStep_Predict(ARKodeMem ark_mem, int istage, N_Vector yguess);
int imexgarkStep_StageSetup(ARKodeMem ark_mem);
int imexgarkStep_NlsInit(ARKodeMem ark_mem);
int imexgarkStep_Nls(ARKodeMem ark_mem, int nflag);
int imexgarkStep_HandleNFlag(ARKodeMem ark_mem, int *nflagPtr, int *ncfPtr);

int imexgarkStep_ComputeSolutions(ARKodeMem ark_mem, realtype *dsm);
int imexgarkStep_DoErrorTest(ARKodeMem ark_mem, int *nflagPtr,
                             int *nefPtr, realtype dsm);
int imexgarkStep_PrepareNextStep(ARKodeMem ark_mem, realtype dsm);

/* private functions passed to nonlinear solver */
int imexgarkStep_NlsResidual(N_Vector yy, N_Vector res, void* arkode_mem);
int imexgarkStep_NlsFPFunction(N_Vector yy, N_Vector res, void* arkode_mem);
int imexgarkStep_NlsLSetup(N_Vector yy, N_Vector res, booleantype jbad,
                           booleantype* jcur, void* arkode_mem);
int imexgarkStep_NlsLSolve(N_Vector yy, N_Vector delta, void* arkode_mem);
int imexgarkStep_NlsConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                             realtype tol, N_Vector ewt, void* arkode_mem);

/* ===============================================================
 * Reusable ARKStep Error Messages
 * ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_IMEXGARKSTEP_NO_MEM    "Time step module memory is NULL."
#define MSG_NLS_INIT_FAIL          "The nonlinear solver's init routine failed."

#ifdef __cplusplus
}
#endif

#endif
