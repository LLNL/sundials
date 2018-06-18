/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * Based on arkode_arkstep_impl.h written by Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2018, Southern Methodist University and 
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
 * -----------------------------------------------------------------------------
 * Implementation header file for ARKode's Implicit-Explicit (IMEX) Generalized
 * Additive Runge-Kutta time stepper module.
 * ---------------------------------------------------------------------------*/

#ifndef _ARKODE_IMEXGARKSTEP_IMPL_H
#define _ARKODE_IMEXGARKSTEP_IMPL_H

#include "arkode/arkode.h"
#include "arkode/arkode_imexgarkstep.h"
#include "arkode_impl.h"

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

  /* IMEX GARK problem specification -- My' = fe(t,y) + fi(t,y) */
  ARKRhsFn    fe;              /* Explicit Rhs function          */
  ARKRhsFn    fi;              /* Implicit Rhs function          */
  booleantype linear;          /* SUNTRUE if fi is linear        */
  booleantype linear_timedep;  /* SUNTRUE if dfi/dy depends on t */

  /* IMEX GARK method storage and parameters */
  N_Vector *Fe;           /* explicit RHS at each stage                    */
  N_Vector *Fi;           /* implicit RHS at each stage                    */
  N_Vector sdata;         /* old stage data in residual                    */
  N_Vector zpred;         /* predicted stage solution                      */
  int q;                  /* method order                                  */
  int p;                  /* embedding order                               */
  int istage;             /* current stage                                 */
  int stages;             /* number of stages                              */
  ARKodeButcherTable Bee; /* Explicit Butcher table                        */
  ARKodeButcherTable Bei; /* Coupling explicit with implicit Butcher table */
  ARKodeButcherTable Bie; /* Coupling implicit with explicit Butcher table */
  ARKodeButcherTable Bii; /* Implicit Butcher table                        */

  /* Time step adaptivity data */
  ARKodeHAdaptMem hadapt_mem;  /* time step adaptivity structure   */
  booleantype     hadapt_pq;   /* choice of using p (0) vs q (1)   */
  int             maxnef;      /* max error test fails in one step */

  /* (Non)Linear solver parameters & data */
  realtype gamma;      /* gamma = h * A(i,i)                       */
  realtype gammap;     /* gamma at the last setup call             */
  realtype gamrat;     /* gamma / gammap                           */
  realtype dgmax;      /* call lsetup if |gamma/gammap-1| >= dgmax */

  int      predictor;  /* implicit prediction method to use        */
  realtype crdown;     /* nonlinear conv rate estimation constant  */
  realtype rdiv;       /* nonlin divergence if del/delp > rdiv     */
  realtype crate;      /* estimated nonlin convergence rate        */
  realtype eRNrm;      /* estimated residual norm, used in nonlin
                          and linear solver convergence tests      */
  realtype nlscoef;    /* coefficient in nonlin. convergence test  */
  int      mnewt;      /* internal Newton iteration counter        */

  int      msbp;       /* positive => max # steps between lsetup 
                          negative => call at each Newton iter     */
  long int nstlp;      /* step number of last setup call           */

  int      maxcor;     /* max num iter. for solving nonlin. eq.    */
  int      maxncf;     /* max num nonlin. conv. fails in one step  */

  booleantype jcur;    /* is Jacobian info for lin solver current? */
  
  /* Fixed-point Solver Data */
  booleantype use_fp;  /* flag for fixed-point solver vs Newton    */
  long int    fp_m;    /* number of vectors to use in acceleration */
  ARKodeFPMem fp_mem;  /* accelerated fixed-point solver structure */

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
  realtype        msetuptime;  /* "t" value at last msetup call                   */
  int             msolve_type; /* interface type: 0=iterative; 1=direct; 2=custom */

  /* Counters */
  long int nst_attempts;  /* num attempted steps                */
  long int nfe;           /* num fe calls                       */
  long int nfi;           /* num fi calls                       */
  long int ncfn;          /* num corrector convergence failures */
  long int netf;          /* num error test failures            */
  long int nni;           /* num Newton iterations performed    */
  long int nsetups;       /* num setup calls                    */

  /* Reusable arrays for fused vector operations */
  realtype *cvals;
  N_Vector *Xvecs;

} *ARKodeIMEXGARKStepMem;


/* =============================================================================
 * ARK time step module private function prototypes
 * ===========================================================================*/

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
int imexgarkStep_Init(void* arkode_mem);
void* imexgarkStep_GetLmem(void* arkode_mem);
void* imexgarkStep_GetMassMem(void* arkode_mem);
ARKRhsFn imexgarkStep_GetImplicitRHS(void* arkode_mem);
int imexgarkStep_GetGammas(void* arkode_mem, realtype *gamma,
                           realtype *gamrat, booleantype **jcur,
                           booleantype *dgamma_fail);
int imexgarkStep_FullRHS(void* arkode_mem, realtype t, 
                         N_Vector y, N_Vector f, int mode);
int imexgarkStep_Step(void* arkode_mem);
int imexgarkStep_Resize(void* arkode_mem, ARKVecResizeFn resize,
                        void *resize_data, sunindextype lrw_diff,
                        sunindextype liw_diff, N_Vector tmpl);
void imexgarkStep_PrintMem(void* arkode_mem, FILE* outfile);
int imexgarkStep_Free(void* arkode_mem);

/* Internal utility routines */  
booleantype imexgarkStep_CheckNVector(N_Vector tmpl);
int imexgarkStep_SetButcherTables(ARKodeMem ark_mem);
int imexgarkStep_CheckButcherTables(ARKodeMem ark_mem);
int imexgarkStep_Predict(ARKodeMem ark_mem, int istage, N_Vector yguess);
int imexgarkStep_StageSetup(ARKodeMem ark_mem);
int imexgarkStep_Nls(ARKodeMem ark_mem, int nflag);
int imexgarkStep_NlsResid(ARKodeMem ark_mem, N_Vector y, 
                          N_Vector fy, N_Vector r);
int imexgarkStep_NlsNewton(ARKodeMem ark_mem, int nflag);
int imexgarkStep_NlsAccelFP(ARKodeMem ark_mem, int nflag);
int imexgarkStep_AndersonAcc(ARKodeMem ark_mem, N_Vector gval, 
                             N_Vector fv, N_Vector x, N_Vector xold, 
                             int iter, realtype *R, realtype *gamma);
int imexgarkStep_Ls(ARKodeMem ark_mem, int nflag);
int imexgarkStep_HandleNFlag(ARKodeMem ark_mem, int *nflagPtr, int *ncfPtr);
  
int imexgarkStep_ComputeSolutions(ARKodeMem ark_mem, realtype *dsm);
int imexgarkStep_DoErrorTest(ARKodeMem ark_mem, int *nflagPtr,
                             int *nefPtr, realtype dsm);
int imexgarkStep_PrepareNextStep(ARKodeMem ark_mem, realtype dsm);

/* =============================================================================
  Reusable ARKStep Error Messages
 * ===========================================================================*/

/* Initialization and I/O error messages */
#define MSG_IMEXGARK_NO_STEP_MEM   "Time step module memory is NULL."
#define MSG_IMEXGARK_NO_ADAPT_MEM  "Adaptivity memory structure not allocated."


#ifdef __cplusplus
}
#endif

#endif
