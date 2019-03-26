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
 * Implementation header file for ARKode's MRI time stepper module.
 * ---------------------------------------------------------------------------*/

#ifndef _ARKODE_MRISTEP_IMPL_H
#define _ARKODE_MRISTEP_IMPL_H

#include "arkode/arkode_mristep.h"
#include "arkode/arkode_arkstep.h"

#include "arkode_impl.h"
#include "arkode_arkstep_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
  MRI function types
  ===============================================================*/

typedef int (*MRIEvolveInner)(void* inner_arkode_mem, realtype t0,
                              N_Vector y0, realtype tout);


/*===============================================================
  MRI time step module data structure
  ===============================================================*/

/*---------------------------------------------------------------
  The type ARKodeMRIStepMem is type pointer to struct
  ARKodeMRIStepMemRec. This structure contains fields to
  perform a MRI time step.
  ---------------------------------------------------------------*/
typedef struct ARKodeMRIStepMemRec {

  /* MRI problem specification */
  ARKRhsFn fs;   /* y' = fs(t,y) + ff1(t,y) + ff2(t,y)  */
  ARKRhsFn ff1;  /* inner forcing is added to this RHS  */
  ARKRhsFn ff2;  /* inner forcing is NOT added this RHS */

  /* Outer RK method storage and parameters */
  N_Vector *F;            /* slow RHS at each stage */
  int q;                  /* method order           */
  int p;                  /* embedding order        */
  int stages;             /* number of stages       */
  ARKodeButcherTable B;   /* MRI Butcher table      */

  /* Inner stepper data */
  void           *inner_mem;          /* inner stepper memory            */
  N_Vector        inner_forcing;      /* RHS forcing vector              */
  int             inner_retval;       /* last inner stepper return value */
  MRISTEP_ID      inner_stepper_id;   /* inner stepper identifier        */
  MRIEvolveInner  inner_evolve;       /* inner stepper evolve function   */

  /* Wrappers for user-supplied inner stepper functions */
  ARKEwtFn             inner_ewtfn;
  ARKPostProcessStepFn inner_postprocessstepfn;
  ARKLsJacFn           inner_jacfn;
  ARKLsJacTimesSetupFn inner_jactimessetupfn;
  ARKLsJacTimesVecFn   inner_jactimesvecfn;
  ARKLsPrecSetupFn     inner_precsetupfn;
  ARKLsPrecSolveFn     inner_precsolvefn;

  /* Counters */
  long int nfs;  /* num fs calls */

  /* Reusable arrays for fused vector operations */
  realtype* cvals;
  N_Vector* Xvecs;

} *ARKodeMRIStepMem;


/*===============================================================
  MRI time step module private function prototypes
  ===============================================================*/

/* Create MRIStep memory structure */
void* mriStep_Create(ARKRhsFn fs, realtype t0, N_Vector y0);

/* Interface routines supplied to ARKode */
int mriStep_Init(void* arkode_mem, int init_type);
int mriStep_FullRHS(void* arkode_mem, realtype t,
                    N_Vector y, N_Vector f, int mode);
int mriStep_TakeStep(void* arkode_mem);

/* Internal utility routines */
int mriStep_AccessStepMem(void* arkode_mem, const char *fname,
                          ARKodeMem *ark_mem, ARKodeMRIStepMem *step_mem);
booleantype mriStep_CheckNVector(N_Vector tmpl);
int mriStep_SetButcherTable(ARKodeMem ark_mem);
int mriStep_CheckButcherTable(ARKodeMem ark_mem);

int mriStep_ComputeErrorEst(ARKodeMem ark_mem, realtype *dsm);
int mriStep_DoErrorTest(ARKodeMem ark_mem, int *nefPtr,
                        realtype dsm);
int mriStep_PrepareNextStep(ARKodeMem ark_mem, realtype dsm);

/* Wrappers for user-supplied functions to the inner stepper */
int mriStep_InnerRhsFn(realtype t, N_Vector y, N_Vector ydot, void *user_data);
int mriStep_InnerRhsFnForcing(realtype t, N_Vector y, N_Vector ydot, void *user_data);
int mriStep_InnerFullRhs(void *arkode_mem, realtype t, N_Vector y,
                         N_Vector ydot);
int mriStep_InnerEwtFn(N_Vector y, N_Vector ewt, void *user_data);
int mriStep_InnerPostProcessStepFn(realtype t, N_Vector y,
                                   void *user_data);
int mriStep_InnerJacFn(realtype t, N_Vector y, N_Vector fy,
                       SUNMatrix Jac, void *user_data,
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int mriStep_InnerJacTimesSetupFn(realtype t, N_Vector y,
                                 N_Vector fy, void *user_data);
int mriStep_InnerJacTimesVecFn(N_Vector v, N_Vector Jv,
                               realtype t, N_Vector y,
                               N_Vector fy, void *user_data,
                               N_Vector tmp);
int mriStep_InnerPrecSetupFn(realtype t, N_Vector y,
                             N_Vector fy, booleantype jok,
                             booleantype *jcurPtr,
                             realtype gamma, void *user_data);
int mriStep_InnerPrecSolveFn(realtype t, N_Vector y,
                             N_Vector fy, N_Vector r,
                             N_Vector z, realtype gamma,
                             realtype delta, int lr,
                             void *user_data);

/* Attach ARKStep inner stepper */
int mriStep_AttachARK(void* arkode_mem, void* inner_arkode_mem);

/* Evolve ARKStep inner stepper */
int mriStep_EvolveInnerARK(void* inner_arkode_mem, realtype tout,
                           N_Vector yout, realtype tspan);

/*===============================================================
  Reusable MRIStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_MRISTEP_NO_MEM    "Time step module memory is NULL."

#ifdef __cplusplus
}
#endif

#endif
