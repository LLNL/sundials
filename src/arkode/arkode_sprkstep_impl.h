/*---------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Implementation header file for ARKODE's ARK time stepper
 * module.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_SPRKSTEP_IMPL_H
#define _ARKODE_SPRKSTEP_IMPL_H

#include <arkode/arkode.h>
#include <arkode/arkode_sprk.h>
#include <arkode/arkode_sprkstep.h>
#include "arkode_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  SPRK time step module constants
  ===============================================================*/


/*===============================================================
  SPRK time step module data structure
  ===============================================================*/

/*---------------------------------------------------------------
  Types : struct ARKodeSPRKStepMemRec, ARKodeSPRKStepMem
  ---------------------------------------------------------------
  The type ARKodeSPRKStepMem is type pointer to struct
  ARKodeSPRKStepMemRec.  This structure contains fields to
  perform an symplectic-partitioned Runge-Kutta time step.
  ---------------------------------------------------------------*/
typedef struct ARKodeSPRKStepMemRec {
  
  /* SPRK method and storage */
  ARKodeSPRKMem method;
  int q;                       /* method order */
  N_Vector sdata;         

  /* SPRK problem specification */
  ARKRhsFn f1;                 /* p' = f1(t,q) = - dV(t,q)/dq  */
  ARKRhsFn f2;                 /* q' = f2(p)   =   dT(p)/dp    */

  /* Counters */
  long int nf1;                /* number of calls to f1        */
  long int nf2;                /* number of calls to f2        */
  int istage;

} *ARKodeSPRKStepMem;


/*===============================================================
  SPRK time step module private function prototypes
  ===============================================================*/

int sprkStep_Init(void* arkode_mem, int init_type);
int sprkStep_FullRHS(void* arkode_mem, realtype t,
                     N_Vector y, N_Vector f, int mode);
int sprkStep_TakeStep(void* arkode_mem, realtype *dsmPtr, int *nflagPtr);
int sprkStep_TakeStep_Compensated(void* arkode_mem, realtype *dsmPtr, int *nflagPtr);

/* Internal utility routines */
int sprkStep_AccessStepMem(void* arkode_mem, const char *fname,
                          ARKodeMem *ark_mem, ARKodeSPRKStepMem *step_mem);
booleantype sprkStep_CheckNVector(N_Vector tmpl);
int sprkStep_f1(ARKodeSPRKStepMem step_mem, sunrealtype tcur, N_Vector ycur, N_Vector f1, void* user_data);
int sprkStep_f2(ARKodeSPRKStepMem step_mem, sunrealtype tcur, N_Vector ycur, N_Vector f2, void* user_data);
int sprkStep_SPRKStage(ARKodeMem ark_mem, ARKodeSPRKStepMem step_mem, N_Vector prev_stage,
                       sunrealtype bi, sunrealtype Bi, N_Vector stage_result);

// /* private functions for interfacing with MRIStep */
// int sprkStep_SetInnerForcing(void* arkode_mem, realtype tshift, realtype tscale,
//                             N_Vector *f, int nvecs);
// int sprkStep_MRIStepInnerEvolve(MRIStepInnerStepper stepper,
//                                realtype t0, realtype tout, N_Vector y);
// int sprkStep_MRIStepInnerFullRhs(MRIStepInnerStepper stepper, realtype t,
//                                 N_Vector y, N_Vector f, int mode);
// int sprkStep_MRIStepInnerReset(MRIStepInnerStepper stepper, realtype tR,
//                               N_Vector yR);


/*===============================================================
  Reusable SPRKStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_SPRKSTEP_NO_MEM    "Time step module memory is NULL."
#define MSG_NLS_INIT_FAIL      "The nonlinear solver's init routine failed."

/* Other error messages */
#define MSG_ARK_MISSING_FE     "Cannot specify that method is explicit without providing a function pointer to fe(t,y)."
#define MSG_ARK_MISSING_FI     "Cannot specify that method is implicit without providing a function pointer to fi(t,y)."
#define MSG_ARK_MISSING_F      "Cannot specify that method is ImEx without providing function pointers to fi(t,y) and fe(t,y)."

#ifdef __cplusplus
}
#endif

#endif
