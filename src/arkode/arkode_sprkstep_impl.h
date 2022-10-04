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

#include <arkode/arkode_sprkstep.h>
#include "arkode_impl.h"

/* access to MRIStepInnerStepper_Create */
#include "arkode/arkode_mristep.h"
#include "sundials/sundials_types.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

struct TakeStepFunctor {
  int (*call)(TakeStep self, void* arkode_mem, sunrealtype *dsmPtr, int *nflagPtr);
  int q;
  int p;
  int stages;
  int num_vecs;
  void* content;
};

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
  TakeStep take_step;
  int q;                      /* method order                 */
  int p;                      /* embedding order              */
  int stages;                 /* number of stages             */
  N_Vector* scratch;          /* vectors that may be used     */

  /* SPRK problem specification */
  ARKRhsFn *fk;               /* array of RHS functions       */
  int num_rhs;                /* number of RHS functions      */

  /* Counters */
  long int* nfk;              /* number of calls to fk        */

  // /* Data for using SPRKStep with external polynomial forcing */
  // booleantype expforcing;  /* add forcing to explicit RHS */
  // booleantype impforcing;  /* add forcing to implicit RHS */
  // realtype    tshift;      /* time normalization shift    */
  // realtype    tscale;      /* time normalization scaling  */
  // N_Vector*   forcing;     /* array of forcing vectors    */
  // int         nforcing;    /* number of forcing vectors   */

} *ARKodeSPRKStepMem;


/*===============================================================
  SPRK time step module private function prototypes
  ===============================================================*/

int sprkStep_Init(void* arkode_mem, int init_type);
int sprkStep_FullRHS(void* arkode_mem, realtype t,
                     N_Vector y, N_Vector f, int mode);
int sprkStep_TakeStep_Sprk(void* arkode_mem, realtype *dsmPtr, int *nflagPtr);
int sprkStep_TakeStep_SprkInc(void* arkode_mem, realtype *dsmPtr, int *nflagPtr);

int sprkStep_Fk(ARKodeSPRKStepMem step_mem, sunrealtype tcur, N_Vector ycur, N_Vector Fk, int k, void* user_data);

/* Internal utility routines */
int sprkStep_AccessStepMem(void* arkode_mem, const char *fname,
                          ARKodeMem *ark_mem, ARKodeSPRKStepMem *step_mem);
booleantype sprkStep_CheckNVector(N_Vector tmpl);

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
