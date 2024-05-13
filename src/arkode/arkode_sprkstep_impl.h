/*---------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Implementation header file for ARKODE's SPRK time stepper
 * module.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_SPRKSTEP_IMPL_H
#define _ARKODE_SPRKSTEP_IMPL_H

#include <arkode/arkode.h>
#include <arkode/arkode_sprk.h>
#include <arkode/arkode_sprkstep.h>

#include "arkode_impl.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
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
typedef struct ARKodeSPRKStepMemRec
{
  /* SPRK method and storage */
  ARKodeSPRKTable method; /* method spec  */
  int q;                  /* method order */
  N_Vector sdata;         /* persisted stage data */
  N_Vector yerr;          /* error vector for compensated summation */

  /* SPRK problem specification */
  ARKRhsFn f1; /* p' = f1(t,q) = - dV(t,q)/dq  */
  ARKRhsFn f2; /* q' = f2(t,p) =   dT(t,p)/dp  */

  /* Counters */
  long int nf1; /* number of calls to f1        */
  long int nf2; /* number of calls to f2        */
  int istage;

}* ARKodeSPRKStepMem;

/*===============================================================
  SPRK time step module private function prototypes
  ===============================================================*/

/* Interface routines supplied to ARKODE */
int sprkStep_Init(ARKodeMem ark_mem, int init_type);
int sprkStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                     int mode);
int sprkStep_TakeStep(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int sprkStep_TakeStep_Compensated(ARKodeMem ark_mem, sunrealtype* dsmPtr,
                                  int* nflagPtr);
int sprkStep_SetUserData(ARKodeMem ark_mem, void* user_data);
int sprkStep_SetDefaults(ARKodeMem ark_mem);
int sprkStep_SetOrder(ARKodeMem ark_mem, int ord);
int sprkStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile, SUNOutputFormat fmt);
int sprkStep_WriteParameters(ARKodeMem ark_mem, FILE* fp);
int sprkStep_Reset(ARKodeMem ark_mem, sunrealtype tR, N_Vector yR);
int sprkStep_Resize(ARKodeMem ark_mem, N_Vector y0, sunrealtype hscale,
                    sunrealtype t0, ARKVecResizeFn resize, void* resize_data);
void sprkStep_Free(ARKodeMem ark_mem);
void sprkStep_PrintMem(ARKodeMem ark_mem, FILE* outfile);

/* Internal utility routines */
int sprkStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                 ARKodeMem* ark_mem, ARKodeSPRKStepMem* step_mem);
int sprkStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                           ARKodeSPRKStepMem* step_mem);
sunbooleantype sprkStep_CheckNVector(N_Vector tmpl);

/* f1 = p' (Force evaluation) */
int sprkStep_f1(ARKodeSPRKStepMem step_mem, sunrealtype tcur, N_Vector ycur,
                N_Vector f1, void* user_data);

/* f2 = q' (Velocity evaluation) */
int sprkStep_f2(ARKodeSPRKStepMem step_mem, sunrealtype tcur, N_Vector ycur,
                N_Vector f2, void* user_data);

/*===============================================================
  Reusable SPRKStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_SPRKSTEP_NO_MEM "Time step module memory is NULL."

#ifdef __cplusplus
}
#endif

#endif
