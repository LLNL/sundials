/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
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
 * Implementation header file for ARKODE's ERK time stepper
 * module.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_ERKSTEP_IMPL_H
#define _ARKODE_ERKSTEP_IMPL_H

#include <arkode/arkode_erkstep.h>

#include "arkode_impl.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ERK time step module constants -- move many items here from
  arkode_impl.h
  ===============================================================*/

/*===============================================================
  ERK time step module data structure
  ===============================================================*/

/*---------------------------------------------------------------
  Types : struct ARKodeERKStepMemRec, ARKodeERKStepMem
  ---------------------------------------------------------------
  The type ARKodeERKStepMem is type pointer to struct
  ARKodeERKStepMemRec.  This structure contains fields to
  perform an explicit Runge-Kutta time step.
  ---------------------------------------------------------------*/
typedef struct ARKodeERKStepMemRec
{
  /* ERK problem specification */
  ARKRhsFn f; /* y' = f(t,y)                */

  /* ARK method storage and parameters */
  N_Vector* F;          /* explicit RHS at each stage */
  int q;                /* method order               */
  int p;                /* embedding order            */
  int stages;           /* number of stages           */
  ARKodeButcherTable B; /* ERK Butcher table          */

  /* Counters */
  long int nfe; /* num fe calls               */

  /* Reusable arrays for fused vector operations */
  sunrealtype* cvals;
  N_Vector* Xvecs;

}* ARKodeERKStepMem;

/*===============================================================
  ERK time step module private function prototypes
  ===============================================================*/

/* Interface routines supplied to ARKODE */
int erkStep_Init(ARKodeMem ark_mem, int init_type);
int erkStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                    int mode);
int erkStep_TakeStep(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int erkStep_SetUserData(ARKodeMem ark_mem, void* user_data);
int erkStep_SetDefaults(ARKodeMem ark_mem);
int erkStep_SetOrder(ARKodeMem ark_mem, int ord);
int erkStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile, SUNOutputFormat fmt);
int erkStep_WriteParameters(ARKodeMem ark_mem, FILE* fp);
int erkStep_Reset(ARKodeMem ark_mem, sunrealtype tR, N_Vector yR);
int erkStep_Resize(ARKodeMem ark_mem, N_Vector y0, sunrealtype hscale,
                   sunrealtype t0, ARKVecResizeFn resize, void* resize_data);
void erkStep_Free(ARKodeMem ark_mem);
void erkStep_PrintMem(ARKodeMem ark_mem, FILE* outfile);
int erkStep_GetEstLocalErrors(ARKodeMem ark_mem, N_Vector ele);

/* Internal utility routines */
int erkStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                ARKodeMem* ark_mem, ARKodeERKStepMem* step_mem);
int erkStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                          ARKodeERKStepMem* step_mem);
sunbooleantype erkStep_CheckNVector(N_Vector tmpl);
int erkStep_SetButcherTable(ARKodeMem ark_mem);
int erkStep_CheckButcherTable(ARKodeMem ark_mem);
int erkStep_ComputeSolutions(ARKodeMem ark_mem, sunrealtype* dsm);

/* private functions for relaxation */
int erkStep_SetRelaxFn(ARKodeMem ark_mem, ARKRelaxFn rfn, ARKRelaxJacFn rjac);
int erkStep_RelaxDeltaE(ARKodeMem ark_mem, ARKRelaxJacFn relax_jac_fn,
                        long int* relax_jac_fn_evals, sunrealtype* delta_e_out);
int erkStep_GetOrder(ARKodeMem ark_mem);

/*===============================================================
  Reusable ERKStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_ERKSTEP_NO_MEM "Time step module memory is NULL."

#ifdef __cplusplus
}
#endif

#endif
