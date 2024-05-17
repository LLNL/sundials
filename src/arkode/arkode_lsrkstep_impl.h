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
 * Implementation header file for ARKODE's LSRK time stepper
 * module.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_LSRKSTEP_IMPL_H
#define _ARKODE_LSRKSTEP_IMPL_H

#include <arkode/arkode_lsrkstep.h>

#include "arkode_impl.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  LSRK time step module data structure
  ===============================================================*/

/*---------------------------------------------------------------
  Types : struct ARKodeLSRKStepMemRec, ARKodeLSRKStepMem
  ---------------------------------------------------------------
  The type ARKodeLSRKStepMem is type pointer to struct
  ARKodeLSRKStepMemRec.  This structure contains fields to
  perform an explicit Runge-Kutta time step.
  ---------------------------------------------------------------*/
typedef struct ARKodeLSRKStepMemRec
{
  /* LSRK problem specification */
  ARKRhsFn fe;

  /* LSRK method storage and parameters */
  N_Vector* F;          /* RHS vector storage */
  int q;                /* method order       */
  int p;                /* embedding order    */
  int stages;           /* number of stages   */

  /* Counters */
  long int nfe;         /* num fe calls       */

  /* Reusable fused vector operation arrays */
  sunrealtype* cvals;
  N_Vector* Xvecs;

}* ARKodeLSRKStepMem;

/*===============================================================
  LSRK time step module private function prototypes
  ===============================================================*/

/* Interface routines supplied to ARKODE */
int lsrkStep_Init(ARKodeMem ark_mem, int init_type);
int lsrkStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                    int mode);
int lsrkStep_TakeStep(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int lsrkStep_SetDefaults(ARKodeMem ark_mem);
int lsrkStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile, SUNOutputFormat fmt);
int lsrkStep_WriteParameters(ARKodeMem ark_mem, FILE* fp);
int lsrkStep_Reset(ARKodeMem ark_mem, sunrealtype tR, N_Vector yR);
int lsrkStep_Resize(ARKodeMem ark_mem, N_Vector y0, sunrealtype hscale,
                   sunrealtype t0, ARKVecResizeFn resize, void* resize_data);
void lsrkStep_Free(ARKodeMem ark_mem);
void lsrkStep_PrintMem(ARKodeMem ark_mem, FILE* outfile);

/* Internal utility routines */
int lsrkStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                ARKodeMem* ark_mem, ARKodeLSRKStepMem* step_mem);
int lsrkStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                          ARKodeLSRKStepMem* step_mem);
sunbooleantype lsrkStep_CheckNVector(N_Vector tmpl);

/*===============================================================
  Reusable LSRKStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_LSRKSTEP_NO_MEM "Time step module memory is NULL."

#ifdef __cplusplus
}
#endif

#endif
