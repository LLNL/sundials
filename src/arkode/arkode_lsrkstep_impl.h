/*---------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
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
  ARKDomEigFn dom_eig_fn;

  N_Vector Fe; /* RHS vector storage */

  int q; /* method order               */
  int p; /* embedding order            */

  int req_stages; /* number of requested stages   */

  ARKODE_LSRKMethodType LSRKmethod;

  /* Counters and stats*/
  long int nfe;               /* num fe calls       */
  long int dom_eig_num_evals; /* num of dom_eig computations   */
  int stage_max;              /* num of max stages used      */
  int stage_max_limit;        /* max allowed num of stages     */
  int dom_eig_nst; /* num of steps that successfully used dom_eig; indicates dom_eig update when 0;  */

  /* Spectral info */
  sunrealtype lambdaR;         /* Real part of the dominated eigenvalue*/
  sunrealtype lambdaI;         /* Imaginary part of the dominated eigenvalue*/
  sunrealtype spectral_radius; /* spectral radius*/
  sunrealtype spectral_radius_max; /* max spectral radius*/
  sunrealtype spectral_radius_min; /* min spectral radius*/
  sunrealtype dom_eig_safety; /* some safety factor for the user provided dom_eig*/
  int dom_eig_freq; /* indicates dom_eig update after dom_eig_freq successful steps*/

  /* Flags */
  sunbooleantype dom_eig_update; /* flag indicating new dom_eig is needed */
  sunbooleantype const_Jac;      /* flag indicating Jacobian is constant */
  sunbooleantype dom_eig_is_current; /* SUNTRUE if dom_eig has been evaluated at tn */
  sunbooleantype is_SSP;             /* flag indicating SSP method*/

  /* Reusable fused vector operation arrays */
  sunrealtype* cvals;
  N_Vector* Xvecs;
  int nfusedopvecs; /* length of cvals and Xvecs arrays */

}* ARKodeLSRKStepMem;

/*===============================================================
  LSRK time step module private function prototypes
  ===============================================================*/

/* Interface routines supplied to ARKODE */
int lsrkStep_Init(ARKodeMem ark_mem, int init_type);
int lsrkStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                     int mode);
int lsrkStep_TakeStepRKC(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int lsrkStep_TakeStepRKL(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int lsrkStep_TakeStepSSPs2(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int lsrkStep_TakeStepSSPs3(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int lsrkStep_TakeStepSSP43(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int lsrkStep_TakeStepSSP104(ARKodeMem ark_mem, sunrealtype* dsmPtr,
                            int* nflagPtr);
int lsrkStep_SetDefaults(ARKodeMem ark_mem);
int lsrkStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile, SUNOutputFormat fmt);
int lsrkStep_WriteParameters(ARKodeMem ark_mem, FILE* fp);
int lsrkStep_Reset(ARKodeMem ark_mem, sunrealtype tR, N_Vector yR);
int lsrkStep_Resize(ARKodeMem ark_mem, N_Vector y0, sunrealtype hscale,
                    sunrealtype t0, ARKVecResizeFn resize, void* resize_data);
void lsrkStep_Free(ARKodeMem ark_mem);
void lsrkStep_PrintMem(ARKodeMem ark_mem, FILE* outfile);
int lsrkStep_GetEstLocalErrors(ARKodeMem ark_mem, N_Vector ele);

/* Internal utility routines */
int lsrkStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                 ARKodeMem* ark_mem, ARKodeLSRKStepMem* step_mem);
int lsrkStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                           ARKodeLSRKStepMem* step_mem);
void lsrkStep_DomEigUpdateLogic(ARKodeMem ark_mem, ARKodeLSRKStepMem step_mem,
                                sunrealtype dsm);
sunbooleantype lsrkStep_CheckNVector(N_Vector tmpl);
int lsrkStep_ComputeNewDomEig(ARKodeMem ark_mem, ARKodeLSRKStepMem step_mem);

/*===============================================================
  Reusable LSRKStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_LSRKSTEP_NO_MEM "Time step module memory is NULL."

#ifdef __cplusplus
}
#endif

#endif
