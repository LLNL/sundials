/*---------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
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
 * TODO
 *--------------------------------------------------------------*/

#ifndef ARKODE_SPLITTINGSTEP_IMPL_H_
#define ARKODE_SPLITTINGSTEP_IMPL_H_

#include <arkode/arkode_splittingstep.h>

typedef struct ARKodeSplittingStepMemRec
{
  SUNStepper *steppers;

  ARKodeSplittingCoeffs coeffs;

  /* Counters */
  long int nfe;       /* num fe calls               */
  long int nfi;       /* num fi calls               */
  long int nsetups;   /* num setup calls            */
  long int nls_iters; /* num nonlinear solver iters */
  long int nls_fails; /* num nonlinear solver fails */

  /* Reusable arrays for fused vector operations */
  sunrealtype* cvals; /* scalar array for fused ops       */
  N_Vector* Xvecs;    /* array of vectors for fused ops   */
  int nfusedopvecs;   /* length of cvals and Xvecs arrays */

}* ARKodeSplittingStepMem;

#endif
