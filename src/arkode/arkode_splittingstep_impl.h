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
  ARKodeSplittingExecutionPolicy policy;
  int partitions;
  int order;
  sunbooleantype own_policy;
}* ARKodeSplittingStepMem;

#endif
