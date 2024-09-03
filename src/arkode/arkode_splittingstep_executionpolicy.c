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
 * This is the implementation file for execution policy.
 *--------------------------------------------------------------*/

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_splittingstep.h>
#include <sundials/sundials_nvector.h>

#include "sundials_macros.h"

/*---------------------------------------------------------------
  This routine frees all the ARKodeSplittingExecutionPolicy
  memory
  ---------------------------------------------------------------*/
void ARKodeSplittingExecutionPolicy_Free(ARKodeSplittingExecutionPolicy* policy)
{
  if (policy != NULL)
  {
    (*policy)->free(*policy);
    free(*policy);
    *policy = NULL;
  }
}
