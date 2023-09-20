/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * ARKUserStability structure and utility routines.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_timestepheuristics.h>
#include "arkode_impl.h"


/* -----------------------------------------------------
 * Function to create an ARKUserStabilityData structure
 */

void* ARKUserStability(void* arkode_mem, ARKExpStabFn EStab, void* estab_data)
{
  /* Return with failure if arkode_mem or EStab are NULL */
  if ((EStab == NULL) || (arkode_mem == NULL)) { return (NULL); }

  /* Create data structure */
  ARKUserStabilityData userstab_data = NULL;
  userstab_data = (ARKUserStabilityData)malloc(sizeof *userstab_data);
  if (userstab_data == NULL) { return(NULL); }

  /* Attach inputs into structure */
  userstab_data->EStab = EStab;
  userstab_data->estab_data = estab_data;
  userstab_data->arkode_mem = arkode_mem;

  return (void*) userstab_data;
}


/* -----------------------------------------------------------------
 * Function to free an ARKUserStabilityData structure
 */

void ARKUserStability_Free(void* sdata)
{
  if (sdata == NULL) { return; }
  ARKUserStabilityData userstab_data = (ARKUserStabilityData) sdata;
  free(userstab_data);
  userstab_data = NULL;
  sdata = NULL;
}

/* -----------------------------------------------------
 * Function to evaluate the user-provided stability fn.
 */

int ARKControlExpStab(realtype *hstab, void* sdata)
{
  int retval;
  if (sdata == NULL) { return(SUNTIMESTEPHEURISTICS_USER_FCN_FAIL); }
  ARKUserStabilityData X = (ARKUserStabilityData) sdata;
  ARKodeMem ark_mem = (ARKodeMem) (X->arkode_mem);
  retval = (X->EStab)(ark_mem->ycur, ark_mem->tn, hstab, X->estab_data);
  if (retval != ARK_SUCCESS) { return(SUNTIMESTEPHEURISTICS_USER_FCN_FAIL); }
  return(SUNTIMESTEPHEURISTICS_SUCCESS);
}
