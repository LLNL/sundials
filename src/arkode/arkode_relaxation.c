/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation file for ARKODE's relaxation (in time)
 * functionality
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "arkode_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>


int arkSetRelaxFn(void* arkode_mem, ARKRelaxFn rfn, ARKRelaxJacFn rjac)
{
  ARKodeMem ark_mem;

  /* Check inputs */
  if (!arkode_mem)
  {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetRelaxFn", MSG_ARK_NO_MEM);
    return ARK_MEM_NULL;
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (!rfn)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkSetRelaxFn",
                    "The relaxation function is NULL.");
    return ARK_ILL_INPUT;
  }

  if (!rjac)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "arkSetRelaxFn",
                    "The relaxation Jacobian function is NULL.");
    return ARK_ILL_INPUT;
  }

  /* Allocate the relaxation memory structure */
  if (!(ark_mem->relax_mem))
  {
    ark_mem->relax_mem = (ARKodeRelaxMem) malloc(sizeof(*(ark_mem->relax_mem)));
    if (!(ark_mem->relax_mem)) return(ARK_MEM_FAIL);

    /* Zero out proj_mem */
    memset(ark_mem->relax_mem, 0, sizeof(struct ARKodeRelaxMemRec));
  }

  ark_mem->relax_mem->rfn  = rfn;
  ark_mem->relax_mem->rjac = rjac;

  return ARK_SUCCESS;
}

/* this could almost be implemented as a user supplied post processing
   if the user had access to the method stages, NOTE this requires the
   stage solutions not the stage derivatives (as posed right now) so
   currently this requires the Z form of take step */
/* Residual and Jacobian are stepper dependent? to account for imex vs erk vs dirk */

static int arkRelaxResidual(realtype gam, realtype* res, ARKodeMem ark_mem)
{
  int flag;

  /* yn + gam * dir = gam ycur + (1 - gam) yn */
  N_VLinearSum(gam, ark_mem->ycur, (ONE - gam), ark_mem->yn,
               ark_mem->tempv1);

  flag = ark_mem->relax_mem->rfn(ark_mem->tempv1, res, ark_mem->user_data);
  if (flag) return flag;

  /* res = r(yn + gam * dir) - r(yn) - gam * est */
  *res = *res - ark_mem->relax_mem->rcur - gam * ark_mem->relax_mem->est;

  return ARK_SUCCESS;
}


int arkRelax(void* arkode_mem, realtype* gam)
{
  int i, flag;
  realtype res;

  ARKodeMem ark_mem;

  /* Check inputs */
  if (!arkode_mem)
  {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE",
                    "arkSetRelaxFn", MSG_ARK_NO_MEM);
    return ARK_MEM_NULL;
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *gam = ONE;
  for (i = 0; i < 100; i++)
  {
    flag = arkRelaxResidual(*gam, &res, ark_mem);
    if (flag) return flag;

    *gam += res;

    if (res < 1.0e-12) break;
  }

  printf("gamma = %g\n", *gam);

  return 0;
}

/*===============================================================
  EOF
  ===============================================================*/
