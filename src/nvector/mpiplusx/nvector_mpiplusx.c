/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * This is the implementation file for the MPIPlusX NVECTOR.
 * -----------------------------------------------------------------*/

#include <nvector/nvector_mpiplusx.h>

#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_errors.h"
#include "sundials_nvector_impl.h"

#define MPIPLUSX_LOCAL_VECTOR(v) (N_VGetSubvector_MPIManyVector(v, 0))

N_Vector N_VMake_MPIPlusX(MPI_Comm comm, N_Vector X, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  N_Vector v;

  SUNAssertNull(X, SUN_ERR_ARG_CORRUPT);

  v = NULL;
  v = N_VMake_MPIManyVector(comm, 1, &X, SUNCTX_);
  SUNCheckLastErrNull();

  /* override certain ops */
  v->ops->nvgetvectorid     = N_VGetVectorID_MPIPlusX;
  v->ops->nvgetarraypointer = N_VGetArrayPointer_MPIPlusX;
  v->ops->nvsetarraypointer = N_VSetArrayPointer_MPIPlusX;
  v->ops->nvgetlocallength  = N_VGetLocalLength_MPIPlusX;

  /* debugging functions */
  if (X->ops->nvprint) { v->ops->nvprint = N_VPrint_MPIPlusX; }

  if (X->ops->nvprintfile) { v->ops->nvprintfile = N_VPrintFile_MPIPlusX; }

  return v;
}

void N_VPrint_MPIPlusX(N_Vector v)
{
  N_Vector x = MPIPLUSX_LOCAL_VECTOR(v);
  if (x->ops->nvprint) { x->ops->nvprint(x); }
}

void N_VPrintFile_MPIPlusX(N_Vector v, FILE* outfile)
{
  N_Vector x = MPIPLUSX_LOCAL_VECTOR(v);
  if (x->ops->nvprintfile) { x->ops->nvprintfile(x, outfile); }
}
