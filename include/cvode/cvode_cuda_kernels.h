/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos and David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for CUDA specific kernels in CVODE.
 * -----------------------------------------------------------------*/

#ifndef _CVODE_CUDA_KERNELS_H
#define _CVODE_CUDA_KERNELS_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

SUNDIALS_EXPORT int CVodeSetUseIntegratorFusedKernels_CUDA(void *cvode_mem,
                                                           booleantype onoff);

SUNDIALS_EXPORT int CVDiagSetUseFusedKernels_CUDA(void *cvode_mem,
                                                  booleantype onoff);

#ifdef __cplusplus
}
#endif

#endif
