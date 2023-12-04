/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_MPI_ERRORS_H
#define _SUNDIALS_MPI_ERRORS_H

#include <sundials/sundials_errors.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_EXPORT
void SUNMPIAbortErrHandlerFn(int line, const char* func, const char* file,
                             const char* msg, SUNErrCode err_code,
                             void* err_user_data, SUNContext sunctx);

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_MPI_ERRORS_H */
