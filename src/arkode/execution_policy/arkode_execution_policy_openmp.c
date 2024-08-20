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

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_splittingstep.h>

#include "sundials_macros.h"

// int execute_openmp(ARKParallelExecuteFn fn, N_Vector *y, sunrealtype *alpha, const int sequential_methods, void *user_data)
// {
//   #pragma omp parallel default(none) private(i) shared(fn, y, alpha, sequential_methods, user_data)
//   {
//     #pragma omp for schedule(static) num_threads(?)
//     for (int i = 1; i < sequential_methods; i++) {
//       // We assume that the "true" version of y_n is stored in y[0], then
//       // distribute it to the other y[i]. It might be possible to remove this if
//       // we can ensure the y[i] are always identical even after a reset,
//       // reinit, etc.
//       N_VScale(1, y[0], y[i]);
//     }

//     #pragma omp for schedule(dynamic) num_threads(?)
//     for (int i = 0; i < sequential_methods; i++) {
//       fn(i, y[i], user_data);
//     }
//   }

//   // This could be done as a parallel reduction. The benefit appears minimal as
//   // it would require at least 4 sequential methods to theoretically get a
//   // speedup.
//   N_VLinearCombination(sequential_methods, alpha, y, y[0]);
// }
