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
#include <sundials/sundials_nvector.h>

static int setup_serial(ARKodeSplittingExecutionPolicy policy, N_Vector y,
                        const int sequential_methods)
{
  // Nothing needed
  return ARK_SUCCESS;
}

static int execute_serial(ARKodeSplittingExecutionPolicy policy,
                          ARKParallelExecuteFn fn, N_Vector yn, N_Vector ycur,
                          N_Vector tmp, sunrealtype* alpha,
                          const int sequential_methods, void* user_data)
{
  N_VScale(1, yn, ycur);
  fn(0, ycur, user_data);

  // For many methods alpha[0] == 1, so this is redundant. TODO: add check or
  // hope compiler optimized this out
  // N_VScale(alpha[0], ycur, ycur);

  // for (int i = 1; i < sequential_methods; i++)
  // {
  //   N_VScale(1, yn, tmp);
  //   int retval = fn(i, tmp, user_data);
  //   if (retval != ARK_SUCCESS)
  //   // TODO: error handling of fn
  //   N_VLinearSum(1, ycur, alpha[i], tmp, yn);
  // }

  return ARK_SUCCESS;
}

static void free_serial(ARKodeSplittingExecutionPolicy policy)
{
  // Nothing needed
}

ARKodeSplittingExecutionPolicy ARKodeSplittingExecutionPolicy_Serial()
{
  ARKodeSplittingExecutionPolicy policy = malloc(sizeof(*policy));
  if (policy == NULL) { return NULL; }
  policy->setup    = setup_serial;
  policy->execute = execute_serial;
  policy->free     = free_serial;
  policy->data     = NULL;
  return policy;
}

void ARKodeSplittingExecutionPolicyFree(ARKodeSplittingExecutionPolicy* policy)
{}

// #include <mpi.h>

// ARKodeSplittingExecutionPolicy SplittingStepExecutionPolicy_MPI(SUNComm comm) {
//   ARKodeSplittingExecutionPolicy policy = malloc(sizeof(*policy));
//   policy->setup = setup_mpi;
//   policy->exectute = execute_mpi;
//   policy->free = free_mpi;
//   SUNComm *data = malloc(sizeof(*data));
//   data[0] = comm;
//   policy->data = data;
// }

// static int mpi_reduce(void *in, void *inout, int *len, MPI_Datatype *datatype) {
//   // Is there any way
// }

// static int setup_mpi(ARKodeSplittingExecutionPolicy policy, N_Vector y, const int sequential_methods) {
//   MPI_Op_create(mpi_reduce, SUNTRUE, )
// }

// static int execute_mpi(ARKodeSplittingExecutionPolicy policy, ARKParallelExecuteFn fn, N_Vector yn, N_Vector ycur, N_Vector tmp, sunrealtype *alpha, const int sequential_methods, void *user_data)
// {
//   SUNComm comm = *((SUNComm*) policy->data);
//   int rank;
//   MPI_Comm_rank(comm, &rank);

//   N_VScale(1, yn, ycur);
//   fn(rank, y, user_data);
//   N_VScale(alpha[rank], ycur, ycur);

//   for (int i = 2 * rank; i < sequential_methods; rank++) {
//     N_VScale(1, yn, tmp);
//     fn(i, tmp, user_data);
//     NN_VLinearSum(1, ycur, alpha[i], tmp);
//   }

//   sunindextype bufferSize;
//   N_VBufSize(ycur, &bufferSize);
// }

// static void free_mpi(ARKodeSplittingExecutionPolicy policy) {
//   MPI_op_free();
//   free(policy->data);
// }

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

// int execute_mpi(ARKParallelExecuteFn fn, N_Vector *y, sunrealtype *alpha, const int sequential_methods, void *user_data)
// {
//   MPI_Bcast()

//   MPI_Barrier(comm);

//   int rank;
//   MPI_Comm_rank(comm, &rank);

//   for (int i = rank; i < sequential_methods; i+=rank) {
//     fn(i, y[i], user_data);
//   }

//   // This could be done as a parallel reduction. The benefit appears minimal as
//   // it would require at least 4 sequential methods to theoretically get a
//   // speedup.
//   N_VLinearCombination(sequential_methods, alpha, y, y[0]);
// }
