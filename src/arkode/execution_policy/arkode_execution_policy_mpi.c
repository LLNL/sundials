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