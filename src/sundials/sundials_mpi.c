/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is implementation of SUNDIALS MPI wrapper functions.
 * -----------------------------------------------------------------*/

#include <sundials/sundials_mpi.h>

int SUNMPI_Comm_size(SUNMPI_Comm comm, int *size)
{
#if SUNDIALS_MPI_ENABLED
  if (comm != SUNMPI_COMM_NULL) {
    return(MPI_Comm_size(comm, size));
  } else {
    *size = 1;
    return(SUNMPI_SUCCESS);
  }
#else
  *size = 1;
  return(SUNMPI_SUCCESS);
#endif
}

int SUNMPI_Allreduce_scalar(realtype d, realtype *result, SUNMPI_Op op, SUNMPI_Comm comm)
{
#if SUNDIALS_MPI_ENABLED
  if (comm != SUNMPI_COMM_NULL) {
    return(MPI_Allreduce(&d, result, 1, PVEC_REAL_MPI_TYPE, op, comm));
  } else {
    *result = d;
    return(SUNMPI_SUCCESS);
  }
#else
  /* If MPI is not enabled, just copy input into output */
  *result = d;
  return(SUNMPI_SUCCESS);
#endif /* ifdef SUNDIALS_MPI_ENABLED */
}


int SUNMPI_Allreduce(realtype *d, int nvec, SUNMPI_Op op, SUNMPI_Comm comm)
{
#if SUNDIALS_MPI_ENABLED
  if (comm != SUNMPI_COMM_NULL) {
    return(MPI_Allreduce(MPI_IN_PLACE, d, nvec, PVEC_REAL_MPI_TYPE, op, comm));
  } else {
    return(SUNMPI_SUCCESS);
  }
#else
  /* If MPI is not enabled don't do reduction */
  return(SUNMPI_SUCCESS);
#endif /* ifdef SUNDIALS_MPI_ENABLED */
}


int SUNMPI_Comm_dup(SUNMPI_Comm comm, SUNMPI_Comm *newcomm)
{
#if SUNDIALS_MPI_ENABLED
  return(MPI_Comm_dup(comm, newcomm));
#else
  *newcomm = comm;
  return(SUNMPI_SUCCESS);
#endif
}


int SUNMPI_Comm_compare(SUNMPI_Comm comm1, SUNMPI_Comm comm2, int *result)
{
#if SUNDIALS_MPI_ENABLED
  return(MPI_Comm_compare(comm1, comm2, result));
#else
  *result = SUNMPI_IDENT;
  return(SUNMPI_SUCCESS);
#endif
}


int SUNMPI_Comm_rank(SUNMPI_Comm comm, int *rank)
{
#if SUNDIALS_MPI_ENABLED
  if (comm != SUNMPI_COMM_NULL) {
    return(MPI_Comm_rank(comm, rank));
  } else {
    *rank = 0;
    return(SUNMPI_SUCCESS);
  }
#else
  *rank = 0;
  return(SUNMPI_SUCCESS);
#endif
}
