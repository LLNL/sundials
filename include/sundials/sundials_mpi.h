/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Slaven Peles @ LLNL
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
 * This header file provides a minimal wrapper for MPI, enabling
 * SUNDIALS to be built either with or without MPI, using a single
 * set of SUNDIALS implementation routines.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_MPI_H
#define _SUNDIALS_MPI_H

#include <sundials/sundials_types.h>
#include <sundials/sundials_mpi_types.h>


#if SUNDIALS_MPI_ENABLED

#include <mpi.h>

/* SUNMPI types and identifiers that point to MPI equivalents */

typedef MPI_Comm SUNMPI_Comm;
typedef MPI_Op   SUNMPI_Op;

#define SUNMPI_COMM_WORLD MPI_COMM_WORLD
#define SUNMPI_COMM_NULL  MPI_COMM_NULL
#define SUNMPI_SUCCESS    MPI_SUCCESS
#define SUNMPI_IDENT      MPI_IDENT
#define SUNMPI_CONGRUENT  MPI_CONGRUENT
#define SUNMPI_SIMILAR    MPI_SIMILAR
#define SUNMPI_UNEQUAL    MPI_UNEQUAL
#define SUNMPI_MAX        MPI_MAX
#define SUNMPI_MIN        MPI_MIN
#define SUNMPI_SUM        MPI_SUM
#define SUNMPI_PROD       MPI_PROD
#define SUNMPI_LAND       MPI_LAND
#define SUNMPI_LOR        MPI_LOR
#define SUNMPI_BAND       MPI_BAND
#define SUNMPI_BOR        MPI_BOR

#else

/* SUNMPI types and identifiers in the absence of MPI */

typedef int SUNMPI_Comm;

#define SUNMPI_COMM_WORLD 0
#define SUNMPI_SUCCESS    0
#define SUNMPI_COMM_NULL -1
#define SUNMPI_IDENT      0
#define SUNMPI_CONGRUENT  1
#define SUNMPI_SIMILAR    2
#define SUNMPI_UNEQUAL   -1

typedef enum {
  SUNMPI_MAX,
  SUNMPI_MIN,
  SUNMPI_SUM,
  SUNMPI_PROD,
  SUNMPI_LAND,
  SUNMPI_LOR,
  SUNMPI_BAND,
  SUNMPI_BOR
} SUNMPI_Op;

#endif


#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  /* SUNMPI functions that wrap MPI equivalents */

  SUNDIALS_EXPORT int SUNMPI_Comm_size(SUNMPI_Comm comm, int *size);
  SUNDIALS_EXPORT int SUNMPI_Allreduce_scalar(realtype d, realtype *result, SUNMPI_Op op, SUNMPI_Comm comm);
  SUNDIALS_EXPORT int SUNMPI_Allreduce(realtype *d, int nvec, SUNMPI_Op op, SUNMPI_Comm comm);
  SUNDIALS_EXPORT int SUNMPI_Comm_dup(SUNMPI_Comm comm, SUNMPI_Comm *newcomm);
  SUNDIALS_EXPORT int SUNMPI_Comm_compare(SUNMPI_Comm comm1, SUNMPI_Comm comm2, int *result);
  SUNDIALS_EXPORT int SUNMPI_Comm_rank(SUNMPI_Comm comm, int *rank);

#ifdef __cplusplus
}
#endif



#endif /* _SUNDIALS_MPI_H */
