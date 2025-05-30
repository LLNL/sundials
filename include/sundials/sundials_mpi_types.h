/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This header file contains definitions of MPI data types, which
 * are used by MPI parallel vector implementations.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_MPI_TYPES_H
#define _SUNDIALS_MPI_TYPES_H

#include <mpi.h>
#include <sundials/sundials_types.h>

/* define MPI data types */

#if defined(SUNDIALS_SINGLE_PRECISION)
#define MPI_SUNREALTYPE    MPI_FLOAT
#define MPI_SUNCOMPLEXTYPE MPI_C_COMPLEX
#elif defined(SUNDIALS_DOUBLE_PRECISION)
#define MPI_SUNREALTYPE    MPI_DOUBLE
#define MPI_SUNCOMPLEXTYPE MPI_C_DOUBLE_COMPLEX
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define MPI_SUNREALTYPE    MPI_LONG_DOUBLE
#define MPI_SUNCOMPLEXTYPE MPI_C_LONG_DOUBLE_COMPLEX
#endif

#if defined(SUNDIALS_SCALAR_TYPE_REAL)
#define MPI_SUNSCALARTYPE MPI_SUNREALTYPE
#else
#define MPI_SUNSCALARTYPE MPI_SUNCOMPLEXTYPE
#endif

#if defined(SUNDIALS_INT64_T)
#define MPI_SUNINDEXTYPE MPI_INT64_T
#elif defined(SUNDIALS_INT32_T)
#define MPI_SUNINDEXTYPE MPI_INT32_T
#endif

#endif /* _SUNDIALS_MPI_TYPES_H */
