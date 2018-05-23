/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This header file contains definitions of MPI data types, which
 * are used by MPI parallel vector implementations.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_MPI_H
#define _SUNDIALS_MPI_H

#include <sundials/sundials_types.h>
#include <sundials/sundials_mpi_types.h>


#ifdef SUNDIALS_MPI_ENABLED

#include <mpi.h>
typedef MPI_Comm SUNDIALS_Comm;
#else
typedef int SUNDIALS_Comm;
  #warning "SUNDIALS_MPI_ENABLED not defined!"
#endif // ifdef SUNDIALS_MPI_ENABLED

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

SUNDIALS_EXPORT realtype SUNDIALS_Reduce(realtype d, int op, SUNDIALS_Comm comm);
SUNDIALS_EXPORT void SUNDIALS_Allreduce(realtype *d, int nvec, int op, SUNDIALS_Comm comm);
SUNDIALS_EXPORT void SUNDIALS_Comm_size(SUNDIALS_Comm comm, int *npes);

#ifdef __cplusplus
}
#endif



#endif /* _SUNDIALS_MPI_H */