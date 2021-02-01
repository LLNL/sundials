/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Vector constructors for the SUNDIALS Matlab interfaces.
 * -----------------------------------------------------------------
 */

#include <stdlib.h>
#include <mpi.h>
#include "nvm.h"
#include <nvector/nvector_serial.h>
#include <nvector/nvector_parallel.h>

int sundials_VecType;
MPI_Comm sundials_comm;

void InitVectors()
{
  const mxArray *mx_comm;
  char *str;

  /* Check if the Matlab global variable sundials_MPI_comm exists
     (mpirun and mpiruns set it) */

  mx_comm = mexGetVariable("global", "sundials_MPI_comm");

  if (mx_comm == NULL) {

    /* If it does not exist, set vector type to 1 (serial) */

    sundials_VecType = 1;

  } else {

    /* If it does exist, set vector type to 2 (parallel) and
       set the MPI communicator */
    
    sundials_VecType = 2;

    str = mxArrayToString(mx_comm);
    if (!strcmp(str,"NULL"  ))      sundials_comm = MPI_COMM_NULL  ;
    else if (!strcmp(str,"WORLD" )) sundials_comm = MPI_COMM_WORLD ;
    else if (!strcmp(str,"SELF"  )) sundials_comm = MPI_COMM_SELF  ;
    else                            sundials_comm = *(MPI_Comm*)mxGetData(mx_comm);
    
  }

}

N_Vector NewVector(sunindextype n)
{
  N_Vector v;
  sunindextype nlocal, nglobal;

  if (sundials_VecType == 1) {

    v = N_VNew_Serial((sunindextype)n);

  } else {

    nlocal = n;
    MPI_Allreduce(&nlocal, &nglobal, 1, MPI_INT, MPI_SUM, sundials_comm);
    v = N_VNew_Parallel(sundials_comm, nlocal, nglobal);

  }

  return(v);
}


