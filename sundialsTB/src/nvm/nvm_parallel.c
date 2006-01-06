/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-01-06 19:00:26 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Vector constructors for the SUNDIALS Matlab interfaces.
 * -----------------------------------------------------------------
 */

#include <stdlib.h>
#include "mpi.h"
#include "nvm.h"
#include "nvector_serial.h"
#include "nvector_parallel.h"

int VecType;
MPI_Comm comm;

void InitVectors(int vec_type, mxArray *mx_comm)
{
  char *str;

  VecType = vec_type;

  if (VecType == 2) {

    str = mxArrayToString(mx_comm);
    if (!strcmp(str,"NULL"  ))      comm = MPI_COMM_NULL  ;
    else if (!strcmp(str,"WORLD" )) comm = MPI_COMM_WORLD ;
    else if (!strcmp(str,"SELF"  )) comm = MPI_COMM_SELF  ;
    else                            comm = *(MPI_Comm*)mxGetData(mx_comm);

  }

}

N_Vector NewVector(int n)
{
  N_Vector v;
  int nlocal, nglobal;

  if (VecType == 1) {

    v = N_VNewEmpty_Serial((long int)n);

  } else {

    nlocal = n;
    MPI_Allreduce(&nlocal, &nglobal, 1, MPI_INT, MPI_SUM, comm);
    v = N_VNewEmpty_Parallel(comm, (long int)nlocal, (long int)nglobal);

  }

  return(v);
}


