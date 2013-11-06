/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2012-03-07 21:41:19 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2013, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
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

N_Vector NewVector(long int n)
{
  N_Vector v;
  long int nlocal, nglobal;

  if (sundials_VecType == 1) {

    v = N_VNew_Serial((long int)n);

  } else {

    nlocal = n;
    MPI_Allreduce(&nlocal, &nglobal, 1, MPI_INT, MPI_SUM, sundials_comm);
    v = N_VNew_Parallel(sundials_comm, nlocal, nglobal);

  }

  return(v);
}


