/*******************************************************************
 *                                                                 *
 * File          : fnvector_parallel.c                             *
 * Programmers   : Radu Serban @ LLNL                              *
 * Version of    : 26 June 2002                                    *
 *                                                                 *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This file, companion of nvector_parallel.c contains the         *
 * implementation of the Fortran interface to M_EnvInit_Parallel   *
 * and M_EnvFree_Parallel.                                         *
 *                                                                 *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h"
#include "nvector_parallel.h"
#include "fnvector_parallel.h"
#include "mpi.h"

/* Define global variable F2C_machEnv */
M_Env F2C_machEnv;

/* Fortran callable interfaces to M_EnvInit_Parallel
   and M_EnvFree_Parallel */

void F_MENVINITP(integertype *nlocal, integertype *nglobal, int *ier)
{
  
  /* Call M_EnvInit_Parallel:
     the first slot is for the communicator. 
     (From Fortran, only MPI_COMM_WORLD is allowed)
     *nlocal  is the local vector length
     *nglobal is the global vector length */

 int dumargc; char **dumargv;

 F2C_machEnv = M_EnvInit_Parallel(MPI_COMM_WORLD, *nlocal, *nglobal,
                                  &dumargc, &dumargv);

 *ier = (F2C_machEnv == NULL) ? -1 : 0 ;
}


void F_MENVFREEP()
{
  M_EnvFree_Parallel(F2C_machEnv);
}

