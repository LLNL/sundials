/*******************************************************************
 *                                                                 *
 * File          : fnvector_parallel.c                             *
 * Programmers   : Radu Serban @ LLNL                              *
 * Version of    : 06 June 2003                                    *
 *                                                                 *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This file, companion of nvector_parallel.c contains the         *
 * implementation of the Fortran interface to NV_SpecInit_Parallel *
 * and NV_SpecFree_Parallel.                                       *
 *                                                                 *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h"
#include "nvector_parallel.h"
#include "fnvector_parallel.h"
#include "mpi.h"

/* Define global variable F2C_machEnv */
NV_Spec F2C_nvspec;

/* Fortran callable interfaces to NV_SpecInit_Parallel
   and NV_SpecFree_Parallel */

void FNV_INITP(integertype *nlocal, integertype *nglobal, int *ier)
{
  
  /* Call NV_SpecInit_Parallel:
     the first slot is for the communicator. 
     (From Fortran, only MPI_COMM_WORLD is allowed)
     *nlocal  is the local vector length
     *nglobal is the global vector length */

 int dumargc; char **dumargv;

 F2C_nvspec = NV_SpecInit_Parallel(MPI_COMM_WORLD, *nlocal, *nglobal,
                                   &dumargc, &dumargv);

 *ier = (F2C_nvspec == NULL) ? -1 : 0 ;
}


void FNV_FREEP()
{
  NV_SpecFree_Parallel(F2C_nvspec);
}

