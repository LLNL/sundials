/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-22 00:12:50 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Fortran/C interface routines for IDA/IDALAPACK.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fida.h"     /* actual function names, prototypes and global vars.*/
#include "ida_impl.h" /* definition of IDAMem type                         */

#include <ida/ida_lapack.h>

/*************************************************/

void FIDA_LAPACKDENSE(int *neq, int *ier)
{

  *ier = 0;

  *ier = IDALapackDense(IDA_idamem, *neq);

  IDA_ls = IDA_LS_LAPACKDENSE;

  return;
}

/*************************************************/

void FIDA_LAPACKBAND(int *neq, int *mupper, int *mlower, int *ier)
{

  *ier = 0;

  *ier = IDALapackBand(IDA_idamem, *neq, *mupper, *mlower);

  IDA_ls = IDA_LS_LAPACKBAND;

  return;
}

/*************************************************/
