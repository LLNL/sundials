/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
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
