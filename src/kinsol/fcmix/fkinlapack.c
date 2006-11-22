/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-22 00:12:51 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Fortran/C interface routines for KINSOL/KINLAPACK
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fkinsol.h"
#include "kinsol_impl.h"

#include <kinsol/kinsol_lapack.h>

/***************************************************************************/

void FKIN_LAPACKDENSE(int *neq, int *ier)
{
  *ier = KINLapackDense(KIN_kinmem, *neq);
  KIN_ls = KIN_LS_LAPACKDENSE;
}

/***************************************************************************/

void FCV_LAPACKBAND(int *neq, int *mupper, int *mlower, int *ier)
{
  *ier = KINLapackBand(KIN_kinmem, *neq, *mupper, *mlower);
  KIN_ls = KIN_LS_LAPACKBAND;
}

/***************************************************************************/

