/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
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

