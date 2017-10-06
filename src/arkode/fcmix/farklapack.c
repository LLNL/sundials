/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * Fortran/C interface routines for ARKODE/ARKLAPACK
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_lapack.h>

/*=============================================================*/

/* Fortran interface to C routine ARKLapackDense; see farkode.h 
   for further information */
void FARK_LAPACKDENSE(int *neq, int *ier)
{
  *ier = ARKLapackDense(ARK_arkodemem, *neq);
  ARK_ls = ARK_LS_LAPACKDENSE;
  return;
}

/* Fortran interface to C routine ARKMassLapackDense; see farkode.h 
   for further information */
void FARK_MASSLAPACKDENSE(int *neq, int *ier)
{
  *ier = ARKMassLapackDense(ARK_arkodemem, *neq, NULL);
  ARK_mass_ls = ARK_LS_LAPACKDENSE;
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKLapackBand; see farkode.h 
   for further information */
void FARK_LAPACKBAND(int *neq, int *mupper, int *mlower, int *ier)
{
  *ier = ARKLapackBand(ARK_arkodemem, *neq, *mupper, *mlower);
  ARK_ls = ARK_LS_LAPACKBAND;
  return;
}

/* Fortran interface to C routine ARKMassLapackBand; see farkode.h 
   for further information */
void FARK_MASSLAPACKBAND(int *neq, int *mupper, int *mlower, int *ier)
{
  *ier = ARKMassLapackBand(ARK_arkodemem, *neq, *mupper, 
			   *mlower, NULL);
  ARK_mass_ls = ARK_LS_LAPACKBAND;
  return;
}

/*===============================================================
   EOF
===============================================================*/

