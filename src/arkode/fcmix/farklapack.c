/*---------------------------------------------------------------
  $Revision: 1.0 $
  $Date: $
 ---------------------------------------------------------------- 
  Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
  Fortran/C interface routines for ARKODE/ARKLAPACK
 --------------------------------------------------------------*/

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

/*=============================================================*/

/* Fortran interface to C routine ARKLapackBand; see farkode.h 
   for further information */
void FARK_LAPACKBAND(int *neq, int *mupper, int *mlower, int *ier)
{
  *ier = ARKLapackBand(ARK_arkodemem, *neq, *mupper, *mlower);
  ARK_ls = ARK_LS_LAPACKBAND;
  return;
}

/*===============================================================
   EOF
===============================================================*/

