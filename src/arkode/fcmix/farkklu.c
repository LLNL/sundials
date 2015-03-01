/*---------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
 Copyright (c) 2014, Southern Methodist University.
 All rights reserved.
 For details, see the LICENSE file.
 ----------------------------------------------------------------
 Fortran/C interface routines for ARKODE/ARKKLU
 --------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_klu.h>

/*=============================================================*/

/* Fortran interface to C routine ARKKLU; see farkode.h for 
   further details */
void FARK_KLU(int *neq, int *nnz, int *ordering, int *ier)
{
  *ier = ARKKLU(ARK_arkodemem, *neq, *nnz);
  ARKKLUSetOrdering(ARK_arkodemem, *ordering);
  ARK_ls = ARK_LS_KLU;
  return;
}

/* Fortran interface to C routine ARKKLUReinit; see farkode.h for 
   further details */
void FARK_KLUReinit(int *neq, int *nnz, *reinit_type, int *ier)
{
  *ier = ARKKLUReinit(ARK_arkodemem, *neq, *nnz, *reinit_type);
}

/*===============================================================
   EOF
===============================================================*/

