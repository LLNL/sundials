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
void FARK_KLU(int *neq, int *nnz, int *ier)
{
  *ier = ARKKLU(ARK_arkodemem, *neq, *nnz);
  ARK_ls = ARK_LS_KLU;
  return;
}

/*===============================================================
   EOF
===============================================================*/

