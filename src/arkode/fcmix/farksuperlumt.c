/*---------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
 Copyright (c) 2014, Southern Methodist University.
 All rights reserved.
 For details, see the LICENSE file.
 ----------------------------------------------------------------
 Fortran/C interface routines for ARKODE/ARKSUPERLUMT
 --------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_superlumt.h>

/*=============================================================*/

/* Fortran interface to C routine ARKSuperLUMT; see farkode.h for 
   further details */
void FARK_SUPERLUMT(int *nthreads, int *neq, int *nnz, int *ordering, int *ier)
{
  *ier = ARKSuperLUMT(ARK_arkodemem, *nthreads, *neq, *nnz);
  ARKSuperLUMTSetOrdering(ARK_arkodemem, *ordering);
  ARK_ls = ARK_LS_SUPERLUMT;
  return;
}

/*===============================================================
   EOF
===============================================================*/

