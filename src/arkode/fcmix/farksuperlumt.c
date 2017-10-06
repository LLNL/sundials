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
 * Fortran/C interface routines for ARKODE/ARKSUPERLUMT
 *--------------------------------------------------------------*/

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


/* Fortran interface to C routine ARKMassSuperLUMT; see farkode.h for 
   further details */
void FARK_MASSSUPERLUMT(int *nthreads, int *neq, int *nnz, int *ordering, int *ier)
{
  *ier = ARKMassSuperLUMT(ARK_arkodemem, *nthreads, *neq, *nnz, NULL);
  ARKMassSuperLUMTSetOrdering(ARK_arkodemem, *ordering);
  ARK_mass_ls = ARK_LS_SUPERLUMT;
  return;
}

/*===============================================================
   EOF
===============================================================*/

