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
 * Fortran/C interface routines for ARKODE/ARKKLU
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_klu.h>

/*=============================================================*/

/* Fortran interface to C routine ARKKLU; see farkode.h for 
   further details */
void FARK_KLU(int *neq, int *nnz, int *sparsetype, int *ordering, int *ier)
{
  *ier = ARKKLU(ARK_arkodemem, *neq, *nnz, *sparsetype);
  ARKKLUSetOrdering(ARK_arkodemem, *ordering);
  ARK_ls = ARK_LS_KLU;
  return;
}

/* Fortran interface to C routine ARKKLUReinit; see farkode.h for 
   further details */
void FARK_KLUREINIT(int *neq, int *nnz, int *reinit_type, int *ier)
{
  *ier = ARKKLUReInit(ARK_arkodemem, *neq, *nnz, *reinit_type);
}

/* Fortran interface to C routine ARKMassKLU; see farkode.h for 
   further details */
void FARK_MASSKLU(int *neq, int *nnz, int *sparsetype, int *ordering, int *ier)
{
  *ier = ARKMassKLU(ARK_arkodemem, *neq, *nnz, *sparsetype, NULL);
  ARKMassKLUSetOrdering(ARK_arkodemem, *ordering);
  ARK_mass_ls = ARK_LS_KLU;
  return;
}

/* Fortran interface to C routine ARKMassKLUReinit; see farkode.h for 
   further details */
void FARK_MASSKLUREINIT(int *neq, int *nnz, int *reinit_type, int *ier)
{
  *ier = ARKMassKLUReInit(ARK_arkodemem, *neq, *nnz, *reinit_type);
}

/*===============================================================
   EOF
===============================================================*/

