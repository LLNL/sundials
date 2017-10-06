/*
 * -----------------------------------------------------------------
 * $Revision: 4654 $
 * $Date: 2016-02-17 20:12:58 -0800 (Wed, 17 Feb 2016) $
 * -----------------------------------------------------------------
 * Programmer(s): Carol Woodward @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2015, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the Fortran interface to
 * the KINKLU solver. See fkinsol.h for usage.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "fkinsol.h"
#include "kinsol_impl.h"
#include <kinsol/kinsol_klu.h>
 
/*
 * ----------------------------------------------------------------
 * Function : FKIN_KLU
 * ----------------------------------------------------------------
 */

void FKIN_KLU(int *neq, int *nnz, int *sparsetype, int *ordering, int *ier)
{
  *ier = KINKLU(KIN_kinmem, *neq, *nnz, *sparsetype);
  KINKLUSetOrdering(KIN_kinmem, *ordering);
  KIN_ls = KIN_LS_KLU;
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_KLUReinit
 * ----------------------------------------------------------------
 */

void FKIN_KLUREINIT(int *neq, int *nnz, int *reinit_type, int *ier)
{
  *ier = KINKLUReInit(KIN_kinmem, *neq, *nnz, *reinit_type);
}

