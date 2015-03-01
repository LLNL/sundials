/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
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
 * the KINSuperLUMT solver. See fkinsol.h for usage.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "fkinsol.h"
#include "kinsol_impl.h"
#include <kinsol/kinsol_superlumt.h>
 
/*
 * ----------------------------------------------------------------
 * Function : FKIN_SUPERLUMT
 * ----------------------------------------------------------------
 */

void FKIN_SUPERLUMT(int *nthreads, int *neq, int *nnz, int *ordering, int *ier)
{
  *ier = KINSuperLUMT(KIN_kinmem, *neq, *nnz);
  KINSuperLUMTSetOrdering(KIN_kinmem, *ordering);
  KIN_ls = KIN_LS_SUPERLUMT;
}


