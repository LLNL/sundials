/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2004-07-22 21:12:43 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban, LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/shared/LICENSE
 * -----------------------------------------------------------------
 * This file (companion of nvector_serial.h) contains the 
 * implementation needed for the Fortran initialization of serial 
 * vector operations
 * -----------------------------------------------------------------
 */
#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h"
#include "nvector_serial.h"
#include "fnvector_serial.h"

/* Define global variable F2C_vec */
N_Vector F2C_vec;

/* Fortran callable interfaces */

void FNV_INITS(long int *neq, int *ier)
{
 F2C_vec = N_VNewEmpty_Serial(*neq);

 *ier = (F2C_vec == NULL) ? -1 : 0 ;
}


void FNV_FREES()
{
  N_VDestroyEmpty_Serial(F2C_vec);
}

