/*******************************************************************
 *                                                                 *
 * File          : fnvector_serial.c                               *
 * Programmers   : Radu Serban @ LLNL                              *
 * Version of    : 07 February 2004                                *
 *                                                                 *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This file, companion of nvector_serial.c contains the           *
 * implementation of the Fortran interface to NV_SpecInit_Serial   *
 * and NV_SpecFree_Serial.                                         *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h"
#include "nvector_serial.h"
#include "fnvector_serial.h"

/* Define global variable F2C_nvspec */
NV_Spec F2C_nvspec;

/* Fortran callable interfaces to NV_SpecInit_Serial
   and NV_SpecFree_Serial */

void FNV_INITS(long int *neq, int *ier)
{
 F2C_nvspec = NV_SpecInit_Serial(*neq);

 *ier = (F2C_nvspec == NULL) ? -1 : 0 ;
}


void FNV_FREES()
{
  NV_SpecFree_Serial(F2C_nvspec);
}

