/*******************************************************************
 *                                                                 *
 * File          : fnvector_serial.c                               *
 * Programmers   : Radu Serban @ LLNL                              *
 * Version of    : 26 June 2002                                    *
 *                                                                 *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This file, companion of nvector_serial.c contains the           *
 * implementation of the Fortran interface to M_EnvInit_Serial     *
 * and M_EnvFree_Serial.                                           *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h"
#include "nvector_serial.h"
#include "fnvector_serial.h"

/* Define global variable F2C_machEnv */
M_Env F2C_machEnv;

/* Fortran callable interfaces to M_EnvInit_Serial
   and M_EnvFree_Serial */

void F_MENVINITS(integertype *neq, int *ier)
{
 F2C_machEnv = M_EnvInit_Serial(*neq);

 *ier = (F2C_machEnv == NULL) ? -1 : 0 ;
}


void F_MENVFREES()
{
  M_EnvFree_Serial(F2C_machEnv);
}

