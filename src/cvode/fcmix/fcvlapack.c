/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-11-10 21:04:11 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Fortran/C interface routines for CVODE/CVLAPACK
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"     /* actual fn. names, prototypes and global vars.*/
#include "cvode_impl.h" /* definition of CVodeMem type                  */

#include <cvode/cvode_lapack.h>

/***************************************************************************/

void FCV_LAPACKDENSE(int *neq, int *ier)
{
  /* neq  is the problem size */

  *ier = CVLapackDense(CV_cvodemem, *neq);

  CV_ls = CV_LS_LAPACKDENSE;
}

/***************************************************************************/

void FCV_LAPACKBAND(int *neq, int *mupper, int *mlower, int *ier)
{
  /* 
     neq        is the problem size
     mupper     is the upper bandwidth
     mlower     is the lower bandwidth 
  */

  *ier = CVLapackBand(CV_cvodemem, *neq, *mupper, *mlower);

  CV_ls = CV_LS_LAPACKBAND;
}

/***************************************************************************/

