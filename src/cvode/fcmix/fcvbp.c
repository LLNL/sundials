/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2007-04-27 18:56:27 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This module contains the routines necessary to interface with the
 * CVBANDPRE module and user-supplied Fortran routines.
 * The routines here call the generically named routines and provide 
 * a standard interface to the C code of the CVBANDPRE package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"                 /* actual fn. names, prototypes and global vars.*/
#include "fcvbp.h"                  /* prototypes of interfaces to CVBANDPRE        */

#include <cvode/cvode_bandpre.h>    /* prototypes of CVBANDPRE functions and macros */
#include <cvode/cvode_sptfqmr.h>    /* prototypes of CVSPTFQMR interface routines   */
#include <cvode/cvode_spbcgs.h>     /* prototypes of CVSPBCG interface routines     */
#include <cvode/cvode_spgmr.h>      /* prototypes of CVSPGMR interface routines     */

/***************************************************************************/

void FCV_BPINIT(long int *N, long int *mu, long int *ml, int *ier)
{
  /* 
     Call CVBandPrecInit to initialize the CVBANDPRE module:
     N      is the vector size
     mu, ml are the half-bandwidths of the retained preconditioner blocks
  */

  *ier = CVBandPrecInit(CV_cvodemem, *N, *mu, *ml);

  return;
}

/***************************************************************************/

/* C function FCVBPOPT to access optional outputs from CVBANDPRE_Data */

void FCV_BPOPT(long int *lenrwbp, long int *leniwbp, long int *nfebp)
{
  CVBandPrecGetWorkSpace(CV_cvodemem, lenrwbp, leniwbp);
  CVBandPrecGetNumRhsEvals(CV_cvodemem, nfebp);
}
