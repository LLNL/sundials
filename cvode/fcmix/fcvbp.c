/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2004-10-21 20:55:05 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * This module contains the routines necessary to interface with the
 * CVBANDPRE module and user-supplied Fortran routines.
 * The routines here call the generically named routines and provide 
 * a standard interface to the C code of the CVBANDPRE package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvbandpre.h"      /* prototypes of CVBANDPRE functions and macros */
#include "cvode.h"          /* CVODE constants and prototypes               */
#include "cvspgmr.h"        /* prototypes of CVSPGMR interface routines     */
#include "fcvbp.h"          /* prototypes of interfaces to CVBANDPRE        */
#include "fcvode.h"         /* actual function names, prototypes and
			       global variables                             */
#include "nvector.h"        /* definition of type N_Vector                  */
#include "sundialstypes.h"  /* definition of type realtype                  */

/***************************************************************************/

void FCV_BPINIT(long int *N, long int *mu, long int *ml, int *ier)
{
  /* 
     Call CVBandPrecAlloc to initialize the CVBANDPRE module:
     N      is the vector size
     mu, ml are the half-bandwidths of the retained preconditioner blocks
  */

  CVBP_Data = CVBandPrecAlloc(CV_cvodemem, *N, *mu, *ml);

  if (CVBP_Data == NULL) *ier = -1; 
  else                   *ier = 0;

  return;
}

/***************************************************************************/

void FCV_BPSPGMR(int *pretype, int *gstype, int *maxl, realtype *delt, int *ier)
{
  /* 
     Call CVBPSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     pretype    is the preconditioner type
     gstype     is the Gram-Schmidt process type
     maxl       is the maximum Krylov dimension
     delt       is the linear convergence tolerance factor 
  */

  *ier = CVBPSpgmr(CV_cvodemem, *pretype, *maxl, CVBP_Data);
  if (*ier != CVSPGMR_SUCCESS) return;

  *ier = CVSpgmrSetGSType(CV_cvodemem, *gstype);
  if (*ier != CVSPGMR_SUCCESS) return;

  *ier = CVSpgmrSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPGMR_SUCCESS) return;

  CV_ls = 4;
}

/***************************************************************************/

/* C function FCVBPOPT to access optional outputs from CVBANDPRE_Data */

void FCV_BPOPT(long int *lenrpw, long int *lenipw, long int *nfe)
{
  CVBandPrecGetWorkSpace(CVBP_Data, lenrpw, lenipw);
  CVBandPrecGetNumRhsEvals(CVBP_Data, nfe);
}

/***************************************************************************/

/* C function FCVBPFREE to interface to CVBandPrecFree, to free memory 
   created by CVBandPrecAlloc */

void FCV_BPFREE()
{
  CVBandPrecFree(CVBP_Data);
}
