/*
 * -----------------------------------------------------------------
 * $Revision: 1.15 $
 * $Date: 2006-06-15 15:38:50 $
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

#include "fcvode.h"           /* actual fn. names, prototypes and global vars.*/
#include "fcvbp.h"            /* prototypes of interfaces to CVBANDPRE        */

#include "cvode.h"            /* CVODE constants and prototypes               */
#include "cvode_bandpre.h"    /* prototypes of CVBANDPRE functions and macros */

#include "cvode_sptfqmr.h"    /* prototypes of CVSPTFQMR interface routines   */
#include "cvode_spbcgs.h"     /* prototypes of CVSPBCG interface routines     */
#include "cvode_spgmr.h"      /* prototypes of CVSPGMR interface routines     */

#include "sundials_nvector.h" /* definition of type N_Vector                  */
#include "sundials_types.h"   /* definition of type realtype                  */

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

void FCV_BPSPTFQMR(int *pretype, int *maxl, realtype *delt, int *ier)
{
  /* 
     Call CVBPSptfqmr to specify the SPTFQMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     pretype    is the preconditioner type
     maxl       is the maximum Krylov dimension
     delt       is the linear convergence tolerance factor 
  */

  *ier = CVBPSptfqmr(CV_cvodemem, *pretype, *maxl, CVBP_Data);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPILS_SUCCESS) return;

  CV_ls = CV_LS_SPTFQMR;
}

/***************************************************************************/

void FCV_BPSPBCG(int *pretype, int *maxl, realtype *delt, int *ier)
{
  /* 
     Call CVBPSpbcg to specify the SPBCG linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     pretype    is the preconditioner type
     maxl       is the maximum Krylov dimension
     delt       is the linear convergence tolerance factor 
  */

  *ier = CVBPSpbcg(CV_cvodemem, *pretype, *maxl, CVBP_Data);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPILS_SUCCESS) return;

  CV_ls = CV_LS_SPBCG;
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
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetGSType(CV_cvodemem, *gstype);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPILS_SUCCESS) return;

  CV_ls = CV_LS_SPGMR;
}

/***************************************************************************/

/* C function FCVBPOPT to access optional outputs from CVBANDPRE_Data */

void FCV_BPOPT(long int *lenrwbp, long int *leniwbp, long int *nfebp)
{
  CVBandPrecGetWorkSpace(CVBP_Data, lenrwbp, leniwbp);
  CVBandPrecGetNumRhsEvals(CVBP_Data, nfebp);
}

/***************************************************************************/

/* C function FCVBPFREE to interface to CVBandPrecFree, to free memory 
   created by CVBandPrecAlloc */

void FCV_BPFREE(void)
{
  CVBandPrecFree(&CVBP_Data);
}
