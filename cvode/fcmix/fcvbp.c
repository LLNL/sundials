/****************************************************************************
 * File         : fcvbp.c                                                   *
 * Programmers  : Radu Serban @ LLNL                                        *
 * Version of   : 19 February 2004                                          *
 ****************************************************************************
 *                                                                          *
 * This module contains the routines necessary to interface with the        *
 * CVBANDPRE module and user-supplied Fortran routines.                     *
 * The routines here call the generically named routines and provide a      *
 * standard interface to the C code of the CVBANDPRE package.               *
 *                                                                          *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of type realtype                  */
#include "nvector.h"       /* definitions of type N_Vector                  */
#include "cvode.h"         /* CVODE constants and prototypes                */
#include "fcvode.h"        /* actual function names, prototypes, global vars*/
#include "fcvbp.h"         /* prototypes of interfaces to CVBANDPRE         */
#include "cvspgmr.h"       /* prototypes of CVSPGMR interface routines      */
#include "cvbandpre.h"     /* prototypes of CVBANDPRE functions, macros     */

/***************************************************************************/

void FCV_BPINIT(long int *N, long int *mu, long int *ml,
                     int *pretype, int *gstype, int *maxl, realtype *delt, 
                     int *ier)
{

  /* 
     First call CVBandPrecAlloc to initialize the CVBANDPRE module:
     *N          is the vector size
     *mu, *ml    are the half-bandwidths of the retained preconditioner blocks
     */

  CVBP_Data = CVBandPrecAlloc(CV_cvodemem, *N, *mu, *ml);
  if (CVBP_Data == NULL) { *ier = -1; return; }

  /* Call CVBPSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     *pretype    is the preconditioner type
     *gstype     is the Gram-Schmidt process type
     *maxl       is the maximum Krylov dimension
     *delt       is the linear convergence tolerance factor */

  *ier = CVBPSpgmr(CV_cvodemem, *pretype, *maxl, CVBP_Data);
  if (*ier != 0) return;

  *ier = CVSpgmrSetGSType(CV_cvodemem, *gstype);
  if (*ier != 0) return;

  *ier = CVSpgmrSetDelt(CV_cvodemem, *delt);
  if (*ier != 0) return;

  CV_ls = 4;

}

/***************************************************************************/

/* C function FCVBPOPT to access optional outputs from CVBANDPRE_Data */

void FCV_BPOPT(long int *lenrpw, long int *lenipw, long int *nfe)
{
  CVBandPrecGetIntWorkSpace(CVBP_Data, lenipw);
  CVBandPrecGetRealWorkSpace(CVBP_Data, lenrpw);
  CVBandPrecGetNumRhsEvals(CVBP_Data, nfe);

}

/***************************************************************************/

/* C function FCVBPFREE to interface to CVBandPrecFree, to free memory 
   created by CVBandPrecAlloc */

void FCV_BPFREE()
{
  CVBandPrecFree(CVBP_Data);
}
