/******************************************************************
 * File          : fcvspgmr21.c                                   *
 * Programmers   : Alan C. Hindmarsh and Radu Serban @ LLNL       *
 * Version of    : 26 June 2002                                   *
 *----------------------------------------------------------------*
 *                                                                *
 * Fortran/C interface routine for CVODE/CVSPGMR, for the case of *
 * preconditioner setup and solve routines, and a Jtimes routine. *
 *                                                                *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of types realtype and integertype   */
#include "nvector.h"       /* definitions of type N_Vector and vector macros  */
#include "fcvode.h"        /* actual function names, prototypes, global vars. */
#include "cvspgmr.h"       /* CVSpgmr prototype                               */


/***************************************************************************/

void FCV_SPGMR21(int *pretype, int *gstype, int *maxl, realtype *delt, int *ier)
{
  /* Call CVSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     *pretype    is the preconditioner type
     *gstype     is the Gram-Schmidt process type
     *maxl       is the maximum Krylov dimension
     *delt       is the linear convergence tolerance factor
     CVPreco     is a pointer to the preconditioner setup routine
     CVPSol      is a pointer to the preconditioner solve routine
     NULL        is the pointer to P_data
     CVJtimes    is a pointer to the Jtimes routine
     NULL        is the pointer to jac_data                             */

  *ier = CVSpgmr(CV_cvodemem, *pretype, *gstype, *maxl, *delt,
                 CVPreco, CVPSol, NULL, CVJtimes, NULL);
}


/***************************************************************************/

void FCV_REINSPGMR21(int *pretype, int *gstype, int *maxl, realtype *delt,
                     int *ier)
{
  /* Call CVReInitSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     *pretype    is the preconditioner type
     *gstype     is the Gram-Schmidt process type
     *maxl       is the maximum Krylov dimension
     *delt       is the linear convergence tolerance factor
     CVPreco     is a pointer to the preconditioner setup routine
     CVPSol      is a pointer to the preconditioner solve routine
     NULL        is the pointer to P_data
     CVJtimes    is a pointer to the Jtimes routine
     NULL        is the pointer to jac_data                             */

  *ier = CVReInitSpgmr(CV_cvodemem, *pretype, *gstype, *maxl, *delt,
                       CVPreco, CVPSol, NULL, CVJtimes, NULL);
}
