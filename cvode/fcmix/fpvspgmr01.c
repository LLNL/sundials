/* File fpvspgmr01.c:  Fortran/C interface routine for PVODE/CVSPGMR, for
   the case of no preconditioner routines, but a user Jtimes routine.
   Version of 11 January 2002 */

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "fpvode.h"   /* actual function names, prototypes, global vars.   */
#include "cvspgmr.h"  /* CVSpgmr prototype                                 */


/***************************************************************************/

void FCV_SPGMR01 (int *pretype, int *gstype, int *maxl, real *delt, int *ier)
{
  /* Call CVSpgmr to specify the SPGMR linear solver:
     CV_cvodemem is the pointer to the CVODE memory block
     *pretype    is the preconditioner type
     *gstype     is the Gram-Schmidt process type
     *maxl       is the maximum Krylov dimension
     *delt       is the linear convergence tolerance factor
     NULL        is a pointer to the preconditioner setup routine
     NULL        is a pointer to the preconditioner solve routine
     NULL        is the pointer to P_data
     CVJtimes    is a pointer to the Jtimes routine
     NULL        is the pointer to jac_data                             */

  *ier = CVSpgmr (CV_cvodemem, *pretype, *gstype, *maxl, *delt, NULL, NULL,
                  NULL, CVJtimes, NULL);
}
