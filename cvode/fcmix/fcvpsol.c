/******************************************************************
 * File          : fcvpsol.c                                      *
 * Programmers   : Alan C. Hindmarsh and Radu Serban @ LLNL       *
 * Version of    : 26 June 2002                                   *
 *----------------------------------------------------------------*
 * This routine interfaces between the user-supplied Fortran      *
 * routine CVPSOL and the various routines that call CVSpgmr.     *
 * See the routines fcvspgmr10.c, fcvspgmr11.c, fcvspgmr20.c,     *
 * and fcvspgmr21.c                                               *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of types realtype and integertype     */
#include "nvector.h"       /* definitions of type N_Vector and vector kernels   */
#include "fcvode.h"        /* actual function names, prototypes, global vars.   */
#include "cvspgmr.h"       /* CVSpgmr prototype                                 */

/********************************************************************/

/* Prototypes of the Fortran routines */
void FCV_PSOL(integertype*, realtype*, realtype*, realtype*, realtype*, 
              realtype*, realtype*, realtype*, 
              long int*, realtype*, int*, realtype*, int*);


/***************************************************************************/

/* C function CVPSol to interface between CVODE and a Fortran subroutine
   CVPSOL for solution of a Krylov preconditioner.
   Addresses of N, t, gamma, delta, lr, y, fy, vtemp, ewt, r, z, and the
   address nfePtr, are passed to CVPSOL, using the routine N_VGetData 
   from NVECTOR.  A return flag ier from CVPSOL is returned by CVPSol.
   Auxiliary data is assumed to be communicated by Common. */

int CVPSol(integertype N, realtype t, N_Vector y, N_Vector fy, N_Vector vtemp,
           realtype gamma, N_Vector ewt, realtype delta, long int *nfePtr,
           N_Vector r, int lr, void *P_data, N_Vector z)
{
  realtype *ydata, *fydata, *vtdata, *ewtdata, *rdata, *zdata;
  int ier = 0;

  ydata = N_VGetData(y);
  fydata = N_VGetData(fy);
  vtdata = N_VGetData(vtemp);
  ewtdata = N_VGetData(ewt);
  rdata = N_VGetData(r);
  zdata = N_VGetData(z);

  FCV_PSOL(&N, &t, ydata, fydata, vtdata, &gamma, ewtdata, &delta,
           nfePtr, rdata, &lr, zdata, &ier);

  N_VSetData(zdata, z);

  return(ier);
}
