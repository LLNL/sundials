/**********************************************************************
 * File          : fcvpsol.c                                          *
 * Programmers   : Alan C. Hindmarsh and Radu Serban @ LLNL           *
 * Version of    : 30 March 2003                                      *
 *--------------------------------------------------------------------*
 * This C function CVPSol is to interface between the CVSPGMR module  *
 * and the user-supplied preconditioner solve routine CVPSOL.         *
 * Note the use of the generic name FCV_PSOL below.                   *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of types realtype and integertype   */
#include "nvector.h"       /* definitions of type N_Vector and vector kernels */
#include "fcvode.h"        /* actual function names, prototypes, global vars. */
#include "cvspgmr.h"       /* CVSpgmr prototype                               */

/********************************************************************/

/* Prototype of the Fortran routine */
void FCV_PSOL(realtype*, realtype*, realtype*, realtype*, 
              realtype*, realtype*, realtype*, 
              long int*, realtype*, int*, realtype*, int*);


/***************************************************************************/

/* C function CVPSol to interface between CVODE and a Fortran subroutine
   CVPSOL for solution of a Krylov preconditioner.
   Addresses of t, gamma, delta, lr, y, fy, vtemp, ewt, r, z, and the
   address nfePtr, are passed to CVPSOL, using the routine N_VGetData 
   from NVECTOR.  A return flag ier from CVPSOL is returned by CVPSol.
   Auxiliary data is assumed to be communicated by Common. */

int CVPSol(realtype t, N_Vector y, N_Vector fy, N_Vector vtemp,
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

  FCV_PSOL(&t, ydata, fydata, vtdata, &gamma, ewtdata, &delta,
           nfePtr, rdata, &lr, zdata, &ier);

  N_VSetData(zdata, z);

  return(ier);
}
