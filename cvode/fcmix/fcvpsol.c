/**********************************************************************
 * File          : fcvpsol.c                                          *
 * Programmers   : Alan C. Hindmarsh and Radu Serban @ LLNL           *
 * Version of    : 27 January 2004                                    *
 *--------------------------------------------------------------------*
 * The C function FCVPSol is to interface between the CVSPGMR module  *
 * and the user-supplied preconditioner solve routine FCVPSOL.        *
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
              realtype*, int*, realtype*, int*);

/***************************************************************************/

void FCV_SPGMRSETPSOL(int *flag, int *ier)
{
  if (*flag == 0) CVSpgmrSetPrecSolveFn(CV_cvodemem, NULL);
  else            CVSpgmrSetPrecSolveFn(CV_cvodemem, FCVPSol);
}

/***************************************************************************/

/* C function FCVPSol to interface between CVODE and a Fortran subroutine
   FCVPSOL for solution of a Krylov preconditioner.
   Addresses of t, gamma, delta, lr, y, fy, vtemp, ewt, r, z, and the
   address nfePtr, are passed to FCVPSOL, using the routine N_VGetData 
   from NVECTOR.  A return flag ier from FCVPSOL is returned by FCVPSol.
   Auxiliary data is assumed to be communicated by common blocks. */

int FCVPSol(realtype t, N_Vector y, N_Vector fy, 
            N_Vector r, N_Vector z,
            realtype gamma, realtype delta,
            int lr, void *P_data, N_Vector vtemp)
{
  N_Vector ewt;
  realtype *ydata, *fydata, *vtdata, *ewtdata, *rdata, *zdata;
  int ier = 0;

  CVodeGetErrWeights(CV_cvodemem, &ewt);

  ydata = N_VGetData(y);
  fydata = N_VGetData(fy);
  vtdata = N_VGetData(vtemp);
  ewtdata = N_VGetData(ewt);
  rdata = N_VGetData(r);
  zdata = N_VGetData(z);

  FCV_PSOL(&t, ydata, fydata, vtdata, &gamma, ewtdata, &delta,
           rdata, &lr, zdata, &ier);

  N_VSetData(zdata, z);

  return(ier);
}
