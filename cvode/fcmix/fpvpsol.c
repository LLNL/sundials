/* File fpvpsol.c:  Fortran/C interface routines for PVODE/CVSPGMR
   Version of 28 December 2001 */

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "fpvode.h"   /* actual function names, prototypes, global vars.   */
#include "cvspgmr.h"  /* CVSpgmr prototype                                 */


/***************************************************************************/

/* C function CVPSol to interface between CVODE and a Fortran subroutine
   PVPSOL for solution of a Krylov preconditioner.
   Addresses of Nlocal, t, gamma, delta, lr, y, fy, vtemp, ewt, r, z, and the
   address nfePtr, are passed to PVPSOL, using the macros N_VLOCLENGTH and
   N_VDATA.  A return flag ier from PVPSOL is returned by CVPSol.
   Auxiliary data is assumed to be communicated by Common. */

int CVPSol(integer N, real t, N_Vector y, N_Vector fy, N_Vector vtemp,
           real gamma, N_Vector ewt, real delta, long int *nfePtr,
           N_Vector r, int lr, void *P_data, N_Vector z)
{
  real *ydata, *fydata, *vtdata, *ewtdata, *rdata, *zdata;
  int ier = 0;
  integer Nlocal;

  Nlocal = N_VLOCLENGTH(y);
  ydata = N_VDATA(y);
  fydata = N_VDATA(fy);
  vtdata = N_VDATA(vtemp);
  ewtdata = N_VDATA(ewt);
  rdata = N_VDATA(r);
  zdata = N_VDATA(z);

  FCV_PSOL (&Nlocal, &t, ydata, fydata, vtdata, &gamma, ewtdata, &delta,
          nfePtr, rdata, &lr, zdata, &ier);

  return(ier);
}
