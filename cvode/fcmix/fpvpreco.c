/* File fpvpreco.c:  Fortran/C interface routine for PVODE/CVSPGMR.
   Version of 28 December 2001 */

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "fpvode.h"   /* actual function names, prototypes, global vars.   */
#include "cvspgmr.h"  /* CVSpgmr prototype                                 */


/***************************************************************************/

/* C function CVPreco to interface between CVODE and a Fortran subroutine
   PVPRECO for setup of a Krylov preconditioner.
   Addresses of Nlocal, t, jok, gamma, h, uround, y, fy, ewt, vtemp1, vtemp2, 
   vtemp3, and the addresses jcurPtr and nfePtr, are passed to PVPRECO, using
   the macros N_VLOCLENGTH and N_VDATA.  A return flag ier from PVPRECO
   is returned by CVPreco.
   Auxiliary data is assumed to be communicated by Common. */

int CVPreco(integer N, real t, N_Vector y, N_Vector fy, boole jok,
            boole *jcurPtr, real gamma, N_Vector ewt, real h,
            real uround, long int *nfePtr, void *P_data,
            N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  real *ydata, *fydata, *ewtdata, *v1data, *v2data, *v3data;
  int ier = 0;
  integer Nlocal;

  Nlocal = N_VLOCLENGTH(y);
  ydata = N_VDATA(y);
  fydata = N_VDATA(fy);
  ewtdata = N_VDATA(ewt);
  v1data = N_VDATA(vtemp1);
  v2data = N_VDATA(vtemp2);
  v3data = N_VDATA(vtemp3);

  FCV_PRECO (&Nlocal, &t, ydata, fydata, &jok, jcurPtr, &gamma, ewtdata,
           &h, &uround, nfePtr, v1data, v2data, v3data, &ier);

  return(ier);
}
