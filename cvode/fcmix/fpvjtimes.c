/* File fpvjtimes.c:  Fortran/C interface routine for PVODE/CVSPGMR.
   Version of 28 December 2001 */

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "fpvode.h"   /* actual function names, prototypes, global vars.   */
#include "cvspgmr.h"  /* CVSpgmr prototype                                 */


/***************************************************************************/

/* C function  CVJtimes to interface between PVODE and  user-supplied
   Fortran routine PVJTIMES for Jacobian*vector product.
   Addresses of Nlocal, v, Jv, t, y, fy, vnrm, ewt, h, uround, ytemp, and
   the address nfePtr, are passed to PVJTIMES, using the macros N_VLOCLENGTH
   and N_VDATA.  A return flag ier from PVJTIMES is returned by CVJtimes.
   Auxiliary data is assumed to be communicated by Common. */

int CVJtimes(integer N, N_Vector v, N_Vector Jv, RhsFn f, 
             void *f_data, real t, N_Vector y, N_Vector fy,
             real vnrm, N_Vector ewt, real h, real uround, 
             void *jac_data, long int *nfePtr, N_Vector work)
{

  real *vdata, *Jvdata, *ydata, *fydata, *ewtdata, *wkdata;
  int ier = 0;
  integer Nlocal;

  Nlocal = N_VLOCLENGTH(y);
  vdata = N_VDATA(v);
  Jvdata = N_VDATA(Jv);
  ydata = N_VDATA(y);
  fydata = N_VDATA(fy);
  ewtdata = N_VDATA(ewt);
  wkdata = N_VDATA(work);

  FCV_JTIMES (&Nlocal, vdata, Jvdata, &t, ydata, fydata, &vnrm,
              ewtdata, &h, &uround, nfePtr, wkdata, &ier);

  return(ier);
}
