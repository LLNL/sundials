/******************************************************************
 * File          : fcvtimes.c                                     *
 * Programmers   : Alan C. Hindmarsh and Radu Serban @ LLNL       *
 * Version of    : 27 March 2002                                  *
 *----------------------------------------------------------------*
 * Routine used to interface between a Fortran main and a user-   *
 * supplied J-times routine FCVJTIMES (Jacobian J times vector v) *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "fcvode.h"   /* actual function names, prototypes, global vars.   */
#include "cvspgmr.h"  /* CVSpgmr prototype                                 */

/* Prototypes of the Fortran routines */
void FCV_JTIMES(integer*, real*, real*, real*, real*, real*, real*,
                real*, real*, real*, long int*, real*, int*);

/***************************************************************************/

/* C function  CVJtimes to interface between CVODE and  user-supplied
   Fortran routine CVJTIMES for Jacobian*vector product.
   Addresses of N, v, Jv, t, y, fy, vnrm, ewt, h, uround, ytemp, and
   the address nfePtr, are passed to CVJTIMES, using the routine
   N_VGetData from NVECTOR. A return flag ier from CVJTIMES is 
   returned by CVJtimes.
   Auxiliary data is assumed to be communicated by Common. */

int CVJtimes(integer N, N_Vector v, N_Vector Jv, RhsFn f, 
             void *f_data, real t, N_Vector y, N_Vector fy,
             real vnrm, N_Vector ewt, real h, real uround, 
             void *jac_data, long int *nfePtr, N_Vector work)
{

  real *vdata, *Jvdata, *ydata, *fydata, *ewtdata, *wkdata;
  int ier = 0;

  vdata = N_VGetData(v);
  Jvdata = N_VGetData(Jv);
  ydata = N_VGetData(y);
  fydata = N_VGetData(fy);
  ewtdata = N_VGetData(ewt);
  wkdata = N_VGetData(work);

  FCV_JTIMES (&N, vdata, Jvdata, &t, ydata, fydata, &vnrm,
              ewtdata, &h, &uround, nfePtr, wkdata, &ier);

  N_VSetData(Jvdata, Jv);

  return(ier);
}
