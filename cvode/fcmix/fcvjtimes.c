/**********************************************************************
 * File          : fcvjtimes.c                                        *
 * Programmers   : Alan C. Hindmarsh and Radu Serban @ LLNL           *
 * Version of    : 30 March 2003                                      *
 *--------------------------------------------------------------------*
 * This C function CVJtimes is to interface between the CVSPGMR module*
 * and the user-supplied Jacobian-times-vector routine CVJTIMES.      *
 * Note the use of the generic name FCV_JTIMES below.                 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of types realtype and integertype   */
#include "nvector.h"       /* definitions of type N_Vector and vector macros  */
#include "fcvode.h"        /* actual function names, prototypes, global vars. */
#include "cvspgmr.h"       /* CVSpgmr prototype                               */

/* Prototype of the Fortran routine */
void FCV_JTIMES(realtype*, realtype*, realtype*, realtype*, 
                realtype*, realtype*, int*);

/***************************************************************************/

/* C function  CVJtimes to interface between CVODE and  user-supplied
   Fortran routine CVJTIMES for Jacobian*vector product.
   Addresses of v, Jv, t, y, fy, vnrm, ewt, h, uround, ytemp, and
   the address nfePtr, are passed to CVJTIMES, using the routine
   N_VGetData from NVECTOR. A return flag ier from CVJTIMES is 
   returned by CVJtimes.
   Auxiliary data is assumed to be communicated by Common. */

int CVJtimes(N_Vector v, N_Vector Jv, realtype t, 
             N_Vector y, N_Vector fy,
             void *jac_data, N_Vector work)
{

  realtype *vdata, *Jvdata, *ydata, *fydata, *ewtdata, *wkdata;
  int ier = 0;

  vdata = N_VGetData(v);
  Jvdata = N_VGetData(Jv);
  ydata = N_VGetData(y);
  fydata = N_VGetData(fy);
  wkdata = N_VGetData(work);

  FCV_JTIMES (vdata, Jvdata, &t, ydata, fydata, wkdata, &ier);

  N_VSetData(Jvdata, Jv);

  return(ier);
}
