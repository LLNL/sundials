/*******************************************************************
 * File          : fcvpreco.c                                      *
 * Programmers   : Alan C. Hindmarsh and Radu Serban @ LLNL        *
 * Version of    : 27 March 2002                                   *
 *-----------------------------------------------------------------*
 * This C function CVPreco is to interface between routines that   *
 * call CVSpgmr(fcvspgmr20, fcvspgmr21) and the user-supplied      *
 * preconditioner setup routine CVPRECO. Note the use of generic   *
 * names below (FCV_PRECO)                                         *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer             */
#include "nvector.h"  /* definitions of type N_Vector and vector macros    */
#include "fcvode.h"   /* actual function names, prototypes, global vars.   */
#include "cvspgmr.h"  /* CVSpgmr prototype                                 */

/*********************************************************************/

/* Prototypes of the Fortran routines */
void FCV_PRECO(integer*, real*, real*, real*, boole*, boole*, real*, real*, 
               real*, real*, long int *, real*, real*, real*, int*);

/***************************************************************************/

/* C function CVPreco to interface between CVODE and a Fortran subroutine
   CVPRECO for setup of a Krylov preconditioner.
   Addresses of Nlocal, t, jok, gamma, h, uround, y, fy, ewt, vtemp1, vtemp2, 
   vtemp3, and the addresses jcurPtr and nfePtr, are passed to CVPRECO, using
   the routine N_VGetData from NVECTOR.  A return flag ier from CVPRECO
   is returned by CVPreco.
   Auxiliary data is assumed to be communicated by Common. */

int CVPreco(integer N, real t, N_Vector y, N_Vector fy, boole jok,
            boole *jcurPtr, real gamma, N_Vector ewt, real h,
            real uround, long int *nfePtr, void *P_data,
            N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  real *ydata, *fydata, *ewtdata, *v1data, *v2data, *v3data;
  int ier = 0;

  ydata   = N_VGetData(y);
  fydata  = N_VGetData(fy);
  ewtdata = N_VGetData(ewt);
  v1data  = N_VGetData(vtemp1);
  v2data  = N_VGetData(vtemp2);
  v3data  = N_VGetData(vtemp3);

  FCV_PRECO (&N, &t, ydata, fydata, &jok, jcurPtr, &gamma, ewtdata,
           &h, &uround, nfePtr, v1data, v2data, v3data, &ier);

  /* Note: there is no need to use N_VSetData since we are not getting back any
    information that should go into an N_Vector */

  return(ier);
}
