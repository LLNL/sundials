/*******************************************************************
 *                                                                 *
 * File          : fkinpreco.c                                     *
 * Programmers   : Allan G Taylor, Alan C. Hindmarsh, and          *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 29 July 2002                                    *
 *-----------------------------------------------------------------*
 * This C function KINPreco is to interface between KINSOL and the *
 * Fortran user-supplied preconditioner setup routine.             *
 * Note the use of generic name K_PRECO below.                     *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of types realtype and integertype */
#include "nvector.h"       /* definitions of type N_Vector                  */
#include "kinsol.h"        /* KINSOL constants and prototypes               */
#include "fkinsol.h"       /* prototypes of interfaces, global variables    */

/*********************************************************************/

/* Prototype of the user-supplied Fortran routine */
void K_PRECO(integertype*, realtype*, realtype*, realtype*, realtype*, 
             realtype*, realtype*, realtype*, long int*, int*);

/*********************************************************************/

/* C function KINPreco to interface between KINSpgmr and KPRECO, 
   the user-supplied Fortran preconditioner setup routine. */

int KINPreco(integertype Neq, N_Vector uu, N_Vector uscale, 
             N_Vector fval, N_Vector fscale,
             N_Vector vtemp1, N_Vector vtemp2,
             SysFn func, realtype u_round,
             long int *nfePtr, void *P_data)
{
  realtype *udata,*uscaledata, *fdata, *fscaledata, *vtemp1data, *vtemp2data;
  int retcode;
  
 udata        = N_VGetData(uu);
 uscaledata   = N_VGetData(uscale);
 fdata        = N_VGetData(fval);
 fscaledata   = N_VGetData(fscale);
 vtemp1data   = N_VGetData(vtemp1);
 vtemp2data   = N_VGetData(vtemp2);
 
 K_PRECO(&Neq, udata, uscaledata, fdata, fscaledata, vtemp1data, vtemp2data,
         &u_round, nfePtr, &retcode);

 /* Note: there is no need to use N_VSetData since we are not getting back any
    information that should go into an N_Vector */

 return(retcode);
}
