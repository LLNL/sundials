/*******************************************************************
 *                                                                 *
 * File          : fkinpreco.c                                     *
 * Programmers   : Allan G Taylor, Alan C. Hindmarsh, and          *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 5 August 2003                                   *
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
#include "kinspgmr.h"      
#include "fkinsol.h"       /* prototypes of interfaces, global variables    */

/*********************************************************************/

/* Prototype of the user-supplied Fortran routine */
void K_PRECO(realtype*, realtype*, realtype*, realtype*, 
             realtype*, realtype*, int*);


/***************************************************************************/

void F_KSPGMRSETPRECO(int *flag, int *ier)
{
  if (*flag == 0) KINSpgmrSetPrecSetupFn(KIN_mem, NULL);
  else            KINSpgmrSetPrecSetupFn(KIN_mem, KINPreco);
}

/*********************************************************************/

/* C function KINPreco to interface between KINSpgmr and KPRECO, 
   the user-supplied Fortran preconditioner setup routine. */

int KINPreco(N_Vector uu, N_Vector uscale,
             N_Vector fval, N_Vector fscale,
             void *P_data,
             N_Vector vtemp1, N_Vector vtemp2)
{
  realtype *udata,*uscaledata, *fdata, *fscaledata, *vtemp1data, *vtemp2data;
  int retcode;
  
 udata        = N_VGetData(uu);
 uscaledata   = N_VGetData(uscale);
 fdata        = N_VGetData(fval);
 fscaledata   = N_VGetData(fscale);
 vtemp1data   = N_VGetData(vtemp1);
 vtemp2data   = N_VGetData(vtemp2);
 
 K_PRECO(udata, uscaledata, fdata, fscaledata, vtemp1data, vtemp2data, &retcode);

 /* Note: there is no need to use N_VSetData since we are not getting back any
    information that should go into an N_Vector */

 return(retcode);
}
