/*******************************************************************
 *                                                                 *
 * File          : fkinpreco.c                                     *
 * Programmers   : Allan G Taylor, Alan C. Hindmarsh, and          *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 8 March 2002                                    *
 *-----------------------------------------------------------------*
 * This C function KINPreco is to interface between routines that  *
 * call KINSpgmr(fkinspgmr20, fkinspgmr21) and the user-supplied   *
 * preconditioner setup routine KPRECO. Note the use of generic    *
 * names below (K_PRECO)                                           *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h"  /* definitions of types real and integer      */
#include "nvector.h"   /* definitions of type N_Vector               */
#include "kinsol.h"    /* KINSOL constants and prototypes            */
#include "fkinsol.h"   /* prototypes of interfaces, global variables */

/*********************************************************************/

/* Prototypes of the Fortran routines */
void K_PRECO(integer*, real*, real*, real*, real*, real*, real*, real*, long int*, int*);

/*********************************************************************/

/* C function KINPreco to interface between KINSpgmr and KPRECO, 
   the user-supplied Fortran preconditioner setup routine. */

int KINPreco(integer Neq, N_Vector uu, N_Vector uscale, 
             N_Vector fval, N_Vector fscale,
             N_Vector vtemp1, N_Vector vtemp2,
             SysFn func, real u_round,
             long int *nfePtr, void *P_data)
{
  real *udata,*uscaledata, *fdata, *fscaledata, *vtemp1data, *vtemp2data;
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
