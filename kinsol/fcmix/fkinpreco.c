/*******************************************************************
 *                                                                 *
 * File          : fkinpreco.c                                     *
 * Programmers   : Allan G Taylor and Alan C. Hindmarsh @ LLNL     *
 * Version of    : 21 December 2001                                *
 *-----------------------------------------------------------------*
 * This C function KINPreco is to interface between routines that  *
 * call KINSpgmr(fkinspgmr20, fkinspgmr21) and the user-supplied   *
 * preconditioner setup routine KPRECO. Note the use of generic    *
 * names below (K_PRECO)                                           *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h"  /* definitions of types real and integer            */
#include "nvector.h"   /* definitions of type N_Vector and related macros  */
#include "kinsol.h"    /* KINSOL constants and prototypes                  */
#include "fkinsol.h"   /* prototypes of interfaces, global variables       */

/***************************************************************************/


/* C function KINPreco to interface between KINSpgmr and KPRECO, the user-
  supplied Fortran preconditioner setup routine. */

int KINPreco(integer Neq, N_Vector uu, N_Vector uscale, 
             N_Vector fval, N_Vector fscale,
	     N_Vector vtemp1, N_Vector vtemp2,
	     SysFn func, real u_round,
	     long int *nfePtr, void *P_data)
{
  real *udata,*uscaledata, *fdata, *fscaledata, *vtemp1data, *vtemp2data;
  int retcode;

 udata        = N_VDATA(uu);
 uscaledata   = N_VDATA(uscale);
 fdata        = N_VDATA(fval);
 fscaledata   = N_VDATA(fscale);
 vtemp1data   = N_VDATA(vtemp1);
 vtemp2data   = N_VDATA(vtemp2);
 
 K_PRECO(&Neq, udata, uscaledata, fdata, fscaledata, vtemp1data, vtemp2data,
       &u_round, nfePtr, &retcode);

 return(retcode);
}
