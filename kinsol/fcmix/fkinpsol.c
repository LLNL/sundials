 /******************************************************************
 *                                                                *
 * File          : fkinpsol.c                                     *
 * Programmers   : Allan G Taylor and Alan C. Hindmarsh @ LLNL    *
 * Version of    : 21 December 2001                               *
 *----------------------------------------------------------------*
 * This routine interfaces between the user-supplied Fortran      *
 * routine KPSOL and the various routines that call KINSpgmr.     *
 * See the routines fkinspgmr10.c, fkinspgmr11.c, fkinspgmr20.c,  *
 * and fkinspgmr21.c                                              *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h" /* definitions of types real and integer            */
#include "nvector.h"  /* definitions of type N_Vector and related macros  */
#include "kinsol.h"   /* KINSOL constants and prototypes                  */
#include "fkinsol.h"  /* prototypes of interfaces, global variables       */

/**************************************************************************/

/* C function KINPSol to interface between KINSpgmr and KPSOL, the user-
  supplied Fortran preconditioner solve routine. */

int KINPSol(integer Neq, N_Vector uu, N_Vector uscale, 
             N_Vector fval, N_Vector fscale,
	     N_Vector vtem, N_Vector ftem,
	     SysFn func, real u_round,
	     long int *nfePtr, void *P_data)
{
  real *udata,*uscaledata, *fdata,*fscaledata, *vtemdata, *ftemdata;
  int retcode;

 udata      = N_VDATA(uu);
 uscaledata = N_VDATA(uscale);
 fdata      = N_VDATA(fval);
 fscaledata = N_VDATA(fscale);
 vtemdata   = N_VDATA(vtem);
 ftemdata   = N_VDATA(ftem);
 
 K_PSOL(&Neq, udata, uscaledata, fdata, fscaledata, vtemdata, ftemdata,
       &u_round, nfePtr,&retcode);

 return(retcode);

}


