 /*****************************************************************
 *                                                                *
 * File          : fkinpsol.c                                     *
 * Programmers   : Allan G Taylor, Alan C. Hindmarsh, and         * 
 *                 Radu Serban @ LLNL                             *
 * Version of    : 27 June 2002                                   *
 *----------------------------------------------------------------*
 * This routine interfaces between the user-supplied Fortran      *
 * routine KPSOL and the various routines that call KINSpgmr.     *
 * See the routines fkinspgmr10.c, fkinspgmr11.c, fkinspgmr20.c,  *
 * and fkinspgmr21.c                                              *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of types realtype and integertype */
#include "nvector.h"       /* definitions of type N_Vector                  */
#include "kinsol.h"        /* KINSOL constants and prototypes               */
#include "fkinsol.h"       /* prototypes of interfaces, global variables    */

/********************************************************************/

/* Prototypes of the Fortran routines */
void K_PSOL(integertype*, realtype*, realtype*, realtype*, realtype*, 
            realtype*, realtype*, realtype*, long int*, int*);

/********************************************************************/

/* C function KINPSol to interface between KINSpgmr and KPSOL, the user-
  supplied Fortran preconditioner solve routine. */

int KINPSol(integertype Neq, N_Vector uu, N_Vector uscale, 
            N_Vector fval, N_Vector fscale,
            N_Vector vtem, N_Vector ftem,
            SysFn func, realtype u_round,
            long int *nfePtr, void *P_data)
{
  realtype *udata,*uscaledata, *fdata,*fscaledata, *vtemdata, *ftemdata;
  int retcode;

  udata      = N_VGetData(uu);
  uscaledata = N_VGetData(uscale);
  fdata      = N_VGetData(fval);
  fscaledata = N_VGetData(fscale);
  vtemdata   = N_VGetData(vtem);
  ftemdata   = N_VGetData(ftem);
  
  K_PSOL(&Neq, udata, uscaledata, fdata, fscaledata, vtemdata, ftemdata,
         &u_round, nfePtr,&retcode);
 
  N_VSetData(vtemdata, vtem);

  return(retcode);

}


