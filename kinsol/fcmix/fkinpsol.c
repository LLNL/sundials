/******************************************************************
 *                                                                *
 * File          : fkinpsol.c                                     *
 * Programmers   : Allan G Taylor, Alan C. Hindmarsh, and         * 
 *                 Radu Serban @ LLNL                             *
 * Version of    : 31 MArch 2003                                  *
 *----------------------------------------------------------------*
 * This routine interfaces between KINSOL and a user-supplied     *
 * Fortran routine KPSOL for the various cases involving KINSpgmr.*
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of types realtype and integertype */
#include "nvector.h"       /* definitions of type N_Vector                  */
#include "kinsol.h"        /* KINSOL constants and prototypes               */
#include "fkinsol.h"       /* prototypes of interfaces, global variables    */

/********************************************************************/

/* Prototype of the user-supplied Fortran routine */
void K_PSOL(realtype*, realtype*, realtype*, realtype*, 
            realtype*, realtype*, realtype*, long int*, int*);

/********************************************************************/

/* C function KINPSol to interface between KINSpgmr and KPSOL, the user-
  supplied Fortran preconditioner solve routine. */

int KINPSol(N_Vector uu, N_Vector uscale, 
            N_Vector fval, N_Vector fscale,
            N_Vector vv, N_Vector ftem,
            SysFn func, realtype u_round,
            long int *nfePtr, void *P_data)
{
  realtype *udata,*uscaledata, *fdata,*fscaledata, *vvdata, *ftemdata;
  int retcode;

  udata      = N_VGetData(uu);
  uscaledata = N_VGetData(uscale);
  fdata      = N_VGetData(fval);
  fscaledata = N_VGetData(fscale);
  vvdata     = N_VGetData(vv);
  ftemdata   = N_VGetData(ftem);
  
  K_PSOL(udata, uscaledata, fdata, fscaledata, vvdata, ftemdata,
         &u_round, nfePtr,&retcode);
 
  N_VSetData(vvdata, vv);

  return(retcode);

}


