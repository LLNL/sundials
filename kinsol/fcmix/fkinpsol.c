/******************************************************************
 * File          : fkinpsol.c                                     *
 * Programmers   : Allan G Taylor, Alan C. Hindmarsh, and         * 
 *                 Radu Serban @ LLNL                             *
 * Version of    : 07 February 2004                               *
 *----------------------------------------------------------------*
 * This routine interfaces between KINSOL and a user-supplied     *
 * Fortran routine FKPSOL.                                        *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h" /* definitions of type realtype               */
#include "nvector.h"       /* definitions of type N_Vector               */
#include "kinsol.h"        /* KINSOL constants and prototypes            */
#include "kinspgmr.h"
#include "fkinsol.h"       /* prototypes of interfaces, global variables */

/********************************************************************/

/* Prototype of the user-supplied Fortran routine */
void FK_PSOL(realtype*, realtype*, realtype*, realtype*, 
             realtype*, realtype*, int*);

/***************************************************************************/

void FKIN_SPGMRSETPSOL(int *flag, int *ier)
{
  if (*flag == 0) KINSpgmrSetPrecSolveFn(KIN_mem, NULL);
  else            KINSpgmrSetPrecSolveFn(KIN_mem, FKINPSol);
}


/********************************************************************/

/* C function FKINPSol to interface between KINSpgmr and FKPSOL, the user-
  supplied Fortran preconditioner solve routine. */

int FKINPSol(N_Vector uu, N_Vector uscale, 
             N_Vector fval, N_Vector fscale, 
             N_Vector vv, void *P_data,
             N_Vector ftem)
{
  realtype *udata,*uscaledata, *fdata,*fscaledata, *vvdata, *ftemdata;
  int retcode;

  udata      = N_VGetData(uu);
  uscaledata = N_VGetData(uscale);
  fdata      = N_VGetData(fval);
  fscaledata = N_VGetData(fscale);
  vvdata     = N_VGetData(vv);
  ftemdata   = N_VGetData(ftem);
  
  FK_PSOL(udata, uscaledata, fdata, fscaledata, vvdata, ftemdata, &retcode);
 
  N_VSetData(vvdata, vv);

  return(retcode);

}


