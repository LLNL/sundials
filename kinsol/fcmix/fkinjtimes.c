/******************************************************************
 *                                                                *
 * File          : fkinjtimes.c                                   *
 * Programmers   : Allan G Taylor, Alan C. Hindmarsh, and         *
 *                 Radu Serban @ LLNL                             *
 * Version of    : 27 January 2004                                *
 *----------------------------------------------------------------*
 * Routine used to interface between KINSOL and a Fortran         *
 * user-supplied routine FKJTIMES (Jacobian J times vector v).    *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h"
#include "nvector.h"
#include "kinsol.h"  /* prototypes of interfaces, global variables  */
#include "kinspgmr.h"
#include "fkinsol.h" /* prototypes of interfaces, global variables  */

/* Prototype of the user-supplied Fortran routine */
void FK_JTIMES(realtype*, realtype*, int*, realtype*, int*);

/***************************************************************************/

void FKIN_SPGMRSETJAC(int *flag, int *ier)
{
  if (*flag == 0) KINSpgmrSetJacTimesVecFn(KIN_mem, NULL);
  else            KINSpgmrSetJacTimesVecFn(KIN_mem, FKINJtimes);
}

/**************************************************************************
 *  C function FKINJtimes is to interface between KINSpgmr/KINSpgmrJTimes  *
 *  and FKJTIMES, the user-supplied Fortran.                               *
 **************************************************************************/

int FKINJtimes(N_Vector v, N_Vector Jv,
               N_Vector uu, booleantype *new_uu, 
               void *J_data)
{
 int retcode;
 realtype *vdata, *Jvdata, *uudata;
 
 vdata  = N_VGetData(v);
 uudata = N_VGetData(uu);
 Jvdata = N_VGetData(Jv);
 
 FK_JTIMES(vdata, Jvdata, (int *)new_uu, uudata, &retcode);

 N_VSetData(Jvdata, Jv);

 return(retcode);
}
