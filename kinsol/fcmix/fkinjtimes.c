/*
 * ----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2004-06-18 21:36:28 $
 * ----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * ----------------------------------------------------------------
 * Routines used to interface between KINSOL and a Fortran
 * user-supplied routine FKJTIMES (Jacobian J times vector v).
 * ----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h"
#include "nvector.h"
#include "kinsol.h"   /* prototypes of interfaces and global variables */
#include "kinspgmr.h"
#include "fkinsol.h"  /* prototypes of interfaces and global variables */

/*
 * ----------------------------------------------------------------
 * prototype of the user-supplied fortran routine
 * ----------------------------------------------------------------
 */

void FK_JTIMES(realtype*, realtype*, int*, realtype*, int*);

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SPGMRSETJAC
 * ----------------------------------------------------------------
 */

void FKIN_SPGMRSETJAC(int *flag, int *ier)
{
  if ((*flag) == 0) KINSpgmrSetJacTimesVecFn(KIN_mem, NULL);
  else KINSpgmrSetJacTimesVecFn(KIN_mem, FKINJtimes);
}

/*
 * ----------------------------------------------------------------
 * Function : FKINJtimes
 * ----------------------------------------------------------------
 * C function FKINJtimes is used to interface between
 * KINSpgmr/KINSpgmrJTimes and FK_JTIMES (user-supplied Fortran
 * routine).
 * ----------------------------------------------------------------
 */

int FKINJtimes(N_Vector v, N_Vector Jv,
               N_Vector uu, booleantype *new_uu, 
               void *J_data)
{
 int retcode;
 realtype *vdata, *Jvdata, *uudata;
 
 vdata  = (realtype *) N_VGetData(v);
 uudata = (realtype *) N_VGetData(uu);
 Jvdata = (realtype *) N_VGetData(Jv);
 
 FK_JTIMES(vdata, Jvdata, (int *)new_uu, uudata, &retcode);

 N_VSetData((void *)Jvdata, Jv);

 return(retcode);
}
