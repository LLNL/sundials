/*
 * ----------------------------------------------------------------
 * $Revision: 1.10 $
 * $Date: 2004-06-18 21:36:28 $
 * ----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * ----------------------------------------------------------------
 * The C function FKINPSet is used to interface between KINSOL and
 * the Fortran user-supplied preconditioner setup routine.
 *
 * Note: The use of the generic name FK_PSET below.
 * ----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "sundialstypes.h"  /* definition of type realtype */
#include "nvector.h"        /* definition of type N_Vector */
#include "kinsol.h"         /* KINSOL constants and prototypes */
#include "kinspgmr.h"      
#include "fkinsol.h"        /* prototypes of interfaces and global variables */

/*
 * ----------------------------------------------------------------
 * prototype of the user-supplied fortran routine
 * ----------------------------------------------------------------
 */

void FK_PSET(realtype*, realtype*, realtype*, realtype*, 
             realtype*, realtype*, int*);

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SPGMRSETPSET
 * ----------------------------------------------------------------
 */

void FKIN_SPGMRSETPSET(int *flag, int *ier)
{
  if ((*flag) == 0) KINSpgmrSetPrecSetupFn(KIN_mem, NULL);
  else KINSpgmrSetPrecSetupFn(KIN_mem, FKINPSet);
}

/*
 * ----------------------------------------------------------------
 * Function : FKINPSet
 * ----------------------------------------------------------------
 * C function FKINPSet is used to interface between FK_PSET and
 * the user-supplied Fortran preconditioner setup routine.
 * ----------------------------------------------------------------
 */

int FKINPSet(N_Vector uu, N_Vector uscale,
             N_Vector fval, N_Vector fscale,
             void *P_data,
             N_Vector vtemp1, N_Vector vtemp2)
{
  realtype *udata,*uscaledata, *fdata, *fscaledata, *vtemp1data, *vtemp2data;
  int retcode;
  
  udata      = (realtype *) N_VGetData(uu);
  uscaledata = (realtype *) N_VGetData(uscale);
  fdata      = (realtype *) N_VGetData(fval);
  fscaledata = (realtype *) N_VGetData(fscale);
  vtemp1data = (realtype *) N_VGetData(vtemp1);
  vtemp2data = (realtype *) N_VGetData(vtemp2);
  
  FK_PSET(udata, uscaledata, fdata, fscaledata, vtemp1data, vtemp2data, &retcode);

 /* Note: There is no need to use N_VSetData since we are not getting back any
    information that should go into an N_Vector */

 return(retcode);
}
