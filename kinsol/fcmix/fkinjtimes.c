/*
 * -----------------------------------------------------------------
 * $Revision: 1.16 $
 * $Date: 2006-06-15 15:39:24 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Routines used to interface between KINSOL and a Fortran
 * user-supplied routine FKJTIMES (Jacobian J times vector v).
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "kinsol_spils.h"
#include "fkinsol.h"
#include "kinsol_impl.h"

#include "sundials_nvector.h"
#include "sundials_types.h"

/*
 * ----------------------------------------------------------------
 * prototype of the user-supplied fortran routine
 * ----------------------------------------------------------------
 */
#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FK_JTIMES(realtype*, realtype*, int*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/*
 * ----------------------------------------------------------------
 * Function : FKIN_SPILSSETJAC
 * ----------------------------------------------------------------
 */

void FKIN_SPILSSETJAC(int *flag, int *ier)
{
  if ((*flag) == 0) KINSpilsSetJacTimesVecFn(KIN_kinmem, NULL, NULL);
  else              KINSpilsSetJacTimesVecFn(KIN_kinmem, FKINJtimes, NULL);

  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKINJtimes
 * ----------------------------------------------------------------
 * C function FKINJtimes is used to interface between
 * KINSp* / KINSp*JTimes and FK_JTIMES (user-supplied Fortran
 * routine).
 * ----------------------------------------------------------------
 */

int FKINJtimes(N_Vector v, N_Vector Jv,
               N_Vector uu, booleantype *new_uu, 
               void *J_data)
{
  int retcode;
  realtype *vdata, *Jvdata, *uudata;

  vdata = Jvdata = uudata = NULL;

  vdata  = N_VGetArrayPointer(v);
  uudata = N_VGetArrayPointer(uu);
  Jvdata = N_VGetArrayPointer(Jv);
 
  FK_JTIMES(vdata, Jvdata, (int *) new_uu, uudata, &retcode);

  return(retcode);
}
