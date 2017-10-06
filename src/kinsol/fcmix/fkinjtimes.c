/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * Routines used to interface between KINSOL and a Fortran
 * user-supplied routine FKJTIMES (Jacobian J times vector v).
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fkinsol.h"
#include "kinsol_impl.h"

#include <kinsol/kinsol_spils.h>

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
  if ((*flag) == 0) KINSpilsSetJacTimesVecFn(KIN_kinmem, NULL);
  else              KINSpilsSetJacTimesVecFn(KIN_kinmem, FKINJtimes);

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
               void *user_data)
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
