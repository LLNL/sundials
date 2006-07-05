/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-07-05 15:32:37 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Fortran/C interface routines for KINSOL/KINDENSE, for the case
 * of a user-supplied Jacobian approximation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fkinsol.h"     /* prototypes of standard interfaces and global vars.*/
#include "kinsol_impl.h" /* definition of KINMem type                         */

#include <kinsol/kinsol_dense.h>

/*
 * ----------------------------------------------------------------
 * prototypes of the user-supplied fortran routines
 * ----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FK_DJAC(long int*, realtype*, realtype*, realtype*,
		    realtype*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/*
 * ----------------------------------------------------------------
 * Function : FKIN_DENSESETJAC
 * ----------------------------------------------------------------
 */

void FKIN_DENSESETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = KINDenseSetJacFn(KIN_kinmem, NULL, NULL);
  }
  else {
    *ier = KINDenseSetJacFn(KIN_kinmem, (KINDenseJacFn) FKINDenseJac, NULL);
  }

  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKINDenseJac
 * ----------------------------------------------------------------
 * C function FKINDenseJac interfaces between KINSOL and a Fortran
 * subroutine FKDJAC for solution of a linear system with dense
 * Jacobian approximation. Addresses are passed to FKDJAC, using
 * the macro DENSE_COL from DENSE and the routine N_VGetArrayPointer
 * from NVECTOR. Auxiliary data is assumed to be communicated by
 * Common.
 * ----------------------------------------------------------------
 */

int FKINDenseJac(long int N, DenseMat J, N_Vector uu, N_Vector fval,
		 void *jac_data, N_Vector vtemp1, N_Vector vtemp2)
{
  realtype *uu_data, *fval_data, *jacdata, *v1_data, *v2_data;
  int ier;

  /* Initialize all pointers to NULL */
  uu_data = fval_data = jacdata = v1_data = v2_data = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  /* Get pointers to vector data */
  uu_data   = N_VGetArrayPointer(uu);
  fval_data = N_VGetArrayPointer(fval);
  v1_data   = N_VGetArrayPointer(vtemp1);
  v2_data   = N_VGetArrayPointer(vtemp2);

  jacdata = DENSE_COL(J,0);

  /* Call user-supplied routine */
  FK_DJAC(&N, uu_data, fval_data, jacdata, v1_data, v2_data, &ier);

  return(ier);
}
