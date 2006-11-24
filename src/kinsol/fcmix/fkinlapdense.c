/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-11-24 19:09:25 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Fortran/C interface routines for KINSOL/KINLAPACKDENSE, for the
 * case of a user-supplied Jacobian approximation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fkinsol.h"     /* prototypes of standard interfaces and global vars.*/
#include "kinsol_impl.h" /* definition of KINMem type                         */

#include <kinsol/kinsol_lapack.h>

/*
 * ----------------------------------------------------------------
 * prototypes of the user-supplied fortran routines
 * ----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FK_DJAC(int*, realtype*, realtype*, realtype*,
                    realtype*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/*
 * ----------------------------------------------------------------
 * Function : FKIN_LAPACKDENSESETJAC
 * ----------------------------------------------------------------
 */

void FKIN_LAPACKDENSESETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = KINDlsSetJacFn(KIN_kinmem, NULL, NULL);
  }
  else {
    *ier = KINDlsSetJacFn(KIN_kinmem, ( void *)FKINLapackDenseJac, NULL);
  }
  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKINLapackDenseJac
 * ----------------------------------------------------------------
 * C function FKINLapackDenseJac interfaces between KINSOL and a
 * Fortran subroutine FKLDJAC for solution of a linear system
 * with dense Jacobian approximation using lapack functinos.
 * Addresses are passed to FKDJAC, using the macro DENSE_COL 
 * and the routine N_VGetArrayPointer from NVECTOR. 
 * Auxiliary data is assumed to be communicated by common blocks.
 * ----------------------------------------------------------------
 */

int FKINLapackDenseJac(int N, N_Vector uu, N_Vector fval,
                       DlsMat J, void *jac_data, N_Vector vtemp1, N_Vector vtemp2)
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
