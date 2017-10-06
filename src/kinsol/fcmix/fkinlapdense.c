/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
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

extern void FK_DJAC(long int*, realtype*, realtype*, realtype*,
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
    *ier = KINDlsSetDenseJacFn(KIN_kinmem, NULL);
  }
  else {
    *ier = KINDlsSetDenseJacFn(KIN_kinmem, FKINLapackDenseJac);
  }
  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKINLapackDenseJac
 * ----------------------------------------------------------------
 * C function FKINLapackDenseJac interfaces between KINSOL and a
 * Fortran subroutine FKDJAC for solution of a linear system
 * with dense Jacobian approximation using lapack functinos.
 * Addresses are passed to FKDJAC, using the macro DENSE_COL 
 * and the routine N_VGetArrayPointer from NVECTOR. 
 * Auxiliary data is assumed to be communicated by common blocks.
 * ----------------------------------------------------------------
 */

int FKINLapackDenseJac(long int N, N_Vector uu, N_Vector fval,
                       DlsMat J, void *user_data, N_Vector vtemp1, N_Vector vtemp2)
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
