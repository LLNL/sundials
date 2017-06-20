/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner, Carol Woodward, Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for a generic SUNMATRIX package.
 * It contains the implementation of the SUNMatrix operations listed
 * in sundials_matrix.h
 * -----------------------------------------------------------------
 */

#include <stdlib.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

/*
 * -----------------------------------------------------------------
 * Functions in the 'ops' structure
 * -----------------------------------------------------------------
 */

SUNMatrix_ID SUNMatrixGetID(SUNMatrix A)
{
  SUNMatrix_ID id;
  id = A->ops->getid(A);
  return(id);
}

SUNMatrix SUNMatrixClone(SUNMatrix A)
{
  SUNMatrix B = NULL;
  B = A->ops->clone(A);
  return(B);
}

void SUNMatrixDestroy(SUNMatrix A)
{
  if (A==NULL) return;
  A->ops->destroy(A);
  return;
}

int SUNMatrixZero(SUNMatrix A)
{
  return((int) A->ops->zero(A));
}

int SUNMatrixScale(realtype c, SUNMatrix A)
{
  return((int) A->ops->scale(c, A));
}

int SUNMatrixCopy(SUNMatrix A, SUNMatrix B)
{
  return((int) A->ops->copy(A, B));
}

int SUNMatrixAddIdentity(SUNMatrix A)
{
  return((int) A->ops->addidentity(A));
}

int SUNMatrixAdd(SUNMatrix A, SUNMatrix B)
{
  return((int) A->ops->add(A,B));
}

int SUNMatrixMatvec(SUNMatrix A, N_Vector x, N_Vector y)
{
  return((int) A->ops->matvec(A,x,y));
}

