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
 * This file (companion of fsunlinsol_superlumt.h) contains the
 * implementation needed for the Fortran initialization of superlumt
 * linear solver operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fsunlinsol_superlumt.h"

/* Define global linsol variables */

SUNLinearSolver F2C_CVODE_linsol;
SUNLinearSolver F2C_IDA_linsol;
SUNLinearSolver F2C_KINSOL_linsol;
SUNLinearSolver F2C_ARKODE_linsol;

/* Declarations of external global variables */

extern SUNMatrix F2C_CVODE_matrix;
extern SUNMatrix F2C_IDA_matrix;
extern SUNMatrix F2C_KINSOL_matrix;
extern SUNMatrix F2C_ARKODE_matrix;

extern N_Vector F2C_CVODE_vec;
extern N_Vector F2C_IDA_vec;
extern N_Vector F2C_KINSOL_vec;
extern N_Vector F2C_ARKODE_vec;

/* Fortran callable interfaces */

void FSUNSUPERLUMT_INIT(int *code, int *ier, int *num_threads)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_linsol = NULL;
    F2C_CVODE_linsol = SUNSuperLUMT(F2C_CVODE_vec,
                                    F2C_CVODE_matrix,
                                    *num_threads);
    if (F2C_CVODE_linsol == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_linsol = NULL;
    F2C_IDA_linsol = SUNSuperLUMT(F2C_IDA_vec,
                                  F2C_IDA_matrix,
                                  *num_threads);
    if (F2C_IDA_linsol == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_linsol = NULL;
    F2C_KINSOL_linsol = SUNSuperLUMT(F2C_KINSOL_vec,
                                     F2C_KINSOL_matrix,
                                     *num_threads);
    if (F2C_KINSOL_linsol == NULL) *ier = -1;
    break;
  case FCMIX_ARKODE:
    F2C_ARKODE_linsol = NULL;
    F2C_ARKODE_linsol = SUNSuperLUMT(F2C_ARKODE_vec,
                                     F2C_ARKODE_matrix,
                                     *num_threads);
    if (F2C_ARKODE_linsol == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}


void FSUNSUPERLUMT_SETORDERING(int *code, int *ordering_choice, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    *ier = SUNSuperLUMTSetOrdering(F2C_CVODE_linsol, *ordering_choice);
    break;
  case FCMIX_IDA:
    *ier = SUNSuperLUMTSetOrdering(F2C_IDA_linsol, *ordering_choice);
    break;
  case FCMIX_KINSOL:
    *ier = SUNSuperLUMTSetOrdering(F2C_KINSOL_linsol, *ordering_choice);
    break;
  case FCMIX_ARKODE:
    *ier = SUNSuperLUMTSetOrdering(F2C_ARKODE_linsol, *ordering_choice);
    break;
  default:
    *ier = -1;
  }
}
