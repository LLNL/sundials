/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This file (companion of fsunlinsol_superlu.c) contains the
 * definitions needed for the initialization of superlu
 * linear solver operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FSUNLINSOL_SUPERLU_H
#define _FSUNLINSOL_SUPERLU_H

#include <sunlinsol/sunlinsol_superlu.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FSUNSUPERLU_INIT            SUNDIALS_F77_FUNC(fsunsuperluinit,        FSUNSUPERLUINIT)
#define FSUNSUPERLU_SETORDERING     SUNDIALS_F77_FUNC(fsunsuperlusetordering, FSUNSUPERLUSETORDERING)
#define FSUNMASSSUPERLU_INIT        SUNDIALS_F77_FUNC(fsunmasssuperluinit,        FSUNMASSSUPERLUINIT)
#define FSUNMASSSUPERLU_SETORDERING SUNDIALS_F77_FUNC(fsunmasssuperlusetordering, FSUNMASSSUPERLUSETORDERING)
#else
#define FSUNSUPERLU_INIT            fsunsuperluinit_
#define FSUNSUPERLU_SETORDERING     fsunsuperlusetordering_
#define FSUNMASSSUPERLU_INIT        fsunmasssuperluinit_
#define FSUNMASSSUPERLU_SETORDERING fsunmasssuperlusetordering_
#endif


/* Declarations of global variables */

extern SUNLinearSolver F2C_CVODE_linsol;
extern SUNLinearSolver F2C_IDA_linsol;
extern SUNLinearSolver F2C_KINSOL_linsol;
extern SUNLinearSolver F2C_ARKODE_linsol;
extern SUNLinearSolver F2C_ARKODE_mass_sol;

/* 
 * Prototypes of exported functions 
 *
 * FSUNSUPERLU_INIT - initializes superlu linear solver for main problem
 * FSUNSUPERLU_SETORDERING - sets the ordering choice used by SUPERLU for main problem
 * FSUNMASSSUPERLU_INIT - initializes superlu linear solver for mass matrix
 * FSUNMASSSUPERLU_SETORDERING - sets the ordering choice used by SUPERLU for mass matrix
 */

void FSUNSUPERLU_INIT(int *code, int *ier);
void FSUNSUPERLU_SETORDERING(int *code, int *ordering, int *ier);
void FSUNMASSSUPERLU_INIT(int *ier);
void FSUNMASSSUPERLU_SETORDERING(int *ordering, int *ier);

#ifdef __cplusplus
}
#endif

#endif
