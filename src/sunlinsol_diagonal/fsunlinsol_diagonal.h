/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds, Ashley Crawford @ SMU
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
 * This file (companion of fsunlinsol_diagonal.c) contains the
 * definitions needed for the initialization of diagonal
 * linear solver operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FSUNLINSOL_DIAGONAL_H
#define _FSUNLINSOL_DIAGONAL_H

#include <sunlinsol/sunlinsol_diagonal.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FSUNDIAGLINSOL_INIT     SUNDIALS_F77_FUNC(fsundiaglinsolinit, FSUNDIAGLINSOLINIT)
#define FSUNMASSDIAGLINSOL_INIT SUNDIALS_F77_FUNC(fsunmassdiaglinsolinit, FSUNMASSDIAGLINSOLINIT)
#else
#define FSUNDIAGLINSOL_INIT     fsundiaglinsolinit_
#define FSUNMASSDIAGLINSOL_INIT fsunmassdiaglinsolinit_
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
 * FSUNDIAGLINSOL_INIT - initializes diagonal linear solver for main problem
 * FSUNMASSDIAGLINSOL_INIT - initializes diagonal linear solver for mass matrix solve
 */

void FSUNDIAGLINSOL_INIT(int *code, int *ier);
void FSUNMASSDIAGLINSOL_INIT(int *ier);

#ifdef __cplusplus
}
#endif

#endif
