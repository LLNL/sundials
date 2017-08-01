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
 * This file (companion of fsunlinsol_spfgmr.c) contains the
 * definitions needed for the initialization of SPFGMR
 * linear solver operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FSUNLINSOL_SPFGMR_H
#define _FSUNLINSOL_SPFGMR_H

#include <sunlinsol/sunlinsol_spfgmr.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FSUNSPFGMR_INIT     SUNDIALS_F77_FUNC(fsunspfgmrinit, FSUNSPFGMRINIT)
#define FSUNMASSSPFGMR_INIT SUNDIALS_F77_FUNC(fsunmassspfgmrinit, FSUNMASSSPFGMRINIT)
#else
#define FSUNSPFGMR_INIT     fsunspfgmrinit_
#define FSUNMASSSPFGMR_INIT fsunmassspfgmrinit_
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
 * FSUNSPFGMR_INIT - initializes SPFGMR linear solver for main problem
 * FSUNMASSSPFGMR_INIT - initializes SPFGMR linear solver for mass matrix solve
 */

void FSUNSPFGMR_INIT(int *code, int *pretype, int *maxl, int *ier);
void FSUNMASSSPFGMR_INIT(int *pretype, int *maxl, int *ier);

#ifdef __cplusplus
}
#endif

#endif
