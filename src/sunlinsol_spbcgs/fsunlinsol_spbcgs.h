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
 * This file (companion of fsunlinsol_spbcgs.c) contains the
 * definitions needed for the initialization of SPBCGS
 * linear solver operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FSUNLINSOL_SPBCGS_H
#define _FSUNLINSOL_SPBCGS_H

#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FSUNSPBCGS_INIT   SUNDIALS_F77_FUNC(fsunspbcgsinit, FSUNSPBCGSINIT)
#else
#define FSUNSPBCGS_INIT   fsunspbcgsinit_
#endif


/* Declarations of global variables */

extern SUNLinearSolver F2C_CVODE_linsol;
extern SUNLinearSolver F2C_IDA_linsol;
extern SUNLinearSolver F2C_KINSOL_linsol;
extern SUNLinearSolver F2C_ARKODE_linsol;

/* 
 * Prototypes of exported functions 
 *
 * FSUNSPBCGS_INIT - initializes SPBCGS linear solver for main problem
 */

void FSUNSPBCGS_INIT(int *code, int *pretype, int *maxl, int *ier);

#ifdef __cplusplus
}
#endif

#endif
