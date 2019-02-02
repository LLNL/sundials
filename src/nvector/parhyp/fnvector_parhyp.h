/*
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban and Alan Hindmarsh @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This file (companion of nvector_parhyp.c) contains the
 * definitions needed for the initialization of pahyp
 * vector operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FNVECTOR_PARHYP_H
#define _FNVECTOR_PARHYP_H

#include <nvector/nvector_parhyp.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FNV_INITPH   SUNDIALS_F77_FUNC(fnvinitph, FNVINITPH)
#else
#define FNV_INITPH   fnvinitph_
#endif

/* Declarations of global variables */

extern N_Vector F2C_CVODE_vec;

extern N_Vector F2C_IDA_vec;

extern N_Vector F2C_KINSOL_vec;

extern N_Vector F2C_ARKODE_vec;

/*
 * Prototype of exported function
 * FNV_INITPH   - initializes parhyp vector operations for main problem
 */

#ifndef SUNDIALS_MPI_COMM_F2C
#define MPI_Fint int
#endif

void FNV_INITPH(MPI_Fint *comm, int *code, long int *L, long int *N, int *ier);

#ifdef __cplusplus
}
#endif

#endif
