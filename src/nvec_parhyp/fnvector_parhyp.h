/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date: 2016-01-13
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban and Alan Hindmarsh @ LLNL
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
