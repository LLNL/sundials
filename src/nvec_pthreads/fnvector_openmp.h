/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2010-12-15 19:40:08 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Steven Smith @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2013, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This file (companion of nvector_openmp.h) contains the
 * definitions needed for the initialization of openmp
 * vector operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FNVECTOR_OPENMP_H
#define _FNVECTOR_OPENMP_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <nvector/nvector_openmp.h>  
#include <sundials/sundials_fnvector.h>

#if defined(SUNDIALS_F77_FUNC)
#define FNV_INITS    SUNDIALS_F77_FUNC(fnvinits, FNVINITS)
#else
#define FNV_INITS    fnvinits_
#endif

#if defined(SUNDIALS_F77_FUNC_)

#define FNV_INITS_Q  SUNDIALS_F77_FUNC_(fnvinits_q, FNVINITS_Q)
#define FNV_INITS_S  SUNDIALS_F77_FUNC_(fnvinits_s, FNVINITS_S)
#define FNV_INITS_B  SUNDIALS_F77_FUNC_(fnvinits_b, FNVINITS_B)
#define FNV_INITS_QB SUNDIALS_F77_FUNC_(fnvinits_qb, FNVINITS_QB)

#else

#define FNV_INITS_Q  fnvinits_q_
#define FNV_INITS_S  fnvinits_s_
#define FNV_INITS_B  fnvinits_b_
#define FNV_INITS_QB fnvinits_qb_

#endif

/* Declarations of global variables */

extern N_Vector F2C_CVODE_vec;
extern N_Vector F2C_CVODE_vecQ;
extern N_Vector *F2C_CVODE_vecS;
extern N_Vector F2C_CVODE_vecB;
extern N_Vector F2C_CVODE_vecQB;

extern N_Vector F2C_IDA_vec;
extern N_Vector F2C_IDA_vecQ;
extern N_Vector *F2C_IDA_vecS;
extern N_Vector F2C_IDA_vecB;
extern N_Vector F2C_IDA_vecQB;

extern N_Vector F2C_KINSOL_vec;

extern N_Vector F2C_ARKODE_vec;

/* 
 * Prototypes of exported functions 
 *
 * FNV_INITS    - initializes openmp vector operations for main problem
 * FNV_INITS_Q  - initializes openmp vector operations for quadratures
 * FNV_INITS_S  - initializes openmp vector operations for sensitivities
 * FNV_INITS_B  - initializes openmp vector operations for adjoint problem
 * FNV_INITS_QB - initializes openmp vector operations for adjoint quadratures
 *
 */

void FNV_INITS(int *code, long int *neq, int *ier);
void FNV_INITS_Q(int *code, long int *Nq, int *ier);
void FNV_INITS_S(int *code, int *Ns, int *ier);
void FNV_INITS_B(int *code, long int *NB, int *ier);
void FNV_INITS_QB(int *code, long int *NqB, int *ier);

#ifdef __cplusplus
}
#endif

#endif
