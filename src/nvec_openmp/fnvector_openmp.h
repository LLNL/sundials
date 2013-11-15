/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2010-12-15 19:40:08 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Steven Smith @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2013, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
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
#define FNV_INITS_OPENMP    SUNDIALS_F77_FUNC(fnvinits_openmp, FNVINITS_OPENMP)
#else
#define FNV_INITS_OPENMP    fnvinits_openmp_
#endif

#if defined(SUNDIALS_F77_FUNC_)

#define FNV_INITS_OPENMP_Q  SUNDIALS_F77_FUNC_(fnvinits_openmp_q, FNVINITS_OPENMP_Q)
#define FNV_INITS_OPENMP_S  SUNDIALS_F77_FUNC_(fnvinits_openmp_s, FNVINITS_OPENMP_S)
#define FNV_INITS_OPENMP_B  SUNDIALS_F77_FUNC_(fnvinits_openmp_b, FNVINITS_OPENMP_B)
#define FNV_INITS_OPENMP_QB SUNDIALS_F77_FUNC_(fnvinits_openmp_qb, FNVINITS_OPENMP_QB)

#else

#define FNV_INITS_OPENMP_Q  fnvinits_openmp_q_
#define FNV_INITS_OPENMP_S  fnvinits_openmp_s_
#define FNV_INITS_OPENMP_B  fnvinits_openmp_b_
#define FNV_INITS_OPENMP_QB fnvinits_openmp_qb_

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
 * FNV_INITS_OPENMP    - initializes openmp vector operations for main problem
 * FNV_INITS_OPENMP_Q  - initializes openmp vector operations for quadratures
 * FNV_INITS_OPENMP_S  - initializes openmp vector operations for sensitivities
 * FNV_INITS_OPENMP_B  - initializes openmp vector operations for adjoint problem
 * FNV_INITS_OPENMP_QB - initializes openmp vector operations for adjoint quadratures
 *
 */

void FNV_INITS_OPENMP(int *code, long int *neq, int *num_threads, int *ier);
void FNV_INITS_OPENMP_Q(int *code, long int *Nq, int *num_threads, int *ier);
void FNV_INITS_OPENMP_S(int *code, int *Ns, int *ier);
void FNV_INITS_OPENMP_B(int *code, long int *NB, int *num_threads, int *ier);
void FNV_INITS_OPENMP_QB(int *code, long int *NqB, int *num_threads, int *ier);

#ifdef __cplusplus
}
#endif

#endif
