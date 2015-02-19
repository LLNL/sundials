/*----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
 Copyright (c) 2013, Southern Methodist University.
 All rights reserved.
 For details, see the LICENSE file.
 -----------------------------------------------------------------
 This file (companion of nvector_parallel_grid.h) contains the
 definitions needed for the initialization of parallel grid
 vector operations in Fortran.
 ---------------------------------------------------------------*/

#ifndef _FNVECTOR_PARALLEL_GRID_H
#define _FNVECTOR_PARALLEL_GRID_H

#include <nvector/nvector_parallel_grid.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FNV_INITP    SUNDIALS_F77_FUNC(fnvinitp, FNVINITP)
#else
#define FNV_INITP    fnvinitp_
#endif

#if defined(SUNDIALS_F77_FUNC_)

#define FNV_INITPG_Q  SUNDIALS_F77_FUNC_(fnvinitpg_q,  FNVINITPG_Q)
#define FNV_INITPG_S  SUNDIALS_F77_FUNC_(fnvinitpg_s,  FNVINITPG_S)
#define FNV_INITPG_B  SUNDIALS_F77_FUNC_(fnvinitpg_b,  FNVINITPG_B)
#define FNV_INITPG_QB SUNDIALS_F77_FUNC_(fnvinitpg_qb, FNVINITPG_QB)

#else

#define FNV_INITPG_Q  fnvinitpg_q_
#define FNV_INITPG_S  fnvinitpg_s_
#define FNV_INITPG_B  fnvinitpg_b_
#define FNV_INITPG_QB fnvinitpg_qb_

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
 * FNV_INITPG    - initializes parallel grid vector operations for main problem
 * FNV_INITPG_Q  - initializes parallel grid vector operations for quadratures
 * FNV_INITPG_S  - initializes parallel grid vector operations for sensitivities
 * FNV_INITPG_B  - initializes parallel grid vector operations for adjoint problem
 * FNV_INITPG_QB - initializes parallel grid vector operations for adjoint quadratures
 *
 */

#ifndef SUNDIALS_MPI_COMM_F2C
#define MPI_Fint int
#endif

void FNV_INITPG(MPI_Fint *comm, int *code, long int *dims, long int *dim_len, 
		long int *dim_alen, long int *dim_off, long int *F_ordering, 
		long int *glob_len, int *ier);
void FNV_INITPG_Q(MPI_Fint *comm, int *code, long int *dims, long int *dim_lenQ, 
		  long int *dim_alenQ, long int *dim_offQ, long int *F_ordering, 
		  long int *glob_lenQ, int *ier);
void FNV_INITPG_B(MPI_Fint *comm, int *code, long int *dims, long int *dim_lenB, 
		  long int *dim_alenB, long int *dim_offB, long int *F_ordering, 
		  long int *glob_lenB, int *ier);
void FNV_INITPG_QB(MPI_Fint *comm, int *code, long int *dims, long int *dim_lenQB, 
		   long int *dim_alenQB, long int *dim_offQB, long int *F_ordering, 
		   long int *glob_lenQB, int *ier);
void FNV_INITPG_S(int *code, long int *Ns, int *ier);

#ifdef __cplusplus
}
#endif

#endif
