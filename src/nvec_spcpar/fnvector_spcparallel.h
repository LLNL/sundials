/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
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
 * This file (companion of nvector_spcparallel.c) contains the
 * definitions needed for the initialization of the spcparallel
 * vector operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FNVECTOR_SPCPARALLEL_H
#define _FNVECTOR_SPCPARALLEL_H

#include <nvector/nvector_spcparallel.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(F77_FUNC)

#define FNV_INITSPCP    F77_FUNC(fnvinitspcp, FNVINITSPCP)
#define FNV_INITSPCP_Q  F77_FUNC_(fnvinitspcp_q, FNVINITSPCP_Q)
#define FNV_INITSPCP_S  F77_FUNC_(fnvinitspcp_s, FNVINITSPCP_S)
#define FNV_INITSPCP_B  F77_FUNC_(fnvinitspcp_b, FNVINITSPCP_B)
#define FNV_INITSPCP_QB F77_FUNC_(fnvinitspcp_qb, FNVINITSPCP_QB)

#else

#define FNV_INITSPCP    fnvinitspcp_
#define FNV_INITSPCP_Q  fnvinitspcp_q_
#define FNV_INITSPCP_S  fnvinitspcp_s_
#define FNV_INITSPCP_B  fnvinitspcp_b_
#define FNV_INITSPCP_QB fnvinitspcp_qb_

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
 * FNV_INITSPCP    - initializes parallel vector operations for main problem
 * FNV_INITSPCP_Q  - initializes parallel vector operations for quadratures
 * FNV_INITSPCP_S  - initializes parallel vector operations for sensitivities
 * FNV_INITSPCP_B  - initializes parallel vector operations for adjoint problem
 * FNV_INITSPCP_QB - initializes parallel vector operations for adjoint quadratures
 *
 */

#ifndef SUNDIALS_MPI_COMM_F2C
#define MPI_Fint int
#endif

void FNV_INITSPCP(MPI_Fint *comm, int *code, 
		  int *Ngrp, int *Nspc,
		  indextype *Nx, indextype *Ny, indextype *Nz,
		  indextype *NGx, indextype *NGy, indextype *NGz,
		  int *ier);
void FNV_INITSPCP_Q(MPI_Fint *comm, int *code, 
		    int *NgrpQ, int *NspcQ,
		    indextype *NxQ, indextype *NyQ, indextype *NzQ,
		    indextype *NGxQ, indextype *NGyQ, indextype *NGzQ,
		    int *ier);
void FNV_INITSPCP_B(MPI_Fint *comm, int *code, 
		    int *NgrpB, int *NspcB,
		    indextype *NxB, indextype *NyB, indextype *NzB,
		    indextype *NGxB, indextype *NGyB, indextype *NGzB,
		    int *ier);
void FNV_INITSPCP_QB(MPI_Fint *comm, int *code, 
		     int *NgrpQB, int *NspcQB,
		     indextype *NxQB, indextype *NyQB, indextype *NzQB,
		     indextype *NGxQB, indextype *NGyQB, indextype *NGzQB,
		     int *ier);
void FNV_INITSPCP_S(int *code, int *Ns, int *ier);

#ifdef __cplusplus
}
#endif

#endif
