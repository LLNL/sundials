/*
 * -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This file (companion of fnvector_parallel_grid.h) contains the 
 * implementation needed for the Fortran initialization of parallel 
 * grid vector operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fnvector_parallel_grid.h"

/* Define global vector variables */

N_Vector F2C_CVODE_vec;
N_Vector F2C_CVODE_vecQ;
N_Vector *F2C_CVODE_vecS;
N_Vector F2C_CVODE_vecB;
N_Vector F2C_CVODE_vecQB;

N_Vector F2C_IDA_vec;
N_Vector F2C_IDA_vecQ;
N_Vector *F2C_IDA_vecS;
N_Vector F2C_IDA_vecB;
N_Vector F2C_IDA_vecQB;

N_Vector F2C_KINSOL_vec;

N_Vector F2C_ARKODE_vec;

#ifndef SUNDIALS_MPI_COMM_F2C
#define MPI_Fint int
#endif

/* Fortran callable interfaces */

void FNV_INITPG(MPI_Fint *comm, int *code, long int *dims, long int *dim_len, 
		long int *dim_alen, long int *dim_off, long int *F_ordering, 
		long int *glob_len, int *ier)
{
  MPI_Comm F2C_comm;

#ifdef SUNDIALS_MPI_COMM_F2C
  F2C_comm = MPI_Comm_f2c(*comm);
#else
  F2C_comm = MPI_COMM_WORLD;
#endif

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vec = NULL;
    F2C_CVODE_vec = N_VNewEmpty_Parallel_Grid(F2C_comm, *dims, dim_len, dim_alen, 
					      dim_off, *F_ordering, *glob_len);
    if (F2C_CVODE_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vec = NULL;
    F2C_IDA_vec = N_VNewEmpty_Parallel_Grid(F2C_comm, *dims, dim_len, dim_alen, 
					    dim_off, *F_ordering, *glob_len);
    if (F2C_IDA_vec == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_vec = NULL;
    F2C_KINSOL_vec = N_VNewEmpty_Parallel_Grid(F2C_comm, *dims, dim_len, dim_alen, 
					       dim_off, *F_ordering, *glob_len);
    if (F2C_KINSOL_vec == NULL) *ier = -1;
    break;
  case FCMIX_ARKODE:
    F2C_ARKODE_vec = NULL;
    F2C_ARKODE_vec = N_VNewEmpty_Parallel_Grid(F2C_comm, *dims, dim_len, dim_alen, 
					       dim_off, *F_ordering, *glob_len);
    if (F2C_ARKODE_vec == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITPG_Q(MPI_Fint *comm, int *code, long int *dims, long int *dim_lenQ, 
		  long int *dim_alenQ, long int *dim_offQ, long int *F_ordering, 
		  long int *glob_lenQ, int *ier)
{
  MPI_Comm F2C_comm;

#ifdef SUNDIALS_MPI_COMM_F2C
  F2C_comm = MPI_Comm_f2c(*comm);
#else
  F2C_comm = MPI_COMM_WORLD;
#endif

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecQ = NULL;
    F2C_CVODE_vecQ = N_VNewEmpty_Parallel_Grid(F2C_comm, *dims, dim_lenQ, dim_alenQ, 
					       dim_offQ, *F_ordering, *glob_lenQ);
    if (F2C_CVODE_vecQ == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecQ = NULL;
    F2C_IDA_vecQ = N_VNewEmpty_Parallel_Grid(F2C_comm, *dims, dim_lenQ, dim_alenQ, 
					     dim_offQ, *F_ordering, *glob_lenQ);
    if (F2C_IDA_vecQ == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITPG_B(MPI_Fint *comm, int *code, long int *dims, long int *dim_lenB, 
		  long int *dim_alenB, long int *dim_offB, long int *F_ordering, 
		  long int *glob_lenB, int *ier)
{
  MPI_Comm F2C_comm;

#ifdef SUNDIALS_MPI_COMM_F2C
  F2C_comm = MPI_Comm_f2c(*comm);
#else
  F2C_comm = MPI_COMM_WORLD;
#endif

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecB = NULL;
    F2C_CVODE_vecB = N_VNewEmpty_Parallel_Grid(F2C_comm, *dims, dim_lenB, dim_alenB, 
					       dim_offB, *F_ordering, *glob_lenB);
    if (F2C_CVODE_vecB == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecB = NULL;
    F2C_IDA_vecB = N_VNewEmpty_Parallel_Grid(F2C_comm, *dims, dim_lenB, dim_alenB, 
					     dim_offB, *F_ordering, *glob_lenB);
    if (F2C_IDA_vecB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITPG_QB(MPI_Fint *comm, int *code, long int *dims, long int *dim_lenQB, 
		   long int *dim_alenQB, long int *dim_offQB, long int *F_ordering, 
		   long int *glob_lenQB, int *ier)
{
  MPI_Comm F2C_comm;

#ifdef SUNDIALS_MPI_COMM_F2C
  F2C_comm = MPI_Comm_f2c(*comm);
#else
  F2C_comm = MPI_COMM_WORLD;
#endif


  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecQB = NULL;
    F2C_CVODE_vecQB = N_VNewEmpty_Parallel_Grid(F2C_comm, *dims, dim_lenQB, dim_alenQB, 
						dim_offQB, *F_ordering, *glob_lenQB);
    if (F2C_CVODE_vecQB == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecQB = NULL;
    F2C_IDA_vecQB = N_VNewEmpty_Parallel_Grid(F2C_comm, *dims, dim_lenQB, dim_alenQB, 
					      dim_offQB, *F_ordering, *glob_lenQB);
    if (F2C_IDA_vecQB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITPG_S(int *code, long int *Ns, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecS = NULL;
    F2C_CVODE_vecS = (N_Vector *) N_VCloneVectorArrayEmpty_Parallel_Grid(*Ns, F2C_CVODE_vec);
    if (F2C_CVODE_vecS == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecS = NULL;
    F2C_IDA_vecS = (N_Vector *) N_VCloneVectorArrayEmpty_Parallel_Grid(*Ns, F2C_IDA_vec);
    if (F2C_IDA_vecS == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}
