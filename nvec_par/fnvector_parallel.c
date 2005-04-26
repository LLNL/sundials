/*
 * -----------------------------------------------------------------
 * $Revision: 1.10 $
 * $Date: 2005-04-26 23:42:05 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/shared/LICENSE.
 * -----------------------------------------------------------------
 * This file (companion of nvector_serial.h) contains the 
 * implementation needed for the Fortran initialization of parallel 
 * vector operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fnvector_parallel.h"
#include "mpi.h"
#include "nvector_parallel.h"
#include "sundialstypes.h"

/* Define global vector variables */

N_Vector F2C_CVODE_vec;

N_Vector F2C_CVODES_vec;
N_Vector F2C_CVODES_vecQ;
N_Vector *F2C_CVODES_vecS;
N_Vector F2C_CVODES_vecB;
N_Vector F2C_CVODES_vecQB;

N_Vector F2C_IDA_vec;

N_Vector F2C_IDAS_vec;
N_Vector F2C_IDAS_vecQ;
N_Vector *F2C_IDAS_vecS;
N_Vector F2C_IDAS_vecB;
N_Vector F2C_IDAS_vecQB;

N_Vector F2C_KINSOL_vec;

/* Fortran callable interfaces */

#ifdef SUNDIALS_MPI_COMM_F2C

void FNV_INITP(MPI_Fint *comm, int *code, long int *L, long int *N, int *ier)
{
  MPI_Comm F2C_comm = MPI_Comm_f2c(*comm);

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vec = N_VNewEmpty_Parallel(F2C_comm, *L, *N);
    if (F2C_CVODE_vec == NULL) *ier = -1;
    break;
  case FCMIX_CVODES:
    F2C_CVODES_vec = N_VNewEmpty_Parallel(F2C_comm, *L, *N);
    if (F2C_CVODES_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vec = N_VNewEmpty_Parallel(F2C_comm, *L, *N);
    if (F2C_IDA_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vec = N_VNewEmpty_Parallel(F2C_comm, *L, *N);
    if (F2C_IDAS_vec == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_vec = N_VNewEmpty_Parallel(F2C_comm, *L, *N);
    if (F2C_KINSOL_vec == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITP_Q(MPI_Fint *comm, int *code, long int *Lq, long int *Nq, int *ier)
{
  MPI_Comm F2C_comm = MPI_Comm_f2c(*comm);

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecQ = N_VNewEmpty_Parallel(F2C_comm, *Lq, *Nq);
    if (F2C_CVODES_vecQ == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecQ = N_VNewEmpty_Parallel(F2C_comm, *Lq, *Nq);
    if (F2C_IDAS_vecQ == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITP_S(MPI_Fint *comm, int *code, int *Ns, long int *L, long int *N, int *ier)
{
  MPI_Comm F2C_comm = MPI_Comm_f2c(*comm);

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecS = N_VNewVectorArrayEmpty_Parallel(*Ns, F2C_comm, *L, *N);
    if (F2C_CVODES_vecS == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecS = N_VNewVectorArrayEmpty_Parallel(*Ns, F2C_comm, *L, *N);
    if (F2C_IDAS_vecS == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}


void FNV_INITP_B(MPI_Fint *comm, int *code, long int *LB, long int *NB, int *ier)
{
  MPI_Comm F2C_comm = MPI_Comm_f2c(*comm);

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecB = N_VNewEmpty_Parallel(F2C_comm, *LB, *NB);
    if (F2C_CVODES_vecB == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecB = N_VNewEmpty_Parallel(F2C_comm, *LB, *NB);
    if (F2C_IDAS_vecB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITP_QB(MPI_Fint *comm, int *code, long int *LqB, long int *NqB, int *ier)
{
  MPI_Comm F2C_comm = MPI_Comm_f2c(*comm);

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecQB = N_VNewEmpty_Parallel(F2C_comm, *LqB, *NqB);
    if (F2C_CVODES_vecQB == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecQB = N_VNewEmpty_Parallel(F2C_comm, *LqB, *NqB);
    if (F2C_IDAS_vecQB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

#else

void FNV_INITP(int *code, long int *L, long int *N, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vec = N_VNewEmpty_Parallel(MPI_COMM_WORLD, *L, *N);
    if (F2C_CVODE_vec == NULL) *ier = -1;
    break;
  case FCMIX_CVODES:
    F2C_CVODES_vec = N_VNewEmpty_Parallel(MPI_COMM_WORLD, *L, *N);
    if (F2C_CVODES_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vec = N_VNewEmpty_Parallel(MPI_COMM_WORLD, *L, *N);
    if (F2C_IDA_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vec = N_VNewEmpty_Parallel(MPI_COMM_WORLD, *L, *N);
    if (F2C_IDAS_vec == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_vec = N_VNewEmpty_Parallel(MPI_COMM_WORLD, *L, *N);
    if (F2C_KINSOL_vec == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITP_Q(int *code, long int *Lq, long int *Nq, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecQ = N_VNewEmpty_Parallel(MPI_COMM_WORLD, *Lq, *Nq);
    if (F2C_CVODES_vecQ == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecQ = N_VNewEmpty_Parallel(MPI_COMM_WORLD, *Lq, *Nq);
    if (F2C_IDAS_vecQ == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITP_S(int *code, int *Ns, long int *L, long int *N, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecS = N_VNewVectorArrayEmpty_Parallel(*Ns, MPI_COMM_WORLD, *L, *N);
    if (F2C_CVODES_vecS == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecS = N_VNewVectorArrayEmpty_Parallel(*Ns, MPI_COMM_WORLD, *L, *N);
    if (F2C_IDAS_vecS == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}


void FNV_INITP_B(int *code, long int *LB, long int *NB, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecB = N_VNewEmpty_Parallel(MPI_COMM_WORLD, *LB, *NB);
    if (F2C_CVODES_vecB == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecB = N_VNewEmpty_Parallel(MPI_COMM_WORLD, *LB, *NB);
    if (F2C_IDAS_vecB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITP_QB(int *code, long int *LqB, long int *NqB, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecQB = N_VNewEmpty_Parallel(MPI_COMM_WORLD, *LqB, *NqB);
    if (F2C_CVODES_vecQB == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecQB = N_VNewEmpty_Parallel(MPI_COMM_WORLD, *LqB, *NqB);
    if (F2C_IDAS_vecQB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

#endif
