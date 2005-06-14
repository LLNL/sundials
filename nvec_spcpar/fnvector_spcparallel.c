/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2005-06-14 19:00:57 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/shared/LICENSE.
 * -----------------------------------------------------------------
 * This file (companion of nvector_sspcparallel.c) contains the 
 * implementation needed for the Fortran initialization of the 
 * spcparallel vector operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fnvector_spcparallel.h"
#include "nvector_spcparallel.h"
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

void FNV_INITSPCP(MPI_Fint *comm, int *code, 
                  int *nspc,
                  long int *Nx, long int *Ny, long int *Nz,
                  long int *NGx, long int *NGy, long int *NGz,
                  int *ier)
{
  MPI_Comm F2C_comm = MPI_Comm_f2c(*comm);

  *ier = 0;

  /* Create the global vector for the desired SUNDIALS solver */

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vec = N_VNewEmpty_SpcParallel(F2C_comm,*nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_CVODE_vec == NULL) *ier = -1;
    break;
  case FCMIX_CVODES:
    F2C_CVODES_vec = N_VNewEmpty_SpcParallel(F2C_comm,*nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_CVODES_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vec = N_VNewEmpty_SpcParallel(F2C_comm,*nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_IDA_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vec = N_VNewEmpty_SpcParallel(F2C_comm,*nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_IDAS_vec == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_vec = N_VNewEmpty_SpcParallel(F2C_comm,*nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_KINSOL_vec == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_Q(MPI_Fint *comm, int *code, 
                    int *nspcQ,
                    long int *NxQ, long int *NyQ, long int *NzQ,
                    long int *NGxQ, long int *NGyQ, long int *NGzQ,
                    int *ier)
{
  MPI_Comm F2C_comm = MPI_Comm_f2c(*comm);

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecQ = N_VNewEmpty_SpcParallel(F2C_comm,*nspcQ,*NxQ,*NyQ,*NzQ,*NGxQ,*NGyQ,*NGzQ);
    if (F2C_CVODES_vecQ == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecQ = N_VNewEmpty_SpcParallel(F2C_comm,*nspcQ,*NxQ,*NyQ,*NzQ,*NGxQ,*NGyQ,*NGzQ);
    if (F2C_IDAS_vecQ == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_S(MPI_Fint *comm, int *code, 
                    int *Ns, 
                    int *nspc,
                    long int *Nx, long int *Ny, long int *Nz,
                    long int *NGx, long int *NGy, long int *NGz,
                    int *ier)
{
  MPI_Comm F2C_comm = MPI_Comm_f2c(*comm);

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecS = N_VNewVectorArrayEmpty_SpcParallel(*Ns, F2C_comm,*nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_CVODES_vecS == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecS = N_VNewVectorArrayEmpty_SpcParallel(*Ns, F2C_comm,*nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_IDAS_vecS == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_B(MPI_Fint *comm, int *code, 
                    int *nspcB,
                    long int *NxB, long int *NyB, long int *NzB,
                    long int *NGxB, long int *NGyB, long int *NGzB,
                    int *ier)
{
  MPI_Comm F2C_comm = MPI_Comm_f2c(*comm);

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecB = N_VNewEmpty_SpcParallel(F2C_comm,*nspcB,*NxB,*NyB,*NzB,*NGxB,*NGyB,*NGzB);
    if (F2C_CVODES_vecB == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecB = N_VNewEmpty_SpcParallel(F2C_comm,*nspcB,*NxB,*NyB,*NzB,*NGxB,*NGyB,*NGzB);
    if (F2C_IDAS_vecB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_QB(MPI_Fint *comm, int *code, 
                     int *nspcQB,
                     long int *NxQB, long int *NyQB, long int *NzQB,
                     long int *NGxQB, long int *NGyQB, long int *NGzQB,
                     int *ier)
{
  MPI_Comm F2C_comm = MPI_Comm_f2c(*comm);

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecQB = N_VNewEmpty_SpcParallel(F2C_comm,*nspcQB,*NxQB,*NyQB,*NzQB,*NGxQB,*NGyQB,*NGzQB);
    if (F2C_CVODES_vecQB == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecQB = N_VNewEmpty_SpcParallel(F2C_comm,*nspcQB,*NxQB,*NyQB,*NzQB,*NGxQB,*NGyQB,*NGzQB);
    if (F2C_IDAS_vecQB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

#else

void FNV_INITSPCP(int *code, 
                  int *nspc,
                  long int *Nx, long int *Ny, long int *Nz,
                  long int *NGx, long int *NGy, long int *NGz,
                  int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vec = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_CVODE_vec == NULL) *ier = -1;
    break;
  case FCMIX_CVODES:
    F2C_CVODES_vec = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_CVODES_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vec = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_IDA_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vec = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_IDAS_vec == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_vec = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_KINSOL_vec == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_Q(int *code, 
                    int *nspcQ,
                    long int *NxQ, long int *NyQ, long int *NzQ,
                    long int *NGxQ, long int *NGyQ, long int *NGzQ,
                    int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecQ = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*nspcQ,*NxQ,*NyQ,*NzQ,*NGxQ,*NGyQ,*NGzQ);
    if (F2C_CVODES_vecQ == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecQ = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*nspcQ,*NxQ,*NyQ,*NzQ,*NGxQ,*NGyQ,*NGzQ);
    if (F2C_IDAS_vecQ == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_S(int *code, 
                    int *Ns, 
                    int *nspc,
                    long int *Nx, long int *Ny, long int *Nz,
                    long int *NGx, long int *NGy, long int *NGz,
                    int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecS = N_VNewVectorArrayEmpty_SpcParallel(*Ns, MPI_COMM_WORLD,*nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_CVODES_vecS == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecS = N_VNewVectorArrayEmpty_SpcParallel(*Ns, MPI_COMM_WORLD,*nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_IDAS_vecS == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_B(int *code, 
                    int *nspcB,
                    long int *NxB, long int *NyB, long int *NzB,
                    long int *NGxB, long int *NGyB, long int *NGzB,
                    int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecB = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*nspcB,*NxB,*NyB,*NzB,*NGxB,*NGyB,*NGzB);
    if (F2C_CVODES_vecB == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecB = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*nspcB,*NxB,*NyB,*NzB,*NGxB,*NGyB,*NGzB);
    if (F2C_IDAS_vecB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_QB(int *code, 
                     int *nspcQB,
                     long int *NxQB, long int *NyQB, long int *NzQB,
                     long int *NGxQB, long int *NGyQB, long int *NGzQB,
                     int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecQB = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*nspcQB,*NxQB,*NyQB,*NzQB,*NGxQB,*NGyQB,*NGzQB);
    if (F2C_CVODES_vecQB == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecQB = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*nspcQB,*NxQB,*NyQB,*NzQB,*NGxQB,*NGyQB,*NGzQB);
    if (F2C_IDAS_vecQB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

#endif
