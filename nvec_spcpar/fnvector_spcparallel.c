/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2005-06-21 19:22:51 $
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

/* Fortran callable interfaces */

#ifdef SUNDIALS_MPI_COMM_F2C

void FNV_INITSPCP(MPI_Fint *comm, int *code, 
                  int *Ngrp, int *Nspc,
                  long int *Nx, long int *Ny, long int *Nz,
                  long int *NGx, long int *NGy, long int *NGz,
                  int *ier)
{
  MPI_Comm F2C_comm = MPI_Comm_f2c(*comm);

  *ier = 0;

  /* Create the global vector for the desired SUNDIALS solver */

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vec = N_VNewEmpty_SpcParallel(F2C_comm,*Ngrp,Nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_CVODE_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vec = N_VNewEmpty_SpcParallel(F2C_comm,*Ngrp,Nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_IDA_vec == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_vec = N_VNewEmpty_SpcParallel(F2C_comm,*Ngrp,Nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_KINSOL_vec == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_Q(MPI_Fint *comm, int *code, 
                    int *NgrpQ, int *NspcQ,
                    long int *NxQ, long int *NyQ, long int *NzQ,
                    long int *NGxQ, long int *NGyQ, long int *NGzQ,
                    int *ier)
{
  MPI_Comm F2C_comm = MPI_Comm_f2c(*comm);

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecQ = N_VNewEmpty_SpcParallel(F2C_comm,*NgrpQ,NspcQ,*NxQ,*NyQ,*NzQ,*NGxQ,*NGyQ,*NGzQ);
    if (F2C_CVODE_vecQ == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecQ = N_VNewEmpty_SpcParallel(F2C_comm,*NgrpQ,NspcQ,*NxQ,*NyQ,*NzQ,*NGxQ,*NGyQ,*NGzQ);
    if (F2C_IDA_vecQ == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_B(MPI_Fint *comm, int *code, 
                    int *NgrpB, int *NspcB,
                    long int *NxB, long int *NyB, long int *NzB,
                    long int *NGxB, long int *NGyB, long int *NGzB,
                    int *ier)
{
  MPI_Comm F2C_comm = MPI_Comm_f2c(*comm);

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecB = N_VNewEmpty_SpcParallel(F2C_comm,*NgrpB,NspcB,*NxB,*NyB,*NzB,*NGxB,*NGyB,*NGzB);
    if (F2C_CVODE_vecB == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecB = N_VNewEmpty_SpcParallel(F2C_comm,*NgrpB,NspcB,*NxB,*NyB,*NzB,*NGxB,*NGyB,*NGzB);
    if (F2C_IDA_vecB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_QB(MPI_Fint *comm, int *code, 
                     int *NgrpQB, int *NspcQB,
                     long int *NxQB, long int *NyQB, long int *NzQB,
                     long int *NGxQB, long int *NGyQB, long int *NGzQB,
                     int *ier)
{
  MPI_Comm F2C_comm = MPI_Comm_f2c(*comm);

  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecQB = N_VNewEmpty_SpcParallel(F2C_comm,*NgrpQB,NspcQB,*NxQB,*NyQB,*NzQB,*NGxQB,*NGyQB,*NGzQB);
    if (F2C_CVODE_vecQB == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecQB = N_VNewEmpty_SpcParallel(F2C_comm,*NgrpQB,NspcQB,*NxQB,*NyQB,*NzQB,*NGxQB,*NGyQB,*NGzQB);
    if (F2C_IDA_vecQB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

#else

void FNV_INITSPCP(int *code, 
                  int *Ngrp, int *Nspc,
                  long int *Nx, long int *Ny, long int *Nz,
                  long int *NGx, long int *NGy, long int *NGz,
                  int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vec = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*Ngrp,Nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_CVODE_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vec = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*Ngrp,Nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_IDA_vec == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_vec = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*Ngrp,Nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_KINSOL_vec == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_Q(int *code, 
                    int *NgrpQ, int *NspcQ,
                    long int *NxQ, long int *NyQ, long int *NzQ,
                    long int *NGxQ, long int *NGyQ, long int *NGzQ,
                    int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecQ = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*NgrpQ,NspcQ,*NxQ,*NyQ,*NzQ,*NGxQ,*NGyQ,*NGzQ);
    if (F2C_CVODE_vecQ == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecQ = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*NgrpQ,NspcQ,*NxQ,*NyQ,*NzQ,*NGxQ,*NGyQ,*NGzQ);
    if (F2C_IDA_vecQ == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_B(int *code, 
                    int *NgrpB, int *NspcB,
                    long int *NxB, long int *NyB, long int *NzB,
                    long int *NGxB, long int *NGyB, long int *NGzB,
                    int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecB = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*NgrpB,NspcB,*NxB,*NyB,*NzB,*NGxB,*NGyB,*NGzB);
    if (F2C_CVODE_vecB == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecB = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*NgrpB,NspcB,*NxB,*NyB,*NzB,*NGxB,*NGyB,*NGzB);
    if (F2C_IDA_vecB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_QB(int *code, 
                     int *NgrpQB, int *NspcQB,
                     long int *NxQB, long int *NyQB, long int *NzQB,
                     long int *NGxQB, long int *NGyQB, long int *NGzQB,
                     int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecQB = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*NgrpQB,NspcQB,*NxQB,*NyQB,*NzQB,*NGxQB,*NGyQB,*NGzQB);
    if (F2C_CVODE_vecQB == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecQB = N_VNewEmpty_SpcParallel(MPI_COMM_WORLD,*NgrpQB,NspcQB,*NxQB,*NyQB,*NzQB,*NGxQB,*NGyQB,*NGzQB);
    if (F2C_IDA_vecQB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

#endif

void FNV_INITSPCP_S(int *code, int *Ns, int *ier) 
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecS = N_VCloneVectorArrayEmpty_SpcParallel(*Ns, F2C_CVODE_vec);
    if (F2C_CVODE_vecS == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecS = N_VCloneVectorArrayEmpty_SpcParallel(*Ns, F2C_IDA_vec);
    if (F2C_IDA_vecS == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

