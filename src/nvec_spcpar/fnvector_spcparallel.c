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
 * This file (companion of nvector_sspcparallel.h) contains the 
 * implementation needed for the Fortran initialization of the 
 * spcparallel vector operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fnvector_spcparallel.h"

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

void FNV_INITSPCP(MPI_Fint *comm, int *code, 
                  int *Ngrp, int *Nspc,
                  long int *Nx, long int *Ny, long int *Nz,
                  long int *NGx, long int *NGy, long int *NGz,
                  int *ier)
{
  MPI_Comm F2C_comm;

#ifdef SUNDIALS_MPI_COMM_F2C
  F2C_comm = MPI_Comm_f2c(*comm);
#else
  F2C_comm = MPI_COMM_WORLD;
#endif

  *ier = 0;

  /* Create the global vector for the desired SUNDIALS solver */

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vec = NULL;
    F2C_CVODE_vec = N_VNewEmpty_SpcParallel(F2C_comm,*Ngrp,Nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_CVODE_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vec = NULL;
    F2C_IDA_vec = N_VNewEmpty_SpcParallel(F2C_comm,*Ngrp,Nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_IDA_vec == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_vec = NULL;
    F2C_KINSOL_vec = N_VNewEmpty_SpcParallel(F2C_comm,*Ngrp,Nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_KINSOL_vec == NULL) *ier = -1;
    break;
  case FCMIX_ARKODE:
    F2C_ARKODE_vec = NULL;
    F2C_ARKODE_vec = N_VNewEmpty_SpcParallel(F2C_comm,*Ngrp,Nspc,*Nx,*Ny,*Nz,*NGx,*NGy,*NGz);
    if (F2C_ARKODE_vec == NULL) *ier = -1;
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
    F2C_CVODE_vecQ = N_VNewEmpty_SpcParallel(F2C_comm,*NgrpQ,NspcQ,*NxQ,*NyQ,*NzQ,*NGxQ,*NGyQ,*NGzQ);
    if (F2C_CVODE_vecQ == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecQ = NULL;
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
    F2C_CVODE_vecB = N_VNewEmpty_SpcParallel(F2C_comm,*NgrpB,NspcB,*NxB,*NyB,*NzB,*NGxB,*NGyB,*NGzB);
    if (F2C_CVODE_vecB == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecB = NULL;
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
    F2C_CVODE_vecQB = N_VNewEmpty_SpcParallel(F2C_comm,*NgrpQB,NspcQB,*NxQB,*NyQB,*NzQB,*NGxQB,*NGyQB,*NGzQB);
    if (F2C_CVODE_vecQB == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecQB = NULL;
    F2C_IDA_vecQB = N_VNewEmpty_SpcParallel(F2C_comm,*NgrpQB,NspcQB,*NxQB,*NyQB,*NzQB,*NGxQB,*NGyQB,*NGzQB);
    if (F2C_IDA_vecQB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITSPCP_S(int *code, int *Ns, int *ier) 
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecS = NULL;
    F2C_CVODE_vecS = N_VCloneVectorArrayEmpty_SpcParallel(*Ns, F2C_CVODE_vec);
    if (F2C_CVODE_vecS == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecS = NULL;
    F2C_IDA_vecS = N_VCloneVectorArrayEmpty_SpcParallel(*Ns, F2C_IDA_vec);
    if (F2C_IDA_vecS == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}
