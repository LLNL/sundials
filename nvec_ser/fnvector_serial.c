/*
 * -----------------------------------------------------------------
 * $Revision: 1.9 $
 * $Date: 2005-04-15 00:43:59 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/shared/LICENSE.
 * -----------------------------------------------------------------
 * This file (companion of nvector_serial.h) contains the
 * implementation needed for the Fortran initialization of serial
 * vector operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fnvector_serial.h"
#include "nvector_serial.h"
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

void FNV_INITS(int *code, long int *N, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vec = N_VNewEmpty_Serial(*N);
    if (F2C_CVODE_vec == NULL) *ier = -1;
    break;
  case FCMIX_CVODES:
    F2C_CVODES_vec = N_VNewEmpty_Serial(*N);
    if (F2C_CVODES_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vec = N_VNewEmpty_Serial(*N);
    if (F2C_IDA_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vec = N_VNewEmpty_Serial(*N);
    if (F2C_IDAS_vec == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_vec = N_VNewEmpty_Serial(*N);
    if (F2C_KINSOL_vec == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITS_Q(int *code, long int *Nq, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecQ = N_VNewEmpty_Serial(*Nq);
    if (F2C_CVODES_vecQ == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecQ = N_VNewEmpty_Serial(*Nq);
    if (F2C_IDAS_vecQ == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITS_S(int *code, int *Ns, long int *N, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecS = N_VNewVectorArrayEmpty_Serial(*Ns, *N);
    if (F2C_CVODES_vecS == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecS = N_VNewVectorArrayEmpty_Serial(*Ns, *N);
    if (F2C_IDAS_vecS == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}


void FNV_INITS_B(int *code, long int *NB, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecB = N_VNewEmpty_Serial(*NB);
    if (F2C_CVODES_vecB == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecB = N_VNewEmpty_Serial(*NB);
    if (F2C_IDAS_vecB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITS_QB(int *code, long int *NqB, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODES:
    F2C_CVODES_vecQB = N_VNewEmpty_Serial(*NqB);
    if (F2C_CVODES_vecQB == NULL) *ier = -1;
    break;
  case FCMIX_IDAS:
    F2C_IDAS_vecQB = N_VNewEmpty_Serial(*NqB);
    if (F2C_IDAS_vecQB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}
