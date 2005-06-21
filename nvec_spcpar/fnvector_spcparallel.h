/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2005-06-21 19:22:51 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/shared/LICENSE.
 * -----------------------------------------------------------------
 * This file (companion of nvector_spcparallel.h) contains the
 * definitions needed for the initialization of the spcparallel
 * vector operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FNVECTOR_SPCPARALLEL_H
#define _FNVECTOR_SPCPARALLEL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "mpi.h"
#include "nvector.h"  
#include "fnvector.h"

#ifndef _SUNDIALS_CONFIG_H
#define _SUNDIALS_CONFIG_H
#include "sundials_config.h"
#endif

#if defined(F77_FUNC)

#define FNV_INITSPCP    F77_FUNC(fnvinitspcp, FNVINITSPCP)
#define FNV_INITSPCP_Q  F77_FUNC_(fnvinitspcp_q, FNVINITSPCP_Q)
#define FNV_INITSPCP_S  F77_FUNC_(fnvinitspcp_s, FNVINITSPCP_S)
#define FNV_INITSPCP_B  F77_FUNC_(fnvinitspcp_b, FNVINITSPCP_B)
#define FNV_INITSPCP_QB F77_FUNC_(fnvinitspcp_qb, FNVINITSPCP_QB)

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)

#define FNV_INITSPCP    fnvinitspcp
#define FNV_INITSPCP_Q  fnvinitspcp_q
#define FNV_INITSPCP_S  fnvinitspcp_s
#define FNV_INITSPCP_B  fnvinitspcp_b
#define FNV_INITSPCP_QB fnvinitspcp_qb

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_UPPER)

#define FNV_INITSPCP    FNVINITSPCP
#define FNV_INITSPCP_Q  FNVINITSPCP_Q
#define FNV_INITSPCP_S  FNVINITSPCP_S
#define FNV_INITSPCP_B  FNVINITSPCP_B
#define FNV_INITSPCP_QB FNVINITSPCP_QB

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_LOWER)

#define FNV_INITSPCP    fnvinitspcp_
#define FNV_INITSPCP_Q  fnvinitspcp_q_
#define FNV_INITSPCP_S  fnvinitspcp_s_
#define FNV_INITSPCP_B  fnvinitspcp_b_
#define FNV_INITSPCP_QB fnvinitspcp_qb_

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_UPPER)

#define FNV_INITSPCP    FNVINITSPCP_
#define FNV_INITSPCP_Q  FNVINITSPCP_Q_
#define FNV_INITSPCP_S  FNVINITSPCP_S_
#define FNV_INITSPCP_B  FNVINITSPCP_B_
#define FNV_INITSPCP_QB FNVINITSPCP_QB_

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_LOWER)

#define FNV_INITSPCP    fnvinitspcp__
#define FNV_INITSPCP_Q  fnvinitspcp_q__
#define FNV_INITSPCP_S  fnvinitspcp_s__
#define FNV_INITSPCP_B  fnvinitspcp_b__
#define FNV_INITSPCP_QB fnvinitspcp_qb__

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_UPPER)

#define FNV_INITSPCP    FNVINITSPCP__
#define FNV_INITSPCP_Q  FNVINITSPCP_Q__
#define FNV_INITSPCP_S  FNVINITSPCP_S__
#define FNV_INITSPCP_B  FNVINITSPCP_B__
#define FNV_INITSPCP_QB FNVINITSPCP_QB__

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

#ifdef SUNDIALS_MPI_COMM_F2C

  void FNV_INITSPCP(MPI_Fint *comm, int *code, 
                    int *Ngrp, int *Nspc,
                    long int *Nx, long int *Ny, long int *Nz,
                    long int *NGx, long int *NGy, long int *NGz,
                    int *ier);
  void FNV_INITSPCP_Q(MPI_Fint *comm, int *code, 
                      int *NgrpQ, int *NspcQ,
                      long int *NxQ, long int *NyQ, long int *NzQ,
                      long int *NGxQ, long int *NGyQ, long int *NGzQ,
                      int *ier);
  void FNV_INITSPCP_B(MPI_Fint *comm, int *code, 
                      int *NgrpB, int *NspcB,
                      long int *NxB, long int *NyB, long int *NzB,
                      long int *NGxB, long int *NGyB, long int *NGzB,
                      int *ier);
  void FNV_INITSPCP_QB(MPI_Fint *comm, int *code, 
                       int *NgrpQB, int *NspcQB,
                       long int *NxQB, long int *NyQB, long int *NzQB,
                       long int *NGxQB, long int *NGyQB, long int *NGzQB,
                       int *ier);

#else

  void FNV_INITSPCP(int *code, 
                    int *Ngrp, int *Nspc,
                    long int *Nx, long int *Ny, long int *Nz,
                    long int *NGx, long int *NGy, long int *NGz,
                    int *ier);
  void FNV_INITSPCP_Q(int *code, 
                      int *NgrpQ, int *NspcQ,
                      long int *NxQ, long int *NyQ, long int *NzQ,
                      long int *NGxQ, long int *NGyQ, long int *NGzQ,
                      int *ier);
  void FNV_INITSPCP_B(int *code, 
                      int *NgrpB, int *NspcB,
                      long int *NxB, long int *NyB, long int *NzB,
                      long int *NGxB, long int *NGyB, long int *NGzB,
                      int *ier);
  void FNV_INITSPCP_QB(int *code, 
                       int *NgrpQB, int *NspcQB,
                       long int *NxQB, long int *NyQB, long int *NzQB,
                       long int *NGxQB, long int *NGyQB, long int *NGzQB,
                       int *ier);

#endif

  void FNV_INITSPCP_S(int *code, int *Ns, int *ier);

#ifdef __cplusplus
}
#endif

#endif
