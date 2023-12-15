/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the CVBBDPRE module, for a
 * band-block-diagonal preconditioner, i.e. a block-diagonal
 * matrix with banded blocks.
 * -----------------------------------------------------------------*/

#ifndef _CVSBBDPRE_H
#define _CVSBBDPRE_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*-----------------
  FORWARD PROBLEMS
  -----------------*/

/* User-supplied function Types */

typedef int (*CVLocalFn)(sunindextype Nlocal, sunrealtype t, N_Vector y,
                         N_Vector g, void* user_data);

typedef int (*CVCommFn)(sunindextype Nlocal, sunrealtype t, N_Vector y,
                        void* user_data);

/* Exported Functions */

SUNDIALS_EXPORT int CVBBDPrecInit(void* cvode_mem, sunindextype Nlocal,
                                  sunindextype mudq, sunindextype mldq,
                                  sunindextype mukeep, sunindextype mlkeep,
                                  sunrealtype dqrely, CVLocalFn gloc,
                                  CVCommFn cfn);

SUNDIALS_EXPORT int CVBBDPrecReInit(void* cvode_mem, sunindextype mudq,
                                    sunindextype mldq, sunrealtype dqrely);

/* Optional output functions */

SUNDIALS_EXPORT int CVBBDPrecGetWorkSpace(void* cvode_mem, long int* lenrwBBDP,
                                          long int* leniwBBDP);

SUNDIALS_EXPORT int CVBBDPrecGetNumGfnEvals(void* cvode_mem,
                                            long int* ngevalsBBDP);

/*------------------
  BACKWARD PROBLEMS
  ------------------*/

/* User-Supplied Function Types */

typedef int (*CVLocalFnB)(sunindextype NlocalB, sunrealtype t, N_Vector y,
                          N_Vector yB, N_Vector gB, void* user_dataB);

typedef int (*CVCommFnB)(sunindextype NlocalB, sunrealtype t, N_Vector y,
                         N_Vector yB, void* user_dataB);

/* Exported Functions */

SUNDIALS_EXPORT int CVBBDPrecInitB(void* cvode_mem, int which,
                                   sunindextype NlocalB, sunindextype mudqB,
                                   sunindextype mldqB, sunindextype mukeepB,
                                   sunindextype mlkeepB, sunrealtype dqrelyB,
                                   CVLocalFnB glocB, CVCommFnB cfnB);

SUNDIALS_EXPORT int CVBBDPrecReInitB(void* cvode_mem, int which,
                                     sunindextype mudqB, sunindextype mldqB,
                                     sunrealtype dqrelyB);

#ifdef __cplusplus
}
#endif

#endif
