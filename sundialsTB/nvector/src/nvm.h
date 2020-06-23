/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2012-03-07 21:41:19 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Header file for the MNVECTOR Matlab interface.
 * -----------------------------------------------------------------
 */

#ifndef _NVM_H
#define _NVM_H

#include <sundials/sundials_nvector.h>
#include "mex.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  /*
   * ------------------
   * Vector functions
   * ------------------
   */
  
  void InitVectors();
  N_Vector NewVector(sunindextype n);

  void PutData(N_Vector v, double *data, sunindextype n);
  void GetData(N_Vector v, double *data, sunindextype n);

#ifdef __cplusplus
}
#endif

#endif
