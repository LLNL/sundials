/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-01-06 19:00:26 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Header file for the MNVECTOR Matlab interface.
 * -----------------------------------------------------------------
 */

#ifndef _NVM_H
#define _NVM_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "nvector.h"
#include "mex.h"

  /*
   * ------------------
   * Vector functions
   * ------------------
   */
  
  void InitVectors(int vec_type, mxArray *mx_comm);
  N_Vector NewVector(int n);

  void PutData(N_Vector v, double *data, int n);
  void GetData(N_Vector v, double *data, int n);

#ifdef __cplusplus
}
#endif

#endif
