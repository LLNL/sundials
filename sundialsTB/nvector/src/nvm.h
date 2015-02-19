/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2012-03-07 21:41:19 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
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
  N_Vector NewVector(long int n);

  void PutData(N_Vector v, double *data, long int n);
  void GetData(N_Vector v, double *data, long int n);

#ifdef __cplusplus
}
#endif

#endif
