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
 * Vector constructors for the SUNDIALS Matlab interfaces.
 * -----------------------------------------------------------------
 */

#include <stdlib.h>
#include "nvm.h"
#include "nvector_serial.h"

void InitVectors(int vec_type, mxArray *mx_comm)
{
  if (vec_type != 1)
    mexPrintf("\n\nNO parallel vector support!!\n\n");
}


N_Vector NewVector(int n)
{
  N_Vector v;

  v = N_VNewEmpty_Serial((long int)n);
  return(v);

}
