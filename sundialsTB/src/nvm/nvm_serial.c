/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-03-07 01:20:06 $
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

void InitVectors()
{}

N_Vector NewVector(int n)
{
  N_Vector v;
  v = N_VNew_Serial((long int)n);
  return(v);
}
