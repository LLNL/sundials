/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2012-03-07 21:41:19 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * begincopyright(llns)
 * Copyright (c) 2013, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * endcopyright(llns)
 * -----------------------------------------------------------------
 * Vector constructors for the SUNDIALS Matlab interfaces.
 * -----------------------------------------------------------------
 */

#include <stdlib.h>
#include "nvm.h"
#include <nvector/nvector_serial.h>

void InitVectors()
{}

N_Vector NewVector(long int n)
{
  N_Vector v;
  v = N_VNew_Serial(n);
  return(v);
}
