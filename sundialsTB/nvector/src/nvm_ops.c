/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2012-03-07 21:41:19 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2013, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * Vector constructors for the SUNDIALS Matlab interfaces.
 * -----------------------------------------------------------------
 */

#include <stdlib.h>
#include "nvm.h"

void PutData(N_Vector v, double *data, long int n)
{
  double *vdata;
  long int i;

  vdata = N_VGetArrayPointer(v);

  for(i=0;i<n;i++) vdata[i] = data[i];

  return;
}


void GetData(N_Vector v, double *data, long int n)
{
  double *vdata;
  long int i;

  vdata = N_VGetArrayPointer(v);

  for(i=0;i<n;i++) data[i] = vdata[i];

  return;
}
