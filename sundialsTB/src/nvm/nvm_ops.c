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

void PutData(N_Vector v, double *data, int n)
{
  double *vdata;
  int i;

  vdata = N_VGetArrayPointer(v);

  for(i=0;i<n;i++) vdata[i] = data[i];

  return;
}


void GetData(N_Vector v, double *data, int n)
{
  double *vdata;
  int i;

  vdata = N_VGetArrayPointer(v);

  for(i=0;i<n;i++) data[i] = vdata[i];

  return;
}
