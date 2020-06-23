/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
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
 * Vector constructors for the SUNDIALS Matlab interfaces.
 * -----------------------------------------------------------------
 */

#include <stdlib.h>
#include "nvm.h"
#include <nvector/nvector_serial.h>

void InitVectors()
{}

N_Vector NewVector(sunindextype n)
{
  N_Vector v;
  v = N_VNew_Serial(n);
  return(v);
}
