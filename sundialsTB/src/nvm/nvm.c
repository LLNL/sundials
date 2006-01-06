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
 * MEX-binding for setting serial/parallel vectors
 * -----------------------------------------------------------------
 */

#include <stdlib.h>
#include "mex.h"

/* Global variables */

int VecType = 1;   /* 1: serial, 2:parallel */
mxArray *mx_comm;  /* MPI communicator      */

/*
 * ---------------------------------------------------------------------------------
 * Main entry point
 * ---------------------------------------------------------------------------------
 */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] )
{
  int mode;
  /* 
   * Modes:
   * 1 - specify parallel vector type
   *     Call: nvm(1,comm)
   * 2 - clear vector type specification
   *     Call: nvm(2)
   * 3 - querry vector type
   *     Call: [vec_type, comm] = nvm(3)
   *
   * nvm is caled with mode=1 and mode=2 by mpirun and mpiruns
   * nvm is called with mode=3 by the initialization functions of
   *     the individual solvers
  */

  mode = (int)mxGetScalar(prhs[0]);
  switch(mode) {
  case 1:
    VecType = 2;
    mx_comm = mxDuplicateArray(prhs[1]);
    mexLock();
    mexMakeArrayPersistent(mx_comm);
    break;
  case 2:
    mexUnlock();
    mxDestroyArray(mx_comm);
    VecType = 1;
    break;
  case 3:
    plhs[0] = mxCreateScalarDouble((double)VecType);
    if (VecType == 1)
      plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL);
    else
      plhs[1] = mxDuplicateArray(mx_comm);
    break;
  }

}

