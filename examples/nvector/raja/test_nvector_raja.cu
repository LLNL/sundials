/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
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
 * This is the testing routine to check the NVECTOR Raja module
 * implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_raja.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

#include <nvector/raja/Vector.hpp>

// Type definitions
typedef sunrajavec::Vector<realtype, sunindextype> vector_type;


/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int           fails = 0;    /* counter for test failures  */
  N_Vector      W, X, Y, Z;   /* test vectors               */
  SUNMPI_Comm comm;
  int           nprocs, myid; /* Number of procs, proc id  */
  const int     print_timing = 0;
  sunindextype  global_length, local_length;
  /*  sunindextype liw, lrw; */


//   /* check input and set vector length */
//   if (argc < 3){
//     printf("ERROR: TWO (2) Inputs required: vector length, print timing \n");
//     return(-1);
//   }
//
//   local_length = atol(argv[1]);
//   if (local_length <= 0) {
//     printf("ERROR: length of vector must be a positive integer \n");
//     return(-1);
//   }
//
//   print_timing = atoi(argv[2]);

  SetTiming(print_timing);

  local_length = 1000;

#ifdef SUNDIALS_MPI_ENABLED
  /* Get processor number and total number of processes */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);
#else
  comm = 0;
  myid = 0;
  nprocs = 1;
#endif
  global_length = nprocs*local_length;

  if (myid == 0)
    printf("\nRunning with vector length %ld \n\n", (long) global_length);

  /* Create vectors */
  W = N_VNewEmpty_Raja(local_length);
  X = N_VNew_Raja(comm, local_length, global_length);

  /* NVector Tests */

  /* Raja specific tests */

  /* Memory allocation tests */
  fails += Test_N_VCloneEmpty(X, myid);
  fails += Test_N_VClone(X, local_length, myid);
  fails += Test_N_VCloneEmptyVectorArray(5, X, myid);
  fails += Test_N_VCloneVectorArray(5, X, local_length, myid);

  Y = N_VClone_Raja(X);
  Z = N_VClone_Raja(X);

  /* Skipped tests */
  /*   fails += Test_N_VSetArrayPointer(W, local_length, myid); */
  /*   fails += Test_N_VGetArrayPointer(X, local_length, myid); */

  /* Vector operation tests */
  fails += Test_N_VConst(X, local_length, myid);
  fails += Test_N_VLinearSum(X, Y, Z, local_length, myid);
  fails += Test_N_VProd(X, Y, Z, local_length, myid);
  fails += Test_N_VDiv(X, Y, Z, local_length, myid);
  fails += Test_N_VScale(X, Z, local_length, myid);
  fails += Test_N_VAbs(X, Z, local_length, myid);
  fails += Test_N_VInv(X, Z, local_length, myid);
  fails += Test_N_VAddConst(X, Z, local_length, myid);
  fails += Test_N_VDotProd(X, Y, local_length, global_length, myid);
  fails += Test_N_VMaxNorm(X, local_length, myid);
  fails += Test_N_VWrmsNorm(X, Y, local_length, myid);
  fails += Test_N_VWrmsNormMask(X, Y, Z, local_length, global_length, myid);
  fails += Test_N_VMin(X, local_length, myid);
  fails += Test_N_VWL2Norm(X, Y, local_length, global_length, myid);
  fails += Test_N_VL1Norm(X, local_length, global_length, myid);
  fails += Test_N_VCompare(X, Z, local_length, myid);
  fails += Test_N_VInvTest(X, Z, local_length, myid);
  fails += Test_N_VConstrMask(X, Y, Z, local_length, myid);
  fails += Test_N_VMinQuotient(X, Y, local_length, myid);

  /* Fused vector operation tests (optional) */
  fails += Test_N_VLinearCombination(X, local_length, myid);
  fails += Test_N_VScaleAddMulti(X, local_length, myid);
  fails += Test_N_VDotProdMulti(X, local_length, global_length, myid);

  /* Vector array operation tests (optional) */
  fails += Test_N_VLinearSumVectorArray(X, local_length, myid);
  fails += Test_N_VScaleVectorArray(X, local_length, myid);
  fails += Test_N_VConstVectorArray(X, local_length, myid);
  fails += Test_N_VWrmsNormVectorArray(X, local_length, myid);
  fails += Test_N_VWrmsNormMaskVectorArray(X, local_length, global_length, myid);
  fails += Test_N_VScaleAddMultiVectorArray(X, local_length, myid);
  fails += Test_N_VLinearCombinationVectorArray(X, local_length, myid);

//   sunindextype lrw, liw;
//   N_VSpace_Raja(X, &lrw, &liw);
//   if (myid == 0)
//     printf("lrw = %ld, liw = %ld\n", lrw, liw);

  /* Free vectors */
  N_VDestroy(W);
  N_VDestroy(X);
  N_VDestroy(Y);
  N_VDestroy(Z);

  /* Print result */
  if (myid == 0) {
    if (fails) {
      printf("FAIL: NVector module failed %i tests \n \n", fails);
    } else {
      printf("SUCCESS: NVector module passed all tests \n \n");
    }
  }

#ifdef SUNDIALS_MPI_ENABLED
  MPI_Finalize();
#endif


  return(fails);
}

/* ----------------------------------------------------------------------
 * Check vector
 * --------------------------------------------------------------------*/
int check_ans(realtype ans, N_Vector X, sunindextype local_length)
{
  int      failure = 0;
  sunindextype i;
  vector_type* xv = sunrajavec::extract<realtype, sunindextype>(X);
  realtype *xdata;

  xv->copyFromDev();

  xdata = xv->host();
  /* check vector data */
  for(i=0; i < local_length; i++){
    failure += FNEQ(xdata[i], ans);
    //printf("%g ?= %g\n", xdata[i], ans);
  }
  return (failure > ZERO) ? (1) : (0);
}

booleantype has_data(N_Vector X)
{
  //sunrajavec::Vector<double, sunindextype>* xv = sunrajavec::extract<realtype, sunindextype>(X);
  auto xv = sunrajavec::extract<realtype, sunindextype>(X);

  return (xv == NULL ? SUNFALSE : SUNTRUE);
}

void set_element(N_Vector X, sunindextype i, realtype val)
{
  vector_type* xv = sunrajavec::extract<realtype, sunindextype>(X);
  xv->copyFromDev();
  (xv->host())[i] = val;
  xv->copyToDev();
}

realtype get_element(N_Vector X, sunindextype i)
{
  vector_type* xv = sunrajavec::extract<realtype, sunindextype>(X);
  xv->copyFromDev();
  return (xv->host())[i];
}
