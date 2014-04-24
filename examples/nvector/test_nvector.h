/*
 * -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date: 2013/08/31 12:25:00 $
 * ----------------------------------------------------------------- 
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the header file contains the prototypes for functions to 
 * test an NVECTOR module implementation. 
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  int Test_N_VCloneVectorArray(int count, N_Vector W, long int local_length, int myid);
  int Test_N_VCloneEmptyVectorArray(int count, N_Vector W, int myid);
  int Test_N_VCloneEmpty(N_Vector W, int myid);
  int Test_N_VClone(N_Vector W, long int local_length, int myid);
  int Test_N_VGetArrayPointer(N_Vector W, long int local_length, int myid);
  int Test_N_VSetArrayPointer(N_Vector W, long int local_length, int myid);
  int Test_N_VLinearSum(N_Vector X, N_Vector Y, N_Vector Z, 
			long int local_length, int myid);
  int Test_N_VConst(N_Vector X, long int local_length, int myid);
  int Test_N_VProd(N_Vector X, N_Vector Y, N_Vector Z, long int local_length, int myid);
  int Test_N_VDiv(N_Vector X, N_Vector Y, N_Vector Z, long int local_length, int myid);
  int Test_N_VScale(N_Vector X, N_Vector Z, long int local_length, int myid);
  int Test_N_VAbs(N_Vector X, N_Vector Z, long int local_length, int myid);
  int Test_N_VInv(N_Vector X, N_Vector Z, long int local_length, int myid);
  int Test_N_VAddConst(N_Vector X, N_Vector Z, long int local_length, int myid);
  int Test_N_VDotProd(N_Vector X, N_Vector Y, 
		      long int local_length, long int global_length, int myid);
  int Test_N_VMaxNorm(N_Vector X, long int local_length, int myid); 
  int Test_N_VWrmsNorm(N_Vector X, N_Vector W, long int local_length, int myid);
  int Test_N_VWrmsNormMask(N_Vector X, N_Vector W, N_Vector ID, 
			   long int local_length, long int global_length, int myid);
  int Test_N_VMin(N_Vector X, long int local_length, int myid); 
  int Test_N_VWL2Norm(N_Vector X, N_Vector W, 
		      long int local_length, long int global_length, int myid);
  int Test_N_VL1Norm(N_Vector X, long int local_length, long int global_length, int myid);
  int Test_N_VCompare(N_Vector X, N_Vector Z, long int local_length, int myid);
  int Test_N_VInvTest(N_Vector X, N_Vector Z, long int local_length, int myid);
  int Test_N_VConstrMask(N_Vector C, N_Vector X, N_Vector M, 
			 long int local_length, int myid);
  int Test_N_VMinQuotient(N_Vector NUM, N_Vector DENOM, long int local_length, int myid);

  void SetTiming(int onoff);
  
#ifdef __cplusplus
}
#endif
