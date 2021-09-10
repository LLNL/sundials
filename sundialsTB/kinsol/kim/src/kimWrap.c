/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * KINSOL wrapper functions.
 * -----------------------------------------------------------------
 */

#include "kim.h"
#include "nvm.h"

static void UpdateUserData(mxArray *new_mtlb_data, kimInterfaceData kimData);

/*
 * ---------------------------------------------------------------------------------
 * Redability replacements
 * ---------------------------------------------------------------------------------
 */

#define N  (kimData->n)
#define ls (kimData->LS)
#define pm (kimData->PM)

/*
 * ---------------------------------------------------------------------------------
 * Wrapper functions
 * ---------------------------------------------------------------------------------
 */

int mxW_KINSys(N_Vector y, N_Vector fy, void *user_data )
{
  kimInterfaceData kimData;
  mxArray *mx_in[3], *mx_out[3];
  int ret;
  
  /* Extract global interface data from user-data */
  kimData = (kimInterfaceData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[1] = kimData->SYSfct;                  /* matlab function handle */ 
  mx_in[2] = kimData->mtlb_data;               /* matlab user data */

  /* Call matlab wrapper */

  GetData(y, mxGetPr(mx_in[0]), N);

  mexCallMATLAB(3,mx_out,3,mx_in,"kim_sys");

  PutData(fy, mxGetPr(mx_out[0]), N);
  ret = (int)*mxGetPr(mx_out[1]);
  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], kimData);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mxW_KINDenseJac(sunindextype Neq,
                    N_Vector y, N_Vector fy, 
                    DlsMat J, void *user_data,
                    N_Vector tmp1, N_Vector tmp2)
{
  kimInterfaceData kimData;
  double *J_data;
  mxArray *mx_in[4], *mx_out[3];
  int i, ret;

  /* Extract global interface data from user-data */
  kimData = (kimInterfaceData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[2] = kimData->JACfct;                   /* matlab function handle */
  mx_in[3] = kimData->mtlb_data;                /* matlab user data */
  
  /* Call matlab wrapper */

  GetData(y, mxGetPr(mx_in[0]), N);
  GetData(fy, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(3,mx_out,4,mx_in,"kim_djac");

  J_data = mxGetPr(mx_out[0]);
  for (i=0;i<N;i++)  memcpy(DENSE_COL(J,i), J_data + i*N, N*sizeof(double));
  ret = (int)*mxGetPr(mx_out[1]);
 
  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], kimData);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}


int mxW_KINBandJac(sunindextype Neq, sunindextype mupper, sunindextype mlower,
                   N_Vector y, N_Vector fy, 
                   DlsMat J, void *user_data,
                   N_Vector tmp1, N_Vector tmp2)
{
  kimInterfaceData kimData;
  double *J_data;
  mxArray *mx_in[4], *mx_out[3];
  sunindextype eband, i;
  int ret;

  /* Extract global interface data from user-data */
  kimData = (kimInterfaceData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[2] = kimData->JACfct;                   /* matlab function handle */
  mx_in[3] = kimData->mtlb_data;                /* matlab user data */
  
  /* Call matlab wrapper */

  GetData(y, mxGetPr(mx_in[0]), N);
  GetData(fy, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(3,mx_out,4,mx_in,"kim_bjac");

  eband =  mupper + mlower + 1;
  J_data = mxGetPr(mx_out[0]);
  for (i=0;i<N;i++) memcpy(BAND_COL(J,i) - mupper, J_data + i*eband, eband*sizeof(double));
  ret = (int)*mxGetPr(mx_out[1]);
 
  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], kimData);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mxW_KINSpilsJac(N_Vector v, N_Vector Jv,
                     N_Vector y, booleantype *new_y, 
                     void *user_data)
{
  kimInterfaceData kimData;
  mxArray *mx_in[5], *mx_out[4];
  int ret;

  /* Extract global interface data from user-data */
  kimData = (kimInterfaceData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* vector v */
  mx_in[2] = mxCreateLogicalScalar(*new_y);     /* */
  mx_in[3] = kimData->JACfct;                   /* matlab function handle */
  mx_in[4] = kimData->mtlb_data;                /* matlab user data */
 
  /* Call matlab wrapper */

  GetData(y, mxGetPr(mx_in[0]), N);
  GetData(v, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(4,mx_out,5,mx_in,"kim_jtv");

  PutData(Jv, mxGetPr(mx_out[0]), N);
  *new_y = mxIsLogicalScalarTrue(mx_out[1]);
  ret = (int)*mxGetPr(mx_out[2]);

  if (!mxIsEmpty(mx_out[3])) {
    UpdateUserData(mx_out[3], kimData);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);
  mxDestroyArray(mx_out[3]);

  return(ret);
}

int mxW_KINSpilsPset(N_Vector y, N_Vector yscale,
                     N_Vector fy, N_Vector fscale,
                     void *user_data, N_Vector vtemp1,
                     N_Vector vtemp2)
{
  kimInterfaceData kimData;
  mxArray *mx_in[6], *mx_out[2];
  int ret;

  /* Extract global interface data from user-data */
  kimData = (kimInterfaceData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yscale */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fscale */
  mx_in[4] = kimData->PSETfct;                  /* matlab function handle */
  mx_in[5] = kimData->mtlb_data;                /* matlab user data */
  
  /* Call matlab wrapper */

  GetData(y,      mxGetPr(mx_in[0]), N);
  GetData(yscale, mxGetPr(mx_in[1]), N);
  GetData(fy,     mxGetPr(mx_in[2]), N);
  GetData(fscale, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(2,mx_out,6,mx_in,"kim_pset");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1], kimData);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);

  return(ret);
}

int mxW_KINSpilsPsol(N_Vector y, N_Vector yscale, 
                     N_Vector fy, N_Vector fscale, 
                     N_Vector v, void *user_data,
                     N_Vector vtemp)
{
  kimInterfaceData kimData;
  mxArray *mx_in[7], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  kimData = (kimInterfaceData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yscale */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current fy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL); /* current fscale */
  mx_in[4] = mxCreateDoubleMatrix(N,1,mxREAL); /* right hand side */
  mx_in[5] = kimData->PSOLfct;                 /* matlab function handle */
  mx_in[6] = kimData->mtlb_data;               /* matlab user data */
  
  /* Call matlab wrapper */

  GetData(y,      mxGetPr(mx_in[0]), N);
  GetData(yscale, mxGetPr(mx_in[1]), N);
  GetData(fy,     mxGetPr(mx_in[2]), N);
  GetData(fscale, mxGetPr(mx_in[3]), N);
  GetData(v,      mxGetPr(mx_in[4]), N);

  mexCallMATLAB(3,mx_out,7,mx_in,"kim_psol");

  PutData(v, mxGetPr(mx_out[0]), N);
  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], kimData);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_in[2]);
  mxDestroyArray(mx_in[3]);
  mxDestroyArray(mx_in[4]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mxW_KINGloc(sunindextype Nlocal, N_Vector y, N_Vector gval, void *user_data)
{
  kimInterfaceData kimData;
  mxArray *mx_in[3], *mx_out[3];
  int ret;

  /* Extract global interface data from user-data */
  kimData = (kimInterfaceData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = kimData->GLOCfct;                  /* matlab function handle */
  mx_in[2] = kimData->mtlb_data;                /* matlab user data */
  
  /* Call matlab wrapper */

  GetData(y, mxGetPr(mx_in[0]), N);
  mexCallMATLAB(3,mx_out,3,mx_in,"kim_gloc");

  PutData(gval, mxGetPr(mx_out[0]), N);

  ret = (int)*mxGetPr(mx_out[1]);

  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2], kimData);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);
}

int mxW_KINGcom(sunindextype Nlocal, N_Vector y, void *user_data)
{
  kimInterfaceData kimData;
  mxArray *mx_in[5], *mx_out[2];
  int ret;

  /* Extract global interface data from user-data */
  kimData = (kimInterfaceData) user_data;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = kimData->GCOMfct;                  /* matlab function handle */
  mx_in[2] = kimData->mtlb_data;                /* matlab user data */
  
  /* Call matlab wrapper */

  GetData(y, mxGetPr(mx_in[0]), N);
  mexCallMATLAB(2,mx_out,3,mx_in,"kim_gcom");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1], kimData);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);

  return(ret);
}

/*
 * ---------------------------------------------------------------------------------
 * Private function to update the user data structure
 * ---------------------------------------------------------------------------------
 */

static void UpdateUserData(mxArray *new_mtlb_data, kimInterfaceData kimData)
{
  mexUnlock();
  mxDestroyArray(kimData->mtlb_data);
  kimData->mtlb_data = mxDuplicateArray(new_mtlb_data);
  mexMakeArrayPersistent(kimData->mtlb_data);
  mexLock();
}
