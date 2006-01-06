/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-01-06 19:00:23 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * KINSOL wrapper functions.
 * -----------------------------------------------------------------
 */

#include "kim.h"

static void UpdateUserData(mxArray *mx_data);

/*
 * ---------------------------------------------------------------------------------
 * Wrapper functions
 * ---------------------------------------------------------------------------------
 */

void mtlb_KINSys(N_Vector y, N_Vector fy, void *f_data )
{
  mxArray *mx_in[3], *mx_out[2];
  
  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[1] = mx_mtlb_SYSfct;                   /* matlab function handle */ 
  mx_in[2] = mx_mtlb_data;                     /* matlab user data */

  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[0]), N);
  mexCallMATLAB(2,mx_out,3,mx_in,"kim_sys");
  PutData(fy, mxGetPr(mx_out[0]), N);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);

}

int mtlb_KINDenseJac(long int N, DenseMat J, 
                     N_Vector y, N_Vector fy, void *jac_data,
                     N_Vector tmp1, N_Vector tmp2)
{
  double *J_data;
  long int i;
  mxArray *mx_in[4], *mx_out[3];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[2] = mx_mtlb_JACfct;                    /* matlab function handle */
  mx_in[3] = mx_mtlb_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[0]), N);
  GetData(fy, mxGetPr(mx_in[1]), N);
  mexCallMATLAB(3,mx_out,4,mx_in,"kim_djac");

  /* Extract data */
  J_data = mxGetPr(mx_out[0]);
  for (i=0;i<N;i++)
    memcpy(DENSE_COL(J,i), J_data + i*N, N*sizeof(double));

  ret = (int)*mxGetPr(mx_out[1]);
 
  if (!mxIsEmpty(mx_out[2])) {
    UpdateUserData(mx_out[2]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_in[1]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
  mxDestroyArray(mx_out[2]);

  return(ret);

}


int mtlb_KINSpilsJac(N_Vector v, N_Vector Jv,
                     N_Vector y, booleantype *new_y, 
                     void *J_data)
{
  mxArray *mx_in[5], *mx_out[4];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* vector v */
  mx_in[2] = mxCreateLogicalScalar(*new_y);     /* */
  mx_in[3] = mx_mtlb_JACfct;                    /* matlab function handle */
  mx_in[4] = mx_mtlb_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[0]), N);
  GetData(v, mxGetPr(mx_in[1]), N);

  mexCallMATLAB(4,mx_out,5,mx_in,"kim_jtv");

  PutData(Jv, mxGetPr(mx_out[0]), N);
  *new_y = mxIsLogicalScalarTrue(mx_out[1]);
  ret = (int)*mxGetPr(mx_out[2]);

  if (!mxIsEmpty(mx_out[3])) {
    UpdateUserData(mx_out[3]);
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

int mtlb_KINSpilsPset(N_Vector y, N_Vector yscale,
                      N_Vector fy, N_Vector fscale,
                      void *P_data, N_Vector vtemp1,
                      N_Vector vtemp2)
{
  mxArray *mx_in[6], *mx_out[2];
  int ret;

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current yscale */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current fscale */
  mx_in[4] = mx_mtlb_PSETfct;                   /* matlab function handle */
  mx_in[5] = mx_mtlb_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y,      mxGetPr(mx_in[0]), N);
  GetData(yscale, mxGetPr(mx_in[1]), N);
  GetData(fy,     mxGetPr(mx_in[2]), N);
  GetData(fscale, mxGetPr(mx_in[3]), N);

  mexCallMATLAB(2,mx_out,6,mx_in,"kim_pset");

  ret = (int)*mxGetPr(mx_out[0]);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1]);
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

int mtlb_KINSpilsPsol(N_Vector y, N_Vector yscale, 
                      N_Vector fy, N_Vector fscale, 
                      N_Vector v, void *P_data,
                      N_Vector vtemp)
{
  mxArray *mx_in[7], *mx_out[3];
  int ret;


  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL); /* current y */
  mx_in[1] = mxCreateDoubleMatrix(N,1,mxREAL); /* current yscale */
  mx_in[2] = mxCreateDoubleMatrix(N,1,mxREAL); /* current fy */
  mx_in[3] = mxCreateDoubleMatrix(N,1,mxREAL); /* current fscale */
  mx_in[4] = mxCreateDoubleMatrix(N,1,mxREAL); /* right hand side */
  mx_in[5] = mx_mtlb_PSOLfct;                  /* matlab function handle */
  mx_in[6] = mx_mtlb_data;                     /* matlab user data */
  
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
    UpdateUserData(mx_out[2]);
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

void mtlb_KINGloc(long int Nlocal, N_Vector y, N_Vector gval, void *f_data)
{
  mxArray *mx_in[3], *mx_out[2];

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mx_mtlb_GLOCfct;                   /* matlab function handle */
  mx_in[2] = mx_mtlb_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[0]), N);
  mexCallMATLAB(2,mx_out,3,mx_in,"kim_gloc");
  PutData(gval, mxGetPr(mx_out[0]), N);

  if (!mxIsEmpty(mx_out[1])) {
    UpdateUserData(mx_out[1]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_out[0]);
  mxDestroyArray(mx_out[1]);
}

void mtlb_KINGcom(long int Nlocal, N_Vector y, void *f_data)
{
  mxArray *mx_in[5], *mx_out[1];

  /* Inputs to the Matlab function */
  mx_in[0] = mxCreateDoubleMatrix(N,1,mxREAL);  /* current y */
  mx_in[1] = mx_mtlb_GCOMfct;                   /* matlab function handle */
  mx_in[2] = mx_mtlb_data;                      /* matlab user data */
  
  /* Call matlab wrapper */
  GetData(y, mxGetPr(mx_in[0]), N);
  mexCallMATLAB(1,mx_out,3,mx_in,"kim_gcom");

  if (!mxIsEmpty(mx_out[0])) {
    UpdateUserData(mx_out[0]);
  }

  /* Free temporary space */
  mxDestroyArray(mx_in[0]);
  mxDestroyArray(mx_out[0]);
}

/*
 * ---------------------------------------------------------------------------------
 * Private function to update the user data structure
 * ---------------------------------------------------------------------------------
 */

static void UpdateUserData(mxArray *mx_data)
{
  mexUnlock();
  mxDestroyArray(mx_mtlb_data);
  mx_mtlb_data = mxDuplicateArray(mx_data);
  mexMakeArrayPersistent(mx_mtlb_data);
  mexLock();
}
