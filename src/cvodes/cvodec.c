/*******************************************************************
 *                                                                 *
 * File          : cvodec.cpp                                      *
 * Programmers   : Radu Serban @ LLNL                              *
 * Version of    : 07 February 2004                                *
 *-----------------------------------------------------------------*
 * Copyright (c) 2006, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see the LICENSE file                        *
 *-----------------------------------------------------------------*
 * This is a set of wrappers for derivative calculations in CVODES *
 * using the complex step method.                                  *
 *                                                                 *
 *******************************************************************/

/*=================================================================*/
/*BEGIN        Import Header Files                                 */
/*=================================================================*/

#include <cvodec.h>

/*=================================================================*/
/*END          Import Header Files                                 */
/*=================================================================*/

#define ZERO   RCONST(0.0)     /* real 0.0 */
#define ONE    RCONST(1.0)     /* real 1.0 */

/*=================================================================*/
/*BEGIN        CVODEC Error Messages                               */
/*=================================================================*/

#define MSG_NO_MEM        "CVode Complex Step-- CVODES memory NULL. \n\n"
#define MSG_MEM_FAIL      "CVodeSetCSDerivs-- memory allocation failed. \n\n"
#define MSG_FCS_NULL      "CVodeSetCSDerivs-- f_cs = NULL illegal. \n\n"
#define MSG_Step_NO_CSMEM "CVodeSetCSStep-- attempt to call before CVodeSetCSDerivs. \n\n"
#define MSG_BAD_DEL       "CVodeSetCSStep-- del_cs = %g is below uround = %g \n\n"
#define MSG_Sens_NO_CSMEM "CVodeSetSensCSRhs-- attempt to call before CVodeSetCSDerivs. \n\n"
#define MSG_PIM_NULL      "CVodeSetSensCSRhs-- p_im = NULL, illeagal. \n\n"
#define MSG_NO_DENSEMEM   "CVDenseSetCSJac-- attempt to call before CVDense. \n\n"
#define MSG_NO_BANDMEM    "CVBandSetCSJac-- attempt to call before CVBand. \n\n"
#define MSG_NO_SPGMRMEM   "CVSpgmrSetCSJac-- attempt to call before CVSpgmr. \n\n"

/*=================================================================*/
/*END          CVODEC Error Messages                               */
/*=================================================================*/

/* Readibility Replacements */
#define errfp  (cv_mem->cv_errfp)
#define uround (cv_mem->cv_uround)
#define lmem   (cv_mem->cv_lmem)
#define csmem  (cv_mem->cv_csmem)
#define nvspec (cv_mem->cv_nvspec)

/*=================================================================*/
/*BEGIN        EXPORTED FUNCTIONS IMPLEMENTATION                   */
/*=================================================================*/

/*------------------  CVodeSetCSDerivs     ------------------------*/
/*
  CVodeSetCSDerivs enables derivative calculations using complex
  step. 
  The user must pass the complexified right hand side routine and
  the pointer to the 'imaginary' user data structure.
  CVodeSetCSDerivs allocates on N_Vector for temporary workspace.
  If successful, it returns SUCCESS. Otherwise it returns an error
  flag and prints an error message.
*/
/*-----------------------------------------------------------------*/

int CVodeSetCSDerivs(void *cvode_mem, void *f_cs, void *f_data_im)
{
  CVodeMem cv_mem;
  CVCSMem  cvcs_mem;

  if (cvode_mem == NULL) {
    fprintf(stderr, MSG_NO_MEM);
    return(CVCS_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  cvcs_mem = (CVCSMem) malloc(sizeof(CVCSMemRec));
  if (cvcs_mem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_MEM_FAIL);
    return(CVCS_MEM_FAIL);
  }

  if (f_cs == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_FCS_NULL);
    return(CVCS_ILL_INPUT);
  }

  if (f_data_im == NULL) {
    cvcs_mem->cvcs_type1  = TRUE;
    cvcs_mem->cvcs_f_cs1  = (RhsCSFn1) f_cs;
    cvcs_mem->cvcs_f_cs2  = NULL;
  } else {
    cvcs_mem->cvcs_type1  = FALSE;
    cvcs_mem->cvcs_f_cs1  = NULL;
    cvcs_mem->cvcs_f_cs2  = (RhsCSFn2) f_cs;
  }

  cvcs_mem->cvcs_f_data_im = f_data_im;

  cvcs_mem->cvcs_del_cs    = uround;

  cvcs_mem->cvcs_p_im      = NULL;
  
  /* attach cs memory to main CVODES memory */
  csmem = cvcs_mem;

  return(SUCCESS);
  
}

/*--------------------  CVodeSetCSstep     ------------------------*/
/*
  Sets the perturbation step used in complex step calculations.
  The default value is uround.
*/
/*-----------------------------------------------------------------*/

int CVodeSetCSStep(void *cvode_mem, realtype del_cs)
{
  CVodeMem cv_mem;
  CVCSMem cvcs_mem;

  if (cvode_mem == NULL) {
    fprintf(stderr, MSG_NO_MEM);
    return(CVCS_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (csmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_Step_NO_CSMEM);
    return(CVCS_NO_CSMEM);
  }

  cvcs_mem = (CVCSMem) csmem;

  if (del_cs < uround) {
    if(errfp!=NULL) fprintf(errfp, MSG_BAD_DEL, del_cs, uround);
    return(CVCS_ILL_INPUT);
  }

  cvcs_mem->cvcs_del_cs = del_cs;

  return(SUCCESS);
}

/*------------------   CVodeSetSensCSRhs   ------------------------*/
/*
  Enables complex step calculation of sensitivity right hand sides
*/
/*-----------------------------------------------------------------*/

int CVodeSetSensCSRhs(void *cvode_mem, realtype *p_im)
{
  CVodeMem cv_mem;
  CVCSMem cvcs_mem;

  if (cvode_mem == NULL) {
    fprintf(stderr, MSG_NO_MEM);
    return(CVCS_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (csmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_Sens_NO_CSMEM);
    return(CVCS_NO_CSMEM);
  }

  cvcs_mem = (CVCSMem) csmem;

  if (p_im == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_PIM_NULL);
    return(CVCS_ILL_INPUT);
  }

  cv_mem->cv_ifS      = ONESENS;
  cv_mem->cv_fSDQ     = FALSE;
  cv_mem->cv_fS1      = CVSensRhsCS;
  cv_mem->cv_fS_data  = cvode_mem;

  cvcs_mem->cvcs_p_im = p_im;

  return(SUCCESS);

}

/*--------------------   CVDenseSetCSJac   ------------------------*/
/*
  Enables complex step calculation of the dense Jacobian
*/
/*-----------------------------------------------------------------*/

int CVDenseSetCSJac(void *cvode_mem)
{
  CVodeMem cv_mem;
  CVDenseMem cvdense_mem;

  if (cvode_mem == NULL) {
    fprintf(stderr, MSG_NO_MEM);
    return(CVCS_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_NO_DENSEMEM);
    return(CVCS_NO_LMEM);
  }

  cvdense_mem = (CVDenseMem) lmem;

  cvdense_mem->d_jac = CVDenseCSJac;
  cvdense_mem->d_J_data = cvode_mem;

  return(SUCCESS);
}

/*---------------------   CVBandSetCSJac   ------------------------*/
/*
  Enables complex step calculation of the band Jacobian
*/
/*-----------------------------------------------------------------*/

int CVBandSetCSJac(void *cvode_mem)
{
  CVodeMem cv_mem;
  CVBandMem cvband_mem;

  if (cvode_mem == NULL) {
    fprintf(stderr, MSG_NO_MEM);
    return(CVCS_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_NO_BANDMEM);
    return(CVCS_NO_LMEM);
  }

  cvband_mem = (CVBandMem) lmem;

  cvband_mem->b_jac = CVBandCSJac;
  cvband_mem->b_J_data = cvode_mem;

  return(SUCCESS);
}

/*-------------   CVSpgmrSetCSJacTimesVec   -----------------------*/
/*
  Enables complex step calculation of Jacobian times Vector
*/
/*-----------------------------------------------------------------*/

int CVSpgmrSetCSJacTimesVec(void *cvode_mem)
{
  CVodeMem cv_mem;
  CVSpgmrMem cvspgmr_mem;

  if (cvode_mem == NULL) {
    fprintf(stderr, MSG_NO_MEM);
    return(CVCS_NO_MEM);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    if(errfp!=NULL) fprintf(errfp, MSG_NO_SPGMRMEM);
    return(CVCS_NO_LMEM);
  }

  cvspgmr_mem = (CVSpgmrMem) lmem;

  cvspgmr_mem->g_jtimes = CVSpgmrCSJtimes;
  cvspgmr_mem->g_j_data = cvode_mem;

  return(SUCCESS);
}

/*=================================================================*/
/*END          EXPORTED FUNCTIONS IMPLEMENTATION                   */
/*=================================================================*/


/*=================================================================*/
/*BEGIN        CS Approximation Routines                           */
/*=================================================================*/

/* Readibility Replacements */
#define f_data_re (cv_mem->cv_f_data)
#define plist     (cv_mem->cv_plist)

#define type1     (cvcs_mem->cvcs_type1)
#define f_cs1     (cvcs_mem->cvcs_f_cs1)
#define f_cs2     (cvcs_mem->cvcs_f_cs2)
#define f_data_im (cvcs_mem->cvcs_f_data_im)
#define del_cs    (cvcs_mem->cvcs_del_cs)
#define p_im      (cvcs_mem->cvcs_p_im)


/*-----------------   CVSensRhsCS   -------------------------------*/
/*
  CVSensRhsCS computes the right hand side of the iS-th sensitivity
  equation using complex step.

  df/dp_j = (1/del) Im[ f(y+i*h*yS_j , p+i*h*e_j)  ] 
*/
/*-----------------------------------------------------------------*/

void CVSensRhsCS(int Ns, realtype t, 
                 N_Vector y, N_Vector ydot, 
                 int iS, N_Vector yS, N_Vector ySdot, 
                 void *fS_data,
                 N_Vector tmp1, N_Vector tmp2)
{
  CVodeMem cv_mem;
  CVCSMem cvcs_mem;
  N_Vector y_re, y_im, f_re, f_im;
  int which;

  /* fS_data points to cvode_mem */
  cv_mem = (CVodeMem) fS_data;
  cvcs_mem = (CVCSMem) csmem;

  if (plist!=NULL)
    which   = abs(plist[iS]) - 1;
  else
    which  = iS;
 
  y_re = y;
  y_im = tmp1;
  f_re = tmp2;
  f_im = ySdot;

  N_VScale(del_cs, yS, y_im);

  p_im[which] = del_cs;

  if (type1)
    f_cs1(t, 0.0, y_re, y_im, f_re, f_im, f_data_re);
  else
    f_cs2(t, 0.0, y_re, y_im, f_re, f_im, f_data_re, f_data_im);

  N_VScale(ONE/del_cs, f_im, ySdot);

  p_im[which] = ZERO;

}

/*---------------------   CVDenseCSJac  ---------------------------*/
/*
  CVDenseCSJac computes the dense Jacobian J using complex step.

  J_j = (1/del) Im[ f(y+i*h*e_j) ]
*/
/*-----------------------------------------------------------------*/

void CVDenseCSJac(long int N, DenseMat J, realtype t, 
                  N_Vector y, N_Vector fy, void *jac_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CVodeMem cv_mem;
  CVCSMem cvcs_mem;
  N_Vector y_re, y_im, f_re, f_im, jthCol;
  realtype *ydata_im;
  long int j;

  /* jac_data points to cvode_mem */
  cv_mem = (CVodeMem) jac_data;
  cvcs_mem = (CVCSMem) csmem;

  y_re = y;
  y_im = tmp1;
  f_re = tmp2;
  f_im = tmp3;

  N_VConst(ZERO, y_im);
  ydata_im = (realtype *) N_VGetData(y_im);

  jthCol = N_VMake((void *)DENSE_COL(J,0), nvspec);

  for (j=0; j < N; j++) {
    N_VSetData((void *)DENSE_COL(J,j), jthCol);
    ydata_im[j] = del_cs;
    if (type1)
      f_cs1(t, 0.0, y_re, y_im, f_re, f_im, f_data_re);
    else
      f_cs2(t, 0.0, y_re, y_im, f_re, f_im, f_data_re, f_data_im);
    N_VScale(ONE/del_cs, f_im, jthCol);
    DENSE_COL(J,j) = (realtype *) N_VGetData(jthCol);
    ydata_im[j] = ZERO;
  }

  N_VDispose(jthCol);

}

/*-------------------   CVBandCSJac   -----------------------------*/
/*
  CVBandCSJac computes the band Jacobian approximation using
  complex step.
*/
/*-----------------------------------------------------------------*/

void CVBandCSJac(long int N, long int mupper, 
                 long int mlower, BandMat J, realtype t, 
                 N_Vector y, N_Vector fy, void *jac_data,
                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CVodeMem cv_mem;
  CVCSMem cvcs_mem;
  N_Vector y_re, y_im, f_re, f_im;
  long int group, i, j, width, ngroups, i1, i2;
  realtype *col_j, *ydata_im, *fdata_im;

  /* jac_data points to cvode_mem */
  cv_mem = (CVodeMem) jac_data;
  cvcs_mem = (CVCSMem) csmem;

  y_re = y;
  y_im = tmp1;
  f_re = tmp2;
  f_im = tmp3;

  N_VConst(ZERO, y_im);
  ydata_im = (realtype *) N_VGetData(y_im);

  /* Set bandwidth and number of column groups */
  width = mlower + mupper + 1;
  ngroups = MIN(width, N);

  for (group=1; group <= ngroups; group++) {

    for(j=group-1; j < N; j+=width)
      ydata_im[j] = del_cs;

    if (type1)
      f_cs1(t, 0.0, y_re, y_im, f_re, f_im, f_data_re);
    else
      f_cs2(t, 0.0, y_re, y_im, f_re, f_im, f_data_re, f_data_im);
    
    fdata_im = (realtype *) N_VGetData(f_im);

    for (j=group-1; j < N; j+=width) {
      ydata_im[j] = ZERO;
      col_j = BAND_COL(J,j);
      i1 = MAX(0, j-mupper);
      i2 = MIN(j+mlower, N-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) = fdata_im[i]/del_cs;
    }
  }

}

/*-----------------   CVSpgmrCSJtimes  ----------------------------*/
/*
  CVSpgmrCSJtimes computes the Jacobian x Vector product using
  complex step.
  
  J*v = (1/del) Im[ f(y+i*h*v) ]
*/
/*-----------------------------------------------------------------*/

int CVSpgmrCSJtimes(N_Vector v, N_Vector Jv, realtype t, 
                    N_Vector y, N_Vector fy,
                    void *jac_data, N_Vector work)
{
  CVodeMem cv_mem;
  CVCSMem cvcs_mem;
  N_Vector y_re, y_im, f_re, f_im;

  /* jac_data points to cvode_mem */
  cv_mem = (CVodeMem) jac_data;
  cvcs_mem = (CVCSMem) csmem;

  y_re = y;
  y_im = v;
  f_re = work;
  f_im = Jv;

  /* Set y_im = h_sc * v */
  N_VScale(del_cs, v, y_im);

  /* Evaluate right hand side in complex arithmetic */
  /* f_im =  f( y+i*h*v) */
  if (type1)
    f_cs1(t, 0.0, y_re, y_im, f_re, f_im, f_data_re);
  else
    f_cs2(t, 0.0, y_re, y_im, f_re, f_im, f_data_re, f_data_im);

  /* Compute Jv = f_im / h_sc */
  N_VScale(ONE/del_cs, f_im, Jv);

  /* Undo scale of v  */
  N_VScale(ONE/del_cs, v, v);

  return(SUCCESS);

}

/*=================================================================*/
/*END          CS Approximation Routines                           */
/*=================================================================*/
