/*
 * -----------------------------------------------------------------
 * $Revision: 1.8 $
 * $Date: 2007-05-11 18:51:32 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials-x.y.z/src/cvodes/LICENSE.
 * -----------------------------------------------------------------
 * Header file for the CVODES Matlab interface.
 * -----------------------------------------------------------------
 */

#ifndef _CVM_H
#define _CVM_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <string.h>
#include <stdlib.h>
#include "mex.h"
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes_bandpre.h>
#include <cvodes/cvodes_bbdpre.h>

/*
 * ---------------------------------------------------------------------------------
 * Constants
 * ---------------------------------------------------------------------------------
 */

/* Tolerance types */

enum {CV_SS, CV_SV, CV_EE};

/* Linear solver types */

enum {LS_NONE, LS_DENSE, LS_DIAG, LS_BAND, LS_SPGMR, LS_SPBCG, LS_SPTFQMR};

/* Preconditioner modules */

enum {PM_NONE, PM_BANDPRE, PM_BBDPRE};

/*
 * ---------------------------------------------------------------------------------
 * Types for global data structures
 * ---------------------------------------------------------------------------------
 */


typedef struct cvmPbData_ {


  int n;             /* problem dimension */
  N_Vector Y;        /* solution vector */

  booleantype Quadr; /* integrate quadratures? */
  int nq;            /* number of quadratures */
  N_Vector YQ;       /* quadratures vector */

  booleantype Fsa;   /* integrate sensitivities? */
  int ns;            /* number of sensitivities */
  N_Vector *YS;      /* sensitivity vectors */

  booleantype RootSet;  /* rootfinding active? */
  int ng;               /* number of root functions */

  booleantype TstopSet; /* tstop active? */

  int LS;            /* linear solver type */
  int PM;            /* preconditioner module */

  booleantype Mon;   /* monitoring? */

  /* Matlab functions and data associated with this problem */

  mxArray *RHSfct;
  mxArray *QUADfct;

  mxArray *JACfct;

  mxArray *PSETfct;
  mxArray *PSOLfct;

  mxArray *GLOCfct;
  mxArray *GCOMfct;
    
  mxArray *Gfct;
    
  mxArray *SRHSfct;

  mxArray *MONfct;
  mxArray *MONdata;

  /* Information for lists of problems (backward problems only) */

  int index;               /* index of this problem */
  struct cvmPbData_ *next; /* pointer to next problem in linked list */

} *cvmPbData;


typedef struct cvmInterfaceData_ {

  void *cvode_mem;    /* CVODES solver memory */

  booleantype asa;    /* Perform ASA? */
  int Nd;             /* number of data points */
  int Nc;             /* number of check points */

  struct cvmPbData_ *fwdPb;
  struct cvmPbData_ *bckPb;

  int NbckPb;         /* Number of backward problems in the linked list bckPb */

  mxArray *mx_data;    /* Matlab user data */

} *cvmInterfaceData;


/*
 * ---------------------------------------------------------------------------------
 * Error handler function
 * ---------------------------------------------------------------------------------
 */

void cvmErrHandler(int error_code, 
                   const char *module, const char *function, 
                   char *msg, void *f_data);

/*
 * ---------------------------------------------------------------------------------
 * Wrapper functions
 * ---------------------------------------------------------------------------------
 */

int mtlb_CVodeRhs(realtype t, N_Vector y, N_Vector yd, void *f_data);

int mtlb_CVodeGfct(realtype t, N_Vector y, double *g, void *g_data);

int mtlb_CVodeQUADfct(realtype t, N_Vector y, N_Vector yQd, void *fQ_data);


int mtlb_CVodeSensRhs1(int Ns, realtype t,
                       N_Vector y, N_Vector ydot,
                       int iS, N_Vector yS, N_Vector ySdot,
                       void *fS_data,
                       N_Vector tmp1, N_Vector tmp2);
int mtlb_CVodeSensRhs(int Ns, realtype t,
                      N_Vector y, N_Vector ydot,
                      N_Vector *yS, N_Vector *ySdot,
                      void *fS_data,
                      N_Vector tmp1, N_Vector tmp2);


int mtlb_CVodeDenseJac(int N, realtype t,
                       N_Vector y, N_Vector fy, 
                       DlsMat J, void *jac_data,
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int mtlb_CVodeBandJac(int N, int mupper, int mlower, realtype t,
                      N_Vector y, N_Vector fy, 
                      DlsMat J, void *jac_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int mtlb_CVodeSpilsJac(N_Vector v, N_Vector Jv, realtype t,
                       N_Vector y, N_Vector fy,
                       void *jac_data, N_Vector tmp);
int mtlb_CVodeSpilsPset(realtype t, N_Vector y, N_Vector fy,
                        booleantype jok, booleantype *jcurPtr,
                        realtype gamma, void *P_data,
                        N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3);
int mtlb_CVodeSpilsPsol(realtype t, N_Vector y, N_Vector fy,
                        N_Vector r, N_Vector z,
                        realtype gamma, realtype delta,
                        int lr, void *P_data, N_Vector tmp);

  
int mtlb_CVodeBBDgloc(int Nlocal, realtype t, N_Vector y, N_Vector g, void *f_data);
int mtlb_CVodeBBDgcom(int Nlocal, realtype t, N_Vector y, void *f_data);



int mtlb_CVodeRhsB(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, void *f_dataB);
int mtlb_CVodeQUADfctB(realtype t, N_Vector y, N_Vector yB, N_Vector qBdot, void *fQ_dataB);
int mtlb_CVodeDenseJacB(int nB, realtype t,
                        N_Vector y, N_Vector yB, N_Vector fyB,
                        DlsMat JB, void *jac_dataB, 
                        N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
int mtlb_CVodeBandJacB(int nB, int mupperB, int mlowerB, realtype t,
                       N_Vector y, N_Vector yB, N_Vector fyB,
                       DlsMat JB, void *jac_dataB, 
                       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
int mtlb_CVodeSpilsJacB(N_Vector vB, N_Vector JvB, realtype t,
                        N_Vector y, N_Vector yB, N_Vector fyB,
                        void *jac_dataB, N_Vector tmpB);
int mtlb_CVodeSpilsPsetB(realtype t, N_Vector y,
                         N_Vector yB, N_Vector fyB,
                         booleantype jokB,
                         booleantype *jcurPtrB, realtype gammaB,
                         void *P_dataB,
                         N_Vector tmp1B, N_Vector tmp2B,
                         N_Vector tmp3B);
int mtlb_CVodeSpilsPsolB(realtype t, N_Vector y,
                         N_Vector yB, N_Vector fyB,
                         N_Vector rB, N_Vector zB,
                         realtype gammaB, realtype deltaB,
                         int lrB, void *P_dataB, N_Vector tmpB);
  
int mtlb_CVodeBBDglocB(int NlocalB, realtype t, N_Vector y, 
                       N_Vector yB, N_Vector gB, void *f_dataB);

int mtlb_CVodeBBDgcomB(int NlocalB, realtype t, N_Vector y, 
                       N_Vector yB, void *f_dataB);

void mtlb_CVodeMonitor(int call, double t, 
                       N_Vector y, N_Vector yQ, N_Vector *yS,
                       cvmInterfaceData cvm_Cdata);
void mtlb_CVodeMonitorB(int call, int idxB, double tB, N_Vector yB, N_Vector yQB,
                        cvmInterfaceData cvm_Cdata);

/*
 * ---------------------------------------------------------------------------------
 * Option handling functions
 * ---------------------------------------------------------------------------------
 */

void get_IntgrOptions(const mxArray *options, cvmPbData thisPb, booleantype fwd,
                      int *lmm, int *iter, int *maxord, booleantype *sld,
                      long int *mxsteps,
                      int *itol, realtype *reltol, double *Sabstol, double **Vabstol,
                      double *hin, double *hmax, double *hmin, 
                      double *tstop);

void get_LinSolvOptions(const mxArray *options, cvmPbData thisPb, booleantype fwd,
                        int *mupper, int *mlower,
                        int *mudq, int *mldq, double *dqrely,
                        int *ptype, int *gstype, int *maxl);

void get_QuadOptions(const mxArray *options, cvmPbData thisPb, booleantype fwd,
                     int Nq,
                     booleantype *errconQ,
                     int *itolQ, double *reltolQ, double *SabstolQ, double **VabstolQ);

void get_FSAOptions(const mxArray *options, cvmPbData thisPb,
                    int *ism,
                    char **pfield_name, int **plist, double **pbar,
                    int *dqtype, double *rho,
                    booleantype *errconS, int *itolS, double *reltolS, 
                    double **SabstolS, double **VabstolS);

#ifdef __cplusplus
}
#endif

#endif
