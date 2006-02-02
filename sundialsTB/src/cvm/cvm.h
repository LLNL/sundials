/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-02-02 00:39:04 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvodes/LICENSE.
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
#include "cvodes.h"
#include "cvodea.h"
#include "cvodes_dense.h"
#include "cvodes_diag.h"
#include "cvodes_band.h"
#include "cvodes_spgmr.h"
#include "cvodes_spbcgs.h"
#include "cvodes_bandpre.h"
#include "cvodes_bbdpre.h"
#include "sundials_nvector.h"

  /*
   * ---------------------------------------------------------------------------------
   * Constants
   * ---------------------------------------------------------------------------------
   */


  /* Linear solver types */

  enum {LS_NONE, LS_DENSE, LS_DIAG, LS_BAND, LS_SPGMR, LS_SPBCG, LS_SPTFQMR};

  /* Preconditioner modules */

  enum {PM_NONE, PM_BANDPRE, PM_BBDPRE};

  /*
   * ---------------------------------------------------------------------------------
   * Declarations for global variables
   * ---------------------------------------------------------------------------------
   */

  /* CVODE data */

  extern void *cvode_mem;   /* CVODES solver memory */
  extern void *bp_data;     /* Preconditioner memory (BandPre or BBDPre) */
  extern N_Vector y;        /* solution vector */
  extern N_Vector yQ;       /* quadratures vector */
  extern N_Vector *yS;      /* sensitivity vectors */
  extern int N;             /* problem dimension */
  extern int Nq;            /* number of quadratures */
  extern int Ng;            /* number of root functions */
  extern int Ns;            /* number of sensitivities */
  extern int Nd;            /* number of data points */
  extern int Nc;            /* number of check points */
  extern int ls;            /* linear solver type */
  extern int pm;            /* preconditioner module */
  extern int ism;           /* sensitivity method */

  extern void *cvadj_mem;   /* CVODES adjoint memory */
  extern int interp;
  extern N_Vector yB;
  extern N_Vector yQB;
  extern int NB;
  extern int NqB;
  extern int lsB;
  extern int pmB;

  /* Matlab data */

  extern mxArray *mx_mtlb_RHSfct;
  extern mxArray *mx_mtlb_QUADfct;
  extern mxArray *mx_mtlb_JACfct;
  extern mxArray *mx_mtlb_PSETfct;
  extern mxArray *mx_mtlb_PSOLfct;
  extern mxArray *mx_mtlb_GLOCfct;
  extern mxArray *mx_mtlb_GCOMfct;

  extern mxArray *mx_mtlb_Gfct;

  extern mxArray *mx_mtlb_SRHSfct;

  extern mxArray *mx_mtlb_RHSfctB;
  extern mxArray *mx_mtlb_QUADfctB;
  extern mxArray *mx_mtlb_JACfctB;
  extern mxArray *mx_mtlb_PSETfctB;
  extern mxArray *mx_mtlb_PSOLfctB;
  extern mxArray *mx_mtlb_GLOCfctB;
  extern mxArray *mx_mtlb_GCOMfctB;

  extern mxArray *mx_mtlb_data;

  /* Monitor */
  
  extern booleantype monitor;
  extern mxArray *mx_mtlb_MONfct;
  extern mxArray *mx_mtlb_MONdata;

  extern booleantype monitorB;
  extern mxArray *mx_mtlb_MONfctB;
  extern mxArray *mx_mtlb_MONdataB;

  /*
   * ---------------------------------------------------------------------------------
   * Wrapper functions
   * ---------------------------------------------------------------------------------
   */

  void mtlb_CVodeErrHandler(int error_code, 
                            const char *module, const char *function, 
                            char *msg, void *eh_data); 
  
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


  int mtlb_CVodeDenseJac(long int N, DenseMat J, realtype t,
                         N_Vector y, N_Vector fy, void *jac_data,
                         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  int mtlb_CVodeBandJac(long int N, long int mupper, long int mlower,
                        BandMat J, realtype t,
                        N_Vector y, N_Vector fy, void *jac_data,
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

  
  int mtlb_CVodeBBDgloc(long int Nlocal, realtype t, N_Vector y,
                        N_Vector g, void *f_data);
  int mtlb_CVodeBBDgcom(long int Nlocal, realtype t, N_Vector y,
                        void *f_data);



  int mtlb_CVodeRhsB(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, void *f_dataB);
  int mtlb_CVodeQUADfctB(realtype t, N_Vector y, N_Vector yB, N_Vector qBdot, void *fQ_dataB);
  int mtlb_CVodeDenseJacB(long int nB, DenseMat JB, realtype t,
                          N_Vector y, N_Vector yB, N_Vector fyB,
                          void *jac_dataB, N_Vector tmp1B,
                          N_Vector tmp2B, N_Vector tmp3B);
  int mtlb_CVodeBandJacB(long int nB, long int mupperB,
                         long int mlowerB, BandMat JB,
                         realtype t, N_Vector y,
                         N_Vector yB, N_Vector fyB,
                         void *jac_dataB, N_Vector tmp1B,
                         N_Vector tmp2B, N_Vector tmp3B);
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
  
  int mtlb_CVodeBBDglocB(long int NlocalB, realtype t, N_Vector y, 
                          N_Vector yB, N_Vector gB, void *f_dataB);

  int mtlb_CVodeBBDgcomB(long int NlocalB, realtype t, N_Vector y, 
                          N_Vector yB, void *f_dataB);

  void mtlb_CVodeMonitor(int call, double t, N_Vector y, N_Vector yQ, N_Vector *yS);
  void mtlb_CVodeMonitorB(int call, double tB, N_Vector yB, N_Vector yQB);

  /*
   * ---------------------------------------------------------------------------------
   * Option handling functions
   * ---------------------------------------------------------------------------------
   */

  int get_IntgrOptions(const mxArray *options,
                       int *lmm, int *iter, int *maxord, booleantype *sld,
                       long int *mxsteps,
                       int *itol, realtype *reltol, double *Sabstol, double **Vabstol,
                       double *hin, double *hmax, double *hmin, double *tstop, booleantype *tstopSet,
                       int *Ng_tmp, mxArray **mx_tmp_Gfct,
                       booleantype *quad, booleantype *fsa, booleantype *asa,
                       booleantype *monitor_tmp, mxArray **mx_tmp_MONfct,
                       mxArray **mx_tmp_MONdata);

  int get_LinSolvOptions(const mxArray *options,
                         int *ls_tmp,
                         int *mupper, int *mlower,
                         int *mudq, int *mldq, double *dqrely,
                         int *ptype, int *gstype, int *maxl, int *pm_tmp,
                         mxArray **mx_tmp_JACfct,
                         mxArray **mx_tmp_PSETfct, mxArray **mx_tmp_PSOLfct,
                         mxArray **mx_tmp_GLOCfct, mxArray **mx_tmp_GCOMfct);
  
  int get_QuadOptions(const mxArray *options,
                      int *Nq_tmp, double **yQ0, mxArray **mx_tmp_QUADfct,
                      booleantype *errconQ,
                      int *itolQ, double *reltolQ, double *SabstolQ, double **VabstolQ);

  int get_FSAOptions(const mxArray *options, 
                     int *Ns_tmp, double **yS0,
                     int *ism_tmp,
                     char **pfield_name, int **plist, double **pbar,
                     int *Srhs, mxArray **mx_tmp_SRHSfct, double *rho,
                     booleantype *errconS, int *itolS, double *reltolS, 
                     double **SabstolS, double **VabstolS);

  int get_ASAOptions(const mxArray *options, int *Nd_tmp, int *interp_tmp);


#ifdef __cplusplus
}
#endif

#endif
